
/*--------------------------------------------------------------------*/
/*--- Cache simulation                                    cg_sim.c ---*/
/*--------------------------------------------------------------------*/

/*
   This file is part of Cachegrind, a high-precision tracing profiler
   built with Valgrind.

   Copyright (C) 2002-2017 Nicholas Nethercote
      njn@valgrind.org

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.

   The GNU General Public License is contained in the file COPYING.
*/

/* Notes:
  - simulates a write-allocate cache
  - (block --> set) hash function uses simple bit selection
  - handling of references straddling two cache blocks:
      - counts as only one cache access (not two)
      - both blocks hit                  --> one hit
      - one block hits, the other misses --> one miss
      - both blocks miss                 --> one miss (not two)
*/

/*------------------------------------------------------------*/
/*--- Types and Data Structures                            ---*/
/*------------------------------------------------------------*/
/* WordSize should be a power of 2, and <=1 WordSize <=64*/
/* Users can set word_size using command opt.*/
/* num_bins is fixed, which should be improved. TODO*/
/* A typical cacheline size is 64*/
/*#define DEFAULT_WORD_SIZE     8
#define MAX_NUM_BINS          8*/
#define DEFAULT_WORD_SIZE     4
#define MAX_NUM_BINS          16

enum {COMP_MISS, CONF_MISS, CAP_MISS, MISS_CATEGORY_COUNT}; //TODO: replace m1_comp, m1_conf, m1_cap
//ULong m1_type[MISS_CATEGORY_COUNT];
//ULong mL_type[MISS_CATEGORY_COUNT];
//int cat = miss_infi ? COMP : (!miss_fa ? CONF : CAP);
//cc->m1[cat]++;           // m1 is now an int m1[3]
//if (cr) cr->m[cat]++;    // m is an int m[3]

typedef
   struct {
      ULong a;  /* total # memory accesses of this kind */
      ULong m1; /* misses in the first level cache */
      ULong mL; /* misses in the second level cache */
      ULong m1_comp, m1_conf, m1_cap;  /* 3 types of cache misses in the first level cache: compulsory, conflict and capacity */
      ULong mL_comp, mL_conf, mL_cap;  /* 3 types of cache misses in the second level cache: compulsory, conflict and capacity */
   }
   CacheCC;

typedef
   struct {
      ULong b;  /* total # branches of this kind */
      ULong mp; /* number of branches mispredicted */
   }
   BranchCC;

//------------------------------------------------------------
// Primary data structure #1: CC table
// - Holds the per-source-line hit/miss stats, grouped by file/function/line.
// - an ordered set of CCs.  CC indexing done by file/function/line (as
//   determined from the instrAddr).
// - Traversed for dumping stats at end in file/func/line hierarchy.

typedef struct {
   HChar* file;
   const HChar* fn;
   Int    line;
}
CodeLoc;

typedef struct {
   CodeLoc  loc; /* Source location that these counts pertain to */
   CacheCC  Ir;  /* Insn read counts */
   CacheCC  Dr;  /* Data read counts */
   CacheCC  Dw;  /* Data write/modify counts */
   BranchCC Bc;  /* Conditional branch counts */
   BranchCC Bi;  /* Indirect branch counts */

/*----------Extension of cache efficiency analysis -----------*/
   ULong num_evicts_D1[MAX_NUM_BINS]; /* The number of cachline evictions with n words used*/
   ULong num_evicts_LL[MAX_NUM_BINS]; /* The number of cachline evictions with n words used*/

   OSet *cr_table_D1; //cacheline replacement list for observing variables in D1, each element is a cacheline_rep_t
   OSet *cr_table_LL; //cacheline replacement list for observing variables in LL, each element is a cacheline_rep_t

   OSet *cu_table_D1; //cacheline spatial usage list for observing variables in D1, each element is a cacheline_usage_t
   OSet *cu_table_LL; //cacheline spatial usage list for observing variables in LL, each element is a cacheline_usage_t
} LineCC;

// First compare file, then fn, then line.
static Word cmp_CodeLoc_LineCC(const void *vloc, const void *vcc)
{
   Word res;
   const CodeLoc* a = (const CodeLoc*)vloc;
   const CodeLoc* b = &(((const LineCC*)vcc)->loc);

   res = VG_(strcmp)(a->file, b->file);
   if (0 != res)
      return res;

   res = VG_(strcmp)(a->fn, b->fn);
   if (0 != res)
      return res;

   return a->line - b->line;
}

/*----------Extension of cache efficiency by JinChao-----------*/
Int CU_DEBUG = 1;

//Setting nth bit in a bitvector on.
static
void bitop_set(ULong* bv, UInt pos)
{
	if(bv == NULL || pos < 0)
		return;

	*bv |= 1 << pos;	
}

//Setting multiple bits in a bitvector on.
__attribute__((always_inline))
static __inline__
void bitop_set_range(ULong* bv, UInt begin, UInt end)
{
	Int i;
	for(i = begin; i <= end; i++)
		*bv |= 1 << i;
}

// Counting non-zero bits in a bit vector using Brian Kernighan’s Algorithm
__attribute__((always_inline))
static __inline__
UInt bitop_count(ULong bv)
{
	UInt count = 0;

	while(bv) {
		bv &= (bv - 1);
		count ++;
	}

	return count;
}

static
Int open_cu_log(void)
{
   const HChar* cu_out_file = "cacheusage.dbg";
//      VG_(expand_file_name)("--cachegrind-out-file", clo_ce_out_file);

   cu_fp = VG_(fopen)(cu_out_file, VKI_O_CREAT|VKI_O_TRUNC|VKI_O_WRONLY,
                                        VKI_S_IRUSR|VKI_S_IWUSR);
   if (cu_fp == NULL) 
      return -1;

   return 0;
}

static
void close_cu_log(void)
{
   if (!cu_fp) 
      return;

   VG_(fclose)(cu_fp);
}

typedef struct {
  UWord        tag;
  ULong        bitvector;   // keep track of spatial usage. bit 0 represents only 1 word been used, bit 1 represents only 2 words been used, and so on ...
  ULong        num_accesses; // keep track of temporary reuse. For every cache hit, it increases by one
  Int          line_num;     // source code line number, for the line that move this cacheline in. It is maily for debugging purpose.
  LineCC       *src_line;    // pointer to LineCC in cg_main.c, for the line that move this cacheline in.
} cacheline_t;

typedef struct {
   Int          size;                   /* bytes */
   Int          assoc;
   Int          line_size;              /* bytes */
   Int          sets;
   Int          sets_min_1;
   Int          line_size_bits;
   Int          tag_shift;
   HChar        desc_line[128];         /* large enough */
//   UWord*       tags;
   UInt         line_mask;
   Int          num_words_per_line;
   Int          word_size_bits;
   cacheline_t  *cachelines;
   UInt         *lru_list;
} cache_t2;


static cache_t2 LL;
static cache_t2 I1;
static cache_t2 D1;
static cache_t2 D1_fa;

static cache_infi INFI;
static cache_fa FA_D1;
static cache_fa FA_LL;

#define MAX_NAME_LEN 32
typedef struct {
   UInt         id;
   Addr         begin, end;
   HChar        *name;
//   HChar        name[MAX_NAME_LEN];
} variable_t;
static UInt var_index = 0;
static OSet *var_table = NULL;
static variable_t **var_array = NULL;

typedef struct {
   UWord pair_id; //id = var_in->id * 1024 + var_out->id; we assume the total number of variables been tracked is less than 1,024.
   ULong m_comp, m_conf, m_cap;  /* 3 types of cache misses: compulsory, conflict and capacity */
} cacheline_rep_t;

typedef struct {
   UWord vid; //variable id
   ULong num_accesses_D1;
   ULong num_accesses_LL;
   ULong num_evicts_D1[MAX_NUM_BINS]; /* The number of cachline evictions with n words used*/
   ULong num_evicts_LL[MAX_NUM_BINS]; /* The number of cachline evictions with n words used*/
} cacheline_usage_t;

static Word cmp_var_range(const void *vleft, const void *vright)
{
   const variable_t* a = (const variable_t*)vleft;
   const variable_t* b = (const variable_t*)vright;

   if(a->end < b->begin)
       return -1;
    
   if(b->end < a->begin)
       return 1;

    return 0; //This means two variables overlap, which is not allowed.
}

static Word comp_addr_var(const void* vkey, const void* velem)
{
   Addr   tag  = *(const Addr*)vkey;
   const variable_t* var = (const variable_t*)velem;
  
   if(tag < var->begin )
       return -1;

   if(tag > var->end )
       return 1;

   return 0;
}
 
static void print_var_table(void)
{
   variable_t *var;

   VG_(OSetGen_ResetIter)(var_table);
   while ( (var = VG_(OSetGen_Next)(var_table)) ) {
      VG_(printf)("%s %u %#lx %#lx\n", var->name, var->id, var->begin, var->end);
   }
}

static void init_var_array(void)
{
   UInt num = VG_(OSetGen_Size)(var_table);
   var_array = VG_(malloc)("cg.sim.var_table.2", sizeof(variable_t*) * num);

   variable_t *var = NULL;
   UInt index = 0;
   VG_(OSetGen_ResetIter)(var_table);
   while ( (var = VG_(OSetGen_Next)(var_table)) ) {
      var_array[index++] = var;
   }

   if(index != var_index)
      VG_(printf)("ERROR: init_var_array(), the number of variables: %u %u \n", index, var_index);

/*   Int i;
   for(i = 0; i < var_index; i++)
      VG_(printf)("%s %u %#lx %#lx\n", var_array[i]->name, var_array[i]->id, var_array[i]->begin, var_array[i]->end);*/
}

static variable_t* find_var_in_table(Addr addr)
{
   variable_t *var;

   VG_(OSetGen_ResetIter)(var_table);
   while ( (var = VG_(OSetGen_Next)(var_table)) ) {
       if(var->begin <= addr && var->end >=addr)
           return var;
   }
   return NULL;
}

/* By this point, the size/assoc/line_size has been checked. */
static void cachesim_initcache(cache_t config, cache_t2* c, UInt word_size)
{
   Int i, j;

   c->size      = config.size;
   c->assoc     = config.assoc;
   c->line_size = config.line_size;

   c->sets           = (c->size / c->line_size) / c->assoc;
   c->sets_min_1     = c->sets - 1;
   c->line_size_bits = VG_(log2)(c->line_size);
   c->tag_shift      = c->line_size_bits + VG_(log2)(c->sets);

   if (c->assoc == 1) {
      VG_(sprintf)(c->desc_line, "%d B, %d B, direct-mapped", 
                                 c->size, c->line_size);
   } else {
      VG_(sprintf)(c->desc_line, "%d B, %d B, %d-way associative",
                                 c->size, c->line_size, c->assoc);
   }

/*   c->tags = VG_(malloc)("cg.sim.ci.1",
                         sizeof(UWord) * c->sets * c->assoc);*/

   c->line_mask = c->line_size - 1;
   c->word_size_bits = VG_(log2)(word_size);
   c->num_words_per_line = c->line_size / word_size; 

   c->cachelines = VG_(malloc)("cg.sim.ci.1",
                         sizeof(cacheline_t) * c->sets * c->assoc);

   for (i = 0; i < c->sets * c->assoc; i++)
   {
//      c->tags[i] = 0;
        c->cachelines[i].tag = 0;
        c->cachelines[i].bitvector = 0;
        c->cachelines[i].num_accesses = 0;
        c->cachelines[i].line_num = 0;
        c->cachelines[i].src_line = NULL;
   }

   c->lru_list = VG_(malloc)("cg.sim.ci.2",
                         sizeof(UInt) * c->sets * c->assoc);

   for (i = 0; i < c->sets; i++)
   {
     for (j = 0; j < c->assoc; j++)
       c->lru_list[i * c->assoc + j] = c->assoc - 1 - j;
   }

   var_table = 
      VG_(OSetGen_Create)(0,
                          cmp_var_range,
                          VG_(malloc), "cg.sim.ci.3",
                          VG_(free));
}

/* This attribute forces GCC to inline the function, getting rid of a
 * lot of indirection around the cache_t2 pointer, if it is known to be
 * constant in the caller (the caller is inlined itself).
 * Without inlining of simulator functions, cachegrind can get 40% slower.
 */
__attribute__((always_inline))
static __inline__
Bool cachesim_setref_is_miss(cache_t2* c, UInt set_no, UWord tag, UInt word_begin, UInt word_end, Int line_num, void* line, cacheline_rep_t **cr_hook)
{
   int i, j;
//   UWord *set;
   cacheline_t *cacheline;
   UInt *id;
   UInt tmp, num_words;
   ULong num_accesses;

//   set = &(c->tags[set_no * c->assoc]);
   cacheline = &(c->cachelines[set_no * c->assoc]);
   id = &(c->lru_list[set_no * c->assoc]);

//   if (CU_DEBUG && cu_fp && c == &D1 && (word_begin > MAX_NUM_BINS || word_end > MAX_NUM_BINS)) 
//      VG_(fprintf)(cu_fp,  "ERROR: set_no: %u, tag: %lx, begin: %u, end: %u\n", set_no, tag, word_begin, word_end);

   /* This loop is unrolled for just the first case, which is the most */
   /* common.  We can't unroll any further because it would screw up   */
   /* if we have a direct-mapped (1-way) cache.                        */
   if (tag == cacheline[id[0]].tag)
   {
      bitop_set_range(&cacheline[id[0]].bitvector, word_begin, word_end);
      cacheline[id[0]].num_accesses++;

      return False;
   }

   /* If the tag is one other than the MRU, move it into the MRU spot  */
   /* and shuffle the rest down.                                       */
   for (i = 1; i < c->assoc; i++) {
      if (tag == cacheline[id[i]].tag) {
         tmp = id[i];
         for (j = i; j > 0; j--) {
            id[j] = id[j - 1];
         }
         id[0] = tmp;

         bitop_set_range(&cacheline[tmp].bitvector, word_begin, word_end); 
         cacheline[tmp].num_accesses++;

         return False;
      }
   }

   /* A miss;  install this tag as MRU, shuffle rest down. */
   UInt evict_id = id[c->assoc - 1];
   cacheline_t evict_line = cacheline[evict_id];
   num_words = bitop_count(evict_line.bitvector);
   num_accesses = evict_line.num_accesses;

//   if (CU_DEBUG && (!num_words || num_words > MAX_NUM_BINS) && evict_line.tag && cu_fp && c == &D1)
//      VG_(fprintf)(cu_fp,  "ERROR: Ev %lx %x, %u, line: %d, %p\n", evict_line.tag, evict_line.bitvector, num_words, evict_line.line_num, evict_line.src_line);

   for (j = c->assoc - 1; j > 0; j--) {
      id[j] = id[j - 1];
   }
   cacheline[evict_id].tag = tag;
   cacheline[evict_id].bitvector = 0;
   cacheline[evict_id].num_accesses = 0;
   cacheline[evict_id].line_num = line_num;
   cacheline[evict_id].src_line = line;
   bitop_set_range(&cacheline[evict_id].bitvector, word_begin, word_end);
   id[0] = evict_id;

   if(num_words && evict_line.tag && evict_line.src_line)
   {
     if(c==&D1)
       evict_line.src_line->num_evicts_D1[num_words-1]++;

     if(c==&LL)
       evict_line.src_line->num_evicts_LL[num_words-1]++;

   }

   LineCC* rline = (LineCC*)line;
   variable_t *var_in = NULL, *var_out = NULL;

   //We only check the start address of each cacheline, which may not be the actually address of data been accessed. 
   //We assume the number of errors generated can be ignored
   Addr out = evict_line.tag << c->line_size_bits;
   var_out = VG_(OSetGen_LookupWithCmp)(var_table, &out, comp_addr_var); 
//   var_out = find_var_in_table(out);
   if(var_out != NULL)
   {
      UWord vid = var_out->id;

      if(c==&D1)
      {
         if(rline->cu_table_D1 == NULL)
            rline->cu_table_D1 = VG_(OSetGen_Create)(offsetof(cacheline_usage_t, vid), NULL, VG_(malloc), "cg.sim.cu.1", VG_(free));
         cacheline_usage_t *cu = VG_(OSetGen_Lookup)(rline->cu_table_D1, &vid);
         if(cu == NULL)
         {
           cu = VG_(OSetGen_AllocNode)(rline->cu_table_D1, sizeof(cacheline_usage_t));
           cu->vid = vid;
           for(i = 0; i < MAX_NUM_BINS; i++)
           {
             cu->num_evicts_D1[i] = 0;
           }
           cu->num_accesses_D1 = 0;
           VG_(OSetGen_Insert)(rline->cu_table_D1, cu);
         }
         cu->num_evicts_D1[num_words-1]++; 
         cu->num_accesses_D1 += num_accesses;
      }

      Addr in = tag << c->line_size_bits;
      var_in = VG_(OSetGen_LookupWithCmp)(var_table, &in, comp_addr_var); //We assume the number of errors can be ignored
//      var_in = find_var_in_table(in);
      if(var_in != NULL)
      {
         UWord pair_id = (var_in->id << 10) + var_out->id;

         if(c==&D1)
         {
            if(rline->cr_table_D1 == NULL)
               rline->cr_table_D1 = VG_(OSetGen_Create)(offsetof(cacheline_rep_t, pair_id), NULL, VG_(malloc), "cg.sim.cr.1", VG_(free));
            cacheline_rep_t *cr = VG_(OSetGen_Lookup)(rline->cr_table_D1, &pair_id);
            if(cr == NULL)
            {
              cr = VG_(OSetGen_AllocNode)(rline->cr_table_D1, sizeof(cacheline_rep_t));
              cr->pair_id = pair_id;
              cr->m_comp = 0;
              cr->m_conf = 0;
              cr->m_cap = 0;
              VG_(OSetGen_Insert)(rline->cr_table_D1, cr);
            }
            *cr_hook = cr;
         }
      }
   }
   return True;
}

__attribute__((always_inline))
static __inline__
Bool cachesim_ref_is_miss(cache_t2* c, Addr a, UChar size, Int line_num, LineCC *line, cacheline_rep_t **cr_hook)
{
   /* A memory block has the size of a cache line */
   UWord block1 =  a         >> c->line_size_bits;
   UWord block2 = (a+size-1) >> c->line_size_bits;
   UInt  set1   = block1 & c->sets_min_1;

   UWord addr_offset = a & c->line_mask; 
   UWord word_begin = addr_offset >> c->word_size_bits;
   UWord word_end1 = (addr_offset + size - 1) >> c->word_size_bits;

   /* Tags used in real caches are minimal to save space.
    * As the last bits of the block number of addresses mapping
    * into one cache set are the same, real caches use as tag
    *   tag = block >> log2(#sets)
    * But using the memory block as more specific tag is fine,
    * and saves instructions.
    */
   UWord tag1   = block1;

   /* Access entirely within line. */
   if (block1 == block2)
      return cachesim_setref_is_miss(c, set1, tag1, word_begin, word_end1, line_num, line, cr_hook);

   /* Access straddles two lines. */
   else if (block1 + 1 == block2) {
      UInt  set2 = block2 & c->sets_min_1;
      UWord tag2 = block2;

      UWord word_end2 = word_end1 - c->num_words_per_line;
      word_end1 = c->num_words_per_line - 1;

      /* always do both, as state is updated as side effect */
      if (cachesim_setref_is_miss(c, set1, tag1, word_begin, word_end1, line_num, line, cr_hook)) {
         cachesim_setref_is_miss(c, set2, tag2, 0, word_end2, line_num, line, cr_hook);
         return True;
      }
      return cachesim_setref_is_miss(c, set2, tag2, 0, word_end2, line_num, line, cr_hook);
   }
   VG_(printf)("addr: %lx  size: %u  blocks: %lu %lu",
               a, size, block1, block2);
   VG_(tool_panic)("item straddles more than two cache sets");
   /* not reached */
   return True;
}

static
void cachesim_collect_undrained_lines(cache_t2* c)
{
   Int i, j;
   UInt id, num_words;
   cacheline_t *cl = c->cachelines;

   for (i = 0; i < c->sets; i++)
   {
     for (j = 0; j < c->assoc; j++)
     {
        id = c->lru_list[i * c->assoc + j];
        if(cl[id].tag && cl[id].src_line) 
        {
           num_words = bitop_count(cl[id].bitvector);

           if(c==&D1)
             cl[id].src_line->num_evicts_D1[num_words-1]++;
           if(c==&LL)
             cl[id].src_line->num_evicts_LL[num_words-1]++;
        }
     }
   }
}

static void cachefa_initcache(cache_t config, cache_fa* c)
{
//   VG_(fprintf)(cu_fp, "cachefa_initcache capacity: %d\n", config.size);
   cachefa_setup(c, (config.size / config.line_size));
}

static void cachesim_initcaches(cache_t I1c, cache_t D1c, cache_t LLc, UInt word_size)
{
   open_cu_log();

   cachesim_initcache(I1c, &I1, word_size);
   cachesim_initcache(D1c, &D1, word_size);
   cachesim_initcache(LLc, &LL, word_size);

   D1c.assoc = D1c.size / D1c.line_size;
   cachesim_initcache(D1c, &D1_fa, word_size);

   //VG_(printf)("cachesim_initcaches word_size: %u, word_size_bits: %d\n", word_size, D1.word_size_bits);

   cachefa_initcache(D1c, &FA_D1);
   cachefa_initcache(LLc, &FA_LL);
}

static void cachesim_finish(void)
{
   cachesim_collect_undrained_lines(&D1);
   cachesim_collect_undrained_lines(&LL);
   close_cu_log();
}

__attribute__((always_inline))
static __inline__
void cachesim_I1_doref_Gen(Addr a, UChar size, ULong* m1, ULong *mL)
{
   cacheline_rep_t *cr = NULL;
   if (cachesim_ref_is_miss(&I1, a, size, 0, NULL, &cr)) {
      (*m1)++;
      if (cachesim_ref_is_miss(&LL, a, size, 0, NULL, &cr))
         (*mL)++;
   }
}

// common special case IrNoX
__attribute__((always_inline))
static __inline__
void cachesim_I1_doref_NoX(Addr a, UChar size, ULong* m1, ULong *mL)
{
   cacheline_rep_t *cr = NULL;

   UWord block  = a >> I1.line_size_bits;
   UInt  I1_set = block & I1.sets_min_1;

   UWord addr_offset = a & I1.line_mask; 
   UWord word_begin = addr_offset >> I1.word_size_bits;
   UWord word_end = (addr_offset + size - 1) >> I1.word_size_bits;

   // use block as tag
   if (cachesim_setref_is_miss(&I1, I1_set, block, word_begin, word_end, 0, NULL, &cr)) {
      UInt  LL_set = block & LL.sets_min_1;
      (*m1)++;
      // can use block as tag as L1I and LL cache line sizes are equal
      if (cachesim_setref_is_miss(&LL, LL_set, block, word_begin, word_end, 0, NULL, &cr))
         (*mL)++;
   }
}

__attribute__((always_inline))
static __inline__
Bool cachesim_D1_doref(Addr a, UChar size, ULong* m1, ULong *mL, int line_num, LineCC* line, CacheCC* cc)
{
   cacheline_rep_t *cr = NULL;

   Bool miss_infi;
   Bool miss_fa;
   Bool miss_fa_LL;

   miss_infi = cacheinfi_ref_is_miss(&INFI, a, size);
   miss_fa = cachefa_ref_is_miss(&FA_D1, a, size);
   miss_fa_LL = cachefa_ref_is_miss(&FA_LL, a, size);

//   Bool miss_fa_D1 = cachesim_ref_is_miss(&D1_fa, a, size, line_num, NULL, NULL);

   if (cachesim_ref_is_miss(&D1, a, size, line_num, line, &cr)) {
      (*m1)++;

     if(miss_infi)
     {
        cc->m1_comp++;
        if(cr!=NULL)
          cr->m_comp++;
      }
      else if(!miss_fa)
//     else if(!miss_fa_D1)
      {
        cc->m1_conf++;
        if(cr!=NULL)
          cr->m_conf++;
      }
      else {
        cc->m1_cap++;
        if(cr!=NULL)
          cr->m_cap++;
      }

      if (cachesim_ref_is_miss(&LL, a, size, line_num, line, &cr)) {
         (*mL)++;

         if(miss_infi)
           cc->mL_comp++;
         else if(miss_fa_LL)
           cc->mL_conf++;
         else
           cc->mL_cap++;
      }

      return True;
   }

   /* messy hack to read variables in the simulated program*/
   /* A file named "varinfo.txt" is generated by the program under simulation */
   /* We monitor its presence and read it to fill var_table */
   static Int ref_counter = 0;
   static Bool done = False;

   struct vg_stat stat_buf;
   VG_(memset)(&stat_buf, 0, sizeof(stat_buf));

   const HChar filename[32] = "varinfo.txt";
   SysRes  var_fd, stat_fd;
   Int     fdno;
   HChar    buffer[512], swap_buffer[256];
   Int      read_len = 256;
   HChar    *read_buffer = buffer, *buffer_end;
   HChar    *line_start, *line_end, *var_name = NULL, *start = NULL;
   Addr    begin, end;
   Int     len;
   Bool    finished;
   variable_t *var;

   if(!done)
   {
     if(ref_counter == 10000)
     {
       ref_counter = 0;

       stat_fd = VG_(stat)(filename, &stat_buf);
       if(!sr_isError(stat_fd))
       {
         var_fd = VG_(open)(filename, VKI_O_RDONLY, VKI_S_IRUSR);

         if (sr_isError(var_fd)) {
           VG_(printf)("Error: can't open %s!\n", filename);
         } else {
           fdno = sr_Res(var_fd);

           while( (len = VG_(read)(fdno, read_buffer, read_len)))
           {
//             VG_(printf)("\nread len: %d\n", len);
             buffer_end = read_buffer + len;
             line_start = buffer;
             line_end = line_start;
             finished = False;
             do
             {
               while(*line_end != '\n' && line_end < buffer_end) //grab a line
               {
                 line_end++;
               }
               
               //find a line. We assume each line follows the format: "name start end". We don't check errors.
               if(*line_end == '\n')
               {
                 var_name = line_start;
                 while(*line_start != ' ') line_start++;
                 *line_start = '\0';
    
                 start = ++line_start;
                 while(*line_start != ' ') line_start++;
                 begin = VG_(strtoll16)(start, &line_start);
    
                 start = ++line_start;
                 end = VG_(strtoll16)(start, &line_end);

                 var = VG_(OSetGen_AllocNode)(var_table, sizeof(variable_t));
                 var->begin =  begin;
                 var->end =  end;
                 var->name = VG_(strdup)("cg.sim.var_table.1", var_name);
                 var->id = var_index++;
                 VG_(OSetGen_Insert)(var_table, var);

                 line_end++;
                 line_start = line_end;

                 if(line_start == buffer_end)
                 {
                    read_buffer = buffer; //The whole buffer is finished.
                    finished = True;
                 }
               }
               else
               {
                  //We need to handle the buffer tail that is not finished
                  VG_(memcpy)(swap_buffer, line_start, buffer_end - line_start);
                  VG_(memcpy)(buffer, swap_buffer, buffer_end - line_start);
                  read_buffer = buffer + (buffer_end-line_start);
                  finished = True;
               }
             } while(!finished);
           }
           VG_(close)(fdno);
           init_var_array();
           /*print_var_table();*/
           VG_(printf)("var_index: %u\n", var_index);
           if(var_index > 0) //Sometime file read returns a NULL buffer. We need to keep trying it.
               done = True;
         }
       }
     }
     ref_counter++;
   }
   return False;
}

/* Check for special case IrNoX. Called at instrumentation time.
 *
 * Does this Ir only touch one cache line, and are L1I/LL cache
 * line sizes the same? This allows to get rid of a runtime check.
 *
 * Returning false is always fine, as this calls the generic case
 */
static Bool cachesim_is_IrNoX(Addr a, UChar size)
{
   UWord block1, block2;

   if (I1.line_size_bits != LL.line_size_bits) return False;
   block1 =  a         >> I1.line_size_bits;
   block2 = (a+size-1) >> I1.line_size_bits;
   if (block1 != block2) return False;

   return True;
}

/*--------------------------------------------------------------------*/
/*--- end                                                 cg_sim.c ---*/
/*--------------------------------------------------------------------*/

