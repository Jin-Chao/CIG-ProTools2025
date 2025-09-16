/*
   This file consists of an infinite cache simulator and a fully associative cache simulator.
   They are used by Cachegrind to identify compulsory and conflict cache misses.
*/

#define HELPER_DEBUG 0

/* The following configuration is only for supporting 64 bytes cacheline */

#define RANGE_SIZE (1024*64*64) /* size of bits*/

#define CB_SIZE 64
#define LOG2CB 6

#define LOG2BM_BLOCK 12  /*each bit stands for a cacheblock, a long is 64 bits*/

#define LOG2RANGE 22

#define CB_MASK (((long)1 << LOG2CB) - 1)
#define BM_BLOCK_MASK ((((long)1 << LOG2BM_BLOCK) - 1) & ~CB_MASK)
#define BITMAP_MASK ((((long)1 << LOG2RANGE) - 1) & ~BM_BLOCK_MASK & ~CB_MASK)
#define RANGE_MASK (~(((long)1 << LOG2RANGE) - 1))

VgFile  *cu_fp = NULL;

typedef struct {
    long      addr;
    char      type;             /* read (r) or write (w)*/
    unsigned short  size;       /* of memory referenced, in bytes */
} MemRef; /* We assume each memory refernce cannot span across more than two cache blocks. */

/*data structure for infinitive cache */
typedef struct {
    Addr     addr;       /* start address of range */
    ULong    *bitmap;     /* one bit stands for a cacheblock*/
} MemRange;

typedef struct {
    int max_num_ranges; /* the number of ranges allocated */
    int cur_num_ranges; /* the number of ranges used */
    MemRange  *ranges;
} cache_infi;

/*data structure for the fully associative cache */
typedef struct cache_block_struct {
    UWord    block_addr;  /* TAG of a chace line*/
    struct cache_block_struct *up;      /* more recently used */
    struct cache_block_struct *down;    /* less recently used */ 
    struct cache_block_struct *bucket;  /* for hash collisions*/
} CacheBlock; 

typedef struct {
  int size;
  CacheBlock **buckets;
} HashTable;

typedef struct {
    int num_blocks;            /* the number of cache lines */
    CacheBlock  *blocks_list;  /*A circular list doubly linked*/
    CacheBlock  *top;
    HashTable *table;
} cache_fa;

/* return 0 if ref has been accessed, 1 if not been accessed */
__attribute__((always_inline))
static __inline__
int cacheinfi_handle_access(cache_infi *cache, Addr a, UChar size )
{
    UInt map_offset, bit_start, bit_end, b;
    int hi, lo, i;
    ULong cb_id;

    cb_id = a >> LOG2CB;
    map_offset = (a & BITMAP_MASK) >> LOG2BM_BLOCK; 
    bit_start = (a & BM_BLOCK_MASK) >> LOG2CB;
    bit_end = ((a + size - 1) & BM_BLOCK_MASK) >> LOG2CB; 

    if(bit_start == bit_end)
        bit_end++;

    if(HELPER_DEBUG)
        VG_(fprintf)(cu_fp, "cb addr: %llx, map_offset: %x, bit:[%x - %x]\n", cb_id, map_offset, bit_start, bit_end);

    hi = cache->cur_num_ranges - 1;
    lo = 0;
    while(lo <= hi)
    {
        i = (hi + lo)/2;
        if(cache->ranges[i].addr + RANGE_SIZE <= a)
           lo = i + 1;
        else if(cache->ranges[i].addr > a)
           hi = i - 1;
        else { /* found the range */
            int sum = 0;

            if(HELPER_DEBUG)
                VG_(fprintf)(cu_fp, "\trange [%d] is found, bitmap: %p %llx\n", i, &(cache->ranges[i].bitmap[map_offset]), cache->ranges[i].bitmap[map_offset]);
            for ( b = bit_start; b < bit_end; b++)
               sum += ((cache->ranges[i].bitmap[map_offset] & ((long)1 << b)) != 0);
            
            if(HELPER_DEBUG)
                VG_(fprintf)(cu_fp, "\tsum: %d\n", sum);

            if (sum == (bit_end - bit_start))
                return 0; /*every byte is accessed*/

            for ( b = bit_start; b < bit_end; b++)  /*set bits*/
               cache->ranges[i].bitmap[map_offset] |= ((long)1 << b);

            if(HELPER_DEBUG)
                VG_(fprintf)(cu_fp, "-=<>=- bitmap: %p %llx\n", &(cache->ranges[i].bitmap[map_offset]), cache->ranges[i].bitmap[map_offset]);
            return 1;
        }
    }

    /* (lo > hi) means a new range should be created*/
    if(cache->cur_num_ranges >= (cache->max_num_ranges)) /* all allocated ranged are used*/
    {
        cache->max_num_ranges = (cache->max_num_ranges + 16) * 2;
        if(cache->ranges == NULL)
            cache->ranges = VG_(calloc)("InfiCache.1", cache->max_num_ranges * sizeof(MemRange), 1);
        else 
            cache->ranges = VG_(realloc)("InfiCache.resize.1", cache->ranges, cache->max_num_ranges * sizeof(MemRange));

        if(cache->ranges == NULL)
        {
            VG_(fprintf)(cu_fp, "Memory allocation error for c->ranges\n");
            VG_(exit)(1);
        }
    }

    for(i = cache->cur_num_ranges - 1; i >= 0; i--)
    {
        if(cache->ranges[i].addr < a)
            break;
        cache->ranges[i+1] = cache->ranges[i];
    }
    cache->cur_num_ranges++;

    cache->ranges[i+1].addr = a & RANGE_MASK;
    cache->ranges[i+1].bitmap = VG_(calloc)("InfiCache.bitmap", RANGE_SIZE / sizeof(long), sizeof(long));

    if(HELPER_DEBUG)
        VG_(fprintf)(cu_fp, "-=<>=- new range [%lx] is added\n", cache->ranges[i+1].addr);

    if(cache->ranges[i+1].bitmap == NULL)
    {
        VG_(fprintf)(cu_fp, "Memory allocation error for c->ranges[i-1].bitmap\n");
        VG_(exit)(1);
    }
    
    for (b = bit_start; b < bit_end; b++)  /*set bits*/
        cache->ranges[i+1].bitmap[map_offset] |= ((long)1 << b);

    if(HELPER_DEBUG)
        VG_(fprintf)(cu_fp, "-=<>=- num_range [%d], bitmap: %p %llx\n", cache->cur_num_ranges, &(cache->ranges[i+1].bitmap[map_offset]), cache->ranges[i+1].bitmap[map_offset]);

    return 1;
}

__attribute__((always_inline))
static __inline__
int cachefa_setup(cache_fa *cache, int size)
{
    //if(cache == NULL || !size || (size & (size - 1)))
    //    return -1;

    cache->num_blocks = size;
    cache->blocks_list = VG_(calloc)("FACache.cachelines", sizeof(CacheBlock), size);

    cache->table = VG_(calloc)("FACache.HashTable", sizeof(HashTable), size);
    cache->table->size = size;
    cache->table->buckets = VG_(calloc)("FACache.HasthTable.Buckets", sizeof(CacheBlock*), size);

    /* build the doubly linked circular list */
    int i;
    for(i = 1; i < size-1; i++)
    {
        cache->blocks_list[i].down = &cache->blocks_list[i+1];
        cache->blocks_list[i].up = &cache->blocks_list[i-1];
    }
    cache->blocks_list[0].up = &cache->blocks_list[size-1];
    cache->blocks_list[0].down = &cache->blocks_list[1];
    cache->blocks_list[size-1].up = &cache->blocks_list[size-2];
    cache->blocks_list[size-1].down = &cache->blocks_list[0];

    cache->top = cache->blocks_list;
    return 0;
}

/*return value: 0 means a cache hit, 1 means a cache miss, -1 means error*/
__attribute__((always_inline))
static __inline__
int cachefa_handle_access(cache_fa *cache, Addr a, UChar size)
{
    ULong block_addr = a >> LOG2CB;

    /* search ref in the hash table */
    int bucket_id = block_addr % cache->table->size;
    CacheBlock *bucket = cache->table->buckets[bucket_id];
    CacheBlock *ptr = NULL;

    while(bucket != NULL && bucket->block_addr != block_addr)
        bucket = bucket->bucket;

    if(HELPER_DEBUG)
        VG_(fprintf)(cu_fp, "cb addr: %llx, bucket_id: %d, bucket: %p\n", block_addr, bucket_id, bucket);

    if(!bucket) /* a cache miss*/
    {
        CacheBlock *bottom = cache->top->up;

        ULong evicted_addr = bottom->block_addr;
        int evicted_bucket_id = evicted_addr % cache->table->size;

        if(HELPER_DEBUG)
            VG_(fprintf)(cu_fp, "evicted cb addr: %llx, bucket_id: %d\n", evicted_addr, evicted_bucket_id);

        /*replace the bottome block, and move it to the top */
        bottom->block_addr = block_addr;
        cache->top = bottom; 

        if(evicted_addr)
        {
            /* remove the evicted block from the hashtable */
            ptr = cache->table->buckets[evicted_bucket_id];

            if(HELPER_DEBUG)
            {
                VG_(fprintf)(cu_fp, "HashTable Buckets[%d]: ", evicted_bucket_id);
                while(ptr != NULL)
                {
                    VG_(fprintf)(cu_fp, "%lx ", ptr->block_addr);
                    ptr = ptr->bucket;
                }
                VG_(fprintf)(cu_fp, "\n");

                ptr = cache->table->buckets[evicted_bucket_id];
            }

            if(ptr == cache->top)
                cache->table->buckets[evicted_bucket_id] = ptr->bucket;
            else if(ptr != NULL)
            {
                while(ptr->bucket != cache->top)
                   ptr = ptr->bucket;
    
                ptr->bucket = cache->top->bucket;
            }

            if(HELPER_DEBUG)
                VG_(fprintf)(cu_fp, "-=<>=- Hash remove bucket[%d]: addr[%llx]\n", evicted_bucket_id, evicted_addr);
        }

        /* add the new block to the hashtable */
        cache->top->bucket = cache->table->buckets[bucket_id];
        cache->table->buckets[bucket_id] = cache->top;
        if(HELPER_DEBUG)
        {
            VG_(fprintf)(cu_fp, "-=<>=- Hash add bucket[%d]: addr[%lx]\n", bucket_id, cache->table->buckets[bucket_id]->block_addr);
            ptr = cache->table->buckets[bucket_id];
            VG_(fprintf)(cu_fp, "HashTable Buckets[%d]: ", bucket_id);
            while(ptr != NULL)
            {
                VG_(fprintf)(cu_fp, "%lx ", ptr->block_addr);
                ptr = ptr->bucket;
            }
            VG_(fprintf)(cu_fp, "\n");
        }

/*        ptr = cache->table->buckets[bucket_id];
        fprintf(stdout, "HashTable Buckets[%d]: ", bucket_id);
        int count = 0;
        while(ptr != NULL)
        {
            ptr = ptr->bucket;
            count++;
        }
        fprintf(stdout, "%d\n", count);*/

        return 1;
    }
    else 
    {
        if(bucket != cache->top) /* a cache hit for a cacheblock which is not the top one*/
        {
            /* move the cache block to the top */
    
            bucket->up->down = bucket->down; /* detach the block from the list*/ 
            bucket->down->up = bucket->up; 
    
            bucket->up = cache->top->up; /* add it between top and bottom */ 
            bucket->down = cache->top;
            cache->top->up->down = bucket;
            cache->top->up = bucket;
    
            cache->top = bucket;  /* update top */
        }

        return 0; 
    }
}

__attribute__((always_inline))
static __inline__
Bool cacheinfi_ref_is_miss(cache_infi* c, Addr a, UChar size)
{
    Addr a1;
    UChar size1;
    UWord block1 = a >> LOG2CB;    //TODO using c->line_size_bits;
    UWord block2 = (a + size - 1) >> LOG2CB; //TODO using c->line_size_bits;

    if(block1 == block2)
        return cacheinfi_handle_access(c, a, size);
    else
    {
        size1 = (1 << LOG2CB) - (CB_MASK & a) - 1; //TODO using c->line_size_bits 
        a1 = a + size1 + 1;
        if(!cacheinfi_handle_access(c, a, size1))
           return cacheinfi_handle_access(c, a1, size - size1);
        else
           return True;
    }
}

__attribute__((always_inline))
static __inline__
Bool cachefa_ref_is_miss(cache_fa* c, Addr a, UChar size)
{
    Addr a1;
    UChar size1;
    UWord block1 = a >> LOG2CB;    //TODO using c->line_size_bits;
    UWord block2 = (a + size - 1) >> LOG2CB; //TODO using c->line_size_bits;

    if(block1 == block2)
        return cachefa_handle_access(c, a, size);
    else
    {
        size1 = (1 << LOG2CB) - (CB_MASK & a) - 1; //TODO using c->line_size_bits 
        a1 = a + size1 + 1;
        if(!cachefa_handle_access(c, a, size1))
           return cachefa_handle_access(c, a1, size - size1);
        else
           return True;
    }
}
