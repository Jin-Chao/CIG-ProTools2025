/*BHEADER****************************************************************
 * (c) 2006   The Regents of the University of California               *
 *                                                                      *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright       *
 * notice and disclaimer.                                               *
 *                                                                      *
 *EHEADER****************************************************************/


//--------------
//  A micro kernel based on IRS
//    http://www.llnl.gov/asci/purple/benchmarks/limited/irs/
//--------------


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include "irsmk.h"
#include "util.h"

int kmin;
int kmax;
int jmin;
int jmax;  
int imin;
int imax;
int kp;
int jp;

int i_lb;
int i_ub;
int x_size;


void allocMem(RadiationData_t *);
void allocMemAos1(RadiationDataAos1_t **rblk);
void allocMemAos2(RadiationDataAos2_t **rblk);
void allocMemAos3(RadiationDataAos3_t **rblk);
void init(Domain_t *domain, RadiationDataAos1_t *rblk_aos1, RadiationDataAos2_t *rblk_aos2, RadiationDataAos3_t *rblk_aos3, RadiationData_t *rblk, double *x, double *b );
void readInput(const char *);
void rmatmult3(Domain_t *domain, RadiationDataAos1_t *rblk_aos1, RadiationDataAos2_t *rblk_aos2, RadiationDataAos3_t *rblk_aos3, RadiationData_t *rblk, double *x, double *b );


void
print_usage(char **argv)
{
    printf("Usage: %s -i input_name [-v varinfo_file]\n", *argv);
}

int main(int argc, char **argv)
{
  Domain_t domain;
  Domain_t *domain_ptr = &domain;

 RadiationData_t rblk;
 RadiationData_t *rblk_ptr = &rblk;

  RadiationDataAos1_t *rblk_aos1_ptr = NULL;
  RadiationDataAos2_t *rblk_aos2_ptr = NULL;
  RadiationDataAos3_t *rblk_aos3_ptr = NULL;

  struct timeval  t0, t1;
  clock_t t0_cpu = 0,
          t1_cpu = 0;

  double *x;
  double *b;

  int i = 0;
#ifdef SMALL_PROBLEM_SIZE
  const int noIter = 250;
#else
  const int noIter = 5000;
#endif

  printf ("\nSequoia Benchmark Version 1.0\n\n");

  char *input_fname = NULL;
  char *varinfo_fname = NULL;
  int cha;
  while ((cha = getopt (argc, argv, "i:v:h")) != -1)
    switch (cha)
    {
      case 'i':
        input_fname = optarg;
        break;
      case 'v':
        varinfo_fname = optarg;
        break;
      case 'h':
      default:
        print_usage(argv);
        exit(0);
    }

  /*if (argc != 2) {
    printf("Usage: %s <input>\n", argv[0]);
    return 1;
  }
  // 
  readInput(argv[1]);
  */

  if(!input_fname){
    print_usage(argv);
    exit(0);
  }

  readInput(input_fname);

  b = (double *)malloc(i_ub*sizeof(double));
  x = (double *)malloc(x_size*sizeof(double));
  
  allocMem(rblk_ptr);
  allocMemAos1(&rblk_aos1_ptr);
  allocMemAos2(&rblk_aos2_ptr);
  allocMemAos3(&rblk_aos3_ptr);

  if(varinfo_fname != NULL)
  {
    int res = varinfo_file_open(varinfo_fname);
    if(res)
      fprintf(stderr, "Cannot open %s for recording varinfo\n", varinfo_fname);

    varinfo_file_print(getName(b), b, b+i_ub);
    varinfo_file_print(getName(x), x, x+x_size);
    varinfo_file_print(getName(rblk_aos1_ptr), rblk_aos1_ptr, &rblk_aos1_ptr[i_ub].dfc);
    varinfo_file_print(getName(rblk_aos2_ptr), rblk_aos2_ptr, &rblk_aos2_ptr[i_ub].cfl);
    varinfo_file_print(getName(rblk_aos3_ptr), rblk_aos3_ptr, &rblk_aos3_ptr[i_ub].ucr);
    varinfo_file_print(getName(rblk_ptr->ufl), rblk_ptr->ufl, rblk_ptr->ufl+i_ub);
    varinfo_file_print(getName(rblk_ptr->ufc), rblk_ptr->ufc, rblk_ptr->ufc+i_ub);
    varinfo_file_print(getName(rblk_ptr->ufr), rblk_ptr->ufr, rblk_ptr->ufr+i_ub);
    varinfo_file_close();
  }
 
  init(domain_ptr, rblk_aos1_ptr, rblk_aos2_ptr, rblk_aos3_ptr, rblk_ptr, x, b);

  for (i=0; i<noIter; ++i) {
     rmatmult3(domain_ptr, rblk_aos1_ptr, rblk_aos2_ptr, rblk_aos3_ptr, rblk_ptr, x, b);
  }

  printf("***** results \n");  
  for (i=0; i<i_ub; i+=i_ub/5) {
    printf("i = %10d      b[i] = %e \n", i, b[i]);
  }

  return(0);
}
