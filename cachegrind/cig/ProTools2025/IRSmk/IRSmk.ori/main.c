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
void init(Domain_t *, RadiationData_t *, double *, double *);
void readInput(const char *);
void rmatmult3(Domain_t *, RadiationData_t *, double *, double *);


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

  if(varinfo_fname != NULL)
  {
    int res = varinfo_file_open(varinfo_fname);
    if(res)
      fprintf(stderr, "Cannot open %s for recording varinfo\n", varinfo_fname);

    varinfo_file_print(getName(b), b, b+i_ub);
    varinfo_file_print(getName(x), x, x+x_size);
    varinfo_file_print(getName(rblk_ptr->dbl), rblk_ptr->dbl, rblk_ptr->dbl+i_ub);
    varinfo_file_print(getName(rblk_ptr->dbc), rblk_ptr->dbc, rblk_ptr->dbc+i_ub);
    varinfo_file_print(getName(rblk_ptr->dbr), rblk_ptr->dbr, rblk_ptr->dbr+i_ub);
    varinfo_file_print(getName(rblk_ptr->dcl), rblk_ptr->dcl, rblk_ptr->dcl+i_ub);
    varinfo_file_print(getName(rblk_ptr->dcc), rblk_ptr->dcc, rblk_ptr->dcc+i_ub);
    varinfo_file_print(getName(rblk_ptr->dcr), rblk_ptr->dcr, rblk_ptr->dcr+i_ub);
    varinfo_file_print(getName(rblk_ptr->dfl), rblk_ptr->dfl, rblk_ptr->dfl+i_ub);
    varinfo_file_print(getName(rblk_ptr->dfc), rblk_ptr->dfc, rblk_ptr->dfc+i_ub);
    varinfo_file_print(getName(rblk_ptr->dfr), rblk_ptr->dfr, rblk_ptr->dfr+i_ub);
    varinfo_file_print(getName(rblk_ptr->cbl), rblk_ptr->cbl, rblk_ptr->cbl+i_ub);
    varinfo_file_print(getName(rblk_ptr->cbc), rblk_ptr->cbc, rblk_ptr->cbc+i_ub);
    varinfo_file_print(getName(rblk_ptr->cbr), rblk_ptr->cbr, rblk_ptr->cbr+i_ub);
    varinfo_file_print(getName(rblk_ptr->ccl), rblk_ptr->ccl, rblk_ptr->ccl+i_ub);
    varinfo_file_print(getName(rblk_ptr->ccc), rblk_ptr->ccc, rblk_ptr->ccc+i_ub);
    varinfo_file_print(getName(rblk_ptr->ccr), rblk_ptr->ccr, rblk_ptr->ccr+i_ub);
    varinfo_file_print(getName(rblk_ptr->cfl), rblk_ptr->cfl, rblk_ptr->cfl+i_ub);
    varinfo_file_print(getName(rblk_ptr->cfc), rblk_ptr->cfc, rblk_ptr->cfc+i_ub);
    varinfo_file_print(getName(rblk_ptr->cfr), rblk_ptr->cfr, rblk_ptr->cfr+i_ub);
    varinfo_file_print(getName(rblk_ptr->ubl), rblk_ptr->ubl, rblk_ptr->ubl+i_ub);
    varinfo_file_print(getName(rblk_ptr->ubc), rblk_ptr->ubc, rblk_ptr->ubc+i_ub);
    varinfo_file_print(getName(rblk_ptr->ubr), rblk_ptr->ubr, rblk_ptr->ubr+i_ub);
    varinfo_file_print(getName(rblk_ptr->ucl), rblk_ptr->ucl, rblk_ptr->ucl+i_ub);
    varinfo_file_print(getName(rblk_ptr->ucc), rblk_ptr->ucc, rblk_ptr->ucc+i_ub);
    varinfo_file_print(getName(rblk_ptr->ucr), rblk_ptr->ucr, rblk_ptr->ucr+i_ub);
    varinfo_file_print(getName(rblk_ptr->ufl), rblk_ptr->ufl, rblk_ptr->ufl+i_ub);
    varinfo_file_print(getName(rblk_ptr->ufc), rblk_ptr->ufc, rblk_ptr->ufc+i_ub);
    varinfo_file_print(getName(rblk_ptr->ufr), rblk_ptr->ufr, rblk_ptr->ufr+i_ub);
    varinfo_file_close();
  }
 
  init(domain_ptr, rblk_ptr, x, b);

  for (i=0; i<noIter; ++i) {
     rmatmult3(domain_ptr, rblk_ptr, x, b);
  }

  printf("***** results \n");  
  for (i=0; i<i_ub; i+=i_ub/5) {
    printf("i = %10d      b[i] = %e \n", i, b[i]);
  }

  return(0);
}
