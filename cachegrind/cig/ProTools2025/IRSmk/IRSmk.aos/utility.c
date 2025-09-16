/*BHEADER****************************************************************
 * (c) 2006   The Regents of the University of California               *
 *                                                                      *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright       *
 * notice and disclaimer.                                               *
 *                                                                      *
 *EHEADER****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "irsmk.h"


//--------------
// read input date
//--------------
void readInput(const char *filename)
{
  FILE *fp;
  if ((fp = fopen(filename, "r"))==NULL) {
    printf("***** ERROR \n");
    printf("      Cannot open input file: %s\n", filename);
    exit(1);
  }

  fscanf(fp, "%d %d %d %d %d %d %d %d",
	 &kmin, &kmax, &jmin, &jmax, &imin, &imax, &kp, &jp);

  printf("***** input  \n");
  printf("kmin = %8d     kmax = %8d \n"  
         "jmin = %8d     jmax = %8d \n"
         "imin = %8d     imax = %8d \n"
         "kp   = %8d     jp   = %8d \n \n \n",
	  kmin, kmax, jmin, jmax, imin, imax, kp, jp);


  //-----
  //  set array size
  //-----
  i_lb = imin + jmin*jp + kmin*kp;   
  i_ub = (imax - 1) + (jmax -1)*jp + (kmax - 1)*kp + 1;

  //-----
  //  sanity check for xdbl
  //-----
  if ((i_lb - kp - jp - 1) < 0) {
    printf("***** ERROR \n");
    printf("      lb of xdbl < 0 \n");
    exit(1);
  }
  
  //-----
  // in theory 
  //  x_size = i_ub           +
  //           kp + jp + 1    +  // for xdbl
  //           kp + jp + 1       // for xufr
  // I just add 10 elements for safety
  //-----
  x_size = i_ub + 2*(kp + jp + 1) + 10;

  printf("***** array bounds  \n");
  printf("i_lb = %10d    i_ub = %10d    x_size = %10d \n \n \n",
	  i_lb, i_ub, x_size);

  fclose(fp);

}


void allocMemAos1(RadiationDataAos1_t **rblk)
{
  int size;

  size = i_ub*sizeof(double);

  *rblk = (RadiationDataAos1_t*)malloc(i_ub*sizeof(RadiationDataAos1_t));
  if ( *rblk == NULL)  printf("*****ERROR: allocMemAos1 out of memory \n");
}

void allocMemAos2(RadiationDataAos2_t **rblk)
{
  int size;

  size = i_ub*sizeof(double);

  *rblk = (RadiationDataAos2_t*)malloc(i_ub*sizeof(RadiationDataAos2_t));
  if ( *rblk == NULL)  printf("*****ERROR: allocMemAos2 out of memory \n");
}

void allocMemAos3(RadiationDataAos3_t **rblk)
{
  int size;

  size = i_ub*sizeof(double);

  *rblk = (RadiationDataAos3_t*)malloc(i_ub*sizeof(RadiationDataAos3_t));
  if ( *rblk == NULL)  printf("*****ERROR: allocMemAos3 out of memory \n");
}

void allocMem(RadiationData_t *rblk)
{
  int size;

  size = i_ub*sizeof(double);

  rblk->ufl = (double *)malloc(size);
  rblk->ufc = (double *)malloc(size);
  rblk->ufr = (double *)malloc(size);

  if ( (rblk->ufr)== NULL)  printf("*****ERROR: allocMem out of memory \n");
}

//--------------
// alloc memory
//--------------
/*
void allocMem(RadiationData_t **rblk)
{
  int size;

  size = i_ub*sizeof(double);

  *rblk = (RadiationData_t*)malloc(i_ub*sizeof(RadiationData_t));
  if ( *rblk == NULL)  printf("*****ERROR: allocMem out of memory \n");
*
  rblk->dbl = (double *)malloc(size);
  rblk->dbc = (double *)malloc(size);
  rblk->dbr = (double *)malloc(size);
  rblk->dcl = (double *)malloc(size);
  rblk->dcc = (double *)malloc(size);
  rblk->dcr = (double *)malloc(size);
  rblk->dfl = (double *)malloc(size);
  rblk->dfc = (double *)malloc(size);
  rblk->dfr = (double *)malloc(size);
  rblk->cbl = (double *)malloc(size);
  rblk->cbc = (double *)malloc(size);
  rblk->cbr = (double *)malloc(size);
  rblk->ccl = (double *)malloc(size);
  rblk->ccc = (double *)malloc(size);
  rblk->ccr = (double *)malloc(size);
  rblk->cfl = (double *)malloc(size);
  rblk->cfc = (double *)malloc(size);
  rblk->cfr = (double *)malloc(size);
  rblk->ubl = (double *)malloc(size);
  rblk->ubc = (double *)malloc(size);
  rblk->ubr = (double *)malloc(size);
  rblk->ucl = (double *)malloc(size);
  rblk->ucc = (double *)malloc(size);
  rblk->ucr = (double *)malloc(size);
  rblk->ufl = (double *)malloc(size);
  rblk->ufc = (double *)malloc(size);
  rblk->ufr = (double *)malloc(size);

  if ( (rblk->ufr)== NULL)  printf("*****ERROR: allocMem out of memory \n");
*

}*/


//-------------- 
//  init
//-------------- 
void init(Domain_t *domain, RadiationDataAos1_t *rblk_aos1, RadiationDataAos2_t *rblk_aos2, RadiationDataAos3_t *rblk_aos3, RadiationData_t *rblk, double *x, double *b )
{
   int i;

   domain->imin = imin;
   domain->imax = imax;
   domain->jmin = jmin;
   domain->jmax = jmax;
   domain->kmin = kmin;   
   domain->kmax = kmax;
   domain->jp   = jp;
   domain->kp   = kp;

   for (i = 0 ; i < i_ub; i++ ) {
     *b++   = 0.0;
     rblk_aos1[i].dbl = i;
     rblk_aos1[i].dbc = i+1;
     rblk_aos1[i].dbr = i+2;
     rblk_aos1[i].dcl = i+3;
     rblk_aos1[i].dcc = i+4;
     rblk_aos1[i].dcr = i+5;
     rblk_aos1[i].dfl = i+6;
     rblk_aos1[i].dfc = i+7; 
  }

   for (i = 0 ; i < i_ub; i++ ) {
     rblk_aos2[i].dfr = i+8; 
     rblk_aos2[i].cbl = i+9;
     rblk_aos2[i].cbc = i+10;
     rblk_aos2[i].cbr = i+11;
     rblk_aos2[i].ccl = i+12;
     rblk_aos2[i].ccc = i+13;
     rblk_aos2[i].ccr = i+14;
     rblk_aos2[i].cfl = i+15; 
   }

   for (i = 0 ; i < i_ub; i++ ) {
     rblk_aos3[i].cfc = i+16; 
     rblk_aos3[i].cfr = i+17; 
     rblk_aos3[i].ubl = i+18; 
     rblk_aos3[i].ubc = i+19;
     rblk_aos3[i].ubr = i+20;
     rblk_aos3[i].ucl = i+21; 
     rblk_aos3[i].ucc = i+22; 
     rblk_aos3[i].ucr = i+23; 
   }

   double *ufl = rblk->ufl ;
   double *ufc = rblk->ufc ;
   double *ufr = rblk->ufr ;

   for (i = 0 ; i < i_ub; i++ ) {
     *ufl++ = i+24; 
     *ufc++ = i+25;
     *ufr++ = i+26;
   }

/*   double *dbl = rblk->dbl ;
   double *dbc = rblk->dbc ;
   double *dbr = rblk->dbr ;
   double *dcl = rblk->dcl ;
   double *dcc = rblk->dcc ;
   double *dcr = rblk->dcr ;
   double *dfl = rblk->dfl ;
   double *dfc = rblk->dfc ;
   double *dfr = rblk->dfr ;
   double *cbl = rblk->cbl ;
   double *cbc = rblk->cbc ;
   double *cbr = rblk->cbr ;
   double *ccl = rblk->ccl ;
   double *ccc = rblk->ccc ;
   double *ccr = rblk->ccr ;
   double *cfl = rblk->cfl ;
   double *cfc = rblk->cfc ;
   double *cfr = rblk->cfr ;
   double *ubl = rblk->ubl ;
   double *ubc = rblk->ubc ;
   double *ubr = rblk->ubr ;
   double *ucl = rblk->ucl ;
   double *ucc = rblk->ucc ;
   double *ucr = rblk->ucr ;
   double *ufl = rblk->ufl ;
   double *ufc = rblk->ufc ;
   double *ufr = rblk->ufr ;


   for (i = 0 ; i < i_ub; i++ ) {
     *b++   = 0.0;
     *dbl++ = i;
     *dbc++ = i+1;
     *dbr++ = i+2;
     *dcl++ = i+3;
     *dcc++ = i+4;
     *dcr++ = i+5;
     *dfl++ = i+6;
     *dfc++ = i+7; 
     *dfr++ = i+8; 
     *cbl++ = i+9;
     *cbc++ = i+10;
     *cbr++ = i+11;
     *ccl++ = i+12;
     *ccc++ = i+13;
     *ccr++ = i+14;
     *cfl++ = i+15; 
     *cfc++ = i+16; 
     *cfr++ = i+17; 
     *ubl++ = i+18; 
     *ubc++ = i+19;
     *ubr++ = i+20;
     *ucl++ = i+21; 
     *ucc++ = i+22; 
     *ucr++ = i+23; 
     *ufl++ = i+24; 
     *ufc++ = i+25;
     *ufr++ = i+26;
   }*/

   for (i=0; i<x_size; ++i ) {
     *x++ = 2.0*i;
   }

}
