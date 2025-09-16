/*BHEADER****************************************************************
 * (c) 2006   The Regents of the University of California               *
 *                                                                      *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright       *
 * notice and disclaimer.                                               *
 *                                                                      *
 *EHEADER****************************************************************/

#include "irsmk.h"

void rmatmult3(Domain_t *domain, RadiationDataAos1_t *rblk_aos1, RadiationDataAos2_t *rblk_aos2, RadiationDataAos3_t *rblk_aos3, RadiationData_t *rblk, double *x, double *b )
{
   char *me = "rmatmult3" ;
   int i, ii, jj, kk ;
   int imin = domain->imin ;
   int imax = domain->imax ;
   int jmin = domain->jmin ;
   int jmax = domain->jmax ;
   int kmin = domain->kmin ;
   int kmax = domain->kmax ;
   int jp   = domain->jp   ;
   int kp   = domain->kp   ;
   /*double *dbl = rblk->dbl ;
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
   double *ucr = rblk->ucr ;*/
   double *ufl = rblk->ufl ;
   double *ufc = rblk->ufc ;
   double *ufr = rblk->ufr ;
   double *xdbl = x - kp - jp - 1 ;
   double *xdbc = x - kp - jp     ;
   double *xdbr = x - kp - jp + 1 ;
   double *xdcl = x - kp      - 1 ;
   double *xdcc = x - kp          ;
   double *xdcr = x - kp      + 1 ;
   double *xdfl = x - kp + jp - 1 ;
   double *xdfc = x - kp + jp     ;
   double *xdfr = x - kp + jp + 1 ;
   double *xcbl = x      - jp - 1 ;
   double *xcbc = x      - jp     ;
   double *xcbr = x      - jp + 1 ;
   double *xccl = x           - 1 ;
   double *xccc = x               ;
   double *xccr = x           + 1 ;
   double *xcfl = x      + jp - 1 ;
   double *xcfc = x      + jp     ;
   double *xcfr = x      + jp + 1 ;
   double *xubl = x + kp - jp - 1 ;
   double *xubc = x + kp - jp     ;
   double *xubr = x + kp - jp + 1 ;
   double *xucl = x + kp      - 1 ;
   double *xucc = x + kp          ;
   double *xucr = x + kp      + 1 ;
   double *xufl = x + kp + jp - 1 ;
   double *xufc = x + kp + jp     ;
   double *xufr = x + kp + jp + 1 ;


   for ( kk = kmin ; kk < kmax ; kk++ ) {
      for ( jj = jmin ; jj < jmax ; jj++ ) {
         for ( ii = imin ; ii < imax ; ii++ ) {
	    i = ii + jj * jp + kk * kp ;
           /* b[i] = dbl[i] * xdbl[i] + dbc[i] * xdbc[i] + dbr[i] * xdbr[i] +
                   dcl[i] * xdcl[i] + dcc[i] * xdcc[i] + dcr[i] * xdcr[i] +
                   dfl[i] * xdfl[i] + dfc[i] * xdfc[i] + dfr[i] * xdfr[i] +
                   cbl[i] * xcbl[i] + cbc[i] * xcbc[i] + cbr[i] * xcbr[i] +
                   ccl[i] * xccl[i] + ccc[i] * xccc[i] + ccr[i] * xccr[i] +
                   cfl[i] * xcfl[i] + cfc[i] * xcfc[i] + cfr[i] * xcfr[i] +
                   ubl[i] * xubl[i] + ubc[i] * xubc[i] + ubr[i] * xubr[i] +
                   ucl[i] * xucl[i] + ucc[i] * xucc[i] + ucr[i] * xucr[i] +
                   ufl[i] * xufl[i] + ufc[i] * xufc[i] + ufr[i] * xufr[i] ;*/
            b[i] = rblk_aos1[i].dbl * xdbl[i] + rblk_aos1[i].dbc * xdbc[i] + rblk_aos1[i].dbr * xdbr[i] +
                   rblk_aos1[i].dcl * xdcl[i] + rblk_aos1[i].dcc * xdcc[i] + rblk_aos1[i].dcr * xdcr[i] +
                   rblk_aos1[i].dfl * xdfl[i] + rblk_aos1[i].dfc * xdfc[i] + rblk_aos2[i].dfr * xdfr[i] +
                   rblk_aos2[i].cbl * xcbl[i] + rblk_aos2[i].cbc * xcbc[i] + rblk_aos2[i].cbr * xcbr[i] +
                   rblk_aos2[i].ccl * xccl[i] + rblk_aos2[i].ccc * xccc[i] + rblk_aos2[i].ccr * xccr[i] +
                   rblk_aos2[i].cfl * xcfl[i] + rblk_aos3[i].cfc * xcfc[i] + rblk_aos3[i].cfr * xcfr[i] +
                   rblk_aos3[i].ubl * xubl[i] + rblk_aos3[i].ubc * xubc[i] + rblk_aos3[i].ubr * xubr[i] +
                   rblk_aos3[i].ucl * xucl[i] + rblk_aos3[i].ucc * xucc[i] + rblk_aos3[i].ucr * xucr[i] +
                   ufl[i] * xufl[i] + ufc[i] * xufc[i] + ufr[i] * xufr[i] ;
	 }
      }
   }

}
