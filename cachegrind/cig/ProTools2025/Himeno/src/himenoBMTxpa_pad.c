/********************************************************************

 This benchmark test program is measuring a cpu performance
 of floating point operation by a Poisson equation solver.

 If you have any question, please ask me via email.
 written by Ryutaro HIMENO, November 26, 2001.
 Version 3.0
 ----------------------------------------------
 Ryutaro Himeno, Dr. of Eng.
 Head of Computer Information Division,
 RIKEN (The Institute of Pysical and Chemical Research)
 Email : himeno@postman.riken.go.jp
 ---------------------------------------------------------------
 You can adjust the size of this benchmark code to fit your target
 computer. In that case, please chose following sets of
 [mimax][mjmax][mkmax]:
 small : 33,33,65
 small : 65,65,129
 midium: 129,129,257
 large : 257,257,513
 ext.large: 513,513,1025
 This program is to measure a computer performance in MFLOPS
 by using a kernel which appears in a linear solver of pressure
 Poisson eq. which appears in an incompressible Navier-Stokes solver.
 A point-Jacobi method is employed in this solver as this method can 
 be easyly vectrized and be parallelized.
 ------------------
 Finite-difference method, curvilinear coodinate system
 Vectorizable and parallelizable on each grid point
 No. of grid points : imax x jmax x kmax including boundaries
 ------------------
 A,B,C:coefficient matrix, wrk1: source term of Poisson equation
 wrk2 : working area, OMEGA : relaxation parameter
 BND:control variable for boundaries and objects ( = 0 or 1)
 P: pressure
********************************************************************/

#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

#include "util.h"

#define NUM_THREADS 1
#define MR(mt,n,r,c,d)  mt->m[(n) * mt->mrows * mt->mcols * mt->mdeps + (r) * mt->mcols* mt->mdeps + (c) * mt->mdeps + (d)]

struct Mat {
  float* buffer;
  float* m;
  int mnums;
  int mrows;
  int mcols;
  int mdeps;
};

/* prototypes */
typedef struct Mat Matrix;

int newMat(Matrix* Mat, int mnums, int mrows, int mcols, int mdeps, int pad);
void clearMat(Matrix* Mat);
bool set_param(int i[],char *size);
void mat_set(Matrix* Mat,int l,float z);
void mat_set_init(Matrix* Mat);
float jacobi(int n,Matrix* M1,Matrix* M2,Matrix* M3,
             Matrix* M4,Matrix* M5,Matrix* M6,Matrix* M7, Matrix* M8, Matrix* M9, Matrix* M10, Matrix* M11, Matrix* M12, Matrix* M13, Matrix* m14);
double second();

float   omega=0.8;
//Matrix  a,b,c,p,bnd,wrk1,wrk2;
Matrix  a0, a1, a2, a3, b0, b1, b2, c0, c1, c2,p,bnd,wrk1,wrk2;

void
print_usage(char **argv)
{
    printf("Usage:  %s -t num_threads -s [grid-size = XS/S/M/L/XL] -l num_loops -v varinfo_file\n", *argv);
    printf(" Grid-size= XS (32x32x64)\n");
    printf("\t    S  (64x64x128)\n");
    printf("\t    M  (128x128x256)\n");
    printf("\t    L  (256x256x512)\n");
    printf("\t    XL (512x512x1024)\n\n");
}

int
main(int argc, char *argv[])
{
  int num_threads = NUM_THREADS;
  int    nn;
  int    imax,jmax,kmax,mimax,mjmax,mkmax,msize[3] = {0, 0, 0};
  float  gosa,target;
  double  cpu0,cpu1,cpu,xmflops2,score,flop;

  char   size[32] = "XS";
  int  num_loops = 60;
  char *varinfo_fname = NULL;

/*  if(argc == 2){
    strcpy(size,argv[1]);
  } else {
    printf("For example: \n");
    printf(" Grid-size= XS (32x32x64)\n");
    printf("\t    S  (64x64x128)\n");
    printf("\t    M  (128x128x256)\n");
    printf("\t    L  (256x256x512)\n");
    printf("\t    XL (512x512x1024)\n\n");
    printf("Grid-size = ");
    scanf("%s",size);
    printf("\n");
  } */

  int cha;
  while ((cha = getopt (argc, argv, "t:s:l:v:h")) != -1)
    switch (cha)
    {
      case 't':
        num_threads    = atol(optarg);
        break;
      case 'l':
        num_loops = atoi(optarg);
        break;
      case 's':
        strcpy(size, optarg);
        break;
      case 'v':
        varinfo_fname = optarg;
        break;
      case 'h':
      default:
        print_usage(argv);
        exit(0);
    }

  if(!set_param(msize,size))
  {
    fprintf(stderr, "ERROR: wrong grid size!\n");
    print_usage(argv);
    exit(0);
  }

  mimax= msize[0];
  mjmax= msize[1];
  mkmax= msize[2];
  imax= mimax-1;
  jmax= mjmax-1;
  kmax= mkmax-1;

  target = 60.0;

  printf("mimax = %d mjmax = %d mkmax = %d\n",mimax,mjmax,mkmax);
  printf("imax = %d jmax = %d kmax =%d\n",imax,jmax,kmax);

  /*
   *    Initializing matrixes
   */
  newMat(&p,1,mimax,mjmax,mkmax, 0);
  newMat(&bnd,1,mimax,mjmax,mkmax, 113);
  newMat(&wrk1,1,mimax,mjmax,mkmax, 129);
  newMat(&wrk2,1,mimax,mjmax,mkmax, 113);
  newMat(&a0,1,mimax,mjmax,mkmax, 147);
  newMat(&a1,1,mimax,mjmax,mkmax, 167);
  newMat(&a2,1,mimax,mjmax,mkmax, 183);
  newMat(&a3,1,mimax,mjmax,mkmax, 17);
  newMat(&b0,1,mimax,mjmax,mkmax, 17);
  newMat(&b1,1,mimax,mjmax,mkmax, 33);
  newMat(&b2,1,mimax,mjmax,mkmax, 51);
  newMat(&c0,1,mimax,mjmax,mkmax, 67);
  newMat(&c1,1,mimax,mjmax,mkmax, 79);
  newMat(&c2,1,mimax,mjmax,mkmax, 97);
/*  newMat(&a,4,mimax,mjmax,mkmax, 17);
  newMat(&b,3,mimax,mjmax,mkmax, 33);
  newMat(&c,3,mimax,mjmax,mkmax, 51);*/

  if(varinfo_fname != NULL)
  {
    int res = varinfo_file_open(varinfo_fname);
    if(res)
      fprintf(stderr, "Cannot open %s for recording varinfo\n", varinfo_fname);

    varinfo_file_print(getName(p), &MR((&p), 0, 0, 0, 0), &MR((&p), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(bnd), &MR((&bnd), 0, 0, 0, 0), &MR((&bnd), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(wrk1), &MR((&wrk1), 0, 0, 0, 0), &MR((&wrk1), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(wrk2), &MR((&wrk2), 0, 0, 0, 0), &MR((&wrk2), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(a0), &MR((&a0), 0, 0, 0, 0), &MR((&a0), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(a1), &MR((&a1), 0, 0, 0, 0), &MR((&a1), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(a2), &MR((&a2), 0, 0, 0, 0), &MR((&a2), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(a3), &MR((&a3), 0, 0, 0, 0), &MR((&a3), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(b0), &MR((&b0), 0, 0, 0, 0), &MR((&b0), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(b1), &MR((&b1), 0, 0, 0, 0), &MR((&b1), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(b2), &MR((&b2), 0, 0, 0, 0), &MR((&b2), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(c0), &MR((&c0), 0, 0, 0, 0), &MR((&c0), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(c1), &MR((&c1), 0, 0, 0, 0), &MR((&c1), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_print(getName(c2), &MR((&c2), 0, 0, 0, 0), &MR((&c2), 0, imax-1, jmax-1, kmax-1));
    varinfo_file_close();
  }
 
  mat_set_init(&p);
  mat_set(&bnd,0,1.0);
  mat_set(&wrk1,0,0.0);
  mat_set(&wrk2,0,0.0);
  mat_set(&a0,0,1.0);
  mat_set(&a1,0,1.0);
  mat_set(&a2,0,1.0);
  mat_set(&a3,0,1.0/6.0);
  mat_set(&b0,0,0.0);
  mat_set(&b1,0,0.0);
  mat_set(&b2,0,0.0);
  mat_set(&c0,0,1.0);
  mat_set(&c1,0,1.0);
  mat_set(&c2,0,1.0);

  omp_set_num_threads(num_threads);

  /*
   *    Start measuring
   */
  nn= 3;
  printf(" Start rehearsal measurement process.\n");
  printf(" Measure the performance in %d times.\n\n",nn);

  cpu0= second();
//  gosa= jacobi(nn,&a,&b,&c,&p,&bnd,&wrk1,&wrk2);
  gosa= jacobi(nn,&a0,&a1,&a2,&a3,&b0,&b1,&b2,&c0,&c1,&c2,&p,&bnd,&wrk1,&wrk2);
  cpu1= second();
  cpu= cpu1 - cpu0;
  flop = (double)(kmax-1)*(double)(jmax-1)*(double)(imax-1)*34.0;

  if(cpu != 0.0)
    xmflops2= flop/cpu*1.e-6*nn;

  printf(" MFLOPS: %f time(s): %f %e\n\n",xmflops2,cpu,gosa);

  //nn= (int)(target/(cpu/3.0));
  nn= num_loops;

  printf(" Now, start the actual measurement process.\n");
  printf(" The loop will be excuted in %d times\n",nn);
  printf(" This will take about one minute.\n");
  printf(" Wait for a while\n\n");

  cpu0 = second();
  gosa= jacobi(nn,&a0,&a1,&a2,&a3,&b0,&b1,&b2,&c0,&c1,&c2,&p,&bnd,&wrk1,&wrk2);
  cpu1 = second();
  cpu = cpu1 - cpu0;

  if(cpu != 0.0)
    xmflops2 = (double)flop/cpu*1.0e-6*nn;

  printf("cpu : %f sec.\n", cpu);
  printf("Loop executed for %d times\n",nn);
  printf("Gosa : %e \n",gosa);
  printf("MFLOPS measured : %f\n",xmflops2);
  score = xmflops2/82.84;
  printf("Score based on Pentium III 600MHz using Fortran 77: %f\n",score);

  /*
   *   Matrix free
   */ 
  clearMat(&p);
  clearMat(&bnd);
  clearMat(&wrk1);
  clearMat(&wrk2);
  clearMat(&a0);
  clearMat(&a1);
  clearMat(&a2);
  clearMat(&a3);
  clearMat(&b0);
  clearMat(&b1);
  clearMat(&b2);
  clearMat(&c0);
  clearMat(&c1);
  clearMat(&c2);
  
  return (0);
}

bool
set_param(int is[],char *size)
{
  if(!strcmp(size,"XS") || !strcmp(size,"xs")){
    is[0]= 32;
    is[1]= 32;
    is[2]= 64;
    return true;
  }
  if(!strcmp(size,"S") || !strcmp(size,"s")){
    is[0]= 64;
    is[1]= 64;
    is[2]= 128;
    return true;
  }
  if(!strcmp(size,"M") || !strcmp(size,"m")){
    is[0]= 128;
    is[1]= 128;
    is[2]= 256;
    return true;
  }
  if(!strcmp(size,"L") || !strcmp(size,"l")){
    is[0]= 256;
    is[1]= 256;
    is[2]= 512;
    return true;
  }
  if(!strcmp(size,"XL") || !strcmp(size,"xl")){
    is[0]= 512;
    is[1]= 512;
    is[2]= 1024;
    return true;
  }
  return false;
}

int
newMat(Matrix* Mat, int mnums,int mrows, int mcols, int mdeps, int pad)
{
  Mat->mnums= mnums;
  Mat->mrows= mrows;
  Mat->mcols= mcols;
  Mat->mdeps= mdeps;
  Mat->buffer= NULL;
  Mat->buffer= (float*) 
    malloc((mnums * mrows * mcols * mdeps + pad) * sizeof(float));
  
  Mat->m = Mat->buffer + pad;

  return(Mat->buffer != NULL) ? 1:0;
}

void
clearMat(Matrix* Mat)
{
  if(Mat->buffer)
    free(Mat->buffer);
  Mat->buffer= NULL;
  Mat->m= NULL;
  Mat->mnums= 0;
  Mat->mcols= 0;
  Mat->mrows= 0;
  Mat->mdeps= 0;
}

void
mat_set(Matrix* Mat, int l, float val)
{
  int i,j,k;

    for(i=0; i<Mat->mrows; i++)
      for(j=0; j<Mat->mcols; j++)
        for(k=0; k<Mat->mdeps; k++)
          MR(Mat,l,i,j,k)= val;
}

void
mat_set_init(Matrix* Mat)
{
  int  i,j,k,l;
  float tt;

  for(i=0; i<Mat->mrows; i++)
    for(j=0; j<Mat->mcols; j++)
      for(k=0; k<Mat->mdeps; k++)
        MR(Mat,0,i,j,k)= (float)(i*i)
          /(float)((Mat->mrows - 1)*(Mat->mrows - 1));
}

float
jacobi(int nn, Matrix* a0, Matrix* a1, Matrix* a2, Matrix* a3, Matrix* b0, Matrix* b1, Matrix*b2, Matrix* c0, Matrix*c1, Matrix*c2,
       Matrix* p,Matrix* bnd,Matrix* wrk1,Matrix* wrk2)
{
  int    i,j,k,n,imax,jmax,kmax;
  float  gosa,gosa1,s0,ss;

  imax= p->mrows-1;
  jmax= p->mcols-1;
  kmax= p->mdeps-1;
#pragma omp parallel shared(a0,a1,a2,a3,p,b0,b1,b2,c0,c1,c2,bnd,wrk1,wrk2,nn,imax,jmax,kmax,omega,gosa) private(i,j,k,s0,ss,gosa1,n)
{
  for(n=0 ; n<nn ; n++){
#pragma omp barrier
#pragma omp master
    {
      gosa = 0.0;
    }
    gosa1= 0.0;
#pragma omp for nowait
    for(i=1 ; i<imax; i++)
      for(j=1 ; j<jmax ; j++)
        for(k=1 ; k<kmax ; k++){
          s0= MR(a0,0,i,j,k)*MR(p,0,i+1,j,  k)
            + MR(a1,0,i,j,k)*MR(p,0,i,  j+1,k)
            + MR(a2,0,i,j,k)*MR(p,0,i,  j,  k+1)
            + MR(b0,0,i,j,k)
             *( MR(p,0,i+1,j+1,k) - MR(p,0,i+1,j-1,k)
              - MR(p,0,i-1,j+1,k) + MR(p,0,i-1,j-1,k) )
            + MR(b1,0,i,j,k)
             *( MR(p,0,i,j+1,k+1) - MR(p,0,i,j-1,k+1)
              - MR(p,0,i,j+1,k-1) + MR(p,0,i,j-1,k-1) )
            + MR(b2,0,i,j,k)
             *( MR(p,0,i+1,j,k+1) - MR(p,0,i-1,j,k+1)
              - MR(p,0,i+1,j,k-1) + MR(p,0,i-1,j,k-1) )
            + MR(c0,0,i,j,k) * MR(p,0,i-1,j,  k)
            + MR(c1,0,i,j,k) * MR(p,0,i,  j-1,k)
            + MR(c2,0,i,j,k) * MR(p,0,i,  j,  k-1)
            + MR(wrk1,0,i,j,k);

          ss= (s0*MR(a3,0,i,j,k) - MR(p,0,i,j,k))*MR(bnd,0,i,j,k);

          gosa1+= ss*ss;

          MR(wrk2,0,i,j,k)= MR(p,0,i,j,k) + omega*ss;
        }
#pragma omp barrier
#pragma omp for nowait
    for(i=1 ; i<imax ; i++)
      for(j=1 ; j<jmax ; j++)
        for(k=1 ; k<kmax ; k++)
          MR(p,0,i,j,k)= MR(wrk2,0,i,j,k);

#pragma omp critical
    {
      gosa+= gosa1;
    }
  } /* end n loop */
}
  return(gosa);
}

double
second()
{

  struct timeval tm;
  double t ;

  static int base_sec = 0,base_usec = 0;

  gettimeofday(&tm, NULL);
  
  if(base_sec == 0 && base_usec == 0)
    {
      base_sec = tm.tv_sec;
      base_usec = tm.tv_usec;
      t = 0.0;
  } else {
    t = (double) (tm.tv_sec-base_sec) + 
      ((double) (tm.tv_usec-base_usec))/1.0e6 ;
  }

  return t ;
}

