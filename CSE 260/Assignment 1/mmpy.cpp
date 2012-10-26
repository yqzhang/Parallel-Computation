// 
// Matrix multiplication benchmark
// written by Scott B. Baden, UCSD, 4/6/11
// Modified with better command line option interface
//
// We compare the naive unblocked algorithm against the blocked version
// As an extra bonus, we also use the fast version of matrix multiplication
// found in the BLAS, e.g ATLAS, AMD ACML, Intel MKL
// (For best perforamnce in the library, we should not block, hence
// set the blocking factor = N)
// MKL will utilize all threads on the core, so this give us
// a not to exceed figure
//
// The benchmark repeats the matrix multiplication computation in order
// to help improve the accuracy of timings
// For values of N=512 or so, 5 repetitions shoudl be sufficient.
// 
// The program takes 3 command line arguments
// N, # repetitions, blocking factor
//

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <sys/time.h>
#include "types.h"
using namespace std;

double getTime();

#ifdef DEBUG
void printMatrix(Grid2D A,int n, char *msg)
#endif

#define fabs(x) ( (x)<0 ? -(x) : (x) )

extern void  cmdLine(int argc, char *argv[], int& n, int& b, int& do_blas, int& do_hand_blocked, int& do_unblocked, int& nreps);


extern void mm_unblocked(Grid2D A, Grid2D B, Grid2D C, int n, int nreps, double& time_un, _DOUBLE_ epsilon);

extern void mm_blas(Grid2D A, Grid2D B, Grid2D C, int n, int BSZ, int nreps, double& time_blas, _DOUBLE_ epsilon);

extern void mm_blocked(Grid2D A, Grid2D B, Grid2D C, Grid2D Aik, Grid2D Bkj, Grid2D Cij, int n, int BSZ, int nreps, double& time_blk, _DOUBLE_ epsilon);


Grid2D Alloc2D(int nx,int ny, const char *mesg){

   Grid2D U = (Grid2D) malloc(sizeof(_DOUBLE_ RESTRICT*)*ny + sizeof(_DOUBLE_ RESTRICT)*nx*ny);
   assert(U);
   if (!U)
      cerr << mesg << endl;
   else{
       int j;
       for(j=0;j<ny;j++)
           U[j] = ( _DOUBLE_  RESTRICT * ) (U+ny) + j*nx;
   }
   return(U);
}

int main(int argc, char **argv)
{

// command line arguments:
// N, # repetitions, blocking factor
int n, nreps, BSZ, do_blas, do_hand_blocked, do_unblocked;
cmdLine(argc, argv,  n, BSZ, do_blas, do_hand_blocked, do_unblocked, nreps);
if ((n%BSZ) != 0){
    fprintf(stderr,"\n *** The block size (%d) doesn't divide n (%d) evenly\n\n",BSZ,n);
    exit(-1);
}

// To avoid problems with round off, we consider two numbers
// to be equal if they differ by not more than epsilon
#ifdef _DOUBLE
_DOUBLE_ epsilon = 1.0e-8;
#else
_DOUBLE_ epsilon = 1.0e-4;
#endif

if (!do_blas && !do_hand_blocked && !do_unblocked){
    printf("No computations enabled.. you must have at least 1 enabled..  exiting.\n");
    exit(-1);
}

printf("\nn: %d\n",n);
printf("nreps: %d\n",nreps);
if (do_hand_blocked)
    printf("BlockSize: %d\n\n",BSZ);
printf("\nMatrix multiply variants:\n");
if (do_unblocked)
    printf("\tRunning unblocked\n");
if (do_blas)
    printf("\tRunning with BLAS to handle blocking\n");
if (do_hand_blocked){
        printf("\tRunning with manual blocking\n");
#ifdef COPY
        printf("\tCopying blocks before multiplying\n");
#endif
}

printf("\n");
#define RESTRICT 
// #define RESTRICT restrict
int i,j,k,r, ii, jj, kk;
// _DOUBLE_ **RESTRICT A=NULL,**RESTRICT B=NULL, **RESTRICT C=NULL;
// _DOUBLE_ **RESTRICT Aik=NULL, **RESTRICT Bkj=NULL, **RESTRICT Cij=NULL;
Grid2D A  = Alloc2D(n,n,"A");
Grid2D B = Alloc2D(n,n,"B");
Grid2D C = Alloc2D(n,n,"C");

#if 0
assert(A = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*n + sizeof(_DOUBLE_)*n*n));
for(i=0;i<n;i++) A[i] = (_DOUBLE_*)(A+n) + i*n;
assert(B = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*n + sizeof(_DOUBLE_)*n*n));
for(i=0;i<n;i++) B[i] = (_DOUBLE_*)(B+n) + i*n;
assert(C = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*n + sizeof(_DOUBLE_)*n*n));
for(i=0;i<n;i++) C[i] = (_DOUBLE_*)(C+n) + i*n;
#endif


Grid2D Aik=Alloc2D(BSZ,BSZ,"Aik");
Grid2D Bkj=Alloc2D(BSZ,BSZ,"Bkj");
Grid2D Cij=Alloc2D(BSZ,BSZ,"Cij");
#if 0
assert(Aik = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*BSZ + sizeof(_DOUBLE_)*BSZ*BSZ));
for(i=0;i<BSZ;i++) Aik[i] = (_DOUBLE_*)(Aik+BSZ) + i*BSZ;
assert(Bkj = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*BSZ + sizeof(_DOUBLE_)*BSZ*BSZ));
for(i=0;i<BSZ;i++) Bkj[i] = (_DOUBLE_*)(Bkj+BSZ) + i*BSZ;
assert(Cij = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*BSZ + sizeof(_DOUBLE_)*BSZ*BSZ));
for(i=0;i<BSZ;i++) Cij[i] = (_DOUBLE_*)(Cij+BSZ) + i*BSZ;
#endif

/* Generates a Hilbert Matrix H(i,j)
  H(I,j) = 1/(i+j+1)
  It's easy to check if the multiplication is correct;
  entry (i,j) of H * H is
  Sum(k) { 1.0/(i+k+1)*(k+j+1) }
 */

for (i=0; i<n; i++)
#pragma ivdep
    for (j=0; j<n; j++){
      B[i][j] = A[i][j] =  1.0 / (_DOUBLE_) (i+j+1);
      C[i][j] = 0.0;
    }
#ifdef DEBUG
// Print out the matrices
  char amsg[]= "--A--------";
  printMatrix("--A--------",A,n);
  char bmsg[]= "--B--------";
  printMatrix(B,n,bmsg);
#endif

_DOUBLE_ alpha = 1.0, d_one=1.0;

double time_blk,  time_un;

_DOUBLE_ sum;
if (do_unblocked){
  mm_unblocked(A, B, C, n, nreps, time_un, epsilon);
}

// Use the blas
// Clear out C so we can use it again
double time_blas;
if (do_blas){
  mm_blas(A, B, C, n, BSZ, nreps, time_blas, epsilon);
}

if (do_hand_blocked){
  mm_blocked(A, B, C, Aik, Bkj, Cij, n, BSZ, nreps, time_blk, epsilon);
}
 
time_un /= (double) nreps;
time_blk /= (double) nreps;

double tg_un = time_un, tg_blk = time_blk;
tg_un /= (double) n; tg_un /= (double) n; tg_un /= (double) n;
tg_blk /= (double) n; tg_blk /= (double) n; tg_blk /= (double) n;
printf("\n");
printf("\nTimes:\n");
if (do_unblocked)
    printf("   Unblocked: %f sec. [tg = %g]\n",time_un,tg_un);
else{
    time_un = tg_un = 0;
}
if (do_hand_blocked)
    printf("Hand Blocked: %f sec. [tg = %g]\n",time_blk, tg_blk );
else{
    time_blk = tg_blk=0;
}

time_blas /= (double) nreps;
double tg_blas = time_blas;
tg_blas /= (double) n; tg_blas /= (double) n; tg_blas /= (double) n;
if (do_blas)
    printf("        BLAS: %f sec. [tg = %g]\n",time_blas, tg_blas);
else{
    time_blas = tg_blas  = 0;
}
double flops = 2*n; flops *= n; flops *= n;
double gflops_un = (flops/time_un)/1.0e9;
double gflops_blk = (flops/time_blk)/1.0e9;
double gflops_blas = (flops/time_blas)/1.0e9;

printf("\nGflop rates:\n");
if (do_blas)
    printf("\tBLAS: %f\n",gflops_blas);
if (do_hand_blocked)
    printf("Hand blocked: %f\n",gflops_blk);
if (do_unblocked)
    printf("   Unblocked: %f\n",gflops_un);
// printf("blas %f, blocked %f, unblocked: %f\n",gflops_blas,gflops_blk, gflops_un);
printf("\n");
if (do_hand_blocked){
    printf("      N   BSZ  Reps    t_un    tg_un    t_blk    tg_blk   t_blas   tg_blas  cpy?\n");
    printf("@ %6d %4d %4d   %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e",n,BSZ,nreps,time_un, tg_un, time_blk, tg_blk, time_blas, tg_blas);
#ifdef COPY
    printf("  Y");
#else
    printf("  N");
#endif
}
else{
    printf("      N   BSZ  Reps    t_un    tg_un    t_blk    tg_blk  t_blas   tg_blas\n");
    printf("@ %6d      %4d   %8.2e %8.2e %8.2e %8.2e %8.2e %8.2e",n,BSZ,nreps,time_un, tg_un, time_blk, tg_blk, time_blas, tg_blas);
}
printf("\n\n");
}
