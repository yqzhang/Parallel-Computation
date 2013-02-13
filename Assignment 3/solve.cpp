/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#ifdef _MPI_
#include "mpi.h"
#endif
#include "types.h"
using namespace std;

// globals
extern int px, py, noComm;

void printT(_DOUBLE_ **E, int m, int n) {
  for(int i = 0; i < m + 2; i++) {
    for(int j = 0; j < n + 2; j++) {
      cout << E[i][j] << " ";
    }
    cout << "\n";
  }
}

void repNorms(ofstream& logfile, _DOUBLE_ **E, _DOUBLE_ dt, int m,int n, int niter, int stats_freq);


// Reports statistics about the computation: the L2 Norm and the Infinity NOrm
// These values should not vary (except to within roundoff)
// when we use different numbers of  processes to solve the problem


// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
// The Linf norm is simply the maximum (absolute) value over
// all points in the array

 _DOUBLE_ stats(_DOUBLE_ **E, int m, int n, _DOUBLE_ *_mx){
     _DOUBLE_ mx = -1;
     _DOUBLE_ l2norm = 0;
     int i, j;
     for (j=1; j<=m+1; j++)
       for (i=1; i<=n+1; i++) {
	   l2norm += E[j][i]*E[j][i];
	   _DOUBLE_ fe = fabs(E[j][i]);
	   if (fe > mx)
	       mx = fe;
      }

     l2norm /= (_DOUBLE_) ((m+1)*(n+1));
     l2norm = sqrt(l2norm);

     *_mx = mx;
     return l2norm;
 }

int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, int niters, _DOUBLE_ alpha, _DOUBLE_ dt, int plot_freq, int stats_freq){

 // Simulated time is different from the integer timestep number
 _DOUBLE_ t = 0.0;

 _DOUBLE_ **E = *_E, **E_prev = *_E_prev;

 _DOUBLE_ *buf1, *buf2, *buf3, *buf4;
 int niter;
 int x, y;
 int nproc, myrank;

#ifdef _MPI_
  // malloc buffer(s)

  buf1 = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * m);
  buf2 = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * m);
  buf3 = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * m);
  buf4 = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * m);
#endif

 // Calculate my position in the grid
#ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&nproc);
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  nproc = 1;
  myrank = 0;
#endif

 y = myrank / px;
 x = myrank % px;

 // We continue to sweep over the mesh until the simulation has reached
 // the desired simulation Time
 // This is different from the number of iterations
  for (niter = 0; niter < niters; niter++){

#ifdef DEBUG
   printMat(E_prev,m,n);
   repNorms(logfile,E_prev,dt,m,n,niter, stats_freq);
   if (plot_freq)
#ifdef PLOT1D
     splot1d(&(E_prev[0][0]),niter,m+1,n+1,WAIT);
#else
    splot(E_prev,niter,m+1,n+1,WAIT);
#endif
#endif

   /* 
    * Copy data from boundary of the computational box to the
    * padding region, set up for differencing computational box's boundary
    *
    * These are physical boundary conditions, and are not to be confused
    * with ghost cells that we would use in an MPI implementation
    *
    * The reason why we copy boundary conditions is to avoid
    * computing single sided differences at the boundaries
    * which increase the running time of solve()
    *
    */
    
   int i,j;
   MPI_Request request1, request2, request3, request4;
   MPI_Request srequest[4];
   int sendCount = 0;

   if (!noComm) {

   if (y != 0) {
     MPI_Irecv(E_prev[0]+1, n, MPI_DOUBLE, myrank - px,
              niter, MPI_COMM_WORLD, &request3);
     MPI_Isend(E_prev[1]+1, n, MPI_DOUBLE, myrank - px, niter, MPI_COMM_WORLD, &srequest[sendCount++]);
   }

   if (y != py-1) {
     MPI_Irecv(E_prev[m+1]+1, n, MPI_DOUBLE, myrank + px,
              niter, MPI_COMM_WORLD, &request4);
     MPI_Isend(E_prev[m]+1, n, MPI_DOUBLE, myrank + px, niter, MPI_COMM_WORLD, &srequest[sendCount++]);
   }

   if (x != 0) {
     MPI_Irecv(buf1, m, MPI_DOUBLE, myrank-1,
              niter, MPI_COMM_WORLD, &request1);
     for (i = 1; i <= m; i++) {
       buf3[i-1] = E_prev[i][1];
     }
     MPI_Isend(buf3, m, MPI_DOUBLE, myrank-1, niter, MPI_COMM_WORLD, &srequest[sendCount++]);
   }

   if (x != px-1) {
     MPI_Irecv(buf2, m, MPI_DOUBLE, myrank+1,
              niter, MPI_COMM_WORLD, &request2);
     for (i = 1; i <= m; i++) {
       buf4[i-1] = E_prev[i][n];
     }
     MPI_Isend(buf4, m, MPI_DOUBLE, myrank+1, niter, MPI_COMM_WORLD, &srequest[sendCount++]);
   }

   }

   if (x == 0 || noComm) {
     // leftmost col
     for (j=1; j<=m; j++) 
       E_prev[j][0] = E_prev[j][2];
   } else {
     MPI_Wait(&request1, MPI_STATUS_IGNORE);
     for (j = 1; j <= m+1; j++) {
       E_prev[j][0] = buf1[j-1];
     }
   }

   if (x == px-1 || noComm) {
     // rightmost col
     for (j=1; j<=m; j++) 
       E_prev[j][n+1] = E_prev[j][n-1];
   } else {
     MPI_Wait(&request2, MPI_STATUS_IGNORE);
     for (j = 1; j <= m+1; j++) {
       E_prev[j][n+1] = buf2[j-1];
     }
   }
 
   if (y == 0 || noComm) {
     // top row
     for (i=1; i<=n; i++) 
       E_prev[0][i] = E_prev[2][i];
   } else {
     MPI_Wait(&request3, MPI_STATUS_IGNORE);
   }

   if (y == py-1 || noComm) {
     // bottom row
     for (i=1; i<=n; i++) 
       E_prev[m+1][i] = E_prev[m-1][i];
   } else {
     MPI_Wait(&request4, MPI_STATUS_IGNORE);
   }

   // Solve for the excitation, a PDE
   for (j=1; j<=m; j++){
     for (i=1; i<=n; i++) {
	E[j][i] = E_prev[j][i]+alpha*(E_prev[j][i+1]+E_prev[j][i-1]-4*E_prev[j][i]+E_prev[j+1][i]+E_prev[j-1][i]);
     }
    }

   /* 
    * Solve the ODE, advancing excitation and recovery variables
    *     to the next timtestep
    */
   for (j=1; j<=m; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (i=1; i<=n; i++, EE++,RR++)
	EE[0] += -dt*(kk*EE[0]*(EE[0]-a)*(EE[0]-1)+EE[0]*RR[0]);
   }

   for (j=1; j<=m; j++){
     _DOUBLE_ *RR = &R[j][1];
     _DOUBLE_ *EE = &E[j][1];
     for (i=1; i<=n; i++, EE++, RR++)
	RR[0] += dt*(epsilon+M1* RR[0]/( EE[0]+M2))*(-RR[0]-kk*EE[0]*(EE[0]-b-1));
   }

   if (stats_freq)
     repNorms(logfile, E,dt,m,n,niter,stats_freq);

   if (plot_freq){
          if (!(niter % plot_freq)){
#ifdef PLOT1D
//            splot1d(&(E[0][0]),t,niter,m+1,n+1,WAIT);
            splot1d(&(E[0][0]),niter,m,n,WAIT);
#else
//            splot(E,t,niter,m+1,n+1,WAIT);
            splot(E,niter,m,n,WAIT);
#endif
        }
    }

  // Swap current and previous
  //MPI_Barrier(MPI_COMM_WORLD);
  if (!noComm)
    MPI_Waitall(sendCount, srequest, MPI_STATUSES_IGNORE);

   _DOUBLE_ **tmp = E; E = E_prev; E_prev = tmp;
 }

  // Store them into the pointers passed in
  *_E = E;
  *_E_prev = E_prev;

  return niter;
}
