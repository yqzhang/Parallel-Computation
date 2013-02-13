/* 
 * Driver for a cardiac elecrophysioly simulatin that uses the
 * Aliev-Panfilov model
 * We use an explicit method
 *
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 */

#include <assert.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <math.h>
#include "time.h"
#include "apf.h"
#ifdef _MPI_
#include "mpi.h"
#endif
#ifdef  _OPENMP
// The "_OPENMP" Macro is defined if you have enabled Openmp via the compiler
#include <omp.h>
#endif
using namespace std;

// Global variables
int px, py, noComm;

#ifndef _MPI_
double getTime();
#endif

// Utilities
// 

// Allocate a 2D array

void print(_DOUBLE_ **E, int m, int n) {
  for(int i = 0; i < m + 2; i++) {
    for(int j = 0; j < n + 2; j++) {
      cout << E[i][j] << " ";
    }
    cout << "\n";
  }
}

_DOUBLE_ **alloc2D(int m,int n){

   _DOUBLE_ **E;
   int nx=n+1, ny=m+1;
   E = (_DOUBLE_**)malloc(sizeof(_DOUBLE_*)*ny + sizeof(_DOUBLE_)*nx*ny);
   assert(E);
   int j;
   for(j=0;j<ny;j++) E[j] = (_DOUBLE_*)(E+ny) + j*nx;
   return(E);
}

void init (_DOUBLE_ **E,_DOUBLE_ **E_prev,_DOUBLE_ **R,int m,int n){
  int nproc, myrank, x, y, myM, myN, startx, starty;

#ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&nproc);
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
 nproc = 1;
 myrank = 0;
#endif

 y = myrank / px;
 x = myrank % px;

 starty = (y < (m + 1) % py) ? (y * ((m + 1 + py - 1) / py) + 1) : (y * ((m + 1) / py) + (m + 1) % py + 1);
 startx = (x < (n + 1) % px) ? (x * ((n + 1 + px - 1) / px) + 1) : (x * ((n + 1) / px) + (n + 1) % px + 1);

 // Calculate the dimensions of the local matrices
 myM = ((m+1 + py - 1) / py) + ((y < (m+1) % py) ? 1 : 0) - 1;
 myN = ((n+1 + px - 1) / px) + ((x < (n+1) % px) ? 1 : 0) - 1;

 //cout << "Rank: " << myrank << " startY:" << starty << " startX:" << startx << " M:" << myM << " N:" << myN << endl;

 int i,j;
 // Initialization
 for (j=1; j<=myM; j++) {
     for (i=1; i<= myN; i++) {
         E_prev[j][i] = R[j][i] = 0;
         if(i + startx - 1 >= n/2 + 2) {
             E_prev[j][i] = 1.0;
         }
         if(j + starty - 1>= m/2 + 2) {
             R[j][i] = 1.0;
         }
     }
 }
}

// External functions
void cmdLine(int argc, char *argv[], int& n, int& stats_freq, int& plot_freq, int& px, int& py, int &noComm, int &niters);
int solve(ofstream& logfile, _DOUBLE_ ***_E, _DOUBLE_ ***_E_prev, _DOUBLE_ **R, int m, int n, int niters, _DOUBLE_ alpha, _DOUBLE_ dt, int plot_freq, int stats_freq);
void printTOD(ofstream& logfile, string mesg);
void ReportStart(ofstream& logfile, _DOUBLE_ dt, int niters, int m, int n, int px, int py, int noComm);
void ReportEnd(ofstream& logfile, int niters, int niter, _DOUBLE_ **E_prev, int m,int n, double t0, int px, int py);

static _DOUBLE_ **res;

// Main program
int main(int argc, char** argv)
{
 /*
  *  Solution arrays
  *   E is the "Excitation" variable, a voltage
  *   R is the "Recovery" variable
  *   E_prev is the Excitation variable for the previous timestep,
  *      and is used in time integration
  */
 _DOUBLE_ **E, **R, **E_prev; // local
 //_DOUBLE_ **res; // global result (for E?)

 // Default values for the command line arguments
 int m=100,n=100;
 int stats_freq = 0;
 int plot_freq = 0;
 //int px = 1, py = 1;
 int niters=100;

 px = 1;
 py = 1;
 noComm = 0;

#ifdef _MPI_
 MPI_Init(&argc,&argv);
#endif

// Parse command line arguments
 cmdLine( argc, argv, n, stats_freq,  plot_freq, px, py, noComm, niters);
 if (n < 26){
    cout << "\n *** N must be larger than 25.  Exiting ... " << endl << endl;
    exit(-1);
 }
 m = n;

 int nproc, myrank, i;
 int x,y;
 int myM, myN; // size of the local matrices
#ifdef _MPI_
 MPI_Comm_size(MPI_COMM_WORLD,&nproc);
 MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  nproc = 1;
  myrank = 0;
#endif
 y = myrank / px;
 x = myrank % px;

#ifdef _MPI_
 if ((px * py) != nproc){
    if (!myrank)
        cout << "\n *** The number of  processors in the geometry (" << px*py << ")  is not the same as the number requested (" << nproc << ")" << endl << endl;
    Stop();
 }
#else
    if (!myrank)
    if ((px * py) > 1){
        cout << "\n *** The number of  processors in the geometry > 1, but you have not enabled MPI\n";
        Stop();
    }
#endif

 // The log file
 // Do not change the file name or remove this call
   ofstream logfile("Log.txt",ios::out);
   printTOD(logfile, "Simulation begins");

 // Calculate the dimensions of the local matrices

 myM = ((m+1 + py - 1) / py) + ((y < (m+1) % py) ? 1 : 0) - 1;
 myN = ((n+1 + px - 1) / px) + ((x < (n+1) % px) ? 1 : 0) - 1;

 // Allocate contiguous memory for solution arrays
 // The computational box is defined on [1:m+1,1:n+1]
 // We pad the arrays in order to facilitate differencing on the 
 // boundaries of the computation box
 E = alloc2D(myM+1,myN+1);
 E_prev = alloc2D(myM+1,myN+1);
 R = alloc2D(myM+1,myN+1);

 if (!myrank) {
   res = alloc2D(m+2, n+2);
 }

 init(E,E_prev,R,m,n);
 //
 // Initization of various simulation variables
 // Do not change the code these assignments statements, as doing so
 // could cause your submission to be graded incorrectly
 //

 // We compute the timestep dt which determines how long 
 // the code will run for

 // We should always use double precision values for the folowing variables:
 //    rp, dte, dtr, ddt
 //
 // This ensures that the computation of dte and especially dt
 // will not lose precision (i.e. if computed as single precision values)

 _DOUBLE_ dx = 1.0/n;
 double rp= kk*(b+1)*(b+1)/4;
 double dte=(dx*dx)/(d*4+((dx*dx))*(rp+kk));
 double dtr=1/(epsilon+((M1/M2)*rp));
 double ddt = (dte<dtr) ? 0.95*dte : 0.95*dtr;
 _DOUBLE_ dt = (_DOUBLE_) ddt;
 _DOUBLE_ alpha = d*dt/(dx*dx);

 // End Initization of various simulation variables

 // Report various information
 // Do not remove this call, it is needed for grading
  ReportStart(logfile, dt, niters, m, n, px, py, noComm);

 // Start the timer
#ifdef _MPI_
 double t0 = -MPI_Wtime();
#else
 double t0 = -getTime();
#endif
 int niter = solve(logfile, &E, &E_prev, R, myM, myN, niters, alpha, dt, plot_freq,stats_freq);

 if (!noComm) {

#ifdef _MPI_
  // Gather info
  MPI_Barrier(MPI_COMM_WORLD); // wait for all nodes to finish solving
#endif

#ifdef _MPI_
  if(myrank == 0) {
    int j;
    for(j=1; j<=myM; j++) {
      for(i=1; i<=myN; i++) {
        res[j][i] = E_prev[j][i];
      }
    }

    _DOUBLE_ *buf = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * (myM * myN));

    int k;
    MPI_Request request;
    for(k=1; k<nproc; k++) {
      int pri_y = k / px;
      int pri_x = k % px;
      int pri_myM = ((m+1+py-1) / py) + ((pri_y < (m+1) % py) ? 1 : 0) - 1;
      int pri_myN = ((n+1+px-1) / px) + ((pri_x < (n+1) % px) ? 1 : 0) - 1;
      int pri_starty = (pri_y < (m + 1) % py) ? (pri_y * ((m + 1 + py - 1) / py) + 1) : (pri_y * ((m + 1) / py) + (m + 1) % py + 1);
      int pri_startx = (pri_x < (n + 1) % px) ? (pri_x * ((n + 1 + px - 1) / px) + 1) : (pri_x * ((n + 1) / px) + (n + 1) % px + 1);
      //cout << "Rank:" << k << " startX:" << pri_startx << " startY:" << pri_starty << " M:" << pri_myM << " N:" << pri_myN << endl;

      MPI_Irecv(buf, pri_myM*pri_myN, MPI_DOUBLE, k,
                k, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, MPI_STATUS_IGNORE);

      for(j=1; j<=pri_myM; j++) {
        for(i=1; i<=pri_myN; i++) {
          res[pri_starty+j-1][pri_startx+i-1] = buf[(j-1)*pri_myN+i-1];
        }
      }
    }
  }
  else {
    _DOUBLE_ *buf = (_DOUBLE_*) malloc(sizeof(_DOUBLE_) * (myM  * myN));
    for(int j=1; j<=myM; j++) {
      for(i=1; i<=myN; i++) {
        buf[(j-1)*myN+i-1] = E_prev[j][i];
      }
    }
    MPI_Send(buf, myN*myM, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
  }
#endif
 }

#ifdef _MPI_
 t0 += MPI_Wtime();
#else
 t0 += getTime();
#endif

if (niter != niters)
   cout << "*** niters should be equal to niters" << endl;
 // Report various information
 // Do not remove this call, it is needed for grading
 ReportEnd(logfile,niter,niter,res,m,n,t0,px, py);

 if (plot_freq){
    printf("\n\nEnter any input to close the program and the plot...");
    int resp;
    scanf("%d",&resp);
  }

 logfile.close();
 free (E);
 free (E_prev);
 free (R);
 if (!myrank)
   free (res);
#ifdef _MPI_
 MPI_Finalize();
#endif

}
