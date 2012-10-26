#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "types.h"

double getTime();
void printMatrix(Grid2D A, int n, char *msg);
void verify(Grid2D A, int n, char *Name, char *msg, _DOUBLE_ epsilon, int nreps);

void mm_unblocked(Grid2D A, Grid2D B, Grid2D C, int n, int nreps, double& time_un, _DOUBLE_ epsilon)
{
  int r, i, j, k;
  _DOUBLE_ sum;
  time_un = -getTime(); // Start the timer
  for (r=0; r<nreps; r++){
    for (i=0; i<n; i++)
      for (j=0; j<n; j++){
	sum = 0;
#pragma vector always
#pragma ivdep
	for (k=0; k<n; k++)
	  sum+= A[i][k] * B[k][j];                 
	C[i][j] += sum;
      }
  } 
  time_un += getTime(); // Stop the timer
  // Verify the results
  char c_msg[] = "C";
  char cu_msg[] = "Unblocked multiplication";
  verify(C,n,c_msg,cu_msg,epsilon,nreps);
#ifdef DEBUG
  char cmsg[]= "--C--------";
  printMatrix(C,n,cmsg);
#endif  
}
