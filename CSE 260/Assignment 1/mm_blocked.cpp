#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "types.h"

using namespace std;

double getTime();
extern void verify(_DOUBLE_ **A, int n, char *Name, char *msg, _DOUBLE_ epsilon, int nreps);
static inline void block_load(Grid2D from, Grid2D to, int x, int y, int BSZ) {
  int i, j, tempX, tempY;

  for (i = 0; i < BSZ; i++) {
    for (j = 0; j < BSZ; j++) {
      tempX = x + i;
      tempY = y + j;
      to[i][j] = from[tempX][tempY];
    }
  }
}

static inline void block_store(Grid2D from, Grid2D to, int x, int y, int BSZ) {
  int i, j, tempX, tempY;

  for (i = 0; i < BSZ; i++) {
    for (j = 0; j < BSZ; j++) {
      tempX = x + i;
      tempY = y + j;
      to[tempX][tempY] = from[i][j];
    }
  }
}

static inline void block_transpose(Grid2D from, Grid2D to, int x, int y, int BSZ) {
  int i, j, tempX, tempY;

  for (i = 0; i < BSZ; i++) {
    for (j = 0; j < BSZ; j++) {
      tempX = x + i;
      tempY = y + j;
      to[j][i] = from[tempX][tempY];
    }
  }
}

static inline void block_multiply(Grid2D Cij, Grid2D Aik, Grid2D Bkj, int BSZ) {
  int i, j, k;
  register _DOUBLE_ RESTRICT sum1, sum2, sum3, sum4;
  register _DOUBLE_ RESTRICT a0, a1, a2, a3;

  for (i = 0; i < BSZ; i++) {
    a0 = Aik[i][0];
    a1 = Aik[i][1];
    a2 = Aik[i][2];
    a3 = Aik[i][3];
    for (j = 0; j < BSZ; j += 4) {
      sum1 = a0 * Bkj[j][0] +
             a1 * Bkj[j][1] +
             a2 * Bkj[j][2] +
             a3 * Bkj[j][3];
      sum2 = a0 * Bkj[j + 1][0] +
             a1 * Bkj[j + 1][1] +
             a2 * Bkj[j + 1][2] +
             a3 * Bkj[j + 1][3];
      sum3 = a0 * Bkj[j + 2][0] +
             a1 * Bkj[j + 2][1] +
             a2 * Bkj[j + 2][2] +
             a3 * Bkj[j + 2][3];
      sum4 = a0 * Bkj[j + 3][0] +
             a1 * Bkj[j + 3][1] +
             a2 * Bkj[j + 3][2] +
             a3 * Bkj[j + 3][3];

      for (k = 4; k < BSZ; k += 4) {
#pragma vector always
#pragma ivdep
        sum1 += Aik[i][k] * Bkj[j][k] +
                Aik[i][k + 1] * Bkj[j][k + 1] +
                Aik[i][k + 2] * Bkj[j][k + 2] +
                Aik[i][k + 3] * Bkj[j][k + 3];
        sum2 += Aik[i][k] * Bkj[j + 1][k] +
                Aik[i][k + 1] * Bkj[j + 1][k + 1] +
                Aik[i][k + 2] * Bkj[j + 1][k + 2] +
                Aik[i][k + 3] * Bkj[j + 1][k + 3];
        sum3 += Aik[i][k] * Bkj[j + 2][k] +
                Aik[i][k + 1] * Bkj[j + 2][k + 1] +
                Aik[i][k + 2] * Bkj[j + 2][k + 2] +
                Aik[i][k + 3] * Bkj[j + 2][k + 3];
        sum4 += Aik[i][k] * Bkj[j + 3][k] +
                Aik[i][k + 1] * Bkj[j + 3][k + 1] +
                Aik[i][k + 2] * Bkj[j + 3][k + 2] +
                Aik[i][k + 3] * Bkj[j + 3][k + 3];
      }
      Cij[i][j] += sum1;
      Cij[i][j + 1] += sum2;
      Cij[i][j + 2] += sum3;
      Cij[i][j + 3] += sum4;
    }
  }
}

void mm_blocked(Grid2D A, Grid2D B, Grid2D C, Grid2D Aik, Grid2D Bkj, Grid2D Cij, int n, int BSZ, int nreps, double& time_blk, _DOUBLE_ epsilon)
{
  int i, j, r, ii, jj;
    // Do the blocked computation
    // First, clear out C, so we can use it again
    // Don't remove this code or matrix multiply will fail to
    // verify for niters > 1
    for (i=0; i<n; i++)
#pragma vector always
#pragma ivdep
      for (j=0; j<n; j++)
        C[i][j] = 0.0;

    time_blk = -getTime(); // Start the timer
    // Your code Goes here
  int k, kk;

  for (r = 0; r < nreps; r++) {
    for (ii = 0; ii < n; ii += BSZ) {
      for (jj = 0; jj < n; jj += BSZ) {
        block_load(C, Cij, ii, jj, BSZ);
        for (kk = 0; kk < n; kk += BSZ) {
          block_load(A, Aik, ii, kk, BSZ);
          block_transpose(B, Bkj, jj, kk, BSZ);
          block_multiply(Cij, Aik, Bkj, BSZ);
        }
        block_store(Cij, C, ii, jj, BSZ);
      }
    }
  }
    // Finish my code
    time_blk += getTime(); // Stop the Timer

    // Verify the results
    char cb_msg[] = "Blocked multiplication";
    char c_msg[] = "C";
    verify(C,n,c_msg,cb_msg,epsilon,nreps);
#ifdef DEBUG
    printMatrix(C,n,cmsg);
#endif
}
