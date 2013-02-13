/* **********************************************************
 *  Author : Urvashi R.V. [04/06/2004]
 * 	Modified by Scott Baden [10/8/06]
 * 	Modified by Pietro Cicoctti [10/8/08]
 * 
 *************************************************************/

/*
 * 
 * Do not change the code in this file, as doing so
 * could cause your submission to be graded incorrectly
 */

#include <stdio.h>

/* Function to plot the 2D array
 * 'gnuplot' is instantiated via a pipe and 
 * the values to be plotted are passed through, along 
 * with gnuplot commands */

#include "types.h"

FILE *gnu=NULL;

void splot(_DOUBLE_ **U, int niter, int m, int n, int WAIT)
{
    int i, j;
    _DOUBLE_ mx, mn;
    int du;

    if(gnu==NULL) gnu = popen("gnuplot","w");
//    sleep(10);
//    printf("hit return when plot is ready");
//    scanf("%d",&du);
    mx= -1; mn = 32768;
    for (j=0; j<m; j++)
       for (i=0; i<n; i++){
       if (U[j][i] > mx)
           mx = U[j][i];
       if (U[j][i] < mn)
           mn = U[j][i];
       }
    fprintf(gnu,"\n");
    fprintf(gnu,"\n");
    fprintf(gnu,"set title \"niter = %d\"\n", niter);
    fprintf(gnu,"set size square\n");
    fprintf(gnu,"set key off\n");
    fprintf(gnu,"set pm3d map\n");
    fprintf(gnu,"set palette defined (-3 \"blue\", 0 \"white\", 1 \"red\")\n");
    /* Various color schemes
     * fprintf(gnu,"set palette rgbformulae 22, 13, 31\n");
     * fprintf(gnu,"set palette rgbformulae 30, 31, 32\n");
    */

    fprintf(gnu,"splot [0:%d] [0:%d][%f:%f] \"-\"\n",m-1,n-1,mn,mx);
    for (j=0; j<m; j++){
       for (i=0; i<n; i++) {
            fprintf(gnu,"%d %d %f\n", i, j, U[i][j]);
       }
       fprintf(gnu,"\n");
    }
    fprintf(gnu,"e\n");
    fflush(gnu);
    return;
}

#define U(i,j) U[i+j*(n+2)]

void splot1d(_DOUBLE_ *U, int niter, int m, int n, int WAIT)
{
    int i, j;
    _DOUBLE_ mx, mn;

    if(gnu==NULL) gnu = popen("gnuplot","w");
    
    mx= -1; mn = 32768;
    for (j=0; j<m; j++)
       for (i=0; i<n; i++){
       if (U(j,i) > mx)
           mx = U(j,i);
       if (U(j,i) < mn)
           mn = U(j,i);
       }
    fprintf(gnu,"set title \"niter = %d\"\n", niter);
    fprintf(gnu,"set size square\n");
    fprintf(gnu,"set key off\n");
    fprintf(gnu,"set pm3d map\n");
    fprintf(gnu,"set palette defined (-3 \"blue\", 0 \"white\", 1 \"red\")\n");
    /* Various color schemes
     * fprintf(gnu,"set palette rgbformulae 22, 13, 31\n");
     * fprintf(gnu,"set palette rgbformulae 30, 31, 32\n");
    */

    fprintf(gnu,"splot [0:%d] [0:%d][%f:%f] \"-\"\n",m-1,n-1,mn,mx);
    for (j=0; j<m; j++){
       for (i=0; i<n; i++) {
            fprintf(gnu,"%d %d %f\n", i, j, U(i,j));
       }
       fprintf(gnu,"\n");
    }
    fprintf(gnu,"e\n");
    fflush(gnu);
    return;
}
