#ifndef apf_h
#define apf_h
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
#include "types.h"
 // Various constants - these definitions shouldn't change
 const _DOUBLE_ a=0.1, b=0.1, kk=8.0, M1= 0.07, M2=0.3, epsilon=0.01, d=5e-5;
 // Wait for response after plotting
#ifdef DEBUG
 const int WAIT = 1;
#else
 const int WAIT = 0;
#endif

#else
#endif

// External functions
extern "C" {
    void splot(_DOUBLE_ **E, int niter, int m, int n, int WAIT);
    void splot1d(_DOUBLE_ *U, int niter, int m, int n, int WAIT);

    void printMat(_DOUBLE_ **U, int m, int n);
    void ABEND();
    void Stop();
}
