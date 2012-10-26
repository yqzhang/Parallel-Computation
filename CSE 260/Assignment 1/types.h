#ifndef types_h
/* Do not change the code in this file, as doing so
 * could cause your submission to be graded incorrectly
 */

/*
 * Include this in every module that uses floating point
 * arithmetic, and declare all floating point values as "_DOUBLE_"
 * With a switch of a command line macro set up in the Makefile
 * we can then change the arithmetic
*/
#define types_h
#ifndef _DOUBLE
#define _DOUBLE_ float
#else
#define _DOUBLE_ double
#endif
#else
#endif

#ifdef __RESTRICT
// Intel uses a different syntax than the Gnu compiler
// #define RESTRICT restrict
// Gnu uses a different syntax than the Intel compiler
#define RESTRICT __restrict
#else
// Turns the restrict keyword into the empty string if not selected
#define RESTRICT
#endif
typedef _DOUBLE_ *RESTRICT *RESTRICT Grid2D;
typedef _DOUBLE_ *RESTRICT *RESTRICT *RESTRICT Grid3D;


#ifdef __INTEL_COMPILER
#include "mkl_cblas.h"
#else

extern "C" 
{
#include "cblas.h"
}
#endif

#include <stdint.h>
