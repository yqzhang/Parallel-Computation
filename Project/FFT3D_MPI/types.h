#ifndef types_h
#define types_h
#ifndef _DOUBLE
#define _DOUBLE_ float
#else
#define _DOUBLE_ double
#endif
#else
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct complexNumber {
  _DOUBLE_ real;
  _DOUBLE_ comp;
} complex;
