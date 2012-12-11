#include "fft_serial.h"

#define PI 3.1415926535897932384626433

void fft3D_serial(complex *** data, int size, int inverse);
void fft1D_serial(complex data[], int length, int isign);
void transposexz_serial(complex *** data, int size);
void transposeyz_serial(complex *** data, int size);
void swap_serial(_DOUBLE_ *a, _DOUBLE_ *b);

void fft3D_serial(complex *** data, int size, int inverse) {
  int x, y;
  for(x = 0; x < size; x++) {
    for(y = 0; y < size; y++) {
      fft1D_serial(data[x][y], size, inverse);
    }
  }
  transposexz_serial(data, size);
  for(x = 0; x < size; x++) {
    for(y = 0; y < size; y++) {
      fft1D_serial(data[x][y], size, inverse);
    }
  }
  transposeyz_serial(data, size);
  for(x = 0; x < size; x++) {
    for(y = 0; y < size; y++) {
      fft1D_serial(data[x][y], size, inverse);
    }
  }
  transposeyz_serial(data, size);
  transposexz_serial(data, size);
}

void fft1D_serial(complex data[], int n, int inverse) {
  int i, j, k, n1, n2, m;
  _DOUBLE_ c, s, exponent, angle, tempImag, tempReal;    
  
  m = log((_DOUBLE_)n) / log(2.0);
  
  j = 0;
  n2 = n/2;
  for (i=1; i < n - 1; i++) {
    n1 = n2;
    while (j >= n1) {
      j = j - n1;
      n1 = n1/2;
    }
    j = j + n1;
               
    if (i < j) {
      swap_serial(&data[i].real, &data[j].real);
      swap_serial(&data[i].comp, &data[j].comp);
    }
  }
                            
  n1 = 0;
  n2 = 1;
                                             
  for (i=0; i < m; i++) {
    n1 = n2;
    n2 = n2 + n2;
    exponent = -2*PI/n2;
    if (inverse) exponent = -exponent;
    angle = 0.0;
                                             
    for (j=0; j < n1; j++) {
      c = cos(angle);
      s = sin(angle);
      angle = angle + exponent;
                                            
      for (k=j; k < n; k=k+n2) {
        tempReal = c*data[k+n1].real - s*data[k+n1].comp;
        tempImag = s*data[k+n1].real + c*data[k+n1].comp;
        data[k+n1].real = data[k].real - tempReal;
        data[k+n1].comp = data[k].comp - tempImag;
        data[k].real = data[k].real + tempReal;
        data[k].comp = data[k].comp + tempImag;
       }
     }
  }

  if (inverse) {
    _DOUBLE_ arg = 1.0 / (_DOUBLE_) n; //sqrt((_DOUBLE_) length);
    for(i = 0; i < n; i++) {
      data[i].real *= arg;
      data[i].comp *= arg;
    }
  }
}

void transposexz_serial(complex *** data, int size) {
  int x, y, z;
 
  for (x = 0; x < size; x++) {
    for (y = 0; y < size; y++) {
      for (z = x+1; z < size; z++) {
        swap_serial(&data[x][y][z].real, &data[z][y][x].real);
        swap_serial(&data[x][y][z].comp, &data[z][y][x].comp);
      }
    }
  }
}

void transposeyz_serial(complex *** data, int size) {
  int x, y, z;

  for (x = 0; x < size; x++) {
    for (y = 0; y < size; y++) {
      for (z = y+1; z < size; z++) {
        swap_serial(&data[x][y][z].real, &data[x][z][y].real);
        swap_serial(&data[x][y][z].comp, &data[x][z][y].comp);
      }
    }
  }
}

void swap_serial(_DOUBLE_ *a, _DOUBLE_ *b) {
  _DOUBLE_ temp;
  temp = *a;
  *a = *b;
  *b = temp;
}

