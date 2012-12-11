#include <string.h>
#include "fft_mpi.h"

#define PI 3.1415926535897932384626433

void fft3D_MPI(complex *** data, int px, int py, int x, int y,
           int sizex, int sizey, int sizez, int noComm, _DOUBLE_ *buf, _DOUBLE_ *recv, int inverse);
void fft1D_MPI(complex data[], int length, int isign);
void transposexz_MPI(complex *** data, int px, int py, int x, int y, int sizex, int sizey, int sizez, _DOUBLE_ *buf, _DOUBLE_ *recv);
void transposeyz_MPI(complex *** data, int px, int py, int x, int y, int sizex, int sizey, int sizez, _DOUBLE_ *buf, _DOUBLE_ *recv);
void swap_MPI(_DOUBLE_ * a, _DOUBLE_ * b);
void gatherInfo_MPI(complex *** global, complex *** data,
                int px, int py, int x, int y,
                int sizex, int sizey, int sizez, _DOUBLE_ *recv);
void gatherBasic_MPI(complex *** global, complex *** data,
                int px, int py, int x, int y,
                int sizex, int sizey, int sizez, _DOUBLE_ *recv);

void gatherInfo_MPI(complex *** global, complex *** data,
                int px, int py, int x, int y,
                int sizex, int sizey, int sizez, _DOUBLE_ *recv) {
  int i, j, k, index, procx, procy;
  _DOUBLE_ * buf;

  int bufSize = 2 * sizex * sizey * sizez;

  buf = (_DOUBLE_*) (((complex **) (data + sizex)) + sizex * sizey);

  MPI_Gather(buf, bufSize, MPI_DOUBLE,
             recv, bufSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(x == 0 && y == 0) {
    index = 0;
    for(procy = 0; procy < py; procy++) {
      for(procx = 0; procx < px; procx++) {
        for(i = 0; i < sizex; i++) {
          for(j = 0; j < sizey; j++) {
            for(k = 0; k < sizez; k++) {
              global[procy * sizey + j][k][procx * sizex + i].real = recv[index++];
              global[procy * sizey + j][k][procx * sizex + i].comp = recv[index++];
            }
          }
        }
      }
    }    
  }
}

void gatherBasic_MPI(complex *** global, complex *** data,
                int px, int py, int x, int y,
                int sizex, int sizey, int sizez, _DOUBLE_ *recv) {
  int i, j, k, index, procx, procy;
  _DOUBLE_ * buf;
  complex ** tmp1;

  int bufSize = 2 * sizex * sizey * sizez;

  tmp1 = (complex **) (data + sizex);
  tmp1 += sizex * sizey;

  buf = (_DOUBLE_*) tmp1;

  MPI_Gather(buf, bufSize, MPI_DOUBLE,
             recv, bufSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(x == 0 && y == 0) {
    index = 0;
    for(procy = 0; procy < py; procy++) {
      for(procx = 0; procx < px; procx++) {
        for(i = 0; i < sizex; i++) {
          memcpy((_DOUBLE_*)(global[procx * sizex + i][procy * sizey]), recv, 2 * sizex * sizez * sizeof(_DOUBLE_));
          recv += 2 * sizey * sizez;
        }
      }
    }
  } 
}

void fft3D_MPI(complex *** data, int px, int py, int x, int y,
           int sizex, int sizey, int sizez, int noComm,
           _DOUBLE_ *buf, _DOUBLE_ *recv, int inverse) {
  int i, j;
  for(i = 0; i < sizex; i++) {
    for(j = 0; j < sizey; j++) {
      fft1D_MPI(data[i][j], sizez, inverse);
    }
  }

  if(!noComm) transposexz_MPI(data, px, py, x, y, sizex, sizey, sizez, buf, recv);
  //transposexz(data, px, py, x, y, sizex, sizey, sizez, buf, recv);

  for(i = 0; i < sizex; i++) {
    for(j = 0; j < sizey; j++) {
      fft1D_MPI(data[i][j], sizez, inverse);
    }
  }

  if(!noComm) transposeyz_MPI(data, px, py, x, y, sizex, sizey, sizez, buf, recv);
  //transposeyz(data, px, py, x, y, sizex, sizey, sizez, buf, recv);

  for(i = 0; i < sizex; i++) {
    for(j = 0; j < sizey; j++) {
      fft1D_MPI(data[i][j], sizez, inverse);
    }
  }

  if(!noComm) transposeyz_MPI(data, px, py, x, y, sizex, sizey, sizez, buf, recv);
  if(!noComm) transposexz_MPI(data, px, py, x, y, sizex, sizey, sizez, buf, recv);
}

void fft1D_MPI(complex data[], int n, int inverse) {
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
      swap_MPI(&data[i].real, &data[j].real);
      swap_MPI(&data[i].comp, &data[j].comp);
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

void transposexz_MPI(complex *** data, int px, int py, int x, int y,
               int sizex, int sizey, int sizez, _DOUBLE_ *buf, _DOUBLE_ *recv) {
 // Second step: init the buffer 
  int i, j, k, l, index;
  MPI_Request * srequest = (MPI_Request *)malloc(sizeof(MPI_Request) * (px - 1));
  MPI_Request * rrequest = (MPI_Request *)malloc(sizeof(MPI_Request) * (px - 1));
  int step = sizex * sizey * sizex * 2;

  index = 0;
  
  for(k = 0; k < px; k++) {
    if(k == x) continue;
    for(l = 0; l < sizex; l++) {
      for(j = 0; j < sizey; j++) {
        for(i = 0; i < sizex; i++) {
          buf[index++] = data[i][j][k * sizex + l].real;
          buf[index++] = data[i][j][k * sizex + l].comp;
        }
      }
    }
  }

  // Third step: MPI_Alltoallv
  for(k = 0, l = 0; k < px; k++) {
    if(k == x) continue;
    //printf("%d receive from: %d\n", y*px+x, y*px+k); 
    MPI_Irecv(&recv[l * step], step, MPI_DOUBLE, y * px + k, y * px + k, MPI_COMM_WORLD, &rrequest[l]);
    MPI_Isend(&buf[l * step], step, MPI_DOUBLE, y * px + k, y * px + x, MPI_COMM_WORLD, &srequest[l]);
    l++;
  }
  MPI_Waitall(px - 1, rrequest, MPI_STATUSES_IGNORE);

  for (i = 0; i < sizex; i++) {
    for (j = 0; j < sizey; j++) {
      for (l = i+1; l < sizex; l++) {
        swap_MPI(&data[i][j][x * sizex + l].real, &data[l][j][x * sizex + i].real);
        swap_MPI(&data[i][j][x * sizex + l].comp, &data[l][j][x * sizex + i].comp);
      }
    }
  }

  index = 0;
  for(k = 0; k < px; k++) {
    if(k == x) continue;
    for(i = 0; i < sizex; i++) {
      for(j = 0; j < sizey; j++) {
        memcpy((_DOUBLE_*)(&data[i][j][k*sizex]), recv, 2 * sizex * sizeof(_DOUBLE_));
        recv += 2 * sizex;   
      }
    }
  }

  MPI_Waitall(px - 1, srequest, MPI_STATUSES_IGNORE);
  free(rrequest);
  free(srequest);
}

void transposeyz_MPI(complex *** data, int px, int py, int x, int y,
               int sizex, int sizey, int sizez, _DOUBLE_ *buf, _DOUBLE_ *recv) {
  // Second step: init the buffer

  int i, j ,k, l, index;
  MPI_Request * srequest = (MPI_Request *)malloc(sizeof(MPI_Request) * (py - 1));
  MPI_Request * rrequest = (MPI_Request *)malloc(sizeof(MPI_Request) * (py - 1));
  int step = sizex * sizey * sizey * 2;

  index = 0;
  for(k = 0; k < py; k++) {
    if(k == y) continue;
    for(i = 0; i < sizex; i++) {
      for(l = 0; l < sizey; l++) {
        for(j = 0; j < sizey; j++) {
          buf[index++] = data[i][j][k * sizey + l].real;
          buf[index++] = data[i][j][k * sizey + l].comp;
        }
      }
    }
  }

  // Third step: send and recv
  for(k = 0, l = 0; k < py; k++) {
    if(k == y) continue;
    //printf("%d receive from: %d\n", y*px+x, k*px+x);
    MPI_Irecv(&recv[l * step], step, MPI_DOUBLE, k * px + x, k * px + x, MPI_COMM_WORLD, &rrequest[l]);
    MPI_Isend(&buf[l * step], step, MPI_DOUBLE, k * px + x, y * px + x, MPI_COMM_WORLD, &srequest[l]);
    l++;
  }
  MPI_Waitall(py - 1, rrequest, MPI_STATUSES_IGNORE);

  for (i = 0; i < sizex; i++) {
    for (j = 0; j < sizey; j++) {
      for(l = j + 1; l < sizey; l++) {
        swap_MPI(&data[i][j][y * sizey + l].real, &data[i][l][y * sizey + j].real);
        swap_MPI(&data[i][j][y * sizey + l].comp, &data[i][l][y * sizey + j].comp);
      }
    }
  }

  index = 0;
  for(k = 0; k < py; k++) {
    if(k == y) continue;
    for(i = 0; i < sizex; i++) {
      for(j = 0; j < sizey; j++) {
        memcpy((_DOUBLE_*)(&data[i][j][k*sizey]), recv, 2 * sizey * sizeof(_DOUBLE_));
        recv += 2 * sizey;
      }
    }
  }

  MPI_Waitall(py - 1, srequest, MPI_STATUSES_IGNORE);
  free(rrequest);
  free(srequest);
}

void swap_MPI(_DOUBLE_ * a, _DOUBLE_ * b) {
  _DOUBLE_ temp;
  temp = * a;
  * a = * b;
  * b = temp;
}
