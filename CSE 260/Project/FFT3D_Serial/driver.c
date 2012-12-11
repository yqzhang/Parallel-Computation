#include "driver.h"

double getTime();

double mean(complex*** data, int size) {
  int i, j, k;
  double mean = 0.0, temp1 = 0.0, temp2 = 0.0;

  for (i = 0; i < size; i++) {
    temp1 = 0.0;
    for (j = 0; j < size; j++) {
      temp2 = 0.0;
      for (k = 0; k < size; k++) {
        temp2 += data[i][j][k].real + data[i][j][k].comp / (2.0 * size);
      }
      temp1 += temp2 / size;
    }
    mean += temp1 / size;
  }
  
  return mean;
}

int isPowOfTwo(int num) {
  return (num != 0) && ((num & (num-1)) == 0);
}

void cmdLine(int argc, char * argv[], int* n) {
  static struct option long_options[] = {
    {"n", required_argument, 0, 'n'},
  };

  int ac;
  for(ac = 1; ac < argc; ac++) {
    int c;
    while((c=getopt_long(argc,argv,"n:",long_options,NULL)) != -1) {
      switch(c) {
        case 'n':
          *n = atoi(optarg);
          break;
        default:
          printf("Usage: driver [-n <problem size>]\n");
          exit(-1);
      }
    }
  }
}

complex *** alloc3D(int sizex, int sizey, int sizez) {
  complex *** E;
  complex ** row;
  complex * col;
  int i, j;

  E = (complex ***)malloc(sizeof(complex **) * sizex + sizeof(complex *) * sizex * sizey + sizeof(complex) * sizex *sizey * sizez);

  for(i = 0; i < sizex; i++) {
    row = (complex **) (E + sizex);
    row += i * sizey;
    E[i] = row;
  }

  for(i = 0; i < sizex; i++) {
    for(j = 0; j < sizey; j++) {
      row = (complex **) (E + sizex);
      row += sizex * sizey;
      col = ((complex *) row) + i * sizey * sizez + j * sizez;
      E[i][j] = col;
    }
  }

  return E;
}

void init(complex *** E, int sizex, int sizey, int sizez) {
  int i, j, k;
  for(i = 0; i < sizex; i++) {
    for(j = 0; j < sizey; j++) {
      for(k = 0; k < sizez; k++) {
        E[i][j][k].real = ((_DOUBLE_)rand()) / RAND_MAX;
        E[i][j][k].comp = ((_DOUBLE_)rand()) / RAND_MAX;
      }
    }
  }
}

int main(int argc, char** argv) {
  complex ***data, ***original;

  int size = 128, i;
  double origMean, newMean;

  cmdLine(argc, argv, &size);

  if (!isPowOfTwo(size)) {
    printf("n must be a power of two!\n");
    return -1;
  }

  data = alloc3D(size, size, size);
  init(data, size, size, size);

  // Start the timer
  double t0 = -getTime();
  fft3D(data, size, 0);

  // Stop the timer
  t0 += getTime();

  double ops = 15 * size * size * size * log(size) / log(2);
  double gflops = 1e-9 * ops / t0;

  printf("Problem size: %d\n", size);
  printf("Running Time: %.6lf sec, GFLOPS: %.4lf\n", t0, gflops);
  printf("\n");

  free(data);

  return 0;
}
