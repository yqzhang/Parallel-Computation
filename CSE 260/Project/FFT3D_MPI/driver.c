#include <time.h>
#include "driver.h"

inline int isPowOfTwo(int num) {
  return (num != 0) && ((num & (num-1)) == 0);
}

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

const double kMicro = 1.0e-6;
double getTime() {
  struct timeval TV;

  const int RC = gettimeofday(&TV, NULL);
  if(RC == -1)
  {
    printf("ERROR: Bad call to gettimeofday\n");
    return(-1);
  }

  return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
}

void cmdLine(int argc, char * argv[], int* n, int* px, int* py, int* noComm, int *verify) {
  static struct option long_options[] = {
    {"n", required_argument, 0, 'n'},
    {"px", required_argument, 0, 'x'},
    {"py", required_argument, 0, 'y'},
    {"nocomm", no_argument, 0, 'k'},
    {"verify", no_argument, 0, 'v'},
  };

  int ac;
  for(ac = 1; ac < argc; ac++) {
    int c;
    while((c=getopt_long(argc,argv,"n:x:y:kp:vp:",long_options,NULL)) != -1) {
      switch(c) {
        case 'n':
          *n = atoi(optarg);
          break;
        case 'x':
          *px = atoi(optarg);
          break;
        case 'y':
          *py = atoi(optarg);
          break;
        case 'k':
          *noComm = 1;
          break;
        case 'v':
          *verify = 1;
          break;
        default:
          printf("Usage: driver [-n <problem size>] [-x <x processor geometry>] [-y <y processor geometry>]\n              [-k <no communication>] [-v <verify>]\n");
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

_DOUBLE_ * alloc1D(int size) {
  return (_DOUBLE_*)malloc(sizeof(_DOUBLE_) * size);
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

void printyz(complex ***data, int sizey, int sizez) {
  int i, j;

  for (i = 0; i < sizey; i++) {
    for (j = 0; j < sizez; j++) {
      printf("%f ", data[0][i][j].real);
    }
    printf("\n");
  }
}

void printxz(complex ***data, int sizex, int sizez) {
  int i, j;

  for (i = 0; i < sizex; i++) {
    for (j = 0; j < sizez; j++) {
      printf("%f ", data[i][0][j].real);
    }
    printf("\n");
  }
}

int main(int argc, char** argv) {
  complex ***global, ***data;
  _DOUBLE_ *buf, *recv, *gath;

  int size = 128;
  int px = 1, py = 1;
  int noComm = 0;
  int verify = 0;

  double serialMean, mpiMean;
  double timeMPI;
  serialMean = mpiMean = 0.0;

  // Init MPI
  MPI_Init(&argc, &argv);

  // command line
  cmdLine(argc, argv, &size, &px, &py, &noComm, &verify);

  // For now, keep these constraints. If there is time, work around them
  // (e.g. pad with zeroes, etc.
  if (!isPowOfTwo(size)) {
    printf("Error: Size of the problem must be a power of 2! Exiting.\n");
    return -1;
  }
  
  if ((size % px != 0) || (size % py != 0)) {
    printf("Error: Processor geometry does not divide the problem size evenly! Exiting.\n");
    return -1;
  }

  int nproc, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // Get the thread ID by x and y 
  int x, y;
  y = myrank / px;
  x = myrank % px;

  int sizex, sizey;
  sizex = size / px;
  sizey = size / py;

  if(!myrank) {
    global = alloc3D(size, size, size);
  }

  data = alloc3D(sizex, sizey, size);
  buf  = alloc1D(2 * sizex * sizey * size);
  recv = alloc1D(2 * sizex * sizey * size);
  gath = alloc1D(2 * size * size * size);
  srand(time(NULL) * (myrank+1));
  init(data, sizex, sizey, size);

  MPI_Barrier(MPI_COMM_WORLD);

  // Do the serial fft if verify = 1
  // Gather into the global matrix
  if(verify) gatherBasic_MPI(global, data, px, py, x, y, sizex, sizey, size, gath);

  MPI_Barrier(MPI_COMM_WORLD);

  if(!myrank) {
    // Do the serial fft3D for verification of the correctness
    if(verify) {
      fft3D_serial(global, size, 0);
      serialMean = mean(global, size);
    }
  }

  // Start the timer
  timeMPI = -MPI_Wtime();

  // Do the MPI fft3D
  fft3D_MPI(data, px, py, x, y, sizex, sizey, size, noComm, buf, recv, 0);
  MPI_Barrier(MPI_COMM_WORLD);

  // Gather info
  if(!noComm) gatherBasic_MPI(global, data, px, py, x, y, sizex, sizey, size, gath);

  // Stop the timer
  timeMPI += MPI_Wtime();

  // Gather info for virification of the correctness
  if(noComm) gatherBasic_MPI(global, data, px, py, x, y, sizex, sizey, size, gath);

  // Get the mean of the MPI version to verify the result
  if(verify) {
    if(!myrank) {
      mpiMean = mean(global, size);
    }
  }

  // Calculate the number of operations and GFlops
  if(!myrank) {
    double ops = 15 * size * size * size * log(size) / log(2);
    double gflops = 1e-9 * ops / timeMPI;
    // Print the result
    printf("Problem size: %d, PX: %d, PY: %d, Communication: ", size, px, py);
    printf("%s\n", (noComm) ? ("OFF") : ("ON"));
    printf("MPI version    - Running Time: %.6lf sec, GFLOPS: %.4lf\n", timeMPI, gflops);

    if(verify) {
      // Verify the result by comparing the mean of the matrix
      printf("Serial Mean: %.6lf, MPI Mean: %.6lf\n", serialMean, mpiMean);
#define EPSILON 0.00001
      if (abs(serialMean - mpiMean) < EPSILON) {
        printf("The calculation seems to be CORRECT!\n");
      } else {
        printf("The calculation is WRONG!\n");
      }
    }

    printf("\n");
  }

  free(data);
  if(!myrank) free(global);

  MPI_Finalize();

  return 0;
}
