#define BLK 16
#define STEP 4
#define VECTOR 16

void setGrid(int n, dim3 &blockDim, dim3 &gridDim)
{
#ifndef USE_CACHE
   // set your block dimensions and grid dimensions here
   blockDim.x = BLK;
   blockDim.y = STEP;
   gridDim.x = (n + BLK - 1) / BLK;
   gridDim.y = (n + BLK - 1) / BLK;
#else
   blockDim.x = VECTOR * 2;
   blockDim.y = 2;
   gridDim.x = (n + 2 * VECTOR - 1) / (2 * VECTOR);
   gridDim.y = (n + 2 * VECTOR - 1) / (2 * VECTOR);
#endif
}

