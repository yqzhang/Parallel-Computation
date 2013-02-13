// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "utils.h"
#include "types.h"
using namespace std;

#define BLK 16
#define STEP 4
#define VECTOR 16 

__global__ void matMul(int N, _DOUBLE_ *C, _DOUBLE_ *A, _DOUBLE_ *B) {
#ifndef USE_CACHE
    const unsigned int bx = blockIdx.x;
    const unsigned int by = blockIdx.y; 
    const unsigned int tx = threadIdx.x; 
    const unsigned int ty = threadIdx.y;

    const unsigned int aBegin = N * (by * BLK); //A(0,by)
    const unsigned int aEnd = aBegin + N;
    const unsigned int aStep = BLK;             //offsetA
    
    const unsigned int bBegin = BLK * bx;       //B(bx,0)
    const unsigned int bStep = BLK * N;         //offsetB
    
    _DOUBLE_ cSub[BLK] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    unsigned int a, b;

    for (a = aBegin,b = bBegin; a < aEnd - BLK; a += aStep,b += bStep) {
        __shared__ _DOUBLE_ As[BLK][BLK];
        __shared__ _DOUBLE_ Bs[BLK][BLK];
#pragma unroll
        for(int i = 0; i < BLK / STEP; i++) {
            As[ty + i * STEP][tx] = A[a + N * (ty + i * STEP) + tx];
            Bs[ty + i * STEP][tx] = B[b + N * (ty + i * STEP) + tx];
        }
        __syncthreads();
#pragma unroll
        for (int i = 0; i < BLK / STEP; i++) {
            for(int k = 0; k < BLK; k++) {
                cSub[i] += As[ty + i * STEP][k] * Bs[k][tx];
            }
        }

        __syncthreads();
    }

    // Maybe out of the matrix
    {
        __shared__ _DOUBLE_ As[BLK][BLK];
        __shared__ _DOUBLE_ Bs[BLK][BLK];
#pragma unroll
        for(int i = 0; i < BLK / STEP; i++) {
            if((ty + i * STEP + by * BLK < N) && ((a - aBegin) + tx < N)) As[ty + i * STEP][tx] = A[a + N * (ty + i * STEP) + tx];
            else As[ty + i * STEP][tx] = 0;
            if((ty + i * STEP + (a - aBegin) < N) && (BLK * bx + tx < N)) Bs[ty + i * STEP][tx] = B[b + N * (ty + i * STEP) + tx];
            else Bs[ty + i * STEP][tx] = 0;
        }
        __syncthreads();
#pragma unroll
        for (int i = 0; i < BLK / STEP; i++) {
            for(int k = 0; k < BLK; k++) {
                cSub[i] += As[ty + i * STEP][k] * Bs[k][tx];
            }
        }

        __syncthreads();
    }

    for(int i = 0; i < BLK / STEP; i++) {
        if(by * BLK + ty  + i * STEP < N && bx * BLK + tx < N) {
            int cIndex = (by * BLK + ty + i * STEP) * N + (bx * BLK + tx);
            C[cIndex] = cSub[i];
        }
    }
#else
    const int I =  blockIdx.x * blockDim.x + threadIdx.x;
    const int J =  (blockIdx.y * blockDim.y + threadIdx.y) * VECTOR;
    _DOUBLE_ _c[VECTOR] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    if(I < N) {
        if(N - J >= VECTOR) {
#pragma unroll
            for (int k = 0; k < N; k++) {
                _DOUBLE_ a = A[I * N + k];
                for(int i = 0; i < VECTOR; i++) {
                    _DOUBLE_ b = B[k * N + J + i];
                    _c[i] += a * b;
                }
            }
#pragma unroll
            for(int i = 0; i < VECTOR; i++) {
                C[I * N + J + i] = _c[i];
            }
        }
        else {
            int upBound =  N - J;
#pragma unroll
            for (int k = 0; k < N; k++) {
                _DOUBLE_ a = A[I * N + k];
                for(int i = 0; i < upBound; i++) {
                    _DOUBLE_ b = B[k * N + J + i];
                    _c[i] += a * b;
                }
            }
#pragma unroll
            for(int i = 0; i < upBound; i++) {
                C[I * N + J + i] = _c[i];
            }
        }
    }
#endif
}
