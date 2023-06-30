#include <cassert>
#include <cuda_runtime.h>
#include "matmul_device.cuh"

#define BLOCK_WIDTH 32
/*
 * Read TODO items below
 */




__global__
void naiveMatmul(float *a, float *b, float *c, int n) {
    int j = blockIdx.x * blockDim.x + threadIdx.x; //column
    int i = blockIdx.y * blockDim.y + threadIdx.y; //row

    float acc = 0;
    for (int k=0; k<n; k++) {
	acc += a[i*n+k] * b[k*n+j]; //each thread access an element of a and b
    //threads access a in the col order and b in the row order at the same time.
    }
    c[i*n+j] = acc;
}


        

__global__
void cacheMatmul(float *a, float *b, float *c, int n) {
    // TODO: replace this function with cache friendly version
    //Devide the three matrix into sub blocks of size BLOCK_WIDTH * BLOCK_WIDTH;
    //The inner loop calculates two sub blocks of a and b of size BLOCK_WIDTH * BLOCK_WIDTH;
    //In the outer loop, each sub block C needs one row of blocks of a and one column of blocks
    //of b.
    
    
    float acc = 0;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int col = blockIdx.x*blockDim.x + threadIdx.x;
 

    for (int i = 0; i < (n/BLOCK_WIDTH); i++) {
        for (int j = 0; j < BLOCK_WIDTH; ++j)
            acc += a[row*n + i*BLOCK_WIDTH + j] * b[(i*BLOCK_WIDTH + j)*n + col];
    }

    c[(blockIdx.y * blockDim.y + threadIdx.y)*n+(blockIdx.x*blockDim.x)+threadIdx.x]= acc;
    
    
    //Threads access a in the row order and b to get c in the row order synchronisedly.
    //The entire a should be access n times in order to get the entire c
    //Each cloumns of thread do multiplication with one a[i, row] * b[row,*] = c[i,*];
    //For all threads within each iteration, it scans onw row of a and entire b to get one row.
    //Although in this way it synchronises access by row, the b is access n times that causes
    //the entire calculation every slow.
    /*
    int row = blockIdx.y * blockDim.y + threadIdx.y; //i
    int col = blockIdx.x * blockDim.x + threadIdx.x; //k
    
    
    for(int i = 0; i < n; i++) {
        float r = b[row * n + col];
        atomicAdd(&c[i * n + col], r * a[i *n + row]);
    }
    __syncthreads();
    */
  
}

__global__
void sharedMatmul(float *a, float *b, float *c, int n) {
    // TODO: replace this function with optimized code using
    // shared memory
    

    __shared__ float sharedMemA[BLOCK_WIDTH][BLOCK_WIDTH];
    __shared__ float sharedMemB[BLOCK_WIDTH][BLOCK_WIDTH];

    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int col = blockIdx.x * blockDim.x + threadIdx.x; 

    float acc = 0;

    for (int i = 0; i < n/BLOCK_WIDTH; i++) {
        sharedMemA[threadIdx.y][threadIdx.x] = a[row * n + (i * BLOCK_WIDTH + threadIdx.x)];
        sharedMemB[threadIdx.y][threadIdx.x] = b[(i * BLOCK_WIDTH + threadIdx.y) * n + col];


        __syncthreads();


        for (int j = 0; j < BLOCK_WIDTH; j++) {
            acc += sharedMemA[threadIdx.y][j] * sharedMemB[j][threadIdx.x];
        }


        __syncthreads();
    }

    c[row * n + col] = acc;
}

void cudaMatmul(float *a, float *b, float *c, int n, MatmulImplementation type)
{
    // TODO: play with the gridSize and blockSize to find the best one
    if (type == NAIVE) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        naiveMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    else if (type == CACHE) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        cacheMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    else if (type == SHARED) {
        dim3 blockSize(32, 32);
        dim3 gridSize(n / 32, n / 32);
        sharedMatmul<<<gridSize, blockSize>>>(a,b,c,n);
    }
    // Unknown type
    else
        assert(false);
}
