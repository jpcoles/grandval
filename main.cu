#include <cuda.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "grandval.h"
#include "nbody.h"
#include "ic.h"

#define WITH_CUDA 1

#include "nbody.cu"
#include "ic.cu"

int cuda_threads_per_block;
int cuda_blocks_x, cuda_blocks_y;

__global__ void cuda_step(struct nbody_potential *phi, struct particle *P, int NP, tyme_t dt)
{
    int i;
    acc_t a[3];

    int pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;

    //printf("blockIdx %i  blockDim %i  threadIdx %i\n", blockIdx.x, blockDim.x, threadIdx.x);
    //printf("particle %i x=%f\n", threadIdx.x, p->x[0]);

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    //phi->accel(phi, p, a);
    accel(phi, p, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

void step(struct potential *phi, struct particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    phi->accel(phi, p, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

void write_positions(struct particle *P, int NP, int first_output)
{
    int i;
    char *mode;
    char fname[256];

    if (first_output)
        mode = "w";
    else
        mode = "a+";

    FILE *fp = stdout;

    for (i=0; i < NP; i++)
    {
        //sprintf(fname, "gv-%05i.out", i);
        //FILE *fp = fopen(fname, mode);
#if WITH_INTEGERS
        fprintf(fp, "%ld %ld %ld\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#else
        fprintf(fp, "%f %f %f\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#endif
        //fclose(fp);
    }
}

void cuda_step_all(struct potential *phi, struct particle *P, int NP, tyme_t dt)
{
    int nblocks  = (int)ceil(NP / (float)cuda_threads_per_block);
    int nthreads = NP < cuda_threads_per_block ? NP : cuda_threads_per_block;
    //fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
    cuda_step<<<nblocks,nthreads>>>((struct nbody_potential *)phi->phi_dev, P, NP, dt);
}

void step_all(struct potential *phi, struct particle *P, int NP, tyme_t dt)
{
    int i;
    #pragma omp parallel for
    for (i=0; i < NP; i++)
        step(phi, P + i, dt);
}

void show_cuda_devices()
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    int device;
    if (deviceCount == 0)
    {
        printf("No CUDA devices found.\n");
    }

    for (device = 0; device < deviceCount; ++device) 
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
        printf("              maxThreadsPerBlock %d.\n", deviceProp.maxThreadsPerBlock);
        printf("              maxThreadsDim %d,%d,%d.\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
        printf("              maxGridSize %d,%d,%d.\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);

        cuda_threads_per_block = deviceProp.maxThreadsPerBlock;
        cuda_threads_per_block = 256;
        cuda_blocks_x = deviceProp.maxGridSize[0];
        cuda_blocks_y = deviceProp.maxGridSize[1];
    }
}

int main(int argc, char **argv)
{
    int NP = 100000;
    int curr_step;

    struct particle *P = (struct particle *)malloc(NP * sizeof(*P));
    assert(P != NULL);

    tyme_t Tmax = 1;
    tyme_t t = 0;
    tyme_t dt = .02;

    struct potential phi;

    show_cuda_devices();

    //ic_random(P, NP, 100);
    ic_circular(P, NP, 100);
    nbody_create_potential(&phi, 1);

    struct particle *Pdev;
    cudaMalloc((void **)&Pdev, NP * sizeof(*P));
    cudaMemcpy(Pdev, P, NP * sizeof(*P), cudaMemcpyHostToDevice);


    for (curr_step = 0, t = dt;
         t < Tmax+dt;
         t += dt, curr_step++)
    {
        //cudaMemcpy(P, Pdev, NP * sizeof(*P), cudaMemcpyDeviceToHost);
        //write_positions(P, NP, curr_step == 0);

        if (t > Tmax) t = Tmax;

        phi.advance(&phi, t);
        cuda_step_all(&phi, Pdev, NP, dt);
    }

    cudaMemcpy(P, Pdev, NP * sizeof(*P), cudaMemcpyDeviceToHost);
    write_positions(P, NP, curr_step == 0);

    return 0;
}

