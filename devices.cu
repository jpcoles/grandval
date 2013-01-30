#include <stdio.h>
#include <cuda.h>
#include "grandval.h"
#include "devices.h"

size_t cuda_threads_per_block;
size_t cuda_blocks_x, cuda_blocks_y;

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
        printf("              totalGlobalMem %ld (approx. %ld particles).\n", deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (sizeof(struct particle)));
        printf("              maxThreadsPerMultiProcessor %d.\n", deviceProp.maxThreadsPerMultiProcessor);
        printf("              maxThreadsPerBlock %d.\n", deviceProp.maxThreadsPerBlock);
        printf("              maxThreadsDim %d,%d,%d.\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
        printf("              maxGridSize %d,%d,%d.\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
        printf("              warpSize %d\n", deviceProp.warpSize);
        printf("              regsPerBlock %d\n", deviceProp.regsPerBlock);

        cuda_threads_per_block = deviceProp.maxThreadsPerBlock;
        cuda_threads_per_block = 512;
        cuda_blocks_x = deviceProp.maxGridSize[0];
        cuda_blocks_y = deviceProp.maxGridSize[1];
    }
}

int select_cuda_device(int device)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0)
    {
        printf("No CUDA devices found.\n");
        return 0;
    }

    if (device >= deviceCount)
        return 0;

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    cudaError_t err = cudaGetLastError();
    if (err)
    {
        eprintf("Kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    printf("Device %d has compute capability %d.%d.\n", device, deviceProp.major, deviceProp.minor);
    printf("              totalGlobalMem %ld (approx. %ld particles).\n", deviceProp.totalGlobalMem, deviceProp.totalGlobalMem / (sizeof(struct particle)));
    printf("              maxThreadsPerMultiProcessor %d.\n", deviceProp.maxThreadsPerMultiProcessor);
    printf("              maxThreadsPerBlock %d.\n", deviceProp.maxThreadsPerBlock);
    printf("              maxThreadsDim %d,%d,%d.\n", deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
    printf("              maxGridSize %d,%d,%d.\n", deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
    printf("              warpSize %d\n", deviceProp.warpSize);
    printf("              regsPerBlock %d\n", deviceProp.regsPerBlock);

    cuda_threads_per_block = deviceProp.maxThreadsPerBlock;
    cuda_threads_per_block = 512;
    cuda_blocks_x = deviceProp.maxGridSize[0];
    cuda_blocks_y = deviceProp.maxGridSize[1];

    return 1;
}

void show_devices()
{
    show_cuda_devices();
}

