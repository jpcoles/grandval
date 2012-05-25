#include <cuda.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "grandval.h"
#include "nbody.h"
#include "ic.h"
#include "io.h"

#define WITH_CUDA 1

size_t cuda_threads_per_block;
size_t cuda_blocks_x, cuda_blocks_y;


#if 0
void step(struct potential *phi, struct particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    phi->accel(phi, p, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}
#endif

void write_positions(struct particle *P, size_t NP, int first_output)
{
    size_t i;
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

#if 0
void step_all(struct potential *phi, struct particle *P, size_t NP, tyme_t dt)
{
    size_t i;
    #pragma omp parallel for
    for (i=0; i < NP; i++)
        step(phi, P + i, dt);
}
#endif

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
        cuda_threads_per_block = 640;
        cuda_blocks_x = deviceProp.maxGridSize[0];
        cuda_blocks_y = deviceProp.maxGridSize[1];
    }
}

int main(int argc, char **argv)
{
    size_t NP = 1000000;
    int curr_step;
    int Ncaptures = 100;

    struct nbody *nbody;

    int red[3] = {255,0,0};
    int grey[3] = {255,255,255};

    struct particle *P = (struct particle *)malloc(NP * sizeof(*P));
    assert(P != NULL);

    struct image image;
    image.nc = 512;
    image.nr = 512;
    image.image  = (unsigned char *)calloc(3 * image.nr * image.nc, sizeof(*image.image));
    image.hist = (int *)calloc(1 * image.nr * image.nc, sizeof(*image.hist));

    tyme_t Tmax = 60;
    tyme_t t = 0;
    tyme_t dt = .02;

    show_cuda_devices();

    //ic_random(P, NP, 100);
    ic_circular(P, NP, 100);

    nbody = nbody_init();
    nbody_set_particles(nbody, P, NP);
    nbody_create_potential(nbody, 2);

    double t_next_capture = Tmax / Ncaptures;
    int curr_capture = 0;
    if (Ncaptures)
    {
        capture(100, P, NP, &image, 1, grey);
        capture_massive(100, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_snapshot(curr_capture, &image);
        curr_capture++;
    }

    for (curr_step = 0, t = dt;
         t < Tmax+dt;
         t += dt, curr_step++)
    {
        if (t > Tmax) t = Tmax;

        nbody_step_all(nbody, dt);
        //phi.advance(&phi, t);
        nbody_advance_potential(nbody, t);

        nbody_get_particles(nbody, &P, &NP);

        if (Ncaptures && t > t_next_capture)
        {
            //write_positions(P, NP, curr_step == 0);
            capture(100, P, NP, &image, 1, grey);
            capture_massive(100, nbody->phi.P, nbody->phi.N, &image, 0, red);
            save_snapshot(curr_capture, &image);
            curr_capture++;
            t_next_capture += Tmax / Ncaptures;
        }
    }

    nbody_get_particles(nbody, &P, &NP);
    if (Ncaptures)
    {
        capture(100, P, NP, &image, 1, grey);
        capture_massive(100, nbody->phi.P, nbody->phi.N, &image, 0, red);
        save_snapshot(curr_capture, &image);
    }
    //write_positions(P, NP, curr_step == 0);

    nbody_free(nbody);

    return 0;
}

