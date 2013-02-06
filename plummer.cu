/* 
 * Implements a plummer potential.
 * 
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "plummer.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void plummer_init(struct potential *phi)
{
    phi->name = "plummer";
    phi->desc = "plummer potential.";

    phi->create         = plummer_create_potential;
    phi->set_particles  = plummer_set_particles;
    phi->get_particles  = plummer_get_particles;
    phi->step_particles = plummer_step_particles;
    phi->advance        = plummer_advance_potential;
    phi->free           = plummer_free;
    phi->energy         = plummer_energy;
}

void plummer_free(void *phi_data)
{
    struct plummer *plummer = (struct plummer *)phi_data;

    if (plummer->P.P != NULL)
    {
        free(plummer->P.P);
        plummer->P.P = NULL;
    }

    if (plummer->P.Pdev != NULL)
    {
        cudaFree(plummer->P.Pdev);
        plummer->P.Pdev = NULL;
    }
}

void plummer_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct plummer *plummer = (struct plummer *)phi_data;
    plummer->P.N = N;
    plummer->P.P = P;
    plummer->P.P_dirty = 1;
    cudaMalloc((void **)&plummer->P.Pdev, N * sizeof(*(plummer->P.Pdev)));
}

int plummer_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct plummer *plummer = (struct plummer *)phi_data;

    if (plummer->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(plummer->P.P, plummer->P.Pdev, plummer->P.N * sizeof(*plummer->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = plummer->P.P;
    *N = plummer->P.N;

    return 1;
}

__device__ void plummer_accel(struct particle *p, mass_t M, dist_t eps2, acc_t *out)
{
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2];

    const dist_t r2 = x*x + y*y + z*z;

    out[0] = -M * x / pow(r2 + eps2, (real)1.5);
    out[1] = -M * y / pow(r2 + eps2, (real)1.5);
    out[2] = -M * z / pow(r2 + eps2, (real)1.5);
}

__global__ void plummer_cuda_step_all(struct particle *P, size_t NP, mass_t M, dist_t eps2, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    plummer_accel(p, M, eps2, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

int plummer_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct plummer *plummer = (struct plummer *)phi_data;

    static int first_time = 1;
    int nblocks  = (plummer->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = plummer->P.N < cuda_threads_per_block ? plummer->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (plummer->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(plummer->P.Pdev, plummer->P.P, plummer->P.N * sizeof(*(plummer->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        plummer->P.P_dirty = 0;
        plummer->P.Pdev_dirty = 0;
    }

    plummer_cuda_step_all<<<nblocks,nthreads>>>(plummer->P.Pdev, plummer->P.N, plummer->phi.M, plummer->phi.eps2, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    plummer->P.Pdev_dirty = 1;

    return 1;
}

void *plummer_create_potential(size_t N)
{
    struct plummer *plummer = (struct plummer *)calloc(1, sizeof(*plummer));
    assert(plummer != NULL);

    plummer->phi.eps2 = 0.05;
    plummer->phi.M = 5;

    return plummer;
}

int plummer_advance_potential(void *phi_data, tyme_t t)
{
    return 1;
}

real plummer_energy(void *phi_data, struct particle *P, size_t N)
{
    size_t i,j;
    real T=0,U=0;

    struct plummer *plummer = (struct plummer *)phi_data;
    const dist_t eps2 = plummer->phi.eps2;
    const mass_t M = plummer->phi.M;

    const tyme_t dt = 0.1;

    for (i=0; i < N; i++)
    {
        vel_t v2 = 0;
        for (j=0; j < 3; j++)
            v2 += pow(P[i].v[j], 2);
        
        T += v2;

        dist_t r2=0;
        for (j=0; j < 3; j++)
            r2 += pow(P[i].x[j] + P[i].v[j] * dt/2, 2);
            //r2 += pow(P[i].x[j], 2);

        U -= 1 / sqrt(r2 + eps2);
    }

    T *= 0.5;
    U *= M;

    return T+U;
}
