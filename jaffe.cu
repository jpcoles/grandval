/* 
 * Implements a jaffe potential.
 * 
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "jaffe.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void jaffe_init(struct potential *phi)
{
    phi->name = "jaffe";
    phi->desc = "Plummer potential.";

    phi->create         = jaffe_create_potential;
    phi->set_particles  = jaffe_set_particles;
    phi->get_particles  = jaffe_get_particles;
    phi->step_particles = jaffe_step_particles;
    phi->advance        = jaffe_advance_potential;
    phi->free           = jaffe_free;
    phi->energy         = jaffe_energy;
}

void jaffe_free(void *phi_data)
{
    if (phi_data == NULL) return;

    struct jaffe *jaffe = (struct jaffe *)phi_data;

    if (jaffe->P.P != NULL)
    {
        free(jaffe->P.P);
        jaffe->P.P = NULL;
    }

    if (jaffe->P.Pdev != NULL)
    {
        cudaFree(jaffe->P.Pdev);
        jaffe->P.Pdev = NULL;
    }
}

void jaffe_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct jaffe *jaffe = (struct jaffe *)phi_data;
    jaffe->P.N = N;
    jaffe->P.P = P;
    jaffe->P.P_dirty = 1;
    cudaMalloc((void **)&jaffe->P.Pdev, N * sizeof(*(jaffe->P.Pdev)));
}

int jaffe_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct jaffe *jaffe = (struct jaffe *)phi_data;

    if (jaffe->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(jaffe->P.P, jaffe->P.Pdev, jaffe->P.N * sizeof(*jaffe->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = jaffe->P.P;
    *N = jaffe->P.N;

    return 1;
}

__device__ void jaffe_accel(struct particle *p, mass_t M, dist_t eps, acc_t *out)
{
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2];

    const dist_t r2 = x*x + y*y + z*z;

    out[0] = -M * x / (r2*(sqrt(r2)+eps));
    out[1] = -M * y / (r2*(sqrt(r2)+eps));
    out[2] = -M * z / (r2*(sqrt(r2)+eps));
}

__global__ void jaffe_cuda_step_all(struct particle *P, size_t NP, mass_t M, dist_t eps, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    jaffe_accel(p, M, eps, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

int jaffe_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct jaffe *jaffe = (struct jaffe *)phi_data;

    static int first_time = 1;
    int nblocks  = (jaffe->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = jaffe->P.N < cuda_threads_per_block ? jaffe->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (jaffe->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(jaffe->P.Pdev, jaffe->P.P, jaffe->P.N * sizeof(*(jaffe->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        jaffe->P.P_dirty = 0;
        jaffe->P.Pdev_dirty = 0;
    }

    jaffe_cuda_step_all<<<nblocks,nthreads>>>(jaffe->P.Pdev, jaffe->P.N, jaffe->phi.M, jaffe->phi.eps, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s", err, cudaGetErrorString(err));
        return 0;
    }

    jaffe->P.Pdev_dirty = 1;

    return 1;
}

void *jaffe_create_potential(size_t N)
{
    struct jaffe *jaffe = (struct jaffe *)calloc(1, sizeof(*jaffe));
    assert(jaffe != NULL);

    jaffe->phi.eps = 0.22;
    jaffe->phi.M = 5;

    return jaffe;
}

int jaffe_advance_potential(void *phi_data, tyme_t t)
{
    return 1;
}

real jaffe_energy(void *phi_data, struct particle *P, size_t N)
{
    size_t i,j;
    real T=0,U=0;

    struct jaffe *jaffe = (struct jaffe *)phi_data;
    const dist_t eps = jaffe->phi.eps;
    const mass_t M = jaffe->phi.M;

    for (i=0; i < N; i++)
    {
        vel_t v2 = 0;
        for (j=0; j < 3; j++)
            v2 += pow(P[i].v[j], 2);
        
        T += v2;

        dist_t r2=0;
        for (j=0; j < 3; j++)
            r2 += pow(P[i].x[j], 2);

        U +=log(sqrt(r2)/(sqrt(r2)+eps)) ;
    }

    T *= 0.5;
    U *= M/eps;

    return T+U;
}
