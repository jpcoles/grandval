/* 
 * Implements a hernquist potential.
 * 
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "hernquist.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void hernquist_init(struct potential *phi)
{
    phi->name = "hernquist";
    phi->desc = "hernquist potential.";

    phi->create         = hernquist_create_potential;
    phi->set_particles  = hernquist_set_particles;
    phi->get_particles  = hernquist_get_particles;
    phi->step_particles = hernquist_step_particles;
    phi->advance        = hernquist_advance_potential;
    phi->free           = hernquist_free;
    phi->energy         = hernquist_energy;
}

void hernquist_free(void *phi_data)
{
    if (phi_data == NULL) return;

    struct hernquist *hernquist = (struct hernquist *)phi_data;

    if (hernquist->P.P != NULL)
    {
        free(hernquist->P.P);
        hernquist->P.P = NULL;
    }

    if (hernquist->P.Pdev != NULL)
    {
        cudaFree(hernquist->P.Pdev);
        hernquist->P.Pdev = NULL;
    }
}

void hernquist_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct hernquist *hernquist = (struct hernquist *)phi_data;
    hernquist->P.N = N;
    hernquist->P.P = P;
    hernquist->P.P_dirty = 1;
    cudaMalloc((void **)&hernquist->P.Pdev, N * sizeof(*(hernquist->P.Pdev)));
}

int hernquist_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct hernquist *hernquist = (struct hernquist *)phi_data;

    if (hernquist->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(hernquist->P.P, hernquist->P.Pdev, hernquist->P.N * sizeof(*hernquist->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = hernquist->P.P;
    *N = hernquist->P.N;

    return 1;
}

__device__ void hernquist_accel(struct particle *p, mass_t M, dist_t eps, acc_t *out)
{
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2];

    const dist_t r = sqrt(x*x + y*y + z*z);

    out[0] = -M * x / (pow(r + eps,(real)2)*r);
    out[1] = -M * y / (pow(r + eps,(real)2)*r);
    out[2] = -M * z / (pow(r + eps,(real)2)*r);
}

__global__ void hernquist_cuda_step_all(struct particle *P, size_t NP, mass_t M, dist_t eps, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    hernquist_accel(p, M, eps, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

int hernquist_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct hernquist *hernquist = (struct hernquist *)phi_data;

    static int first_time = 1;
    int nblocks  = (hernquist->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = hernquist->P.N < cuda_threads_per_block ? hernquist->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (hernquist->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(hernquist->P.Pdev, hernquist->P.P, hernquist->P.N * sizeof(*(hernquist->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        hernquist->P.P_dirty = 0;
        hernquist->P.Pdev_dirty = 0;
    }

    hernquist_cuda_step_all<<<nblocks,nthreads>>>(hernquist->P.Pdev, hernquist->P.N, hernquist->phi.M, hernquist->phi.eps, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    hernquist->P.Pdev_dirty = 1;

    return 1;
}

void *hernquist_create_potential(size_t N)
{
    struct hernquist *hernquist = (struct hernquist*)calloc(1, sizeof(*hernquist));
    assert(hernquist != NULL);

    hernquist->phi.eps = 0.22;
    hernquist->phi.M = 5;

    return hernquist;
}

int hernquist_advance_potential(void *phi_data, tyme_t t)
{
    return 1;
}

real hernquist_energy(void *phi_data, struct particle *P, size_t N)
{
    size_t i,j;
    real T=0,U=0;

    struct hernquist *hernquist = (struct hernquist *)phi_data;
    const dist_t eps = hernquist->phi.eps;
    const mass_t M = hernquist->phi.M;

    const tyme_t dt = 0.1;

    for (i=0; i < N; i++)
    {
        vel_t v2 = 0;
        for (j=0; j < 3; j++)
            v2 += pow(P[i].v[j], 2);
        
        T += v2;

        dist_t r2=0;
        for (j=0; j < 3; j++)
            r2 += pow(P[i].x[j], 2);
            

        U -= 1 / (sqrt(r2) + eps);
    }

    T *= 0.5;
    U *= M;
printf("%f %f\n", T, U); 
    return T+U;

}
