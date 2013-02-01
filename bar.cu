/* Implements a bar potential of the form
 *
 * \phi(r, \theta) = \frac12 \log(1 + r^2(1 + \epsilon \cos^2(\theta - \Omega t)))
 *
 * The derivatives are
 *
 * \frac{\partial}{\partial r} = \frac{r (\epsilon \cos^2(\Omega t - \theta) + 1)}{\epsilon r^2 (\cos^2(\Omega t - \theta) + 1) + 1}
 * \frac{\partial}{\partial \theta} = \frac{\epsilon r^2 \sin(\Omega t - \theta) \cos(\Omega t - \theta)}{r^2 (\epsilon \cos^2(\Omega t - \theta) + 1) + 1}
 *
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "bar.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void bar_init(struct potential *phi)
{
    phi->name = "bar";
    phi->desc = "Bar potential.";

    phi->create         = bar_create_potential;
    phi->set_particles  = bar_set_particles;
    phi->get_particles  = bar_get_particles;
    phi->step_particles = bar_step_particles;
    phi->advance        = bar_advance_potential;
    phi->free           = bar_free;
}

void bar_free(void *phi_data)
{
    struct bar *bar = (struct bar *)phi_data;

    if (bar->P.P != NULL)
    {
        free(bar->P.P);
        bar->P.P = NULL;
    }

    if (bar->P.Pdev != NULL)
    {
        cudaFree(bar->P.Pdev);
        bar->P.Pdev = NULL;
    }
}

void bar_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct bar *bar = (struct bar *)phi_data;
    bar->P.N = N;
    bar->P.P = P;
    bar->P.P_dirty = 1;
    cudaMalloc((void **)&bar->P.Pdev, N * sizeof(*(bar->P.Pdev)));
}

int bar_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct bar *bar = (struct bar *)phi_data;

    if (bar->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(bar->P.P, bar->P.Pdev, bar->P.N * sizeof(*bar->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = bar->P.P;
    *N = bar->P.N;

    return 1;
}

__device__ void bar_accel(struct particle *p, tyme_t t, mass_t M, freq_t omega, acc_t *out)
{
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2];

    const dist_t r2 = x*x + y*y + z*z;
    const dist_t r = sqrt(r2);
    const dist_t theta = atan2(y, x);

    const dist_t dr_dx = x / r;
    const dist_t dr_dy = y / r;

    const dist_t dtheta_dx = -y / pow(r,2);
    const dist_t dtheta_dy =  x / pow(r,2);

    const real cs = cos(omega*t - theta);
    const real sn = sin(omega*t - theta);

    const dist_t q0 = M * pow(cs,2) + 1;
    const dist_t q1 = r2 * q0 + 1;

    const dist_t d_dr     = r * q0 / q1;
    const dist_t d_dtheta = M * r2 * sn * cs / q1;

    out[0] = -0.5 * (d_dr * dr_dx + d_dtheta * dtheta_dx);
    out[1] = -0.5 * (d_dr * dr_dy + d_dtheta * dtheta_dy);
    out[2] = 0;
}

__global__ void bar_cuda_step_all(struct particle *P, size_t NP, tyme_t t, mass_t M, freq_t omega, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    bar_accel(p, t, M, omega, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

int bar_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct bar *bar = (struct bar *)phi_data;

    static int first_time = 1;
    int nblocks  = (bar->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = bar->P.N < cuda_threads_per_block ? bar->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (bar->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(bar->P.Pdev, bar->P.P, bar->P.N * sizeof(*(bar->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        bar->P.P_dirty = 0;
        bar->P.Pdev_dirty = 0;
    }

    bar_cuda_step_all<<<nblocks,nthreads>>>(bar->P.Pdev, bar->P.N, bar->phi.t, bar->phi.M, bar->phi.omega, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    bar->P.Pdev_dirty = 1;

    return 1;
}

void *bar_create_potential(size_t N)
{
    struct bar *bar = (struct bar *)calloc(1, sizeof(*bar));
    assert(bar != NULL);

    bar->phi.t = 0;
    bar->phi.omega = 0.1;
    bar->phi.M = 15;

    return bar;
}

int bar_advance_potential(void *phi_data, tyme_t t)
{
    struct bar *bar = (struct bar *)phi_data;

    assert(t >= bar->phi.t);

    if (t != bar->phi.t)
    {
        fprintf(stderr, "Advancing potential to time %f\n", t);
        bar->phi.t = t;
    }

    return 1;
}


