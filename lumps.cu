#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "lumps.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void lumps_init(struct potential *phi)
{
    phi->name = "lumps";
    phi->desc = "Two body potential.";

    phi->create         = lumps_create_potential;
    phi->set_particles  = lumps_set_particles;
    phi->get_particles  = lumps_get_particles;
    phi->step_particles = lumps_step_particles;
    phi->advance        = lumps_advance_potential;
    phi->free           = lumps_free;
}

void lumps_free(void *phi_data)
{
    if (phi_data == NULL) return;

    struct lumps *lumps = (struct lumps *)phi_data;

    if (lumps->P.P != NULL)
    {
        free(lumps->P.P);
        lumps->P.P = NULL;
    }

    if (lumps->P.Pdev != NULL)
    {
        cudaFree(lumps->P.Pdev);
        lumps->P.Pdev = NULL;
    }

    if (lumps->phi.P != NULL)
    {
        cudaFree(lumps->P.P);
        lumps->phi.P = NULL;
    }

    if (lumps->phi.Pdev != NULL)
    {
        cudaFree(lumps->P.Pdev);
        lumps->phi.Pdev = NULL;
    }
}

void lumps_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct lumps *lumps = (struct lumps *)phi_data;
    lumps->P.N = N;
    lumps->P.P = P;
    lumps->P.P_dirty = 1;
    cudaMalloc((void **)&lumps->P.Pdev, N * sizeof(*(lumps->P.Pdev)));
}

int lumps_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct lumps *lumps = (struct lumps *)phi_data;

    if (lumps->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(lumps->P.P, lumps->P.Pdev, lumps->P.N * sizeof(*lumps->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = lumps->P.P;
    *N = lumps->P.N;

    return 1;
}

__device__ void lumps_accel(struct particle *p, struct massive_particle *Pm, size_t NPm, const mass_t M, const dist_t eps2, const dist_t eps2_g,  acc_t *out)
{
    size_t i;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];
    const dist_t r2_g = x*x + y*y + z*z;

    for (i=0; i < NPm; i++)
    {
        const dist_t dx = Pm[i].x[0] - x;
        const dist_t dy = Pm[i].x[1] - y;
        const dist_t dz = Pm[i].x[2] - z;
        const dist_t r2 = dx*dx + dy*dy + dz*dz;

        //double rinv = 1.0 / sqrt(r2);

        //double re = sqrt(r2) * sqrt(e2);

        #define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
        #define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
        //double gradphi = GRADPHI(re);

        // round'ing is important. Conserves momentum perfectly.
        //#define Fhat(x) round(((Pm[i].m*e2)*(((x)*rinv)*gradphi)))
#if WITH_INTEGERS
        #define Fhat(x) round((Pm[i].m*(x) * pow(eps2 + r2, (real)-3./2.)))
#else
        #define Fhat(x) ((Pm[i].m*(x) * pow(eps2 + r2, (real)-1.5)))
        #define F_global_hat(x) ((-M*(x) * pow(eps2_g + r2_g, (real)-1.5)))
#endif
        out[0] += Fhat(dx)+F_global_hat(x);
        out[1] += Fhat(dy)+F_global_hat(y);
        out[2] += Fhat(dz)+F_global_hat(z);

        //int j; for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    //fprintf(stderr, "\n");
}

__global__ void lumps_cuda_step_all(struct particle *P, size_t NP, struct massive_particle *Pm, size_t NPm, mass_t M, dist_t eps2, dist_t eps2_g,  tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle p = P[pi];
    for (i=0; i < 3; i++) p.x[i] += p.v[i] * dt/2;
    lumps_accel(&p, Pm, NPm, M, eps2, eps2_g, a);
    for (i=0; i < 3; i++) p.v[i] +=   a[i] * dt;
    for (i=0; i < 3; i++) p.x[i] += p.v[i] * dt/2;
    P[pi] = p;
}

int lumps_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct lumps *lumps = (struct lumps *)phi_data;

    static int first_time = 1;
    int nblocks  = (lumps->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = lumps->P.N < cuda_threads_per_block ? lumps->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (lumps->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(lumps->P.Pdev, lumps->P.P, lumps->P.N * sizeof(*(lumps->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        lumps->P.P_dirty = 0;
        lumps->P.Pdev_dirty = 0;
    }

    /* If the local potential particles are dirty, copy them to the device */
    if (lumps->phi.P_dirty)
    {
        //fprintf(stderr, "Updating device potential.\n");
        cudaMemcpy(lumps->phi.Pdev, lumps->phi.P, lumps->phi.N * sizeof(*lumps->phi.Pdev), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }
        lumps->phi.P_dirty = 0;
        lumps->phi.Pdev_dirty = 0;
    }

    lumps_cuda_step_all<<<nblocks,nthreads>>>(lumps->P.Pdev, lumps->P.N, lumps->phi.Pdev, lumps->phi.N, lumps->phi.M,  lumps->phi.eps2, lumps->phi.eps2_g, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    lumps->P.Pdev_dirty = 1;

    return 1;
}


static void lumps_accel_massive(struct lumps *lumps, struct massive_particle *p, acc_t *out)
{
    size_t i;
    struct massive_particle *mp = lumps->phi.P;
    const dist_t e2 = lumps->phi.eps2;
    const dist_t e2_g = lumps->phi.eps2_g;
    const mass_t M = lumps->phi.M;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];
    const dist_t r2_g = x*x + y*y + z*z;

    for (i=0; i < lumps->phi.N; i++)
    {
        if (mp+i == p) continue;

        const dist_t dx = mp[i].x[0] - x;
        const dist_t dy = mp[i].x[1] - y;
        const dist_t dz = mp[i].x[2] - z;
        const dist_t r2 = dx*dx + dy*dy + dz*dz;

        dist_t rinv = 1.0 / sqrt(r2);

        dist_t re = sqrt(r2); // * sqrt(e2);

        #define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
        #define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
        double gradphi = GRADPHI(re);

        // round'ing is important. Conserves momentum perfectly.
        //#define Fhat(x) round(((mp[i].m*e2)*(((x)*rinv)*gradphi)))
#if WITH_INTEGERS
        #define Fhat(x) round((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#else
        #define FhatMP(x) ((mp[i].m*(x) * pow(e2 + r2, -1.5F)))
        #define FMP_global_hat(x) ((-M*(x) * pow(e2_g + r2_g, (real)-1.5)))

#endif
        out[0] += FhatMP(dx) + FMP_global_hat(x);
        out[1] += FhatMP(dy) + FMP_global_hat(y);
        out[2] += FhatMP(dz) + FMP_global_hat(z);

        //fprintf(stderr, "%i %e %e %e\n", i, FhatMP(dx), FhatMP(dy), FhatMP(dz));
    }
    //fprintf(stderr, "\n");
}

static void lumps_drift(struct lumps *lumps, struct massive_particle *p, tyme_t dt)
{
    int i;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt;
}

static void lumps_kick(struct lumps *lumps, struct massive_particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    lumps_accel_massive(lumps, p, a);
    for (i=0; i < 3; i++) p->v[i] += a[i] * dt;
}

static void lumps_step_potential(struct lumps *lumps, tyme_t dt)
{
    size_t i;
    for (i=0; i < lumps->phi.N; i++) lumps_drift(lumps, lumps->phi.P+i, dt/2);
    for (i=0; i < lumps->phi.N; i++) lumps_kick (lumps, lumps->phi.P+i, dt);
    for (i=0; i < lumps->phi.N; i++) lumps_drift(lumps, lumps->phi.P+i, dt/2);
    lumps->phi.P_dirty = 1;
}

void *lumps_create_potential(size_t N)
{
    struct lumps *lumps = (struct lumps *)calloc(1, sizeof(*lumps));
    assert(lumps != NULL);

    assert(N == 2);

    lumps->phi.N = N;
    lumps->phi.M = 2e2;
    lumps->phi.eps2 = 3.5;
    lumps->phi.eps2_g= 5;
    lumps->phi.dt = 0.1;
    lumps->phi.t = 0;

    lumps->phi.P = (struct massive_particle *)malloc(N * sizeof(*lumps->phi.P));

    mass_t m  = 1e2;    // Mass
    dist_t x  = 3;     // Position rel. to 0
    dist_t d = 2*x;     // Relative particle distance

    d = sqrt(d*d + lumps->phi.eps2);    // Softened distance
    vel_t v = sqrt(m / (2*d));         // Velocity


    lumps->phi.P[0].x[0] = x;
    lumps->phi.P[0].x[1] = 0;
    lumps->phi.P[0].x[2] = 0;

    lumps->phi.P[0].v[0] = 0;
    lumps->phi.P[0].v[1] = v;
    lumps->phi.P[0].v[2] = 0;

    lumps->phi.P[1].x[0] = -x;
    lumps->phi.P[1].x[1] = 0;
    lumps->phi.P[1].x[2] = 0;

    lumps->phi.P[1].v[0] = 0;
    lumps->phi.P[1].v[1] = -v;
    lumps->phi.P[1].v[2] = 0;

    lumps->phi.P[0].m = m;
    lumps->phi.P[1].m = m;

    int i; for (i=0; i < lumps->phi.N; i++) lumps_kick (lumps, lumps->phi.P+i, lumps->phi.dt/2);

    cudaMalloc((void **)&lumps->phi.Pdev,     N * sizeof(*(lumps->phi.Pdev)));
    cudaMemcpy(lumps->phi.Pdev, lumps->phi.P, N * sizeof(*(lumps->phi.Pdev)), cudaMemcpyHostToDevice);

    return lumps;
}

int lumps_advance_potential(void *phi_data, tyme_t t)
{
    struct lumps *lumps = (struct lumps *)phi_data;

    assert(t >= lumps->phi.t);

    if (t != lumps->phi.t)
    {
        fprintf(stderr, "Advancing potential to time %f\n", t);

        for (; lumps->phi.t < t-lumps->phi.dt; lumps->phi.t += lumps->phi.dt)
            lumps_step_potential(lumps, lumps->phi.dt);

        if (t > lumps->phi.t)
        {
            lumps_step_potential(lumps, t - lumps->phi.t);
            lumps->phi.t += t - lumps->phi.t;
        }
    }

    return 1;
}


