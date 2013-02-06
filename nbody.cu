#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "nbody.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void nbody_init(struct potential *phi)
{
    phi->name = "nbody";
    phi->desc = "Two body potential.";

    phi->create         = nbody_create_potential;
    phi->set_particles  = nbody_set_particles;
    phi->get_particles  = nbody_get_particles;
    phi->step_particles = nbody_step_particles;
    phi->advance        = nbody_advance_potential;
    phi->free           = nbody_free;
}

void nbody_free(void *phi_data)
{
    struct nbody *nbody = (struct nbody *)phi_data;

    if (nbody->P.P != NULL)
    {
        free(nbody->P.P);
        nbody->P.P = NULL;
    }

    if (nbody->P.Pdev != NULL)
    {
        cudaFree(nbody->P.Pdev);
        nbody->P.Pdev = NULL;
    }

    if (nbody->phi.P != NULL)
    {
        cudaFree(nbody->P.P);
        nbody->phi.P = NULL;
    }

    if (nbody->phi.Pdev != NULL)
    {
        cudaFree(nbody->P.Pdev);
        nbody->phi.Pdev = NULL;
    }
}

void nbody_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct nbody *nbody = (struct nbody *)phi_data;
    nbody->P.N = N;
    nbody->P.P = P;
    nbody->P.P_dirty = 1;
    cudaMalloc((void **)&nbody->P.Pdev, N * sizeof(*(nbody->P.Pdev)));
}

int nbody_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct nbody *nbody = (struct nbody *)phi_data;

    if (nbody->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(nbody->P.P, nbody->P.Pdev, nbody->P.N * sizeof(*nbody->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s\n", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = nbody->P.P;
    *N = nbody->P.N;

    return 1;
}

__device__ void nbody_accel(struct particle *p, struct massive_particle *Pm, size_t NPm, const dist_t eps2, acc_t *out)
{
    size_t i;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

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
        #define Fhat(x) round((Pm[i].m*(x) * pow(eps2 + r2, -3./2.)))
#else
        #define Fhat(x) ((Pm[i].m*(x) * pow(eps2 + r2, -1.5F)))
#endif
        out[0] += Fhat(dx);
        out[1] += Fhat(dy);
        out[2] += Fhat(dz);

        //int j; for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    //fprintf(stderr, "\n");
}

__global__ void nbody_cuda_step_all(struct particle *P, size_t NP, struct massive_particle *Pm, size_t NPm, dist_t eps2, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle p = P[pi];
    for (i=0; i < 3; i++) p.x[i] += p.v[i] * dt/2;
    nbody_accel(&p, Pm, NPm, eps2, a);
    for (i=0; i < 3; i++) p.v[i] +=   a[i] * dt;
    for (i=0; i < 3; i++) p.x[i] += p.v[i] * dt/2;
    P[pi] = p;
}

int nbody_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct nbody *nbody = (struct nbody *)phi_data;

    static int first_time = 1;
    int nblocks  = (nbody->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = nbody->P.N < cuda_threads_per_block ? nbody->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (nbody->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(nbody->P.Pdev, nbody->P.P, nbody->P.N * sizeof(*(nbody->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        nbody->P.P_dirty = 0;
        nbody->P.Pdev_dirty = 0;
    }

    /* If the local potential particles are dirty, copy them to the device */
    if (nbody->phi.P_dirty)
    {
        //fprintf(stderr, "Updating device potential.\n");
        cudaMemcpy(nbody->phi.Pdev, nbody->phi.P, nbody->phi.N * sizeof(*nbody->phi.Pdev), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }
        nbody->phi.P_dirty = 0;
        nbody->phi.Pdev_dirty = 0;
    }

    nbody_cuda_step_all<<<nblocks,nthreads>>>(nbody->P.Pdev, nbody->P.N, nbody->phi.Pdev, nbody->phi.N, nbody->phi.eps2, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        return 0;
    }

    nbody->P.Pdev_dirty = 1;

    return 1;
}


static void nbody_accel_massive(struct nbody *nbody, struct massive_particle *p, acc_t *out)
{
    size_t i;
    struct massive_particle *mp = nbody->phi.P;
    const dist_t e2 = nbody->phi.eps2;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

    for (i=0; i < nbody->phi.N; i++)
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
#endif
        out[0] += FhatMP(dx);
        out[1] += FhatMP(dy);
        out[2] += FhatMP(dz);

        //fprintf(stderr, "%i %e %e %e\n", i, FhatMP(dx), FhatMP(dy), FhatMP(dz));
    }
    //fprintf(stderr, "\n");
}

static void nbody_drift(struct nbody *nbody, struct massive_particle *p, tyme_t dt)
{
    int i;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt;
}

static void nbody_kick(struct nbody *nbody, struct massive_particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    nbody_accel_massive(nbody, p, a);
    for (i=0; i < 3; i++) p->v[i] += a[i] * dt;
}

static void nbody_step_potential(struct nbody *nbody, tyme_t dt)
{
    size_t i;
    for (i=0; i < nbody->phi.N; i++) nbody_drift(nbody, nbody->phi.P+i, dt/2);
    for (i=0; i < nbody->phi.N; i++) nbody_kick (nbody, nbody->phi.P+i, dt);
    for (i=0; i < nbody->phi.N; i++) nbody_drift(nbody, nbody->phi.P+i, dt/2);
    nbody->phi.P_dirty = 1;
}

void *nbody_create_potential(size_t N)
{
    struct nbody *nbody = (struct nbody *)calloc(1, sizeof(*nbody));
    assert(nbody != NULL);

    assert(N == 2);

    nbody->phi.N = N;
    nbody->phi.eps2 = 3.5;
    nbody->phi.dt = 0.1;
    nbody->phi.t = 0;

    nbody->phi.P = (struct massive_particle *)malloc(N * sizeof(*nbody->phi.P));

    mass_t m  = 1e2;    // Mass
    dist_t x  = 3;     // Position rel. to 0
    dist_t d = 2*x;     // Relative particle distance

    d = sqrt(d*d + nbody->phi.eps2);    // Softened distance
    vel_t v = sqrt(m / (2*d));         // Velocity


    nbody->phi.P[0].x[0] = x;
    nbody->phi.P[0].x[1] = 0;
    nbody->phi.P[0].x[2] = 0;

    nbody->phi.P[0].v[0] = 0;
    nbody->phi.P[0].v[1] = v;
    nbody->phi.P[0].v[2] = 0;

    nbody->phi.P[1].x[0] = -x;
    nbody->phi.P[1].x[1] = 0;
    nbody->phi.P[1].x[2] = 0;

    nbody->phi.P[1].v[0] = 0;
    nbody->phi.P[1].v[1] = -v;
    nbody->phi.P[1].v[2] = 0;

    nbody->phi.P[0].m = m;
    nbody->phi.P[1].m = m;

    int i; for (i=0; i < nbody->phi.N; i++) nbody_kick (nbody, nbody->phi.P+i, nbody->phi.dt/2);

    cudaMalloc((void **)&nbody->phi.Pdev,     N * sizeof(*(nbody->phi.Pdev)));
    cudaMemcpy(nbody->phi.Pdev, nbody->phi.P, N * sizeof(*(nbody->phi.Pdev)), cudaMemcpyHostToDevice);

    return nbody;
}

int nbody_advance_potential(void *phi_data, tyme_t t)
{
    struct nbody *nbody = (struct nbody *)phi_data;

    assert(t >= nbody->phi.t);

    if (t != nbody->phi.t)
    {
        fprintf(stderr, "Advancing potential to time %f\n", t);

        for (; nbody->phi.t < t-nbody->phi.dt; nbody->phi.t += nbody->phi.dt)
            nbody_step_potential(nbody, nbody->phi.dt);

        if (t > nbody->phi.t)
        {
            nbody_step_potential(nbody, t - nbody->phi.t);
            nbody->phi.t += t - nbody->phi.t;
        }
    }

    return 1;
}


