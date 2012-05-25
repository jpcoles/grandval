#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "nbody.h"

static const size_t cuda_threads_per_block = 640;

struct nbody *nbody_init()
{
    struct nbody *nb = (struct nbody *)calloc(1, sizeof(*nb));
    assert(nb != NULL);
    return nb;
}

void nbody_free(struct nbody *nbody)
{
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

void nbody_set_particles(struct nbody *nbody, struct particle *P, size_t N)
{
    nbody->P.N = N;
    nbody->P.P = P;
    nbody->P.P_dirty = 1;
    cudaMalloc((void **)&nbody->P.Pdev, N * sizeof(*(nbody->P.Pdev)));
}

void nbody_get_particles(struct nbody *nbody, struct particle **P, size_t *N)
{
    if (nbody->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            fprintf(stderr, "Kernel call had error status %i: %s\n", err, cudaGetErrorString(err));
        }
        cudaMemcpy(nbody->P.P, nbody->P.Pdev, nbody->P.N * sizeof(*nbody->P.P), cudaMemcpyDeviceToHost);
    }

    *P = nbody->P.P;
    *N = nbody->P.N;
}

__device__ void nbody_accel(struct particle *p, struct massive_particle *Pm, size_t NPm, const double eps2, acc_t *out)
{
    size_t i;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

    for (i=0; i < NPm; i++)
    {
        const double dx = Pm[i].x[0] - x;
        const double dy = Pm[i].x[1] - y;
        const double dz = Pm[i].x[2] - z;
        const double r2 = dx*dx + dy*dy + dz*dz;

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
        #define Fhat(x) ((Pm[i].m*(x) * pow(eps2 + r2, -1.5)))
#endif
        out[0] += Fhat(dx);
        out[1] += Fhat(dy);
        out[2] += Fhat(dz);

        //int j; for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    //fprintf(stderr, "\n");
}

__global__ void nbody_cuda_step_all(struct particle *P, size_t NP, struct massive_particle *Pm, size_t NPm, double eps2, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    nbody_accel(p, Pm, NPm, eps2, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

void nbody_step_all(struct nbody *nbody, tyme_t dt)
{
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
        cudaMemcpy(nbody->P.Pdev, nbody->P.P, nbody->P.N * sizeof(*(nbody->P.Pdev)), cudaMemcpyHostToDevice);
        nbody->P.P_dirty = 0;
        nbody->P.Pdev_dirty = 0;
    }

    /* If the local potential particles are dirty, copy them to the device */
    if (nbody->phi.P_dirty)
    {
        cudaMemcpy(nbody->phi.Pdev, nbody->phi.P, nbody->phi.N * sizeof(*nbody->phi.Pdev), cudaMemcpyHostToDevice);
        nbody->phi.P_dirty = 0;
        nbody->phi.Pdev_dirty = 0;
    }

    cudaGetLastError();
    nbody->P.Pdev_dirty = 1;
    nbody_cuda_step_all<<<nblocks,nthreads>>>(nbody->P.Pdev, nbody->P.N, nbody->phi.Pdev, nbody->phi.N, nbody->phi.eps2, dt);
}


static void nbody_accel_massive(struct nbody *nbody, struct massive_particle *p, acc_t *out)
{
    size_t i;
    struct massive_particle *mp = nbody->phi.P;
    const double e2 = nbody->phi.eps2;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

    for (i=0; i < nbody->phi.N; i++)
    {
        if (mp+i == p) continue;

        const double dx = mp[i].x[0] - x;
        const double dy = mp[i].x[1] - y;
        const double dz = mp[i].x[2] - z;
        const double r2 = dx*dx + dy*dy + dz*dz;

        double rinv = 1.0 / sqrt(r2);

        double re = sqrt(r2) * sqrt(e2);

        #define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
        #define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
        double gradphi = GRADPHI(re);

        // round'ing is important. Conserves momentum perfectly.
        //#define Fhat(x) round(((mp[i].m*e2)*(((x)*rinv)*gradphi)))
#if WITH_INTEGERS
        #define Fhat(x) round((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#else
        #define FhatMP(x) ((mp[i].m*(x) * pow(e2 + r2, -1.5)))
#endif
        out[0] += FhatMP(dx);
        out[1] += FhatMP(dy);
        out[2] += FhatMP(dz);

        //int j; for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    //fprintf(stderr, "\n");
}

static void nbody_drift(struct nbody *nbody, struct massive_particle *p, tyme_t dt)
{
    int i;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
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
    for (i=0; i < nbody->phi.N; i++) nbody_drift(nbody, nbody->phi.P+i, dt);
    for (i=0; i < nbody->phi.N; i++) nbody_kick (nbody, nbody->phi.P+i, dt);
    for (i=0; i < nbody->phi.N; i++) nbody_drift(nbody, nbody->phi.P+i, dt);
    nbody->phi.P_dirty = 1;
}

void nbody_create_potential(struct nbody *nbody, int N)
{
    assert(N == 2);

    nbody->phi.N = N;
    nbody->phi.eps2 = 0; //.01;
    nbody->phi.dt = 0.0001;
    nbody->phi.t = 0;

    nbody->phi.P = (struct massive_particle *)malloc(N * sizeof(*nbody->phi.P));

    double m  = 1e2;
    double x  = 10;
    double mu = m/2;
    double v  = sqrt(mu / fabs(2*x));

    //double d = 2*x;
    //v = sqrt(m * d*x/pow(d*d + nbody->phi.eps2, 1.5));


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

    cudaMalloc((void **)&nbody->phi.Pdev,     N * sizeof(*(nbody->phi.Pdev)));
    cudaMemcpy(nbody->phi.Pdev, nbody->phi.P, N * sizeof(*(nbody->phi.Pdev)), cudaMemcpyHostToDevice);
}

void nbody_advance_potential(struct nbody *nbody, tyme_t t)
{
    assert(t >= nbody->phi.t);

    if (t != nbody->phi.t)
    {
        fprintf(stderr, "Advancing potential to time %f\n", t);

        for (; nbody->phi.t < t-nbody->phi.dt; nbody->phi.t += nbody->phi.dt)
            nbody_step_potential(nbody, nbody->phi.dt);

        if (t > nbody->phi.t)
            nbody_step_potential(nbody, t - nbody->phi.dt);
    }
}


