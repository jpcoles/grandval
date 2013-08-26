/* 
 * Implements a plummer potential + lumps.
 * 
 */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "plummer_lumps.h"
#include "util.h"

static const size_t cuda_threads_per_block = 512;

void plummer_lumps_init(struct potential *phi)
{
    phi->name = "plummer_lumps";
    phi->desc = "Plummer potential.";

    phi->create         = plummer_lumps_create_potential;
    phi->set_particles  = plummer_lumps_set_particles;
    phi->get_particles  = plummer_lumps_get_particles;
    phi->step_particles = plummer_lumps_step_particles;
    phi->advance        = plummer_lumps_advance_potential;
    phi->free           = plummer_lumps_free;
    phi->energy         = plummer_lumps_energy;
}

void plummer_lumps_free(void *phi_data)
{
    if (phi_data == NULL) return;

    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)phi_data;

    if (plummer_lumps->P.P != NULL)
    {
        free(plummer_lumps->P.P);
        plummer_lumps->P.P = NULL;
    }

    if (plummer_lumps->P.Pdev != NULL)
    {
        cudaFree(plummer_lumps->P.Pdev);
        plummer_lumps->P.Pdev = NULL;
    }
}

void plummer_lumps_set_particles(void *phi_data, struct particle *P, size_t N)
{
    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)phi_data;
    plummer_lumps->P.N = N;
    plummer_lumps->P.P = P;
    plummer_lumps->P.P_dirty = 1;
    cudaMalloc((void **)&plummer_lumps->P.Pdev, N * sizeof(*(plummer_lumps->P.Pdev)));
}

int plummer_lumps_get_particles(void *phi_data, struct particle **P, size_t *N)
{
    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)phi_data;

    if (plummer_lumps->P.Pdev_dirty)
    {
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA synchronize failed %i: %s", err, cudaGetErrorString(err));
            return 0;
        }
        //fprintf(stderr, "Updating host particles.\n");
        cudaMemcpy(plummer_lumps->P.P, plummer_lumps->P.Pdev, plummer_lumps->P.N * sizeof(*plummer_lumps->P.P), cudaMemcpyDeviceToHost);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy from device failed %i: %s", err, cudaGetErrorString(err));
            return 0;
        }
    }

    *P = plummer_lumps->P.P;
    *N = plummer_lumps->P.N;

    return 1;
}
//next function compute the acceleration due exclusively to lumps
__host__ __device__ real *lump_accel(struct particle *p, double lump_coord[3], mass_t m, dist_t eps2)
{   
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2]; 
    
    const dist_t x_0 = lump_coord[0];//lump coordinates
    const dist_t y_0 = lump_coord[1];
    const dist_t z_0 = lump_coord[2]; 

    const dist_t r2_0= pow(x-x_0, 2) + pow(y-y_0, 2) + pow(z-z_0, 2);
    
    real *out;
    out = new real(3);
    out[0] = - m * (x-x_0) / pow(r2_0 + eps2, (real)1.5);
    out[1] = - m * (y-y_0) / pow(r2_0 + eps2, (real)1.5);
    out[2] = - m * (y-y_0) / pow(r2_0 + eps2, (real)1.5);
    return out;
}

//next function computes the overall acceleration, both plummer+lumps.
//lumps contribution is considerated calling the previous function.
//the parts commented out are a previous rough way to introduce the effect of one and 2 lumps.
//One thing to be fixed for sure is the fact you have to set lump coordinate here in the code.

__device__ void plummer_lumps_accel(struct particle *p, mass_t M, dist_t eps2, acc_t *out)
{
    const dist_t x = p->x[0];
    const dist_t y = p->x[1];
    const dist_t z = p->x[2];

    const dist_t r2 = x*x + y*y + z*z;
/*
    const dist_t x_0=0.5;//lump coordinates
    const dist_t y_0=0.0;
    const dist_t z_0=0.0; 

    const dist_t r2_0= pow(x-x_0, 2) + pow(y-y_0, 2) + pow(z-z_0, 2);

    const dist_t x_1=-0.5;//lump coordinates
    const dist_t y_1=0.0;
    const dist_t z_1=0.0; 

    const dist_t r2_1= pow(x-x_1, 2) + pow(y-y_1, 2) + pow(z-z_1, 2);
*/    
    const mass_t m = 1; // lump mass
    double lump[3];
    lump[0] = 0.5;
    lump[1] = 0;
    lump[2] = 0;
    real *ac_l;
    ac_l = new real(3);
    ac_l = lump_accel(p, lump, m, eps2);    

    out[0] = -M * x / pow(r2 + eps2, (real)1.5) + ac_l[0];//- m * (x-x_0) / pow(r2_0 + eps2, (real)1.5);// - m* (x-x_1) / pow(r2_1 + eps2, (real)1.5);
    out[1] = -M * y / pow(r2 + eps2, (real)1.5) + ac_l[1];//- m * (y-y_0) / pow(r2_0 + eps2, (real)1.5);// - m* (y-y_1) / pow(r2_1 + eps2, (real)1.5);
    out[2] = -M * z / pow(r2 + eps2, (real)1.5) + ac_l[2];//- m * (y-y_0) / pow(r2_0 + eps2, (real)1.5);// - m* (z-z_1) / pow(r2_1 + eps2, (real)1.5);
}

__global__ void plummer_lumps_cuda_step_all(struct particle *P, size_t NP, mass_t M, dist_t eps2, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    plummer_lumps_accel(p, M, eps2, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

__global__ void plummer_lumps_cuda_step_all_first_time(struct particle *P, size_t NP, mass_t M, dist_t eps2, tyme_t dt)
{
    int i;
    acc_t a[3];

    size_t pi = blockIdx.x * blockDim.x + threadIdx.x;

    if (pi >= NP) return;

    struct particle *p = P + pi;
    plummer_lumps_accel(p, M, eps2, a);
    for (i=0; i < 3; i++) p->x[i] += a[i]*dt*dt/4;
}

int plummer_lumps_step_particles(void *phi_data, tyme_t dt)
{
    cudaError_t err;
    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)phi_data;

    static int first_time = 1;
    int nblocks  = (plummer_lumps->P.N + cuda_threads_per_block - 1) / cuda_threads_per_block;
    int nthreads = plummer_lumps->P.N < cuda_threads_per_block ? plummer_lumps->P.N : cuda_threads_per_block;

    if (first_time)
    {
        fprintf(stderr, "Using %i blocks, %i threads\n", nblocks, nthreads);
        first_time = 0;
    }

    /* If the local particles are dirty, copy them to the device */
    if (plummer_lumps->P.P_dirty)
    {
        //fprintf(stderr, "Updating device particles.\n");
        cudaMemcpy(plummer_lumps->P.Pdev, plummer_lumps->P.P, plummer_lumps->P.N * sizeof(*(plummer_lumps->P.Pdev)), cudaMemcpyHostToDevice);
        err = cudaGetLastError();
        if (err)
        {
            errmsg("CUDA memory copy to device failed (%i): %s", err, cudaGetErrorString(err));
            return 0;
        }

        plummer_lumps->P.P_dirty = 0;
        plummer_lumps->P.Pdev_dirty = 0;
    }

    if (first_time)
        plummer_lumps_cuda_step_all_first_time<<<nblocks,nthreads>>>(plummer_lumps->P.Pdev, plummer_lumps->P.N, plummer_lumps->phi.M, plummer_lumps->phi.eps2, dt);

    plummer_lumps_cuda_step_all<<<nblocks,nthreads>>>(plummer_lumps->P.Pdev, plummer_lumps->P.N, plummer_lumps->phi.M, plummer_lumps->phi.eps2, dt);
    err = cudaGetLastError();
    if (err)
    {
        errmsg("CUDA kernel call had error status %i: %s", err, cudaGetErrorString(err));
        return 0;
    }

    plummer_lumps->P.Pdev_dirty = 1;

    return 1;
}

void *plummer_lumps_create_potential(size_t N)
{
    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)calloc(1, sizeof(*plummer_lumps));
    assert(plummer_lumps != NULL);

    plummer_lumps->phi.eps2 = 0.05;
    plummer_lumps->phi.M = 2;

    return plummer_lumps;
}

int plummer_lumps_advance_potential(void *phi_data, tyme_t t)
{
    return 1;
}

//don't look at this it is wrong
real plummer_lumps_energy(void *phi_data, struct particle *P, size_t N)
{
    size_t i,j;
    real T=0,U=0;

    struct plummer_lumps *plummer_lumps = (struct plummer_lumps *)phi_data;
    const dist_t eps2 = plummer_lumps->phi.eps2;
    const mass_t M = plummer_lumps->phi.M;
    
    const dist_t x_0=0.5;//lump coordinates
    const dist_t y_0=0.0;
    const dist_t z_0=0.0;

    const dist_t x_1=-0.5;//lump coordinates
    const dist_t y_1=0.0;
    const dist_t z_1=0.0; 
     
    const mass_t m = 1; // lump mass

    for (i=0; i < N; i++)
    {
        vel_t v2 = 0;
        for (j=0; j < 3; j++)
            v2 += pow(P[i].v[j], 2);
        
        T += 0.5*v2;

        dist_t r2=0;
        dist_t r2_0=0;
        dist_t r2_1=0;

        for (j=0; j < 3; j++)
            r2 += pow(P[i].x[j], 2);
        
        r2_0= pow(P[i].x[0]-x_0, 2) + pow(P[i].x[1]-y_0, 2) + pow(P[i].x[2]-z_0, 2);
        r2_1= pow(P[i].x[0]-x_1, 2) + pow(P[i].x[1]-y_1, 2) + pow(P[i].x[2]-z_1, 2);
        

        U -= M / sqrt(r2 + eps2);
        U -= m / sqrt(r2_0 + eps2) + m / sqrt(r2_1 + eps2);
        P[i].energy_pp=0.5*v2-M/sqrt(r2+eps2)-m/sqrt(r2_0+eps2);//-m/sqrt(r2_1+eps2);
        P[i].id = i;
    }


    return T+U;
}
