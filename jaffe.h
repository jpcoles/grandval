#ifndef jaffe_H
#define jaffe_H

#include "grandval.h"
#include "potential.h"

struct jaffe_potential
{
    dist_t eps;
    mass_t M;
};

struct jaffe_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct jaffe
{
    struct jaffe_potential phi;
    struct jaffe_particles P;
};

void jaffe_init(struct potential *phi);
void jaffe_free(void *phi_data);
void jaffe_set_particles(void *phi_data, struct particle *P, size_t N);
real jaffe_energy(void *phi_data, struct particle *P, size_t N);
int jaffe_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void jaffe_accel(struct particle *p, mass_t M, dist_t eps, acc_t *out);
int jaffe_step_particles(void *phi_data, tyme_t dt);
void *jaffe_create_potential(size_t N);
int jaffe_advance_potential(void *phi_data, tyme_t t);

#endif
