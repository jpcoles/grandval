#ifndef hernquist_H
#define hernquist_H

#include "grandval.h"
#include "potential.h"

struct hernquist_potential
{
    dist_t eps;
    mass_t M;
};

struct hernquist_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct hernquist
{
    struct hernquist_potential phi;
    struct hernquist_particles P;
};

void hernquist_init(struct potential *phi);
void hernquist_free(void *phi_data);
void hernquist_set_particles(void *phi_data, struct particle *P, size_t N);
real hernquist_energy(void *phi_data, struct particle *P, size_t N);
int hernquist_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void hernquist_accel(struct particle *p, mass_t M, dist_t eps, acc_t *out);
int hernquist_step_particles(void *phi_data, tyme_t dt);
void *hernquist_create_potential(size_t N);
int hernquist_advance_potential(void *phi_data, tyme_t t);

#endif
