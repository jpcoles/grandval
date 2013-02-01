#ifndef BAR_H
#define BAR_H

#include "grandval.h"
#include "potential.h"

struct bar_potential
{
    tyme_t t;
    freq_t omega;
    mass_t M;
};

struct bar_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct bar
{
    struct bar_potential phi;
    struct bar_particles P;
};

void bar_init(struct potential *phi);
void bar_free(void *phi_data);
void bar_set_particles(void *phi_data, struct particle *P, size_t N);
int bar_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void bar_accel(struct particle *p, tyme_t t, mass_t M, freq_t omega, acc_t *out);
int bar_step_particles(void *phi_data, tyme_t dt);
void *bar_create_potential(size_t N);
int bar_advance_potential(void *phi_data, tyme_t t);

#endif
