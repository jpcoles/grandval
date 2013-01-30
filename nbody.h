#ifndef NBODY_H
#define NBODY_H

#include "grandval.h"
#include "potential.h"

struct nbody_potential
{
    size_t N;
    struct massive_particle *P;
    struct massive_particle *Pdev;
    dist_t eps2;
    tyme_t t;
    tyme_t dt;
    int P_dirty;
    int Pdev_dirty;
};

struct nbody_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct nbody
{
    struct nbody_potential phi;
    struct nbody_particles P;
};

void nbody_init(struct potential *phi);
void nbody_free(void *phi_data);
void nbody_set_particles(void *phi_data, struct particle *P, size_t N);
int nbody_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void nbody_accel(struct particle *p, struct particle *Pm, size_t NPm, const dist_t eps2, acc_t *out);
int nbody_step_particles(void *phi_data, tyme_t dt);
void *nbody_create_potential(size_t N);
int nbody_advance_potential(void *phi_data, tyme_t t);

#endif
