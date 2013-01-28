#ifndef NBODY_H
#define NBODY_H

#include "grandval.h"

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

struct nbody *nbody_init();
void nbody_free(struct nbody *nbody);
void nbody_set_particles(struct nbody *nbody, struct particle *P, size_t N);
void nbody_get_particles(struct nbody *nbody, struct particle **P, size_t *N);
__device__ void nbody_accel(struct particle *p, struct particle *Pm, size_t NPm, const dist_t eps2, acc_t *out);
void nbody_step_all(struct nbody *nbody, tyme_t dt);
void nbody_create_potential(struct nbody *nbody, int N);
void nbody_advance_potential(struct nbody *nbody, tyme_t t);

#endif
