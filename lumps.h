#ifndef LUMPS_H
#define LUMPS_H

#include "grandval.h"
#include "potential.h"

struct lumps_potential
{
    size_t N;
    struct massive_particle *P;
    struct massive_particle *Pdev;
    dist_t eps2;
    tyme_t t;
    tyme_t dt;
    int P_dirty;
    int Pdev_dirty;
    mass_t M;
    dist_t eps2_g;
};

struct lumps_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct lumps
{
    struct lumps_potential phi;
    struct lumps_particles P;
};

void lumps_init(struct potential *phi);
void lumps_free(void *phi_data);
void lumps_set_particles(void *phi_data, struct particle *P, size_t N);
int lumps_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void lumps_accel(struct particle *p, struct particle *Pm, size_t NPm, const dist_t eps2, const dist_t eps2_g, const mass_t M,  acc_t *out);
int lumps_step_particles(void *phi_data, tyme_t dt);
void *lumps_create_potential(size_t N);
int lumps_advance_potential(void *phi_data, tyme_t t);

#endif
