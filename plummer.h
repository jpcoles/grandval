#ifndef plummer_H
#define plummer_H

#include "grandval.h"
#include "potential.h"

struct plummer_potential
{
    dist_t eps2;
    mass_t M;
};

struct plummer_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct plummer
{
    struct plummer_potential phi;
    struct plummer_particles P;
};

void plummer_init(struct potential *phi);
void plummer_free(void *phi_data);
void plummer_set_particles(void *phi_data, struct particle *P, size_t N);
real plummer_energy(void *phi_data, struct particle *P, size_t N);
int plummer_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void plummer_accel(struct particle *p, mass_t M, dist_t eps2, acc_t *out);
int plummer_step_particles(void *phi_data, tyme_t dt);
void *plummer_create_potential(size_t N);
int plummer_advance_potential(void *phi_data, tyme_t t);

#endif
