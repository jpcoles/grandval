#ifndef plummeri_lumps_H
#define plummer_lumps_H

#include "grandval.h"
#include "potential.h"

struct plummer_lumps_potential
{
    dist_t eps2;
    mass_t M;
};

struct plummer_lumps_particles
{
    size_t N;
    struct particle *P;
    struct particle *Pdev;
    int P_dirty;
    int Pdev_dirty;
};

struct plummer_lumps
{
    struct plummer_lumps_potential phi;
    struct plummer_lumps_particles P;
};

void plummer_lumps_init(struct potential *phi);
void plummer_lumps_free(void *phi_data);
void plummer_lumps_set_particles(void *phi_data, struct particle *P, size_t N);
real plummer_lumps_energy(void *phi_data, struct particle *P, size_t N);
int plummer_lumps_get_particles(void *phi_data, struct particle **P, size_t *N);
__device__ void plummer_lumps_accel(struct particle *p, mass_t M, dist_t eps2, acc_t *out);
int plummer_lumps_step_particles(void *phi_data, tyme_t dt);
void *plummer_lumps_create_potential(size_t N);
int plummer_lumps_advance_potential(void *phi_data, tyme_t t);

#endif
