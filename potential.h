#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "grandval.h"

struct potential
{
    void (*set_particles)(void *phi_data, struct particle *p, size_t N);
    void (*get_particles)(void *phi_data, struct particle **p, size_t *N);
    void (*step_particles)(void *phi_data, tyme_t dt);
    void (*advance)(void *phi_data, tyme_t t);
    void (*free)(void *phi_data);
};

#endif
