#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "grandval.h"

struct potential
{
    char *name;
    char *desc;
    void * (*create)(size_t N);
    void (*set_particles)(void *phi_data, struct particle *p, size_t N);
    int (*get_particles)(void *phi_data, struct particle **p, size_t *N);
    int (*step_particles)(void *phi_data, tyme_t dt);
    int (*advance)(void *phi_data, tyme_t t);
    void (*free)(void *phi_data);
};

void show_potentials();
int find_potential(char *name, struct potential *p);
void add_potential(struct potential *p);


#endif
