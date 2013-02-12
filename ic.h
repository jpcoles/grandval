#ifndef IC_H
#define IC_H

#include "grandval.h"

typedef void (*ic_function)(struct particle *p, size_t N, pos_t R);

struct iclist
{
    ic_function f;
    char *name;
    char *desc;
};

void show_initial_conditions();
int find_ic(char *name, ic_function *f);
void ic_random(struct particle *p, size_t N, pos_t R);
void ic_line(struct particle *p, size_t N, pos_t R);
void ic_circular_plummer(struct particle *p, size_t N, pos_t R);
void ic_circular_hernquist(struct particle *p, size_t N, pos_t R);
void ic_droplet(struct particle *p, size_t N, pos_t R);
void ic_psdroplet(struct particle *p, size_t N, pos_t R);
void ic_pscube(struct particle *p, size_t N, pos_t R);
void ic_disk(struct particle *p, size_t N, pos_t R);

extern struct iclist ics[];

#endif

