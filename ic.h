#ifndef IC_H
#define IC_H

#include "grandval.h"

typedef void (*ic_function)(struct particle *p, int N, pos_t R);

struct iclist
{
    ic_function f;
    char *name;
    char *desc;
};

void show_initial_conditions();
int find_ic(char *name, ic_function *f);
void ic_random(struct particle *p, int N, pos_t R);
void ic_line(struct particle *p, int N, pos_t R);
void ic_droplet(struct particle *p, int N, pos_t R);

extern struct iclist ics[];

#endif

