#ifndef IC_H
#define IC_H

#include "grandval.h"

struct iclist
{
    void (*f)(struct particle *p, int N, pos_t R);
    char *name;
    char *desc;
};

void ic_random(struct particle *p, int N, pos_t R);
void ic_circular(struct particle *p, int N, pos_t R);

#endif

