#include <stdio.h>
#include <stdlib.h>
#include "grandval.h"
#include "ic.h"

struct iclist ics[] = 
{ {ic_random,              "random",  ""},
  {ic_line,                "line",  ""},
  {ic_droplet,             "droplet",  ""},
  {NULL, NULL, NULL}
};

int find_ic(char *name, ic_function *f)
{
    int i;
    for (i=0; ics[i].name != NULL; i++)
    {
        if (!strcmp(name, ics[i].name))
        {
            *f = ics[i].f;
            return 1;
        }
    }

    return 0;
}

void show_initial_conditions()
{
    int i;
    for (i=0; ics[i].name != NULL; i++)
        eprintf("%15s - %s\n", ics[i].name, ics[i].desc);
}

//==============================================================================
//                                 ic_random
//==============================================================================
void ic_random(struct particle *p, int N, pos_t R)
{
    int i;
    for (i=0; i < N; i++)
    {
        p[i].x[0] = (pos_t)(R * (2*drand48()-1));
        p[i].x[1] = (pos_t)(R * (2*drand48()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*drand48()-1));
        p[i].v[0] = 0;
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }
}

//==============================================================================
//                                 ic_circular
//==============================================================================
void ic_line(struct particle *p, int N, pos_t R)
{
    int i;
    for (i=0; i < N; i++)
    {
        p[i].x[0] = 0; //(pos_t)(R * (2*drand48()-1));
        p[i].x[1] = (pos_t)(R * (2*drand48()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*drand48()-1));
        p[i].v[0] = 0;
        p[i].v[1] = 0; //sqrt(2*1e1 / fabs(p[i].x[0]));
        p[i].v[2] = 0;
    }
}

//==============================================================================
//                                 ic_droplet
//==============================================================================
void ic_droplet(struct particle *p, int N, pos_t R)
{
    dist_t r = 0.01;

    int i;
    for (i=0; i < N; i++)
    {
        pos_t z = 2.0 * drand48() - 1.0;
        pos_t t = 2.0 * M_PI * drand48();
        pos_t x = sqrt(1-(z*z)) * cos(t);
        pos_t y = sqrt(1-(z*z)) * sin(t);

        x = 0 + x*r;
        y = R + y*r;
        z = 0 + z*r;

        p[i].x[0] = x;
        p[i].x[1] = y;
        p[i].x[2] = z;
        p[i].v[0] = 0;
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }
}
