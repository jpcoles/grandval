#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "grandval.h"
#include "nbody.h"
#include "ic.h"


void step(struct potential *phi, struct particle *p, tyme_t dt)
{
    int i;
    acc_t a[3];

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    phi->accel(phi, p, a);
    //int j;
    //for (j=0; j < 3; j++) fprintf(stderr, "%i %ld\n", j, a[j]);

    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

void write_positions(struct particle *P, int NP, int first_output)
{
    int i;
    char *mode;
    char fname[256];

    if (first_output)
        mode = "w";
    else
        mode = "a+";

    for (i=0; i < NP; i++)
    {
        sprintf(fname, "gv-%05i.out", i);
        FILE *fp = fopen(fname, mode);
#if WITH_INTEGERS
        fprintf(fp, "%ld %ld %ld\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#else
        fprintf(fp, "%f %f %f\n", P[i].x[0], P[i].x[1], P[i].x[2]);
#endif
        fclose(fp);
    }
}

int main(int argc, char **argv)
{
    int i;
    int NP = 1;

    struct particle *P = malloc(NP * sizeof(*P));
    assert(P != NULL);

    tyme_t Tmax = 100;
    tyme_t t = 0;
    tyme_t dt = .02;

    struct potential phi;

    ic_random(P, NP, 100);

    nbody_create_potential(&phi, 1);

    int curr_step=0;

    for (t=dt; t < Tmax+dt; t += dt, curr_step++)
    {
        write_positions(P, NP, curr_step == 0);

        if (t > Tmax) t = Tmax;

        phi.advance(&phi, t);

        for (i=0; i < NP; i++)
            step(&phi, &P[i], dt);

    }

    write_positions(P, NP, curr_step == 0);

    return 0;
}

