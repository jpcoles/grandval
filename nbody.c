#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "grandval.h"
#include "nbody.h"

void accel(struct potential *phi0, struct particle *p, acc_t *out)
{
    int i;
    struct nbody_potential *phi = phi0->phi;
    struct massive_particle *mp = phi->p;
    const double e2 = phi->eps2;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

    for (i=0; i < phi->N; i++)
    {
        const double dx = mp[i].x[0] - x;
        const double dy = mp[i].x[1] - y;
        const double dz = mp[i].x[2] - z;
        const double r2 = dx*dx + dy*dy + dz*dz;

        double rinv = 1.0 / sqrt(r2);

        double re = sqrt(r2) * sqrt(e2);

        #define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
        #define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
        double gradphi = GRADPHI(re);

        // round'ing is important. Conserves momentum perfectly.
        //#define Fhat(x) round(((mp[i].m*e2)*(((x)*rinv)*gradphi)))
#if WITH_INTEGERS
        #define Fhat(x) round((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#else
        #define Fhat(x) ((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#endif
        out[0] += Fhat(dx);
        out[1] += Fhat(dy);
        out[2] += Fhat(dz);

        int j;
        for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    fprintf(stderr, "\n");
}

void nbody_accel(struct nbody_potential *phi, struct massive_particle *p, acc_t *out)
{
    int i;
    struct massive_particle *mp = phi->p;
    const double e2 = phi->eps2;

    for (i=0; i < 3; i++)
        out[i] = 0;

    const pos_t x = p->x[0];
    const pos_t y = p->x[1];
    const pos_t z = p->x[2];

    for (i=0; i < phi->N; i++)
    {
        const double dx = mp[i].x[0] - x;
        const double dy = mp[i].x[1] - y;
        const double dz = mp[i].x[2] - z;
        const double r2 = dx*dx + dy*dy + dz*dz;

        double rinv = 1.0 / sqrt(r2);

        double re = sqrt(r2) * sqrt(e2);

        #define PHI(r)     ( 1./sqrt(1.+pow((r),2)))
        #define GRADPHI(r) (pow(sqrt(1.+pow((r),2)), -3) * (r))
        double gradphi = GRADPHI(re);

        // round'ing is important. Conserves momentum perfectly.
        //#define Fhat(x) round(((mp[i].m*e2)*(((x)*rinv)*gradphi)))
#if WITH_INTEGERS
        #define Fhat(x) round((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#else
        #define Fhat(x) ((mp[i].m*(x) * pow(e2 + r2, -3./2.)))
#endif
        out[0] += Fhat(dx);
        out[1] += Fhat(dy);
        out[2] += Fhat(dz);

        int j;
        for (j=0; j < 3; j++) fprintf(stderr, "%i %e %e %e\n", j, Fhat(dx), Fhat(dy), Fhat(dz));
    }
    fprintf(stderr, "\n");
}

void nbody_step(struct potential *phi, struct particle *p, tyme_t t, tyme_t dt)
{
    int i;
    acc_t a[3];

    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
    phi->accel(phi, p, a);
    for (i=0; i < 3; i++) p->v[i] +=    a[i] * dt;
    for (i=0; i < 3; i++) p->x[i] += p->v[i] * dt/2;
}

void nbody_create_potential(struct potential *phi0, int N)
{
    struct nbody_potential *phi = malloc(sizeof(*phi)); 

    phi->N = N;
    phi->eps2 = 1;
    phi->p = malloc(N * sizeof(*phi->p));

    phi->p[0].x[0] =
    phi->p[0].x[1] =
    phi->p[0].x[2] = 0;

    phi->p[0].v[0] =
    phi->p[0].v[1] =
    phi->p[0].v[2] = 0;

    phi->p[0].m = 1e5;

    phi0->accel = accel;
    phi0->advance = nbody_advance_phi;
    phi0->phi = phi;
}

void nbody_advance_phi(struct potential *phi0, tyme_t t)
{
    struct nbody_potential *phi = phi0->phi;
    assert(t >= phi->t);

    if (t != phi->t)
    {
        for (; phi->t < t-phi->dt; phi->t += phi->dt)
        {
            for (i=0; i < phi->N; i++)
                nbody_step(&phi, &phi->p[i], phi->dt);
        }

        for (i=0; i < phi->N; i++)
            step(&phi, &phi->p[i], t-phi->t);
    }
}


