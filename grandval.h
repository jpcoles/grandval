#ifndef GRANDVAL_H
#define GRANDVAL_H

#include <inttypes.h>

#define WITH_INTEGERS 0
#define WITH_FLOATS 1

#if WITH_INTEGERS
typedef int64_t tyme_t;
typedef int64_t pos_t;
typedef int64_t dist_t;
typedef int64_t vel_t;
typedef int64_t acc_t;
typedef double  force_t;
#define TIMET  I_TIMET 
#define POST   I_POST  
#define DISTT  I_DISTT  
#define VELT   I_VELT  
#define ACCT   I_ACCT  
#define RHOT   I_RHOT  
#define MASST  I_MASST 
#define ENGYT  I_ENGYT 
#define FORCET I_FORCET
#define SOFTT  I_SOFTT
#endif

#if WITH_FLOATS
typedef double pos_t;
typedef double dist_t;
typedef double vel_t;
typedef double acc_t;
typedef double tyme_t;
typedef double force_t;
#define TIMET  F_TIMET 
#define DISTT  F_DISTT 
#define POST   F_POST  
#define VELT   F_VELT  
#define ACCT   F_ACCT  
#define RHOT   F_RHOT  
#define MASST  F_MASST 
#define ENGYT  F_ENGYT 
#define FORCET F_FORCET
#define SOFTT  F_SOFTT
#endif

typedef double mass_t;             /* Msun             */
typedef double energy_t;           /* */
typedef double soft_t;

struct particle
{
    pos_t x[3];
    vel_t v[3];
};

struct massive_particle
{
    pos_t x[3];
    vel_t v[3];
    mass_t m;
};

struct potential;

typedef void (*accel_fn)(struct potential *phi, struct particle *p, acc_t *out);
typedef void (*advance_fn)(struct potential *phi, tyme_t t);

struct potential
{
    accel_fn accel;
    advance_fn advance;
    void *phi;
};

struct nbody_potential
{
    int N;
    double eps2;
    tyme_t t;
    tyme_t dt;
    struct massive_particle *p;
};

#endif
