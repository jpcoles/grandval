#ifndef GRANDVAL_H
#define GRANDVAL_H

#include <inttypes.h>

#define WITH_INTEGERS 0
#define WITH_FLOATS 1

typedef float real;

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
typedef real pos_t;
typedef real dist_t;
typedef real vel_t;
typedef real acc_t;
typedef real tyme_t;
typedef real force_t;
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

typedef real mass_t;             /* Msun             */
typedef real energy_t;           /* */
typedef real soft_t;

struct image
{
    int nr, nc;
    unsigned char *image;
    int *hist;
};

struct __align__(16) particle
{
    pos_t x[3];
    vel_t v[3];
};

struct __align__(16) massive_particle
{
    pos_t x[3];
    vel_t v[3];
    mass_t m;
};


#endif
