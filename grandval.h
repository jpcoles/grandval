#ifndef GRANDVAL_H
#define GRANDVAL_H

#include <inttypes.h>

#define GRANDVAL_FULL_PROGRAM_NAME "grandval v0.1"


#define WITH_INTEGERS 0
#define WITH_FLOATS 1

#define SET_OPTION(o,v) do { o = v; o ## _set = 1; } while (0)
#define OPTSTR(s) (!strcmp(s, long_options[option_index].name))
#define eprintf(...) do { fprintf (stderr, __VA_ARGS__); } while (0)

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
    char *name;
    char *format;
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

struct program_options
{
    size_t Nparticles;      int Nparticles_set;
    tyme_t dt;              int dt_set;
    tyme_t Tmax;            int Tmax_set;
    int Nimages;            int Nimages_set;
    dist_t Rimages;         int Rimages_set;
    dist_t R;               int R_set;
    char *ic_name;          int ic_name_set;
    char *potential_name;   int potential_name_set;
    int cuda_device;        int cuda_device_set;

    size_t Nsnapshots;      int Nsnapshots_set;
    char *snapshot_name;    int snapshot_name_set;
    char *snapshot_format;  int snapshot_format_set;

    char *image_name;       int image_name_set;
    char *image_format;     int image_format_set;
};


#endif
