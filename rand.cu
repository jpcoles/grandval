#include <math.h>
#include "rand.h"

void rand_seed(long v)
{
    srand48(v);
}

void normal(double *x, double *y)
{
    double U,V;
   
    while ((U=drand48()) == 0);
    while ((V=drand48()) == 0);

    *x = sqrt(-2 * log(U)) * cos(2*M_PI*V);
    *y = sqrt(-2 * log(U)) * sin(2*M_PI*V);
}

double randU()
{
    return drand48();
}
