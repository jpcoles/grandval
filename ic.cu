#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "grandval.h"
#include "rand.h"
#include "ic.h"

struct iclist ics[] = 
{ {ic_random,              "random",  ""},
  {ic_line,                "line",  ""},
  {ic_circular_plummer,    "circular_plummer",  ""},
  {ic_circular_hernquist,  "circular_hernquist",""},  
  {ic_circular_jaffe,      "circular_jaffe",  ""},
  {ic_droplet,             "droplet",  ""},
  {ic_psdroplet,           "psdroplet",  ""},
  {ic_disk,                "disk",  ""},
  {ic_pscube,              "pscube",  ""},
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
void ic_random(struct particle *p, size_t N, pos_t R)
{
    size_t i;
    for (i=0; i < N; i++)
    {
        p[i].x[0] = (pos_t)(R * (2*randU()-1));
        p[i].x[1] = (pos_t)(R * (2*randU()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = 0;
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }
}

//==============================================================================
//                                 ic_circular_plummer
//==============================================================================
void ic_circular_plummer(struct particle *p, size_t N, pos_t R)
{
    mass_t M = 2;
    dist_t eps2 = 0.05;
    //dist_t Rmin = 2*sqrt(eps2);
    dist_t Rmin = 0.1*R;
    size_t i;
    for (i=0; i < N; i++)
    {
        dist_t r = (pos_t)((R-Rmin) * (2*randU()-1));
        r += Rmin * (2*(r>0)-1);
        p[i].x[0] = 0; //(pos_t)(R * (2*randU()-1));
        p[i].x[1] = r;
        //p[i].x[1] = R;
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = sqrt(M * r*r / pow(r*r + eps2, 1.5F)) * (2*(r > 0)-1);
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }
   
} 
//=============================================================================
//                                 ic_circular_hernquist
//============================================================================
void ic_circular_hernquist (struct particle *p, size_t N, pos_t R)
 { 
   mass_t M = 5;
   dist_t eps = 0.22;
   dist_t Rmin = 2*eps;   
   size_t i;
   for (i=0; i < N; i++)
   {  
       dist_t r = (pos_t)((R-Rmin) * (2*randU()-1));
       r += Rmin * (2*(r>0)-1);
       p[i].x[0] = 0; 
       p[i].x[1] = r;
       p[i].x[2] = 0;
       p[i].v[0] = sqrt(M * abs(r) / pow(abs(r)+eps, 2.0F)) * (2*(r > 0)-1);
       p[i].v[1] = 0;
       p[i].v[2] = 0;
    }  
 }

//==============================================================================
//                                 ic_circular_jaffe
//==============================================================================
void ic_circular_jaffe(struct particle *p, size_t N, pos_t R)
{
    mass_t M = 5;
    dist_t eps = 0.22;
    dist_t Rmin = 2*eps;
    size_t i;
    for (i=0; i < N; i++)
    {
        dist_t r = (pos_t)((R-Rmin) * (2*randU()-1));
        r += Rmin * (2*(r>0)-1);
        p[i].x[0] = 0; //(pos_t)(R * (2*randU()-1));
        p[i].x[1] = r;
        //p[i].x[1] = R;
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = sqrt(M / (abs(r) + eps)) * (2*(r > 0)-1);
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }

}
//============================================================================
//                                  ic_line
//============================================================================

//next ic code was used to calculate phase difference between 3-velcoity initial lines
/*
void ic_line(struct particle *p, size_t N, pos_t R)
{
    size_t i;
    dist_t Rmin = 0.1 * R;
    for (i=0; i < N; i+=3)
    {
        p[i].x[0] = (pos_t)(R* (2.0*i/N -1));//(pos_t)(R*(2*randU()-1));//(pos_t)(R * ((0+2*i/N)-1));
        p[i].x[1] = 0;//(pos_t)(R * (2*randU()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = 0; //2*randU()-1;
        p[i].v[1] = 0; //sqrt(2*1e1 / fabs(p[i].x[0]));
        p[i].v[2] = 0;
        
        p[i].energy_pp=0.5*pow(p[i].v[0], 2)-2.0/sqrt(pow(p[i].x[0], 2)+0.05);
    }
    
    for (i=1; i<N; i+=3)
    {
        p[i].x[1] = 0;//(pos_t)(R * (2*randU()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = 1; //2*randU()-1;
        p[i].v[1] = 0; //sqrt(2*1e1 / fabs(p[i].x[0]));
        p[i].v[2] = 0;
        
        if (p[i-1].x[0] >= 0)
        {p[i].x[0] = (pos_t)(sqrt(pow((2/(0.5*p[i].v[0]*p[i].v[0]-p[i-1].energy_pp)),2)-0.05));}
        else
        {p[i].x[0] = (pos_t)(-1*sqrt(pow((2/(0.5*p[i].v[0]*p[i].v[0]-p[i-1].energy_pp)),2)-0.05));}


    }
 for (i=2; i<N; i+=3)
    {
        p[i].x[1] = 0;//(pos_t)(R * (2*randU()-1));
        p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
        p[i].v[0] = -1; //2*randU()-1;
        p[i].v[1] = 0; //sqrt(2*1e1 / fabs(p[i].x[0]));
        p[i].v[2] = 0;
        
        if (p[i-2].x[0] >= 0)
        {p[i].x[0] = (pos_t)(sqrt(pow((2/(0.5*p[i].v[0]*p[i].v[0]-p[i-2].energy_pp)),2)-0.05));}
        else
        {p[i].x[0] = (pos_t)(-1*sqrt(pow((2/(0.5*p[i].v[0]*p[i].v[0]-p[i-2].energy_pp)),2)-0.05));}

    }


}

*/

void ic_line(struct particle *p, size_t N, pos_t R)
{
    size_t i;
    dist_t Rmin = 0.1 * R;
    for (i=0; i < N; i++)
    {
       p[i].x[0] = (pos_t)(R * (2*randU()-1));
       p[i].x[1] = 0;//(pos_t)(R * (2*randU()-1));
       p[i].x[2] = 0; //(pos_t)(env->radius * (2*randU()-1));
       p[i].v[0] = 0; //2*randU()-1;
       p[i].v[1] = 0; //sqrt(2*1e1 / fabs(p[i].x[0]));
       p[i].v[2] = 0;
    }
}


//==============================================================================
//                                 ic_droplet
//==============================================================================
void ic_droplet(struct particle *p, size_t N, pos_t R)
{
    size_t i;
    dist_t rmax = 1.01;

    p[0].x[0] = 0;
    p[0].x[1] = R;
    p[0].x[2] = 0;
    p[0].v[0] = 0;
    p[0].v[1] = 0;
    p[0].v[2] = 0;

    for (i=1; i < N; i++)
    {
        pos_t z = 2.0 * randU() - 1.0;
        pos_t t = 2.0 * M_PI * randU();
        pos_t x = sqrt(1-(z*z)) * cos(t);
        pos_t y = sqrt(1-(z*z)) * sin(t);
        dist_t r = rmax * sqrt(randU());

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

//==============================================================================
//                                 ic_psdroplet
//==============================================================================
void ic_psdroplet(struct particle *p, size_t N, pos_t R)
{
    int j;
    size_t i;
    dist_t rmax = 1e-6;

    const int D = 1;

    srand48(0);

    //p[0].x[0] = R;
    //p[0].x[1] = 0;
    //p[0].x[2] = 0;
    //p[0].v[0] = 0;
    //p[0].v[1] = 0;
    //p[0].v[2] = 0;

    for (i=0; i < N; i++)
    {
        double x[3];
        double v[3];

        memset(x, 0, sizeof(*x) * 3);
        memset(v, 0, sizeof(*v) * 3);

        for (j=0; j < D; j++)
            normal(&x[j],&v[j]);

        dist_t r2=0;
        for (j=0; j < D; j++)
        {
            r2 += pow(x[j], 2);
            r2 += pow(v[j], 2);
        }

        const dist_t r = sqrt(r2);
        const dist_t radius = rmax * pow(1-randU(), 1./(2*D));

        for (j=0; j < D; j++)
        {
            x[j] *= radius / r;
            v[j] *= radius / r;
        }

        p[i].x[0] = x[0] + R;
        p[i].x[1] = x[1];
        p[i].x[2] = x[2];
        p[i].v[0] = v[0];
        p[i].v[1] = v[1];
        p[i].v[2] = v[2];
    }
}

//==============================================================================
//                                 ic_droplet
//==============================================================================
void ic_disk(struct particle *p, size_t N, pos_t R)
{
    dist_t rmax = 1.01;

    size_t i;
    for (i=0; i < N; i++)
    {
        pos_t t = 2.0 * M_PI * randU();
        pos_t x = cos(t);
        pos_t y = sin(t);
        dist_t r = rmax * sqrt(randU());

        x = 0 + x*r;
        y = R + y*r;

        p[i].x[0] = x;
        p[i].x[1] = y;
        p[i].x[2] = 0;
        p[i].v[0] = 0;
        p[i].v[1] = 0;
        p[i].v[2] = 0;
    }
}

void ic_pscube(struct particle *p, size_t N, pos_t R)
{
    size_t i;
    assert(N == 8);

    struct
    {
        pos_t x,y,z;
        vel_t vx,vy,vz;
    } q[8] = { {-1, -1, -1, 0,0,0},
               {-1, -1, +1, 0,0,0},
               {-1, +1, -1, 0,0,0},
               {-1, +1, +1, 0,0,0},
               {+1, -1, -1, 0,0,0},
               {+1, -1, +1, 0,0,0},
               {+1, +1, -1, 0,0,0},
               {+1, +1, +1, 0,0,0} };

    for (i=0; i < N; i++)
    {
        p[i].x[0] = q[i].x;
        p[i].x[1] = q[i].y;
        p[i].x[2] = q[i].z;
        p[i].v[0] = q[i].x * 0.1;
        p[i].v[1] = q[i].y * 0.1;
        p[i].v[2] = q[i].z * 0.1;
    }
}

