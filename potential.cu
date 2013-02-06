#include <stdio.h>
#include "grandval.h"
#include "potential.h"

static int Npotentials;
static struct potential *potentials; 

void show_potentials()
{
    int i;
    for (i=0; i < Npotentials; i++)
        eprintf("%15s - %s\n", potentials[i].name, potentials[i].desc);
}

int find_potential(char *name, struct potential *p)
{
    int i;
    for (i=0; i < Npotentials; i++)
    {
        if (!strcmp(name, potentials[i].name))
        {
            *p = potentials[i];
            return 1;
        }
    }

    return 0;
}

void add_potential(struct potential *p)
{
    Npotentials++;
    potentials = (struct potential *)realloc(potentials, Npotentials);
    potentials[Npotentials-1] = *p;
}

