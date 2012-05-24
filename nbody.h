#ifndef NBODY_H
#define NBODY_H

void nbody_create_potential(struct potential *phi, int N);
void nbody_advance_phi(struct potential *phi, tyme_t t);

#endif
