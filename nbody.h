#ifndef NBODY_H
#define NBODY_H

void nbody_create_potential(struct potential *phi, int N);
void nbody_advance_phi(struct potential *phi, tyme_t t);
__device__ void accel(struct nbody_potential *phi0, struct particle *p, acc_t *out);

#endif
