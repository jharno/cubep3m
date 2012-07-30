//#include "../util/pca_utils.h"

#define NDIM 3
#define BLOCK 32
#define EPS 0.1



typedef struct nbody_struct_s {
  float x[NDIM];    /*the particle positions*/
  float v[NDIM];    /*the particle velocities*/
  float f[NDIM];    /*the forces on the particles*/
  float mass;
  float PE;        /* potential energy */
} NBody;

void calculate_forces_gpu(NBody *data, int n);
