#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <getopt.h>
//#include "../util/pca_utils.h"
#include "nbody.h"
#ifdef PGPLOT
#include <cpgplot.h>
#endif

#include <cuda.h>
//#include <nbody-cuda.h>

//#define SOFT 0.0001

#define GRAVCONST 1.         /* sets timescale */
#define SIM_RANDOM  1        /* which ICs are we going to use? */
#define SIM_2GAL 2
#define SIM_ROT2GAL 3

//#define DO_ASSERTS
//#define NO_MEMCPY

/*--------------------------------------------------------------------------------*/
float *vector(int n)
{
  float *vec=(float *)malloc(n*sizeof(float));
#ifdef DO_ASSERT
  assert(vec!=NULL);
#endif
  return vec;
}

/*--------------------------------------------------------------------------------*/
float *cuda_vector(int n)
{
  float *vec;  
#ifdef DO_ASSERTS
  assert(cudaMalloc ((void **) &vec, sizeof (float) * n)==cudaSuccess);
#else
  cudaMalloc ((void **) &vec, sizeof (float) * n);
#endif
  cudaMemset(vec,0,sizeof(float)*n);
  return vec;
}
/*--------------------------------------------------------------------------------*/
float *cuda_copy(float *vec, int n)
//Routine to take in input float vector, allocate space on the GPU for it, and copy the contents over.
//Returns a pointed to the device vector.
{
  float *dvec=cuda_vector(n);
#ifndef NO_MEMCPY
#ifdef DO_ASSERTS
  assert(cudaMemcpy(dvec,vec,sizeof(float)*n,cudaMemcpyHostToDevice)==cudaSuccess);
#else
  cudaMemcpy(dvec,vec,sizeof(float)*n,cudaMemcpyHostToDevice);
#endif
#endif

  return dvec;

}
/*--------------------------------------------------------------------------------*/
float *cuda_blocked_copy(float *vec, int n, int nn)
//Routine to take in input float vector, allocate space on the GPU for it, and copy the contents over.
//Returns a pointed to the device vector.
{
  float *dvec=cuda_vector(nn);
#ifndef NO_MEMCPY
#ifdef DO_ASSERTS
  assert(cudaMemset (dvec, 0, sizeof(float)*nn)==cudaSuccess);
  assert(cudaMemcpy(dvec,vec,sizeof(float)*n,cudaMemcpyHostToDevice)==cudaSuccess);
#else
  cudaMemset (dvec, 0, sizeof(float)*nn);
  cudaMemcpy(dvec,vec,sizeof(float)*n,cudaMemcpyHostToDevice);
#endif
#endif
  return dvec;

}

/*--------------------------------------------------------------------------------*/
__global__ void cuda_set_float_vec(int n, float *x, float val)
{
  int myind=blockIdx.x*blockDim.x+threadIdx.x;
  if (myind<n)
    x[myind]=val;
}

/*--------------------------------------------------------------------------------*/
float *get_blocked_masses( int n, int nn)  
//get a vector that is one up to n, and zero up to nn
{
  if (1) {
    float *dvec=cuda_vector(nn);
    //should do a memset here, but there's one in cuda_vector
    int nb=n/BLOCK+1;
    cuda_set_float_vec<<<nb,BLOCK>>>(n,dvec,1.0);
    return dvec;
  }
  else {
    float *vec=vector(nn);
    memset(vec,0,nn*sizeof(float));
    int i;
    for (i=0;i<n;i++)
      vec[i]=1;
    float *dvec=cuda_blocked_copy(vec,n,nn);
    free(vec);
    return dvec;
  }
  return NULL;
}

/*--------------------------------------------------------------------------------*/
void fill_random_vec(float *x, int n)
{
  int i;
  for (i=0;i<n;i++)
    x[i]=((float)rand())/((float)RAND_MAX);


}

/*--------------------------------------------------------------------------------*/

//cubepm call will look like:

__global__ void get_forces_cubepm(int n1, float *x, float *y, float *z, float *fx, float *fy, float *fz, int n2, float *x2, float *y2, float *z2, float *m2)
{
  
  //The __shared__ memory is a small block of cache shared between a block of threads.
  //To run fast, each block of threads copies global thread records into a local fast-memory
  //area.  
  float lx,ly,lz,lfx,lfy,lfz;
  __shared__ float lx2[BLOCK],ly2[BLOCK],lz2[BLOCK],lm2[BLOCK];

  //which particle do I own?
  int myind=blockIdx.x*blockDim.x+threadIdx.x;

  //These are now vector calls copying in contiguous chunks of memory into the local shared memory.
  lx=x[myind];
  ly=y[myind];
  lz=z[myind];

  //clear the forces on the particle I own.
  lfx=0;
  lfy=0;
  lfz=0;
  
  //make sure everybody is done before progressing.
  //syncthreads applies to a single block of threads, so 
  //it is extremely fast - something like 4 clock cycles.
  __syncthreads();
  float dx,dy,dz,r,r3;
  int i,blck,nblock;
  nblock=n2/BLOCK;
  int myind2;

  //goal is to pull blocks of particles into __shared__ memory, then loop over all
  //of these particles while they're around.  I have yet to figure out a good way to take 
  //advantage of the force from my particle on other particles, so going to be doing
  //double work.
  for (blck=0;blck<nblock;blck++)  {
    
    //which particle in the new block do I own?
    myind2=blck*blockDim.x+threadIdx.x;
    //load the particle block into cache.
    lx2[threadIdx.x]=x2[myind2];
    ly2[threadIdx.x]=y2[myind2];
    lz2[threadIdx.x]=z2[myind2];
    lm2[threadIdx.x]=m2[myind2];

    //make sure all the data is loaded before going on.
    __syncthreads();
    for (i=0;i<BLOCK;i++) {
      //everybody calculates the force from particle i in the block on 
      //themselves.
      dx=lx-lx2[i];
      dy=ly-ly2[i];
      dz=lz-lz2[i];
      r=dx*dx+dy*dy+dz*dz;// +EPS*EPS;
      if(sqrt(r) < EPS) {
         lfx+=0.0;
         lfy+=0.0;
         lfz+=0.0;
      }
      else{

         //rsqrt is the reciprocal square root.
         r=rsqrt(r);
         r3=r*r*r*lm2[i];  //lm2 is the mass of the second particle
         lfx+=r3*dx;
         lfy+=r3*dy;
         lfz+=r3*dz;
      }
      //end loop over particles in the temporary block
    }
    //why do I need a syncthreads here?  I get wrong answers without it.
    __syncthreads();    
    //end loop over blocks.
  }
  //now I have all the forces on the particles owned by this block of threads.  Send 'em back.
  fx[myind]=-lfx;
  fy[myind]=-lfy;
  fz[myind]=-lfz;
}

/*--------------------------------------------------------------------------------*/
void initialize_positions(float *x, float *y, float *z, int n)
{
  int i;
  for (i=0;i<n;i++) {
    x[i]=2.0*rand()/RAND_MAX - 1;
    y[i]=2.0*rand()/RAND_MAX - 1;
    z[i]=2.0*rand()/RAND_MAX - 1;
    
  }
}
/*--------------------------------------------------------------------------------*/
int get_blocked_size(int n)
{
  int nb=n/BLOCK;
  int nn=nb*BLOCK;
  if (nn<n)
    nn+=BLOCK;
  return nn;
}

/*--------------------------------------------------------------------------------*/
void calculate_pp_forces_cpu(float *x1, float *y1, float *z1, float *fx1, float *fy1, float *fz1, int n1,
			 float *x2, float *y2, float *z2, float *fx2, float *fy2, float *fz2, int n2)
{
  int i1,i2;

  memset(fx1,0,sizeof(float)*n1);
  memset(fy1,0,sizeof(float)*n1);
  memset(fz1,0,sizeof(float)*n1);
  memset(fx2,0,sizeof(float)*n2);
  memset(fy2,0,sizeof(float)*n2);
  memset(fz2,0,sizeof(float)*n2);

  for (i1=0;i1<n1;i1++) 
    for (i2=0;i2<n2;i2++) {
      float dx=x1[i1]-x2[i2];
      float dy=y1[i1]-y2[i2];
      float dz=z1[i1]-z2[i2];
      float rsqr=dx*dx+dy*dy+dz*dz; //+EPS*EPS;
      // Apply hard cutoff here
      if (sqrt(rsqr) < EPS) {

         fx1[i1]+=0.0;
         fy1[i1]+=0.0;
         fz1[i1]+=0.0;

         fx2[i2]+=0.0;
         fy2[i2]+=0.0;
         fz2[i2]+=0.0;

      } else {

         float rr=1/sqrt(rsqr);
         float r3=rr=rr*rr*rr;
         fx1[i1]-=r3*dx;
         fy1[i1]-=r3*dy;
         fz1[i1]-=r3*dz;
      
         fx2[i2]+=r3*dx;
         fy2[i2]+=r3*dy;
         fz2[i2]+=r3*dz;
      }
    }
  
}


/*--------------------------------------------------------------------------------*/
void calculate_pp_forces(float *x1, float *y1, float *z1, float *fx1, float *fy1, float *fz1, int n1,
			 float *x2, float *y2, float *z2, float *fx2, float *fy2, float *fz2, int n2)
{

  if (n1*n2<500*500*1) {
  //  printf("Working pp on CPU\n");
    calculate_pp_forces_cpu(x1,y1,z1,fx1,fy1,fz1,n1,x2,y2,z2,fx2,fy2,fz2,n2);
    return;
  }
  //printf("Working pp on GPU\n");

  int do_forces=1;
  if (n1<0) {
    n1*=-1;
    do_forces=0;
  }
  int nn1=get_blocked_size(n1);
  int nn2=get_blocked_size(n2);

  float *dx1=cuda_blocked_copy(x1,n1,nn1);
  float *dy1=cuda_blocked_copy(y1,n1,nn1);
  float *dz1=cuda_blocked_copy(z1,n1,nn1);
  float *dfx1=cuda_vector(nn1);
  float *dfy1=cuda_vector(nn1);
  float *dfz1=cuda_vector(nn1);

  float *dx2=cuda_blocked_copy(x2,n2,nn2);
  float *dy2=cuda_blocked_copy(y2,n2,nn2);
  float *dz2=cuda_blocked_copy(z2,n2,nn2);
  float *dfx2=cuda_vector(nn2);
  float *dfy2=cuda_vector(nn2);
  float *dfz2=cuda_vector(nn2);

  float *dm1=get_blocked_masses(n1,nn1);
  float *dm2=get_blocked_masses(n2,nn2);

  int nb1=nn1/BLOCK;
  int nb2=nn2/BLOCK;

  
  if (do_forces) {
    get_forces_cubepm<<<nb1,BLOCK>>>( nn1, dx1, dy1, dz1, dfx1, dfy1, dfz1, nn2, dx2,dy2,dz2, dm2);
    get_forces_cubepm<<<nb2,BLOCK>>>( nn2, dx2, dy2, dz2, dfx2, dfy2, dfz2, nn1, dx1,dy1,dz1, dm1);
  }

#ifndef NO_MEMCPY
  
#ifdef DO_ASSERTS
  assert(cudaMemcpy(fx1,dfx1,sizeof(float)*n1,cudaMemcpyDeviceToHost)==cudaSuccess);
  assert(cudaMemcpy(fy1,dfy1,sizeof(float)*n1,cudaMemcpyDeviceToHost)==cudaSuccess);
  assert(cudaMemcpy(fz1,dfz1,sizeof(float)*n1,cudaMemcpyDeviceToHost)==cudaSuccess);

  assert(cudaMemcpy(fx2,dfx2,sizeof(float)*n2,cudaMemcpyDeviceToHost)==cudaSuccess);
  assert(cudaMemcpy(fy2,dfy2,sizeof(float)*n2,cudaMemcpyDeviceToHost)==cudaSuccess);
  assert(cudaMemcpy(fz2,dfz2,sizeof(float)*n2,cudaMemcpyDeviceToHost)==cudaSuccess);
#else
  cudaMemcpy(fx1,dfx1,sizeof(float)*n1,cudaMemcpyDeviceToHost);
  cudaMemcpy(fy1,dfy1,sizeof(float)*n1,cudaMemcpyDeviceToHost);
  cudaMemcpy(fz1,dfz1,sizeof(float)*n1,cudaMemcpyDeviceToHost);

  cudaMemcpy(fx2,dfx2,sizeof(float)*n2,cudaMemcpyDeviceToHost);
  cudaMemcpy(fy2,dfy2,sizeof(float)*n2,cudaMemcpyDeviceToHost);
  cudaMemcpy(fz2,dfz2,sizeof(float)*n2,cudaMemcpyDeviceToHost);
#endif
#endif
  
  cudaFree(dx1);
  cudaFree(dy1);
  cudaFree(dz1);
  cudaFree(dfx1);
  cudaFree(dfy1);
  cudaFree(dfz1);
  cudaFree(dm1);

  cudaFree(dx2);
  cudaFree(dy2);
  cudaFree(dz2);
  cudaFree(dfx2);
  cudaFree(dfy2);
  cudaFree(dfz2);
  cudaFree(dm2);
}

/*--------------------------------------------------------------------------------*/

extern "C"
{
  
  void pp_force_c_(int *n1, float *x, float *y, float *z, float *fx, float *fy, float *fz, int *n2, float *x2, float *y2, float *z2)
  {
  float *fx2=(float *)malloc(sizeof(float)*(*n2));
  float *fy2=(float *)malloc(sizeof(float)*(*n2));
  float *fz2=(float *)malloc(sizeof(float)*(*n2));
  
  calculate_pp_forces(x,y,z,fx,fy,fz,*n1,x2,y2,z2,fx2,fy2,fz2,*n2);
  //printf("position of particle 1 is %17.7e %17.7e %17.7e\n",x[0],y[0],z[0]);
  //printf("position of particle 2 is %17.7e %17.7e %17.7e\n",x2[1],y2[1],z2[1]);
  //printf("force on particle 1 is %17.7e %17.7e %17.7e\n",fx[0],fy[0],fz[0]);
  //printf("force on particle 2 is %17.7e %17.7e %17.7e\n",fx2[1],fy2[1],fz2[1]);
  free(fx2);
  free(fy2);
  free(fz2);
  
  }
}
#if 0
/*================================================================================*/
int main(int argc, char *argv[])
{
  printf("hello!\n");
  int n1,n2;
  if (argc>1) {
    n1=atoi(argv[1]);
  }
  else
  n1=1000;
  if (argc>2) 
    n2=atoi(argv[2]);
  else
    n2=n1;
  
  printf("starting particles are %d %d\n",n1,n2);
  
  float *x1=vector(n1);
  float *y1=vector(n1);
  float *z1=vector(n1);

  float *fx1=vector(n1);
  float *fy1=vector(n1);
  float *fz1=vector(n1);

  float *x2=vector(n2);
  float *y2=vector(n2);
  float *z2=vector(n2);

  float *fx2=vector(n2);
  float *fy2=vector(n2);
  float *fz2=vector(n2);
  
  memset(fx1,0,sizeof(float)*n1);
  memset(fy1,0,sizeof(float)*n1);
  memset(fz1,0,sizeof(float)*n1);

  memset(fx2,0,sizeof(float)*n2);
  memset(fy2,0,sizeof(float)*n2);
  memset(fz2,0,sizeof(float)*n2);


  initialize_positions(x1,y1,z1,n1);
  initialize_positions(x2,y2,z2,n2);



  //do a first call as takes a while to get up and running
  calculate_pp_forces(x1,y1,z1,fx1,fy1,fz1,n1,x2,y2,z2,fx2,fy2,fz2,n2);
  

  float frep=1e10/(n1+2000.0)/(n2+2000.0);
  int nrep=frep;  
  printf("going to do %d repetitions\n",nrep);

  double t1=omp_get_wtime();
  for (int i=0;i<nrep;i++)
    calculate_pp_forces(x2,y2,z2,fx2,fy2,fz2,n2,x1,y1,z1,fx1,fy1,fz1,n1);
  double gpu_time=omp_get_wtime();
  gpu_time-=t1;
  printf("force on particle 1 is %17.7e %17.7e, took %14.4g seconds\n",fx1[0],fx2[0],gpu_time/nrep);

  t1=omp_get_wtime();
  for (int i=0;i<nrep;i++)
    calculate_pp_forces(x2,y2,z2,fx2,fy2,fz2,-1*n2,x1,y1,z1,fx1,fy1,fz1,n1);
  double gpu_overhead_time=omp_get_wtime();
  gpu_overhead_time-=t1;
  printf("force on particle 1 is %17.7e %17.7e, overhead took %14.4g seconds\n",fx1[0],fx2[0],gpu_overhead_time/nrep);



  t1=omp_get_wtime();

    calculate_pp_forces_cpu(x2,y2,z2,fx2,fy2,fz2,n2,x1,y1,z1,fx1,fy1,fz1,n1);
  double cpu_time=omp_get_wtime();
  cpu_time-=t1;
  printf("force on particle 1 is %17.7e %17.7e, took %14.4g seconds\n",fx1[0],fx2[0],cpu_time);
  printf("ratio is %14.4f\n",cpu_time/gpu_time*nrep);
  printf("pure compute ratio is %14.4f\n",cpu_time/(gpu_time-gpu_overhead_time)*nrep);
  

  return 0;
}
#endif
