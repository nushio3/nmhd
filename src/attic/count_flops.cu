//typedef float Real;
//#define USE_SINGLE_PRECISION
typedef double Real;
#define USE_DOUBLE_PRECISION

const Real kGamma = 5.0f/3.0f;

#include "../library/math.h"
#include "hlld.inl"


__global__
void hlld_flux_kernel
(Real *Bx, 
	       
 Real *dens_L, Real *dens_R, 
 Real *velx_L, Real *velx_R, 
 Real *vely_L, Real *vely_R, 
 Real *velz_L, Real *velz_R, 
 Real *pres_L, Real *pres_R, 
 Real *By_L,   Real *By_R,   
 Real *Bz_L,   Real *Bz_R,   
 
 Real *Fdens, 
 Real *Fmomx, Real *Fmomy, Real *Fmomz,
 Real *Fetot, 
 Real *F_By,  Real *F_Bz) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  hlld_flux(Bx[tid], 
	       
            dens_L[tid], dens_R[tid], 
            velx_L[tid], velx_R[tid], 
            vely_L[tid], vely_R[tid], 
            velz_L[tid], velz_R[tid], 
            pres_L[tid], pres_R[tid], 
            By_L[tid],   By_R[tid],   
            Bz_L[tid],   Bz_R[tid],   
            
            Fdens[tid], 
            Fmomx[tid], Fmomy[tid], Fmomz[tid],
            Fetot[tid], 
            F_By[tid],  F_Bz[tid]); 
}

int main () {}
