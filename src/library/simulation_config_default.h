#pragma once
#include "library/direction.h"


#undef USE_MPI

#define USE_SINGLE_PRECISION
#undef  USE_DOUBLE_PRECISION


#ifdef USE_SINGLE_PRECISION
typedef float Real;
#endif
#ifdef USE_DOUBLE_PRECISION
typedef double Real;
#endif

// floating point type used for calculations with high precision, such as accumlations
typedef double PreciseReal;

// mathematical functions defined
#include "library/math.h"

// physical global constants
#define kGamma  Real(5.0f/3.0f) // unfortunate nvcc bug, const cannot be used
const Real kCflFactor = 0.5f;

// type used to sum up the crystal address
typedef int addr_t;

// crystal addressing mode periodicity
#undef CRYSTAL_PERIODIC_X
#undef CRYSTAL_PERIODIC_Y 
#undef CRYSTAL_PERIODIC_Z

// crystal size, computational and physical
const int kGlobalSizeX = 10;
const int kGlobalSizeY = 10;
const int kGlobalSizeZ = 10;

const int kParallelSizeX = 1;				
const int kParallelSizeY = 1;				
const int kParallelSizeZ = 1;
const int kParallelSize = kParallelSizeX * kParallelSizeY * kParallelSizeZ;

const int kMarginSizeX = 0;
const int kMarginSizeY = kParallelSizeY > 1 ? 4 : 0;
const int kMarginSizeZ = kParallelSizeZ > 1 ? 4 : 0;

__device__ const int gMarginSizeX = kMarginSizeX;
__device__ const int gMarginSizeY = kMarginSizeY;
__device__ const int gMarginSizeZ = kMarginSizeZ;



const int kLocalSizeX = kGlobalSizeX / kParallelSizeX;
const int kLocalSizeY = kGlobalSizeY / kParallelSizeY;
const int kLocalSizeZ = kGlobalSizeZ / kParallelSizeZ;

const int kSizeX = kLocalSizeX + 2*kMarginSizeX;	      
const int kSizeY = kLocalSizeY + 2*kMarginSizeY;	      
const int kSizeZ = kLocalSizeZ + 2*kMarginSizeZ;	      
const size_t kSize = size_t(kSizeX) * size_t(kSizeY) * size_t(kSizeZ);

__device__ const int gSizeX = kSizeX;
__device__ const int gSizeY = kSizeY;
__device__ const int gSizeZ = kSizeZ;
__device__ const int gSize = kSize;


const Real kGlobalExtentX = 1.0;
const Real kGlobalExtentY = 1.0;
const Real kGlobalExtentZ = 1.0;

__device__ const Real gGlobalExtentX = kGlobalExtentX;
__device__ const Real gGlobalExtentY = kGlobalExtentY;
__device__ const Real gGlobalExtentZ = kGlobalExtentZ;

const Real kExtentX = kGlobalExtentX / Real(kGlobalSizeX) * Real(kSizeX);
const Real kExtentY = kGlobalExtentY / Real(kGlobalSizeY) * Real(kSizeY);
const Real kExtentZ = kGlobalExtentZ / Real(kGlobalSizeZ) * Real(kSizeZ);

__device__ const Real gExtentX = kExtentX;
__device__ const Real gExtentY = kExtentY;
__device__ const Real gExtentZ = kExtentZ;

const Real kLocalExtentX = kGlobalExtentX / Real(kGlobalSizeX) * Real(kLocalSizeX);
const Real kLocalExtentY = kGlobalExtentY / Real(kGlobalSizeY) * Real(kLocalSizeY);
const Real kLocalExtentZ = kGlobalExtentZ / Real(kGlobalSizeZ) * Real(kLocalSizeZ);

__device__ const Real gLocalExtentX = kLocalExtentX;
__device__ const Real gLocalExtentY = kLocalExtentY;
__device__ const Real gLocalExtentZ = kLocalExtentZ;



const Real kGlobalCornerX = 0.0;
const Real kGlobalCornerY = 0.0;
const Real kGlobalCornerZ = 0.0;

const Real kDX = kExtentX / kSizeX;
const Real kDY = kExtentY / kSizeY;
const Real kDZ = kExtentZ / kSizeZ;

const Real kMinDR = min(kDX, min(kDY,kDZ));

namespace {
namespace simulation_config {
__device__ __host__ int size (XAxis) { return kSizeX; }
__device__ __host__ int size (YAxis) { return kSizeY; }
__device__ __host__ int size (ZAxis) { return kSizeZ; }

__device__ __host__ Real extent (XAxis) { return kExtentX; }
__device__ __host__ Real extent (YAxis) { return kExtentY; }
__device__ __host__ Real extent (ZAxis) { return kExtentZ; }
}
}

namespace {
__device__ __host__ Real dR(XAxis x) { return kDX; }
__device__ __host__ Real dR(YAxis y) { return kDY; }
__device__ __host__ Real dR(ZAxis z) { return kDZ; }
}
