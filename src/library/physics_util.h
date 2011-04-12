#pragma once
#include "library/crystal.h"

typedef Mesh<7,Real> BodyScalar;
typedef Triplet<Mesh<7,Real> > BodyVector;
typedef Triad<Mesh<6,Real>, Mesh<5,Real>, Mesh<3,Real> > FaceVector;
typedef Triad<Mesh<1,Real>, Mesh<2,Real>, Mesh<4,Real> > EdgeVector;

#define CRYSTAL_MAP(addr)						\
  for (int addr = blockIdx.x * blockDim.x + threadIdx.x; addr < kSize;	\
       addr += blockDim.x * gridDim.x) 

#define CUSTOM_CRYSTAL_MAP(addr, size)					\
  for (int addr = blockIdx.x * blockDim.x + threadIdx.x; addr < size;	\
       addr += blockDim.x * gridDim.x) 
