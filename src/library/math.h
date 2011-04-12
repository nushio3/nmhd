#pragma once

const Real kPi = 3.14159265358979323846;

#ifdef USE_SINGLE_PRECISION

#define absR fabs
#define sqrtR sqrtf
#define sinR sinf
#define cosR cosf
#define expR expf

#endif
#ifdef USE_DOUBLE_PRECISION

#define absR abs
#define sqrtR sqrt
#define sinR sin
#define cosR cos
#define expR exp

#endif

template <class T>
__host__ __device__ T sq (const T & x) { return x*x; }

template <class T>
__host__ __device__ T average (const T & a,const T & b) { return T(0.5f)*(a+b); }

template <class T>
__host__ __device__ T enpack (const T & s1, const T & s2, const T & x1, const T & x2, const T & x3) {
  return ((x3 * s2) + x2) * s1 + x1;
}

template <class T>
__host__ __device__ void depack (const T & input, const T & s1, const T & s2, T & x1, T & x2, T & x3) {
  x1 = input % s1;
  T tmp = input / s1;
  x2 = tmp % s2;
  x3 = tmp / s2;
}
