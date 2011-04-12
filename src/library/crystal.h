#pragma once

#include "library/direction.h"
#include "simulation_config.h"


struct PositioningSetup {
  PositioningSetup(Triplet<Real> corner0, int y0) : corner(corner0), mpi_rank_y(y0) {}
  Triplet<Real> corner; // Physical coordinate of the computational corner.
  int mpi_rank_y;
};

template<int Bit> struct RCo{
  int x_,y_,z_; addr_t addr_;
  __device__ __host__ RCo (int x, int y, int z, addr_t addr0) : x_(x), y_(y), z_(z), addr_(addr0) {}
  __device__ __host__ RCo (int x, int y, int z) : x_(x), y_(y), z_(z) {
    calc_addr();
  }
  __device__ __host__ explicit RCo (addr_t addr0) : addr_(addr0) {
    x_ = addr0 % kSizeX;
    y_ = (addr0 / kSizeX) % kSizeY;
    z_ = (addr0 / kSizeX) / kSizeY;
  }
  __device__ __host__
  int bit() const { return Bit; }
  __device__ __host__
  int addr() const {
    return addr_;
  }
  __device__ __host__
  int calc_addr() {
    return addr_ = (
#ifdef CRYSTAL_PERIODIC_Z
  ((z_ + kSizeZ) % kSizeZ)
#else
  (z_)
#endif
 * kSizeY + 
#ifdef CRYSTAL_PERIODIC_Y
  ((y_ + kSizeY) % kSizeY)
#else
  (y_)
#endif
) * kSizeX + 
#ifdef CRYSTAL_PERIODIC_X
  ((x_ + kSizeX) % kSizeX)
#else
  (x_)
#endif
 ;
  }
  __device__ __host__
  int operator[](XAxis) const { return 
#ifdef CRYSTAL_PERIODIC_X
  ((x_ + kSizeX) % kSizeX)
#else
  (x_)
#endif
; }
  __device__ __host__
  int operator[](YAxis) const { return 
#ifdef CRYSTAL_PERIODIC_Y
  ((y_ + kSizeY) % kSizeY)
#else
  (y_)
#endif
; }
  __device__ __host__
  int operator[](ZAxis) const { return 
#ifdef CRYSTAL_PERIODIC_Z
  ((z_ + kSizeZ) % kSizeZ)
#else
  (z_)
#endif
; }
  
  __device__ __host__
  Triplet<Real> position(const Triplet<Real> &computational_corner) const {
 ;
    return Triplet<Real> (    computational_corner[X] + kDX * (Real(
#ifdef CRYSTAL_PERIODIC_X
  ((x_ + kSizeX) % kSizeX)
#else
  (x_)
#endif
) + (Bit & 1 == 0 ? Real(0) : Real(0.5f)))

,
    computational_corner[Y] + kDY * (Real(
#ifdef CRYSTAL_PERIODIC_Y
  ((y_ + kSizeY) % kSizeY)
#else
  (y_)
#endif
) + (Bit & 2 == 0 ? Real(0) : Real(0.5f)))

,
    computational_corner[Z] + kDZ * (Real(
#ifdef CRYSTAL_PERIODIC_Z
  ((z_ + kSizeZ) % kSizeZ)
#else
  (z_)
#endif
) + (Bit & 4 == 0 ? Real(0) : Real(0.5f)))
);
  }
  __device__ __host__
  Triplet<Real> position(const PositioningSetup &pos_setup) const {
    return position(pos_setup.corner);
  }
};


template<int Bit, class T> struct Mesh {
  T *ptr_;
  __device__ __host__ Mesh(T *ptr) : ptr_(ptr) {}
  __device__ __host__ T& operator[](RCo<Bit> coord) { return ptr_[coord.addr()]; }
  template<int Bit2> struct Shift {
    typedef Mesh<Bit ^ Bit2, T> t;
  };
  template<int Bit2>
  Mesh<Bit ^ Bit2, T> shift () const {
    return Mesh<Bit ^ Bit2, T>(ptr_);
  }
  template<class Axis>
  Mesh<Bit ^ Axis::mask, T> shift (Axis a) const {
    return Mesh<Bit ^ Axis::mask, T>(ptr_);
  }
  T *ptr() const { return ptr_; }
};
template <int Bit>
__device__ __host__
RCo<Bit> operator+(const RCo<Bit> rco, XAxis) {
  return RCo<Bit>
    (
     rco.x_ + 1
,     rco.y_
,     rco.z_
     );
}
template <int Bit>
__device__ __host__
RCo<Bit> operator+(const RCo<Bit> rco, YAxis) {
  return RCo<Bit>
    (
     rco.x_
,     rco.y_ + 1
,     rco.z_
     );
}
template <int Bit>
__device__ __host__
RCo<Bit> operator+(const RCo<Bit> rco, ZAxis) {
  return RCo<Bit>
    (
     rco.x_
,     rco.y_
,     rco.z_ + 1
     );
}
template <int Bit>
__device__ __host__
RCo<Bit> operator-(const RCo<Bit> rco, XAxis) {
  return RCo<Bit>
    (
     rco.x_ - 1
,     rco.y_
,     rco.z_
     );
}
template <int Bit>
__device__ __host__
RCo<Bit> operator-(const RCo<Bit> rco, YAxis) {
  return RCo<Bit>
    (
     rco.x_
,     rco.y_ - 1
,     rco.z_
     );
}
template <int Bit>
__device__ __host__
RCo<Bit> operator-(const RCo<Bit> rco, ZAxis) {
  return RCo<Bit>
    (
     rco.x_
,     rco.y_
,     rco.z_ - 1
     );
}
namespace {
__device__ __host__  
RCo<1> operator+(const RCo<0> rco, XHalf) {
;
  return RCo<1>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<0> operator+(const RCo<1> rco, XHalf) {
;
  return RCo<0>(  rco.x_ + 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<3> operator+(const RCo<2> rco, XHalf) {
;
  return RCo<3>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<2> operator+(const RCo<3> rco, XHalf) {
;
  return RCo<2>(  rco.x_ + 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<5> operator+(const RCo<4> rco, XHalf) {
;
  return RCo<5>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<4> operator+(const RCo<5> rco, XHalf) {
;
  return RCo<4>(  rco.x_ + 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<7> operator+(const RCo<6> rco, XHalf) {
;
  return RCo<7>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<6> operator+(const RCo<7> rco, XHalf) {
;
  return RCo<6>(  rco.x_ + 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<2> operator+(const RCo<0> rco, YHalf) {
;
  return RCo<2>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<3> operator+(const RCo<1> rco, YHalf) {
;
  return RCo<3>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<0> operator+(const RCo<2> rco, YHalf) {
;
  return RCo<0>(  rco.x_ 
,  rco.y_ + 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<1> operator+(const RCo<3> rco, YHalf) {
;
  return RCo<1>(  rco.x_ 
,  rco.y_ + 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<6> operator+(const RCo<4> rco, YHalf) {
;
  return RCo<6>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<7> operator+(const RCo<5> rco, YHalf) {
;
  return RCo<7>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<4> operator+(const RCo<6> rco, YHalf) {
;
  return RCo<4>(  rco.x_ 
,  rco.y_ + 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<5> operator+(const RCo<7> rco, YHalf) {
;
  return RCo<5>(  rco.x_ 
,  rco.y_ + 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<4> operator+(const RCo<0> rco, ZHalf) {
;
  return RCo<4>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<5> operator+(const RCo<1> rco, ZHalf) {
;
  return RCo<5>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<6> operator+(const RCo<2> rco, ZHalf) {
;
  return RCo<6>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<7> operator+(const RCo<3> rco, ZHalf) {
;
  return RCo<7>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<0> operator+(const RCo<4> rco, ZHalf) {
;
  return RCo<0>(  rco.x_ 
,  rco.y_ 
,  rco.z_ + 1
);
}
}
namespace {
__device__ __host__  
RCo<1> operator+(const RCo<5> rco, ZHalf) {
;
  return RCo<1>(  rco.x_ 
,  rco.y_ 
,  rco.z_ + 1
);
}
}
namespace {
__device__ __host__  
RCo<2> operator+(const RCo<6> rco, ZHalf) {
;
  return RCo<2>(  rco.x_ 
,  rco.y_ 
,  rco.z_ + 1
);
}
}
namespace {
__device__ __host__  
RCo<3> operator+(const RCo<7> rco, ZHalf) {
;
  return RCo<3>(  rco.x_ 
,  rco.y_ 
,  rco.z_ + 1
);
}
}
namespace {
__device__ __host__  
RCo<1> operator-(const RCo<0> rco, XHalf) {
;
  return RCo<1>(  rco.x_ - 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<0> operator-(const RCo<1> rco, XHalf) {
;
  return RCo<0>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<3> operator-(const RCo<2> rco, XHalf) {
;
  return RCo<3>(  rco.x_ - 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<2> operator-(const RCo<3> rco, XHalf) {
;
  return RCo<2>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<5> operator-(const RCo<4> rco, XHalf) {
;
  return RCo<5>(  rco.x_ - 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<4> operator-(const RCo<5> rco, XHalf) {
;
  return RCo<4>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<7> operator-(const RCo<6> rco, XHalf) {
;
  return RCo<7>(  rco.x_ - 1
,  rco.y_ 
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<6> operator-(const RCo<7> rco, XHalf) {
;
  return RCo<6>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<2> operator-(const RCo<0> rco, YHalf) {
;
  return RCo<2>(  rco.x_ 
,  rco.y_ - 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<3> operator-(const RCo<1> rco, YHalf) {
;
  return RCo<3>(  rco.x_ 
,  rco.y_ - 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<0> operator-(const RCo<2> rco, YHalf) {
;
  return RCo<0>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<1> operator-(const RCo<3> rco, YHalf) {
;
  return RCo<1>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<6> operator-(const RCo<4> rco, YHalf) {
;
  return RCo<6>(  rco.x_ 
,  rco.y_ - 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<7> operator-(const RCo<5> rco, YHalf) {
;
  return RCo<7>(  rco.x_ 
,  rco.y_ - 1
,  rco.z_ 
);
}
}
namespace {
__device__ __host__  
RCo<4> operator-(const RCo<6> rco, YHalf) {
;
  return RCo<4>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<5> operator-(const RCo<7> rco, YHalf) {
;
  return RCo<5>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<4> operator-(const RCo<0> rco, ZHalf) {
;
  return RCo<4>(  rco.x_ 
,  rco.y_ 
,  rco.z_ - 1
);
}
}
namespace {
__device__ __host__  
RCo<5> operator-(const RCo<1> rco, ZHalf) {
;
  return RCo<5>(  rco.x_ 
,  rco.y_ 
,  rco.z_ - 1
);
}
}
namespace {
__device__ __host__  
RCo<6> operator-(const RCo<2> rco, ZHalf) {
;
  return RCo<6>(  rco.x_ 
,  rco.y_ 
,  rco.z_ - 1
);
}
}
namespace {
__device__ __host__  
RCo<7> operator-(const RCo<3> rco, ZHalf) {
;
  return RCo<7>(  rco.x_ 
,  rco.y_ 
,  rco.z_ - 1
);
}
}
namespace {
__device__ __host__  
RCo<0> operator-(const RCo<4> rco, ZHalf) {
;
  return RCo<0>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<1> operator-(const RCo<5> rco, ZHalf) {
;
  return RCo<1>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<2> operator-(const RCo<6> rco, ZHalf) {
;
  return RCo<2>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}
namespace {
__device__ __host__  
RCo<3> operator-(const RCo<7> rco, ZHalf) {
;
  return RCo<3>(  rco.x_ 
,  rco.y_ 
,  rco.z_ 
,rco.addr_);
}
}

