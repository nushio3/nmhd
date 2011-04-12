
#include "library/mhd_solver.h"
#include "library/math.h"

#include <thrust/extrema.h>
#include <thrust/transform.h>

const int kThreadDim = 32, kBlockDim = 32;

MHDSolver::MHDSolver () :
  density_(raw(buffer_[0])), velocity_(raw(buffer_[1]),raw(buffer_[2]),raw(buffer_[3])),
  pressure_(raw(buffer_[4])), magnet_(raw(buffer_[5]),raw(buffer_[6]),raw(buffer_[7])),
  energy_(raw(buffer_[8])), momentum_(raw(buffer_[9]),raw(buffer_[10]),raw(buffer_[11])) ,

  cfl_time_(raw(buffer_[12])),

  density_pace_(raw(buffer_[13])), momentum_pace_(raw(buffer_[14]),raw(buffer_[15]),raw(buffer_[16])),
  energy_pace_(raw(buffer_[17])), magnet_pace_(raw(buffer_[18]),raw(buffer_[19]),raw(buffer_[20])),
  
  density_flux_(raw(buffer_[21])), momentum_flux_(raw(buffer_[22]),raw(buffer_[23]),raw(buffer_[24])),
  energy_flux_(raw(buffer_[25])), magnet_flux_ey_(raw(buffer_[26])), magnet_flux_ez_(raw(buffer_[27])),
  
  // buffer_ will be constructed *earlier* than above elements
  buffer_(28, thrust::device_vector<Real>(kSize)),

  current_time_(0), generation_(0)
{
  for (int i = 0; i < buffer_.size(); ++i) {
    Real *p = raw(buffer_[i]);
    buffer_dict_[p] = i;
  }
  initial_condition();
}

__global__ void initial_condition_kernel
(BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    set_initial_condition
      (o, density, velocity, pressure, magnet);
  }
}
void MHDSolver::initial_condition () {
  initial_condition_kernel<<<kBlockDim, kThreadDim>>>
    (density_, velocity_, pressure_, magnet_);
}


__global__ void boundary_condition_kernel
(BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    set_boundary_condition
      (o, density, velocity, pressure, magnet);
  }
}
void MHDSolver::boundary_condition () {
  boundary_condition_kernel<<<kBlockDim, kThreadDim>>>
    (density_, velocity_, pressure_, magnet_);  
}



template <class Direction>
__device__ Real cfl_condition_inner
(const Direction ex, const RCo<7> &o,
 BodyScalar density, BodyVector velocity, BodyScalar pressure, Real b2) {
  const Real len = dR(ex);
  const Real speed = absR(velocity[ex][o] + sqrtR((kGamma * pressure[o] + b2) / density[o]));
  return kCflFactor * len / speed;
}
__global__ void cfl_condition_kernel
(BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
 const Real max_time,  BodyScalar cfl_time) {
  CRYSTAL_MAP(addr) {
    Real min_time = max_time;
    RCo<7> o(addr);
    const Real b2 =
      max(sq(magnet[X][o-HX]), sq(magnet[X][o+HX]))
+      max(sq(magnet[Y][o-HY]), sq(magnet[Y][o+HY]))
+      max(sq(magnet[Z][o-HZ]), sq(magnet[Z][o+HZ]))
;
      
    min_time = min(min_time, cfl_condition_inner(X, o, density, velocity, pressure, b2));
    min_time = min(min_time, cfl_condition_inner(Y, o, density, velocity, pressure, b2));
    min_time = min(min_time, cfl_condition_inner(Z, o, density, velocity, pressure, b2));
    cfl_time[o] = min_time;
  }
}

__global__ void clear_pace_kernel
(BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace) {
  CRYSTAL_MAP(addr) {
    const RCo<7> o(addr);
    density_pace[o] = 0;
    momentum_pace[X][o] = 0;
    momentum_pace[Y][o] = 0;
    momentum_pace[Z][o] = 0;
    energy_pace[o] = 0;
    magnet_pace[X][o-HX] = 0;
    magnet_pace[Y][o-HY] = 0;
    magnet_pace[Z][o-HZ] = 0;
  }
}

#include "library/hlld.inl"

template<class EX>
__global__
void calc_flux_kernel
(EX ex,
 BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
 typename BodyScalar::Shift<EX::mask>::t density_flux ,
 Triplet<typename BodyScalar::Shift<EX::mask>::t> momentum_flux ,
 typename BodyScalar::Shift<EX::mask>::t energy_flux,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ey,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ez
 ) {
  typename EX::Next ey; typename EX::Prev ez;
  typename EX::Half hex; typename EX::Half::Next hey; typename EX::Half::Prev hez;
  CRYSTAL_MAP(addr) {
    const RCo<7> right(addr);
    const RCo<7^EX::mask> center = right - hex;
    const RCo<7> left = right - ex;
    hlld_flux(magnet[ex][center],
	      density[left], 
	      velocity[ex][left], velocity[ey][left], velocity[ez][left],
	      pressure[left],
	      average(magnet[ey][left-hey], magnet[ey][left+hey]),
	      average(magnet[ez][left-hez], magnet[ez][left+hez]),
	      
	      density[right], 
	      velocity[ex][right], velocity[ey][right], velocity[ez][right],
	      pressure[right],
	      average(magnet[ey][right-hey], magnet[ey][right+hey]),
	      average(magnet[ez][right-hez], magnet[ez][right+hez]),

	      density_flux[center],
	      momentum_flux[ex][center], momentum_flux[ey][center], momentum_flux[ez][center],
	      energy_flux[center],
	      magnet_flux_ey[center], magnet_flux_ez[center]);
  }
}

/*
  the EX direction flux of the EY component of the magnetic field at coordinate <corner>.
 */
template<class EX, class EY, class Coord> __device__ Real corner_magnet_flux
(EX ex, EY ey, Coord corner, BodyVector wind,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ey) {
  typename EX::Half hex; typename EY::Half hey;
  Real sum_wind = 
 wind[ey][corner+hex+hey] 
+ wind[ey][corner+hex-hey] 
+ wind[ey][corner-hex+hey] 
+ wind[ey][corner-hex-hey] 
;
  const Real f0 = magnet_flux_ey[corner - hey];
  const Real f1 = magnet_flux_ey[corner + hey];
  return sum_wind > 0
    ? f0
    : sum_wind < 0
    ? f1
    : average(f0, f1);
}

/*
  the EX direction flux of the EY component of the magnetic field at coordinate <face>,
  differenciated in HEX1 direction.
 */
template<class EX, class EY, class HEX1, class Coord> __device__ Real magnet_flux_difference
(EX ex, EY ey, HEX1 hex1, Coord face,
 BodyVector wind, typename BodyScalar::Shift<EX::mask>::t magnet_flux_ey
 ) {
  return
    + corner_magnet_flux(ex, ey, face + hex1, wind, magnet_flux_ey)
    - corner_magnet_flux(ex, ey, face - hex1, wind, magnet_flux_ey);
}

template<class EX> __global__ void add_pace_kernel
(EX ex,
 BodyVector wind,
 typename BodyScalar::Shift<EX::mask>::t density_flux ,
 Triplet<typename BodyScalar::Shift<EX::mask>::t> momentum_flux ,
 typename BodyScalar::Shift<EX::mask>::t energy_flux,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ey,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ez,
 BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace
 ) {
  typedef typename EX::Next EY; typedef typename EX::Prev EZ;
  EY ey; EZ ez;
  typename EX::Half hex; typename EY::Half hey; typename EZ::Half hez;
  CRYSTAL_MAP(addr) {
    const RCo<7> center(addr);
    const RCo<7^EX::mask> left  = center - hex;
    const RCo<7^EX::mask> right = center + hex;
    
    density_pace[center] += (density_flux[left] - density_flux[right]) / dR(ex);
    momentum_pace[X][center] += (momentum_flux[X][left] - momentum_flux[X][right]) / dR(ex);
    momentum_pace[Y][center] += (momentum_flux[Y][left] - momentum_flux[Y][right]) / dR(ex);
    momentum_pace[Z][center] += (momentum_flux[Z][left] - momentum_flux[Z][right]) / dR(ex);
    energy_pace[center] += (energy_flux[left] - energy_flux[right]) / dR(ex);

    const RCo<7^EX::mask> ox = center - hex;
    const RCo<7^EY::mask> oy = center - hey;
    const RCo<7^EZ::mask> oz = center - hez;

    const Real magic_factor = 0.5f;
    
    magnet_pace[ey][oy] -= magic_factor *
      magnet_flux_difference(ex, ey, hex, oy, wind, magnet_flux_ey) / dR(ex);
    magnet_pace[ez][oz] -= magic_factor *
      magnet_flux_difference(ex, ez, hex, oz, wind, magnet_flux_ez) / dR(ex);
    magnet_pace[ex][ox] += magic_factor *
      magnet_flux_difference(ex, ey, hey, ox, wind, magnet_flux_ey) / dR(ey) + 
      magnet_flux_difference(ex, ez, hez, ox, wind, magnet_flux_ez) / dR(ez);
  }
}



Real max_diff(thrust::device_vector<Real> xs, thrust::device_vector<Real> ys) {
  Real ans = 0;
  for (int i = 0; i < xs.size(); ++i) {
    ans = max(ans, absR(xs[i] - ys[i]));
  }
  return ans;
}

void MHDSolver::add_pace_kernel_caller
(BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
 BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace) {
  {
    Mesh<7^XAxis::mask, Real> density_flux = density_flux_.shift(X) ;
    Triplet<Mesh<7^XAxis::mask,Real> > momentum_flux
			 (
 momentum_flux_[X].shift(X) 
, momentum_flux_[Y].shift(X) 
, momentum_flux_[Z].shift(X) 
);
    Mesh<7^XAxis::mask, Real> energy_flux = energy_flux_.shift(X) ;
    Mesh<7^XAxis::mask, Real> magnet_flux_ey =  magnet_flux_ey_.shift(X);
    Mesh<7^XAxis::mask, Real> magnet_flux_ez =  magnet_flux_ez_.shift(X);

    calc_flux_kernel<<<kBlockDim, kThreadDim>>>
      (X, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_pace_kernel<<<kBlockDim, kThreadDim>>>
      (X,
       velocity, 
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez,
       density_pace, momentum_pace, energy_pace, magnet_pace);
  }
  {
    Mesh<7^YAxis::mask, Real> density_flux = density_flux_.shift(Y) ;
    Triplet<Mesh<7^YAxis::mask,Real> > momentum_flux
			 (
 momentum_flux_[X].shift(Y) 
, momentum_flux_[Y].shift(Y) 
, momentum_flux_[Z].shift(Y) 
);
    Mesh<7^YAxis::mask, Real> energy_flux = energy_flux_.shift(Y) ;
    Mesh<7^YAxis::mask, Real> magnet_flux_ey =  magnet_flux_ey_.shift(Y);
    Mesh<7^YAxis::mask, Real> magnet_flux_ez =  magnet_flux_ez_.shift(Y);

    calc_flux_kernel<<<kBlockDim, kThreadDim>>>
      (Y, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_pace_kernel<<<kBlockDim, kThreadDim>>>
      (Y,
       velocity, 
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez,
       density_pace, momentum_pace, energy_pace, magnet_pace);
  }
  {
    Mesh<7^ZAxis::mask, Real> density_flux = density_flux_.shift(Z) ;
    Triplet<Mesh<7^ZAxis::mask,Real> > momentum_flux
			 (
 momentum_flux_[X].shift(Z) 
, momentum_flux_[Y].shift(Z) 
, momentum_flux_[Z].shift(Z) 
);
    Mesh<7^ZAxis::mask, Real> energy_flux = energy_flux_.shift(Z) ;
    Mesh<7^ZAxis::mask, Real> magnet_flux_ey =  magnet_flux_ey_.shift(Z);
    Mesh<7^ZAxis::mask, Real> magnet_flux_ez =  magnet_flux_ez_.shift(Z);

    calc_flux_kernel<<<kBlockDim, kThreadDim>>>
      (Z, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_pace_kernel<<<kBlockDim, kThreadDim>>>
      (Z,
       velocity, 
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez,
       density_pace, momentum_pace, energy_pace, magnet_pace);
  }
}


__global__ void to_conserved_kernel
(BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
 BodyVector momentum, BodyScalar energy) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    momentum[X][o] = density[o] * velocity[X][o];
    momentum[Y][o] = density[o] * velocity[Y][o];
    momentum[Z][o] = density[o] * velocity[Z][o];
    energy[o] =  pressure[o] / (kGamma-1)
      + Real(0.5f) * density[o] * (sq(velocity[X][o]) + sq(velocity[Y][o]) + sq(velocity[Z][o]))
      + Real(0.25f) * (sq(magnet[X][o-HX]) + sq(magnet[X][o+HX]) +
		       sq(magnet[Y][o-HY]) + sq(magnet[Y][o+HY]) +
		       sq(magnet[Z][o-HZ]) + sq(magnet[Z][o+HZ]));
  }
}

__global__ void to_primitive_kernel
(BodyScalar density, BodyVector momentum, BodyScalar energy, FaceVector magnet,
 BodyVector velocity, BodyScalar pressure) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    velocity[X][o] = momentum[X][o] / density[o];
    velocity[Y][o] = momentum[Y][o] / density[o];
    velocity[Z][o] = momentum[Z][o] / density[o];
    const Real internal_energy = energy[o]
      - Real(0.5f) * (sq(momentum[X][o]) + sq(momentum[Y][o]) + sq(momentum[Z][o])) / density[o]
      - Real(0.25f) * (sq(magnet[X][o-HX]) + sq(magnet[X][o+HX]) +
		       sq(magnet[Y][o-HY]) + sq(magnet[Y][o+HY]) +
		       sq(magnet[Z][o-HZ]) + sq(magnet[Z][o+HZ]));
    pressure[o] = internal_energy * (kGamma-1);
  }
}

__global__ void update_kernel
(Real dt,
 BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace,
 BodyScalar density, BodyVector momentum, BodyScalar energy, FaceVector magnet) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    density[o] += dt * density_pace[o];
    momentum[X][o] += dt * momentum_pace[X][o];
    momentum[Y][o] += dt * momentum_pace[Y][o];
    momentum[Z][o] += dt * momentum_pace[Z][o];
    energy[o] += dt * energy_pace[o];
    magnet[X][o-HX] += dt * magnet_pace[X][o-HX];
    magnet[Y][o-HY] += dt * magnet_pace[Y][o-HY];
    magnet[Z][o-HZ] += dt * magnet_pace[Z][o-HZ];
  }
}


void MHDSolver::proceed (const Real integral_time) {
  bool break_flag = false;
  Real elapsed_time = 0;
  const Real start_time = current_time();

  do {
    if(monitor()) break;
    std::cerr << "elapsed: " << elapsed_time << std::flush;
    
    boundary_condition();
    cfl_condition_kernel<<<kBlockDim, kThreadDim>>>
      (density_, velocity_, pressure_, magnet_, integral_time, cfl_time_);
    double dt = 0.0f;
    {
      static thrust::device_vector<Real>::iterator
	cfl_begin = buffer_dict(cfl_time_).begin(),
	cfl_end = buffer_dict(cfl_time_).end();
      thrust::device_vector<Real>::iterator min_it = thrust::min_element(cfl_begin, cfl_end);
      dt = *min_it;
    }
    std::cerr << "\tdt = " << dt << std::endl;
    
    const double next_time = elapsed_time + dt;
    if (next_time > integral_time) {
      dt = integral_time - elapsed_time;
      break_flag = true;
    }
    
    //// Zero-Initialize Pace
    clear_pace_kernel<<<kBlockDim, kThreadDim>>>
      (density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Calculate Flux and add it to Pace
    add_pace_kernel_caller
      (density_, velocity_, pressure_, magnet_,
       density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Convert to Conserved Variables
    to_conserved_kernel<<<kBlockDim, kThreadDim>>>
      (density_, velocity_, pressure_, magnet_,
       momentum_, energy_);
    
    //// Update Conserved Variables
    update_kernel<<<kBlockDim, kThreadDim>>>
      (dt,
       density_pace_, momentum_pace_, energy_pace_, magnet_pace_,
       density_, momentum_, energy_, magnet_);
    
    //// Convert to Primitive Variables
    to_primitive_kernel<<<kBlockDim, kThreadDim>>>
      (density_, momentum_, energy_, magnet_,
       velocity_, pressure_);

    elapsed_time += dt;
    ++generation_;
    current_time_ = start_time + elapsed_time;
  } while(!break_flag);
}

