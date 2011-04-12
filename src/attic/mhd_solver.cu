

#include "library/direction.h"
#include "simulation_config.h"
#include "library/physics_util.h"


__device__ void set_initial_condition
(const RCo<7> &o,
 BodyScalar &density, BodyVector &velocity, BodyScalar &pressure, FaceVector &magnet) {
  ShockNormalAxis ex;
  typename ShockNormalAxis::Next ey;
  typename ShockNormalAxis::Prev ez;
  
  const Triplet<Real> r = o.position();
  if (r[ex] < 0) {
    density[o] = 1.0f;
    pressure[o] = 1.0f;
    magnet[ey][o-direction::half(ey)] = 1.0f;
  } else {
    density[o] = 0.125f;
    pressure[o] = 0.1f;
    magnet[ey][o-direction::half(ey)] = -1.0f;
  }
  magnet[ex][o-direction::half(ex)] = 0.75f;
  magnet[ez][o-direction::half(ez)] = 0.0f;
 velocity[X][o] = 0.0f; 
 velocity[Y][o] = 0.0f; 
 velocity[Z][o] = 0.0f; 
}

__device__ void set_boundary_condition
(const RCo<7> &o,
 BodyScalar &density, BodyVector &velocity, BodyScalar &pressure, FaceVector &magnet) {
  /*
  const int kMargin = 4;
  if(o[X] < kMargin || o[X] > kSizeX - kMargin - 1) {
    set_initial_condition(o, density, velocity, pressure, magnet);
    }*/
}

#include "library/mhd_solver.inl"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

void MHDSolver::dump (string tag) {
  typename ShockNormalAxis::Next ey;
  typename ShockNormalAxis::Prev ez;

  string dir = "stat/" + simulation_tag;
  
  system(("mkdir -p " + dir).c_str());
  {
    ofstream ofs((dir + "/dump_"+tag+"_hydro.txt").c_str());
    for (int addr = 0; addr < kSize; ++addr) {
      RCo<7> o(addr);

	
      if (      (o[ey] == simulation_config::size(ey)/2) && (o[ez] == simulation_config::size(ez)/2)
) {
	ofs << o.position()
	    << " " << buffer_dict(density_)[addr]
	    << " " << buffer_dict(velocity_[X])[addr]
	    << " " << buffer_dict(velocity_[Y])[addr]
	    << " " << buffer_dict(velocity_[Z])[addr]
	    << " " << buffer_dict(pressure_)[addr]
	    <<endl;
      }
    }
  }
  {
    ofstream ofs((dir + "/dump_"+tag+"_bx.txt").c_str());
    for (int addr = 0; addr < kSize; ++addr) {
      RCo<7 ^ XAxis::mask> o(addr);
      if (      (o[ey] == simulation_config::size(ey)/2) && (o[ez] == simulation_config::size(ez)/2)
) {
	ofs << o.position()
	    << " " << buffer_dict(magnet_[X])[addr]
	    <<endl;
      }
    }
  }
  {
    ofstream ofs((dir + "/dump_"+tag+"_by.txt").c_str());
    for (int addr = 0; addr < kSize; ++addr) {
      RCo<7 ^ YAxis::mask> o(addr);
      if (      (o[ey] == simulation_config::size(ey)/2) && (o[ez] == simulation_config::size(ez)/2)
) {
	ofs << o.position()
	    << " " << buffer_dict(magnet_[Y])[addr]
	    <<endl;
      }
    }
  }
  {
    ofstream ofs((dir + "/dump_"+tag+"_bz.txt").c_str());
    for (int addr = 0; addr < kSize; ++addr) {
      RCo<7 ^ ZAxis::mask> o(addr);
      if (      (o[ey] == simulation_config::size(ey)/2) && (o[ez] == simulation_config::size(ez)/2)
) {
	ofs << o.position()
	    << " " << buffer_dict(magnet_[Z])[addr]
	    <<endl;
      }
    }
  }
}

int MHDSolver::monitor () {
  return 0;
}

