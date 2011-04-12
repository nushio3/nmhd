#include "simulation_config.h"


#include "library/mhd_solver.h"
#include "library/get_time.h"
#include "library/math.h"



#include <thrust/extrema.h>
#include <thrust/transform.h>


MHDSolver::MHDSolver (LocalSetup &local_setup) :
  local_setup_(local_setup),
  
  density_(raw(buffer_[0])), velocity_(raw(buffer_[1]),raw(buffer_[2]),raw(buffer_[3])),
  pressure_(raw(buffer_[4])), magnet_(raw(buffer_[5]),raw(buffer_[6]),raw(buffer_[7])),
  energy_(raw(buffer_[8])), momentum_(raw(buffer_[9]),raw(buffer_[10]),raw(buffer_[11])) ,

  cfl_time_(raw(buffer_[12])),

  density_pace_(raw(buffer_[13])), momentum_pace_(raw(buffer_[14]),raw(buffer_[15]),raw(buffer_[16])),
  energy_pace_(raw(buffer_[17])), magnet_pace_(raw(buffer_[18]),raw(buffer_[19]),raw(buffer_[20])),
  
  density_flux_(raw(buffer_[21])), momentum_flux_(raw(buffer_[22]),raw(buffer_[23]),raw(buffer_[24])),
  energy_flux_(raw(buffer_[25])), magnet_flux_ey_(raw(buffer_[26])), magnet_flux_ez_(raw(buffer_[27])),

  density_predictor_(raw(buffer_[28])), velocity_predictor_(raw(buffer_[29]),raw(buffer_[30]),raw(buffer_[31])),
  pressure_predictor_(raw(buffer_[32])), magnet_predictor_(raw(buffer_[33]),raw(buffer_[34]),raw(buffer_[35])),
  
  // buffer_ will be constructed *earlier* than above elements
  buffer_(36, thrust::device_vector<Real>(kSize)),
  host_buffer_(kSize),
  mpi_buffer_initialized(false)  
{
  for (int i = 0; i < buffer_.size(); ++i) {
    Real *p = raw(buffer_[i]);
    buffer_dict_[p] = i;
  }
}

__global__ void initial_condition_kernel
(PositioningSetup pos_setup, BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    set_initial_condition
      (pos_setup, o, density, velocity, pressure, magnet);
  }
}
void MHDSolver::initial_condition () {
  initial_condition_kernel <<<32, 32>>>
  (local_setup().positioning_setup() , density_, velocity_, pressure_, magnet_);
}

__global__ void boundary_condition_kernel
(PositioningSetup pos_setup, BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    set_boundary_condition
      (pos_setup, o, density, velocity, pressure, magnet);
  }
}
void MHDSolver::boundary_condition () {
  boundary_condition_kernel <<<32, 32>>>
    (local_setup().positioning_setup(), density_, velocity_, pressure_, magnet_);  
}
__device__ Real absolute_smaller(Real a, Real b) {
  return a*b<=0
    ? 0
    : absR(a) < absR(b)
    ? a : b;
}

__device__ void interpolate(const Real x0, const Real x1, const Real x2, const Real x3, Real &left, Real &right) {
  const Real d01 = x1-x0;
  const Real d12 = x2-x1;
  const Real d23 = x3-x2;
  const Real d1 = absolute_smaller(d01,d12);
  const Real d2 = absolute_smaller(d12,d23);
  left  = x1 + d1/2;
  right = x2 - d2/2;
}




#ifdef USE_MPI

__global__ void communicate_gather_kernel_y
(int displacement_int_inc, Real displacement_real_inc, Real relative_velocity_inc,
 int displacement_int_dec, Real displacement_real_dec, Real relative_velocity_dec,

 Real *buf_inc, Real *buf_dec, Real *density, Real *velocity_x, Real *velocity_y, Real *velocity_z, Real *pressure, Real *magnet_x, Real *magnet_y, Real *magnet_z ) {
  const int kUnitSizeY = gSizeX * gMarginSizeY * gSizeZ;

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = density[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = density[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = density[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = density[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[0 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[0 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = velocity_x[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = velocity_x[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = velocity_x[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = velocity_x[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[1 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
      -relative_velocity_inc ;
    buf_dec[1 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
      +relative_velocity_dec ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = velocity_y[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = velocity_y[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = velocity_y[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = velocity_y[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[2 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[2 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = velocity_z[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = velocity_z[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = velocity_z[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = velocity_z[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[3 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[3 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = pressure[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = pressure[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = pressure[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = pressure[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[4 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[4 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = magnet_x[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = magnet_x[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = magnet_x[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = magnet_x[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[5 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[5 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = magnet_y[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = magnet_y[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = magnet_y[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = magnet_y[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[6 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[6 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }

  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY , sx, sy ,sz);
    int inc_x0 = (sx + displacement_int_inc    ) % gSizeX;
    int inc_x1 = (sx + displacement_int_inc + 1) % gSizeX;
    int dec_x0 = (sx - displacement_int_dec - 1  + gSizeX) % gSizeX;
    int dec_x1 = (sx - displacement_int_dec      + gSizeX) % gSizeX;
    Real val_inc0 = magnet_z[ enpack(gSizeX, gSizeY, inc_x0, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_inc1 = magnet_z[ enpack(gSizeX, gSizeY, inc_x1, gSizeY - 2 * gMarginSizeY + sy, sz) ];
    Real val_dec0 = magnet_z[ enpack(gSizeX, gSizeY, dec_x0, gMarginSizeY + sy, sz) ];
    Real val_dec1 = magnet_z[ enpack(gSizeX, gSizeY, dec_x1, gMarginSizeY + sy, sz) ];
    buf_inc[7 * kUnitSizeY + addr] = (Real(1)-displacement_real_inc) * val_inc0 + displacement_real_inc * val_inc0 
       ;
    buf_dec[7 * kUnitSizeY + addr] = displacement_real_dec * val_dec0 + (Real(1)-displacement_real_dec) * val_dec0 
       ;
  }
}

__global__ void communicate_scatter_kernel_y
(Real *buf_inc, Real *buf_dec, Real *density, Real *velocity_x, Real *velocity_y, Real *velocity_z, Real *pressure, Real *magnet_x, Real *magnet_y, Real *magnet_z ) {
  const int kUnitSizeY = gSizeX * gMarginSizeY * gSizeZ;
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    density[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[0 * kUnitSizeY + addr];
    density[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[0 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    velocity_x[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[1 * kUnitSizeY + addr];
    velocity_x[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[1 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    velocity_y[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[2 * kUnitSizeY + addr];
    velocity_y[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[2 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    velocity_z[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[3 * kUnitSizeY + addr];
    velocity_z[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[3 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    pressure[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[4 * kUnitSizeY + addr];
    pressure[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[4 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    magnet_x[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[5 * kUnitSizeY + addr];
    magnet_x[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[5 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    magnet_y[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[6 * kUnitSizeY + addr];
if (sy > 0)
    magnet_y[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[6 * kUnitSizeY + addr];
  }
  CUSTOM_CRYSTAL_MAP(addr, kUnitSizeY) {
    int sx, sy, sz;
    depack(addr, gSizeX, gMarginSizeY, sx, sy ,sz);
    magnet_z[ enpack(gSizeX, gSizeY, sx, sy, sz) ] = buf_inc[7 * kUnitSizeY + addr];
    magnet_z[ enpack(gSizeX, gSizeY, sx, gSizeY - gMarginSizeY + sy, sz) ] = buf_dec[7 * kUnitSizeY + addr];
  }
}

void MHDSolver::relative_displacement (PreciseReal t, int index, Real &normalized_displacement, int &int_displacement, Real &minor_displacement) {
  PreciseReal displacement0 = initial_displacement(index) + t * relative_velocity(index);
  int folding_number = int(floor(displacement0 / kExtentX));
  PreciseReal displacement1 = displacement0 - folding_number * kExtentX;
  normalized_displacement = displacement1;
  int_displacement = int(floor(displacement1 / kDX));
  minor_displacement = displacement1 / kDX - int_displacement;
}

void MHDSolver::update_displacement () {
  int displacement_int; Real normalized_displacement, displacement_real;
  relative_displacement(global_setup().current_time(), local_setup().mpi_rank_3d()[Y],
			  normalized_displacement, displacement_int, displacement_real);
  local_setup().relative_displacement() = normalized_displacement;
  local_setup().relative_velocity() = relative_velocity(local_setup().mpi_rank_3d()[Y]);
}



void MHDSolver::communicate () {
  
  
  static int message_tag = 0;
  
  const int kUnitSizeY = kSizeX * kMarginSizeY * kSizeZ;
  const int kUnitSizeZ = kSizeX * kSizeY * kMarginSizeZ;
  const int kTotalSizeY = kUnitSizeY * 8;
  const int kTotalSizeZ = kUnitSizeZ * 8;
  if (!mpi_buffer_initialized) {
    
    mpi_buffer_y_inc_send.resize(kTotalSizeY);
    mpi_buffer_y_inc_recv.resize(kTotalSizeY);
    mpi_buffer_y_inc_device.resize(kTotalSizeY);
    mpi_buffer_y_dec_send.resize(kTotalSizeY);
    mpi_buffer_y_dec_recv.resize(kTotalSizeY);
    mpi_buffer_y_dec_device.resize(kTotalSizeY);

    mpi_buffer_z_inc_send.resize(kTotalSizeZ);
    mpi_buffer_z_inc_recv.resize(kTotalSizeZ);
    mpi_buffer_z_dec_send.resize(kTotalSizeZ);
    mpi_buffer_z_dec_recv.resize(kTotalSizeZ);
    mpi_buffer_initialized = true;
    
  }

  const int my_rank_x = local_setup().mpi_rank_3d()[X];
  const int my_rank_y = local_setup().mpi_rank_3d()[Y];
  const int my_rank_z = local_setup().mpi_rank_3d()[Z];
  const int my_rank_y_inc = (my_rank_y + 1) % kParallelSizeY;
  const int my_rank_y_dec = (my_rank_y + kParallelSizeY - 1) % kParallelSizeY;
  const int my_rank_z_inc = (my_rank_z + 1) % kParallelSizeZ;
  const int my_rank_z_dec = (my_rank_z + kParallelSizeZ - 1) % kParallelSizeZ;
  const int rank_y_inc = (my_rank_z * kParallelSizeY + my_rank_y_inc) * kParallelSizeX + my_rank_x;
  const int rank_y_dec = (my_rank_z * kParallelSizeY + my_rank_y_dec) * kParallelSizeX + my_rank_x;
  const int rank_z_inc = (my_rank_z_inc * kParallelSizeY + my_rank_y) * kParallelSizeX + my_rank_x;
  const int rank_z_dec = (my_rank_z_dec * kParallelSizeY + my_rank_y) * kParallelSizeX + my_rank_x;

  if (global_setup().generation() <= 0)
    local_setup().ofs_log() << "communicating neighbours of ("
                          << local_setup().mpi_rank() << ") = "
                          << rank_y_dec << " "
                          << rank_y_inc << " "
                          << rank_z_dec << " "
                          << rank_z_inc << std::endl;
  
  if (kMarginSizeY > 0) {
    //// //// Begin Communications in Y Direction //// ////

    
    

    {
      int displacement_int_inc, displacement_int_dec;
      Real displacement_real_inc, displacement_real_dec;
      Real normalized_displacement_inc, normalized_displacement_dec;
      
      relative_displacement(global_setup().current_time(), my_rank_y,  normalized_displacement_inc, displacement_int_inc, displacement_real_inc);
      relative_displacement(global_setup().current_time(), rank_y_dec, normalized_displacement_dec, displacement_int_dec, displacement_real_dec);
      
      communicate_gather_kernel_y  <<<32, 32>>>
      (
       displacement_int_inc, displacement_real_inc, relative_velocity(my_rank_y), 
       displacement_int_dec, displacement_real_dec, relative_velocity(rank_y_dec),
       
       raw(mpi_buffer_y_inc_device),
       raw(mpi_buffer_y_dec_device),
       density_.ptr(), velocity_[X].ptr(), velocity_[Y].ptr(), velocity_[Z].ptr(), pressure_.ptr(), magnet_[X].ptr(), magnet_[Y].ptr(), magnet_[Z].ptr()
       );
    }
    
    cudaMemcpy(raw(mpi_buffer_y_inc_send),
	       raw(mpi_buffer_y_inc_device),
	       sizeof(Real) * kTotalSizeY,
	       cudaMemcpyDeviceToHost);
    cudaThreadSynchronize(); //Redundant?

    cudaMemcpy(raw(mpi_buffer_y_dec_send),
	       raw(mpi_buffer_y_dec_device),
	       sizeof(Real) * kTotalSizeY,
	       cudaMemcpyDeviceToHost);
    cudaThreadSynchronize(); //Redundant?

  cudaThreadSynchronize();
  
  MPI_Barrier(MPI_COMM_WORLD); // wait until all the nodes are ready to communicate
  {
    //// initiate rank-increasing communication for Y direction
    MPI_Request send_req_inc, recv_req_inc;
    MPI_Status send_stat_inc, recv_stat_inc;

    
    MPI_Irecv(raw(mpi_buffer_y_inc_recv), kTotalSizeY, kMPIDatatypeReal,
	      rank_y_dec, message_tag, MPI_COMM_WORLD, &send_req_inc);
    MPI_Isend(raw(mpi_buffer_y_inc_send), kTotalSizeY, kMPIDatatypeReal,
	      rank_y_inc, message_tag, MPI_COMM_WORLD, &recv_req_inc);
    ++message_tag;
    //// initiate rank-decreasing communication for Y direction
    MPI_Request send_req_dec, recv_req_dec;
    MPI_Status send_stat_dec, recv_stat_dec;

    
    MPI_Irecv(raw(mpi_buffer_y_dec_recv), kTotalSizeY, kMPIDatatypeReal,
	      rank_y_inc, message_tag, MPI_COMM_WORLD, &send_req_dec);
    MPI_Isend(raw(mpi_buffer_y_dec_send), kTotalSizeY, kMPIDatatypeReal,
	      rank_y_dec, message_tag, MPI_COMM_WORLD, &recv_req_dec);
    ++message_tag;
    //// wait until completion of the access
    
    MPI_Wait(&recv_req_inc, &recv_stat_inc);
    MPI_Wait(&send_req_inc, &send_stat_inc);
    
    //// wait until completion of the access
    
    MPI_Wait(&recv_req_dec, &recv_stat_dec);
    MPI_Wait(&send_req_dec, &send_stat_dec);
    
  }
  cudaMemcpy(raw(mpi_buffer_y_inc_device),
	     raw(mpi_buffer_y_inc_recv),
	     sizeof(Real) * kTotalSizeY,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_y_dec_device),
	     raw(mpi_buffer_y_dec_recv),
	     sizeof(Real) * kTotalSizeY,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  
  communicate_scatter_kernel_y  <<<32, 32>>>
  (raw(mpi_buffer_y_inc_device),
   raw(mpi_buffer_y_dec_device),
   density_.ptr(), velocity_[X].ptr(), velocity_[Y].ptr(), velocity_[Z].ptr(), pressure_.ptr(), magnet_[X].ptr(), magnet_[Y].ptr(), magnet_[Z].ptr()
   );
  //// //// End Communications in Y(y) Direction //// ////
  }
  if (kMarginSizeZ > 0) {
    //// //// Begin Communications in Z Direction //// ////

    
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 0 * kUnitSizeZ,
	     density_.ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 0 * kUnitSizeZ,
	     density_.ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 1 * kUnitSizeZ,
	     velocity_[X].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 1 * kUnitSizeZ,
	     velocity_[X].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 2 * kUnitSizeZ,
	     velocity_[Y].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 2 * kUnitSizeZ,
	     velocity_[Y].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 3 * kUnitSizeZ,
	     velocity_[Z].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 3 * kUnitSizeZ,
	     velocity_[Z].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 4 * kUnitSizeZ,
	     pressure_.ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 4 * kUnitSizeZ,
	     pressure_.ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 5 * kUnitSizeZ,
	     magnet_[X].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 5 * kUnitSizeZ,
	     magnet_[X].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 6 * kUnitSizeZ,
	     magnet_[Y].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 6 * kUnitSizeZ,
	     magnet_[Y].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  // bring data to be communicated in Z direction from GPU to CPU
  
  cudaMemcpy(raw(mpi_buffer_z_inc_send) + 7 * kUnitSizeZ,
	     magnet_[Z].ptr() + (kSize - 2 * kUnitSizeZ),
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(raw(mpi_buffer_z_dec_send) + 7 * kUnitSizeZ,
	     magnet_[Z].ptr() + kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyDeviceToHost);
  cudaThreadSynchronize(); //Redundant?
 
  cudaThreadSynchronize();
  
  MPI_Barrier(MPI_COMM_WORLD); // wait until all the nodes are ready to communicate
  {
    //// initiate rank-increasing communication for Z direction
    MPI_Request send_req_inc, recv_req_inc;
    MPI_Status send_stat_inc, recv_stat_inc;

    
    MPI_Irecv(raw(mpi_buffer_z_inc_recv), kTotalSizeZ, kMPIDatatypeReal,
	      rank_z_dec, message_tag, MPI_COMM_WORLD, &send_req_inc);
    MPI_Isend(raw(mpi_buffer_z_inc_send), kTotalSizeZ, kMPIDatatypeReal,
	      rank_z_inc, message_tag, MPI_COMM_WORLD, &recv_req_inc);
    ++message_tag;
    //// initiate rank-decreasing communication for Z direction
    MPI_Request send_req_dec, recv_req_dec;
    MPI_Status send_stat_dec, recv_stat_dec;

    
    MPI_Irecv(raw(mpi_buffer_z_dec_recv), kTotalSizeZ, kMPIDatatypeReal,
	      rank_z_inc, message_tag, MPI_COMM_WORLD, &send_req_dec);
    MPI_Isend(raw(mpi_buffer_z_dec_send), kTotalSizeZ, kMPIDatatypeReal,
	      rank_z_dec, message_tag, MPI_COMM_WORLD, &recv_req_dec);
    ++message_tag;
    //// wait until completion of the access
    
    MPI_Wait(&recv_req_inc, &recv_stat_inc);
    MPI_Wait(&send_req_inc, &send_stat_inc);
    
    //// wait until completion of the access
    
    MPI_Wait(&recv_req_dec, &recv_stat_dec);
    MPI_Wait(&send_req_dec, &send_stat_dec);
    
  }
  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(density_.ptr(),
	     raw(mpi_buffer_z_inc_recv) + 0 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(density_.ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 0 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(velocity_[X].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 1 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(velocity_[X].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 1 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(velocity_[Y].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 2 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(velocity_[Y].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 2 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(velocity_[Z].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 3 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(velocity_[Z].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 3 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(pressure_.ptr(),
	     raw(mpi_buffer_z_inc_recv) + 4 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(pressure_.ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 4 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(magnet_[X].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 5 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(magnet_[X].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 5 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(magnet_[Y].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 6 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(magnet_[Y].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 6 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  // bring data that was communicated in Z direction from CPU to GPU
  cudaMemcpy(magnet_[Z].ptr(),
	     raw(mpi_buffer_z_inc_recv) + 7 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  cudaMemcpy(magnet_[Z].ptr() + (kSize - 1 * kUnitSizeZ),
	     raw(mpi_buffer_z_dec_recv) + 7 * kUnitSizeZ,
	     sizeof(Real) * kUnitSizeZ,
	     cudaMemcpyHostToDevice);
  cudaThreadSynchronize(); //Redundant?

  //// //// End Communications in Z(z) Direction //// ////
  }
  
}
#endif //#ifdef USE_MPI

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
(EX ex, Real dt,
 BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
 typename BodyScalar::Shift<EX::mask>::t density_flux ,
 Triplet<typename BodyScalar::Shift<EX::mask>::t> momentum_flux ,
 typename BodyScalar::Shift<EX::mask>::t energy_flux,
 typename BodyScalar::Shift<EX::mask>::t magnet_flux_ey,
 typename BodyScalar::Shift<

 EX::mask>::t magnet_flux_ez
 ) {
  typename EX::Next ey; typename EX::Prev ez;
  typename EX::Half hex; typename EX::Half::Next hey; typename EX::Half::Prev hez;
  CRYSTAL_MAP(addr) {
    const RCo<7^EX::mask> center(addr);
    const RCo<7> c0 = center - hex - ex;
    const RCo<7> c1 = center - hex;
    const RCo<7> c2 = center + hex;
    const RCo<7> c3 = center + hex + ex;
 Real density_left, density_right;
    interpolate(density[c0],density[c1],density[c2],density[c3], density_left, density_right);
 Real velocity_ex_left, velocity_ex_right;
    interpolate(velocity[ex][c0],velocity[ex][c1],velocity[ex][c2],velocity[ex][c3], velocity_ex_left, velocity_ex_right);
 Real velocity_ey_left, velocity_ey_right;
    interpolate(velocity[ey][c0],velocity[ey][c1],velocity[ey][c2],velocity[ey][c3], velocity_ey_left, velocity_ey_right);
 Real velocity_ez_left, velocity_ez_right;
    interpolate(velocity[ez][c0],velocity[ez][c1],velocity[ez][c2],velocity[ez][c3], velocity_ez_left, velocity_ez_right);
 Real pressure_left, pressure_right;
    interpolate(pressure[c0],pressure[c1],pressure[c2],pressure[c3], pressure_left, pressure_right);
 Real magnet_ey_left, magnet_ey_right;
    interpolate(average(magnet[ey][c0-hey],magnet[ey][c0+hey]),average(magnet[ey][c1-hey],magnet[ey][c1+hey]),average(magnet[ey][c2-hey],magnet[ey][c2+hey]),average(magnet[ey][c3-hey],magnet[ey][c3+hey]), magnet_ey_left, magnet_ey_right);
 Real magnet_ez_left, magnet_ez_right;
    interpolate(average(magnet[ez][c0-hez],magnet[ez][c0+hez]),average(magnet[ez][c1-hez],magnet[ez][c1+hez]),average(magnet[ez][c2-hez],magnet[ez][c2+hez]),average(magnet[ez][c3-hez],magnet[ez][c3+hez]), magnet_ez_left, magnet_ez_right);
    
    Real mesh_boundary_velocity = Real(0.5) * dt *
      external_acceleration(density, velocity, pressure, magnet, center)[ex];
    velocity_ex_left  += mesh_boundary_velocity;
    velocity_ex_right += mesh_boundary_velocity;
    


    
    hlld_flux
      (magnet[ex][center],
       
       density_left, density_right, velocity_ex_left, velocity_ex_right, velocity_ey_left, velocity_ey_right, velocity_ez_left, velocity_ez_right, pressure_left, pressure_right, magnet_ey_left, magnet_ey_right, magnet_ez_left, magnet_ez_right,
       
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

template<class EX> __global__ void add_flux_kernel
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
  const Real magic_factor = 0.5f;
  const Real black_magic_factor = 1.0f;

  CRYSTAL_MAP(addr) {
    const RCo<7> center(addr);
    const RCo<7^EX::mask> left  = center - hex;
    const RCo<7^EX::mask> right = center + hex;
    
    density_pace[center] += (density_flux[left] - density_flux[right]) / dR(ex);
    momentum_pace[X][center] += (momentum_flux[X][left] - momentum_flux[X][right]) / dR(ex);
    momentum_pace[Y][center] += (momentum_flux[Y][left] - momentum_flux[Y][right]) / dR(ex);
    momentum_pace[Z][center] += (momentum_flux[Z][left] - momentum_flux[Z][right]) / dR(ex);
    energy_pace[center] += (energy_flux[left] - energy_flux[right]) / dR(ex);

    
    
    const RCo<7^EY::mask> oy(addr);
    magnet_pace[ey][oy] -= magic_factor *
      magnet_flux_difference(ex, ey, hex, oy, wind, magnet_flux_ey) / dR(ex);

    
    
    const RCo<7^EZ::mask> oz(addr);
    magnet_pace[ez][oz] -= magic_factor *
      magnet_flux_difference(ex, ez, hex, oz, wind, magnet_flux_ez) / dR(ex);

    
    
    const RCo<7^EX::mask> ox(addr);
    magnet_pace[ex][ox] += magic_factor * black_magic_factor * 
      (magnet_flux_difference(ex, ey, hey, ox, wind, magnet_flux_ey) / dR(ey) +
       magnet_flux_difference(ex, ez, hez, ox, wind, magnet_flux_ez) / dR(ez) );
  }
}


// debug
Real debug_max_diff(thrust::device_vector<Real> xs, thrust::device_vector<Real> ys) {
  Real ans = 0;
  for (int i = 0; i < xs.size(); ++i) {
    ans = max(ans, absR(xs[i] - ys[i]));
  }
  return ans;
}

Real debug_min(thrust::device_vector<Real> xs) {
  Real ans = xs[0];
  for (int i = 0; i < xs.size(); ++i) {
    ans = min(ans, xs[i]);
  }
  return ans;
}


void MHDSolver::add_flux_kernel_caller
(Real dt, BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
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

    calc_flux_kernel <<<32, 32>>>
      (X, dt, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_flux_kernel <<<32, 32>>>
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

    calc_flux_kernel <<<32, 32>>>
      (Y, dt, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_flux_kernel <<<32, 32>>>
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

    calc_flux_kernel <<<32, 32>>>
      (Z, dt, density, velocity, pressure, magnet,
       density_flux, momentum_flux, energy_flux, magnet_flux_ey, magnet_flux_ez);

    add_flux_kernel <<<32, 32>>>
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

__global__ void add_pace_kernel
(Real dt,
 BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace,
 BodyScalar density_old,  FaceVector magnet_old,
 BodyScalar density_new, BodyVector momentum, BodyScalar energy, FaceVector magnet_new
 ) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    density_new[o] = density_old[o] + dt * density_pace[o];
    
    momentum[X][o] += dt * momentum_pace[X][o];
    
    momentum[Y][o] += dt * momentum_pace[Y][o];
    
    momentum[Z][o] += dt * momentum_pace[Z][o];
    
    energy[o] += dt * energy_pace[o];
    
    magnet_new[X][o-HX] = magnet_old[X][o-HX] + dt * magnet_pace[X][o-HX];
    
    magnet_new[Y][o-HY] = magnet_old[Y][o-HY] + dt * magnet_pace[Y][o-HY];
    
    magnet_new[Z][o-HZ] = magnet_old[Z][o-HZ] + dt * magnet_pace[Z][o-HZ];
    
  }
}

__global__ void external_acceleration_kernel
(Real dt,
 BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet
 ) {
  CRYSTAL_MAP(addr) {
    RCo<7> o(addr);
    Triplet<Real> a = external_acceleration(density, velocity, pressure, magnet, o);
    velocity[X][o] += dt * a[X];
    velocity[Y][o] += dt * a[Y];
    velocity[Z][o] += dt * a[Z];
  }
}

bool MHDSolver::proceed (const Real goal_time, const PreciseReal wallclock_time_limit) {
  bool break_flag = false;
  bool timeout = false;
  bool goal_reached = false;
  PreciseReal wallclock_begin = get_time<PreciseReal>();
  
  //std::cerr << "kMPIDatatypeReal = " << kMPIDatatypeReal << std::endl;
  //std::cerr << "MPI_REAL = " << MPI_DOUBLE_PRECISION << " " << MPI_REAL << std::endl;

  do {
    if(monitor()) break;
    if (wallclock_time_limit > 0 && get_time<PreciseReal>() > wallclock_begin + wallclock_time_limit) {
      timeout = true;
    }
    
#ifdef USE_MPI
    communicate();
#endif
    Real dt;
    boundary_condition();
    cfl_condition_kernel
      <<<32, 32>>>
    (density_, velocity_, pressure_, magnet_, goal_time, cfl_time_);
    dt = 0.0f;
    {
      static thrust::device_vector<Real>::iterator
	cfl_begin = buffer_dict(cfl_time_).begin(),
	cfl_end = buffer_dict(cfl_time_).end();
      thrust::device_vector<Real>::iterator min_it = thrust::min_element(cfl_begin, cfl_end);
      Real my_dt = *min_it;
      if (timeout) my_dt = -1;
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD); // wait until all the nodes are ready
      MPI_Allreduce(&my_dt, &dt, 1, kMPIDatatypeReal,
		    MPI_MIN, MPI_COMM_WORLD);
#else
      dt = my_dt;
#endif
    }
    if (dt < 0) {
      timeout = true;
      goal_reached = false;
      break;
    }

    if (local_setup().mpi_rank()==0) std::cerr << current_time() << " " << dt << std::endl;
    
    const PreciseReal next_time = current_time() + dt;
    if (next_time > goal_time) {
      dt = (goal_time - current_time());
      break_flag = true; goal_reached = true;
    }



    //// Zero-Initialize Pace
    clear_pace_kernel <<<32, 32>>>
    (density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Calculate Flux and add it to Pace
    add_flux_kernel_caller
      (dt/2,
       density_, velocity_, pressure_, magnet_,
       density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Convert to Conserved Variables
    to_conserved_kernel <<<32, 32>>>
    (density_, velocity_, pressure_, magnet_,
     momentum_, energy_);
    
    //// Add to Conserved Variables
    add_pace_kernel <<<32, 32>>>
    (dt/2,
     density_pace_, momentum_pace_, energy_pace_, magnet_pace_,
     density_, magnet_,
     density_predictor_, momentum_, energy_, magnet_predictor_);
    //// Convert to Primitive Variables
    to_primitive_kernel <<<32, 32>>>
    (density_predictor_, momentum_, energy_, magnet_predictor_,
     velocity_predictor_, pressure_predictor_);
    external_acceleration_kernel <<<32, 32>>>
    (dt/2,
     density_predictor_, velocity_predictor_, pressure_predictor_, magnet_predictor_);

    //// Zero-Initialize Pace
    clear_pace_kernel <<<32, 32>>>
    (density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Calculate Flux and add it to Pace
    add_flux_kernel_caller
      (dt,
       density_predictor_, velocity_predictor_, pressure_predictor_, magnet_predictor_,
       density_pace_, momentum_pace_, energy_pace_, magnet_pace_);
    
    //// Convert to Conserved Variables
    to_conserved_kernel <<<32, 32>>>
    (density_, velocity_, pressure_, magnet_,
     momentum_, energy_);
    
    //// Overwrite Conserved Variables
    update_kernel <<<32, 32>>>
    (dt,
     density_pace_, momentum_pace_, energy_pace_, magnet_pace_,
     density_, momentum_, energy_, magnet_);
    //// Convert to Primitive Variables
    to_primitive_kernel <<<32, 32>>>
    (density_, momentum_, energy_, magnet_,
     velocity_, pressure_);
    external_acceleration_kernel <<<32, 32>>>
    (dt,
     density_, velocity_, pressure_, magnet_);
      
    global_setup().proceed(1, dt);

    update_displacement();    

    
  } while(!break_flag);

  return goal_reached;
}
bool MHDSolver::read(std::string filename)  {
  FILE *fp=fopen(filename.c_str(), "r");
  if (!fp) return false;
  
  if (!fread_singlet(fp, kSizeX)) return false;
  if (!fread_singlet(fp, kSizeY)) return false;
  if (!fread_singlet(fp, kSizeZ)) return false;
  
  fread_device(fp,buffer_dict(density_));
  fread_device(fp,buffer_dict(velocity_[X]));
  fread_device(fp,buffer_dict(velocity_[Y]));
  fread_device(fp,buffer_dict(velocity_[Z]));
  fread_device(fp,buffer_dict(pressure_));
  fread_device(fp,buffer_dict(magnet_[X]));
  fread_device(fp,buffer_dict(magnet_[Y]));
  fread_device(fp,buffer_dict(magnet_[Z]));
  fclose(fp);
  return true;
}
bool MHDSolver::write(std::string filename)  {
  FILE *fp=fopen(filename.c_str(), "w");
  if (!fp) return false;
  
  if (!fwrite_singlet(fp, kSizeX)) return false;
  if (!fwrite_singlet(fp, kSizeY)) return false;
  if (!fwrite_singlet(fp, kSizeZ)) return false;
  
  fwrite_device(fp,buffer_dict(density_));
  fwrite_device(fp,buffer_dict(velocity_[X]));
  fwrite_device(fp,buffer_dict(velocity_[Y]));
  fwrite_device(fp,buffer_dict(velocity_[Z]));
  fwrite_device(fp,buffer_dict(pressure_));
  fwrite_device(fp,buffer_dict(magnet_[X]));
  fwrite_device(fp,buffer_dict(magnet_[Y]));
  fwrite_device(fp,buffer_dict(magnet_[Z]));
  fclose(fp);
  return true;
}

void MHDSolver::initialize_state () {
  initial_condition();
}

void MHDSolver::save_state (std::string tag) {
  std::string prefix = local_setup().mkdir_prefix(tag);
  write(prefix + ".bin");

  update_displacement();
  std::ofstream ofs((prefix + "_setup.txt").c_str());
  ofs << global_setup() << std::endl
      << "################" << std::endl
      << local_setup() << std::endl;
}

void MHDSolver::load_state (std::string tag) {
  std::string prefix = local_setup().mkdir_prefix(tag);

  std::ifstream ifs((prefix + "_setup.txt").c_str());
  ifs >> global_setup_;
  std::string spacer; ifs >> spacer;
  ifs >> local_setup_;

  if(!read(prefix + ".bin")) {
    local_setup().ofs_log() << "error reading file: " << prefix + ".bin" << std::endl;
  }
}
