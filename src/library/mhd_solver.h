
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>


#pragma once
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "simulation_config.h"
#include "library/thrust_vector.h"
#include "library/device_mesh.h"
#include "library/physics_util.h"

#ifdef USE_DOUBLE_PRECISION
const std::string kPrecisionTag = "double";
#else
const std::string kPrecisionTag = "float";
#endif

#ifdef USE_MPI
#include <mpi.h>
#ifdef USE_DOUBLE_PRECISION
const MPI_Datatype kMPIDatatypeReal = MPI_DOUBLE_PRECISION;
#else
const MPI_Datatype kMPIDatatypeReal = MPI_REAL;
#endif
#endif
class GlobalSetup {
 protected:
int generation_;
int real_size_;
std::string real_precision_;
PreciseReal current_time_;
int parallel_size_;
int global_size_x_;
int parallel_size_x_;
Real global_extent_x_;
Real global_corner_x_;
Real dX_;
int global_size_y_;
int parallel_size_y_;
Real global_extent_y_;
Real global_corner_y_;
Real dY_;
int global_size_z_;
int parallel_size_z_;
Real global_extent_z_;
Real global_corner_z_;
Real dZ_;
 public:
  const int& generation () const { return generation_; }
  int& generation () { return generation_; }
  const int& real_size () const { return real_size_; }
  const std::string& real_precision () const { return real_precision_; }
  const PreciseReal& current_time () const { return current_time_; }
  PreciseReal& current_time () { return current_time_; }
  const int& parallel_size () const { return parallel_size_; }
  const int& global_size_x () const { return global_size_x_; }
  const int& parallel_size_x () const { return parallel_size_x_; }
  const Real& global_extent_x () const { return global_extent_x_; }
  const Real& global_corner_x () const { return global_corner_x_; }
  const Real& dX () const { return dX_; }
  const int& global_size_y () const { return global_size_y_; }
  const int& parallel_size_y () const { return parallel_size_y_; }
  const Real& global_extent_y () const { return global_extent_y_; }
  const Real& global_corner_y () const { return global_corner_y_; }
  const Real& dY () const { return dY_; }
  const int& global_size_z () const { return global_size_z_; }
  const int& parallel_size_z () const { return parallel_size_z_; }
  const Real& global_extent_z () const { return global_extent_z_; }
  const Real& global_corner_z () const { return global_corner_z_; }
  const Real& dZ () const { return dZ_; }
  

  
  GlobalSetup () :
generation_(0),real_size_(sizeof(Real)),real_precision_(kPrecisionTag),current_time_(0),parallel_size_(kParallelSize),global_size_x_(kGlobalSizeX),parallel_size_x_(kParallelSizeX),global_extent_x_(kGlobalExtentX),global_corner_x_(kGlobalCornerX),dX_(kDX),global_size_y_(kGlobalSizeY),parallel_size_y_(kParallelSizeY),global_extent_y_(kGlobalExtentY),global_corner_y_(kGlobalCornerY),dY_(kDY),global_size_z_(kGlobalSizeZ),parallel_size_z_(kParallelSizeZ),global_extent_z_(kGlobalExtentZ),global_corner_z_(kGlobalCornerZ),dZ_(kDZ)  {}

  
  bool operator==(const GlobalSetup &other) const {
    if (generation() != other.generation()) return false;
    if (real_size() != other.real_size()) return false;
    if (real_precision() != other.real_precision()) return false;
    if (current_time() != other.current_time()) return false;
    if (parallel_size() != other.parallel_size()) return false;
    if (global_size_x() != other.global_size_x()) return false;
    if (parallel_size_x() != other.parallel_size_x()) return false;
    if (global_extent_x() != other.global_extent_x()) return false;
    if (global_corner_x() != other.global_corner_x()) return false;
    if (dX() != other.dX()) return false;
    if (global_size_y() != other.global_size_y()) return false;
    if (parallel_size_y() != other.parallel_size_y()) return false;
    if (global_extent_y() != other.global_extent_y()) return false;
    if (global_corner_y() != other.global_corner_y()) return false;
    if (dY() != other.dY()) return false;
    if (global_size_z() != other.global_size_z()) return false;
    if (parallel_size_z() != other.parallel_size_z()) return false;
    if (global_extent_z() != other.global_extent_z()) return false;
    if (global_corner_z() != other.global_corner_z()) return false;
    if (dZ() != other.dZ()) return false;
    return true;
  }


  std::string encode_member (const int &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, int &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const double &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, double &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const float &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, float &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const size_t &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, size_t &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  
  std::string encode_member (const std::string &x) const {
    std::ostringstream ret;
    char buf[3];
    for (int i = 0; i < x.size(); ++i) {
      unsigned char c = x[i];
      if (33<=c && c <= 126 && c!='%') {
	ret << c;
      } else {
	sprintf(buf, "%02x", c);
	ret << '%' << buf; 
      }
    }
    return ret.str();
  }
  void decode_member (const std::string &x, std::string &ret0) const {
    std::ostringstream ret;
    char buf[3]; buf[2] = 0;

    for (int i = 0; i < x.size(); ++i) {
      char c = x[i];
      if (c!='%') {
	ret << c;
      } else {
	int ic;
	buf[0] = x[++i]; buf[1] = x[++i];
	sscanf(buf, "%02x", &ic);
	ret << char(ic);
      }
    }
    ret0 = ret.str();
  }

  
void proceed (int generation_increase, PreciseReal time_increase) {
  generation() += generation_increase;
  current_time() += time_increase;
}

};


namespace {
std::ostream& operator<<(std::ostream& ostr, const GlobalSetup & x) {

  ostr << 20 << std::endl;
  ostr << "generation" << " " << x.encode_member(x.generation()) << std::endl;
  ostr << "real_size" << " " << x.encode_member(x.real_size()) << std::endl;
  ostr << "real_precision" << " " << x.encode_member(x.real_precision()) << std::endl;
  ostr << "current_time" << " " << x.encode_member(x.current_time()) << std::endl;
  ostr << "parallel_size" << " " << x.encode_member(x.parallel_size()) << std::endl;
  ostr << "global_size_x" << " " << x.encode_member(x.global_size_x()) << std::endl;
  ostr << "parallel_size_x" << " " << x.encode_member(x.parallel_size_x()) << std::endl;
  ostr << "global_extent_x" << " " << x.encode_member(x.global_extent_x()) << std::endl;
  ostr << "global_corner_x" << " " << x.encode_member(x.global_corner_x()) << std::endl;
  ostr << "dX" << " " << x.encode_member(x.dX()) << std::endl;
  ostr << "global_size_y" << " " << x.encode_member(x.global_size_y()) << std::endl;
  ostr << "parallel_size_y" << " " << x.encode_member(x.parallel_size_y()) << std::endl;
  ostr << "global_extent_y" << " " << x.encode_member(x.global_extent_y()) << std::endl;
  ostr << "global_corner_y" << " " << x.encode_member(x.global_corner_y()) << std::endl;
  ostr << "dY" << " " << x.encode_member(x.dY()) << std::endl;
  ostr << "global_size_z" << " " << x.encode_member(x.global_size_z()) << std::endl;
  ostr << "parallel_size_z" << " " << x.encode_member(x.parallel_size_z()) << std::endl;
  ostr << "global_extent_z" << " " << x.encode_member(x.global_extent_z()) << std::endl;
  ostr << "global_corner_z" << " " << x.encode_member(x.global_corner_z()) << std::endl;
  ostr << "dZ" << " " << x.encode_member(x.dZ()) << std::endl;
  return ostr;
}

std::istream& operator>>(std::istream& istr, GlobalSetup & x) {
  int n; istr >> n;
  std::map<std::string, std::string> buf;

  for (int i = 0; i < n; ++i) {
    std::string key, val;
    istr >> key >> val;
    buf[key] = val;
  }
    
  if (buf.count("generation")) {
    std::istringstream iss(buf["generation"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.generation());
  }
  
  if (buf.count("current_time")) {
    std::istringstream iss(buf["current_time"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.current_time());
  }
  
  return istr;
}
}

class LocalSetup {
 protected:
int mpi_rank_;
std::string output_directory_;
std::string mpi_prefix_;
Triplet<int> mpi_rank_3d_;
Triplet<Real> computational_corner_;
Triplet<Real> computational_extent_;
PreciseReal relative_displacement_;
PreciseReal relative_velocity_;
std::ofstream ofs_log_;
std::ofstream ofs_bench_;
size_t computational_size_;
int computational_size_x_;
int margin_x_;
Real computational_corner_x_;
Real computational_extent_x_;
Real dX_;
int computational_size_y_;
int margin_y_;
Real computational_corner_y_;
Real computational_extent_y_;
Real dY_;
int computational_size_z_;
int margin_z_;
Real computational_corner_z_;
Real computational_extent_z_;
Real dZ_;
 public:
  const int& mpi_rank () const { return mpi_rank_; }
  int& mpi_rank () { return mpi_rank_; }
  const std::string& output_directory () const { return output_directory_; }
  std::string& output_directory () { return output_directory_; }
  const std::string& mpi_prefix () const { return mpi_prefix_; }
  std::string& mpi_prefix () { return mpi_prefix_; }
  const Triplet<int>& mpi_rank_3d () const { return mpi_rank_3d_; }
  Triplet<int>& mpi_rank_3d () { return mpi_rank_3d_; }
  const Triplet<Real>& computational_corner () const { return computational_corner_; }
  Triplet<Real>& computational_corner () { return computational_corner_; }
  const Triplet<Real>& computational_extent () const { return computational_extent_; }
  Triplet<Real>& computational_extent () { return computational_extent_; }
  const PreciseReal& relative_displacement () const { return relative_displacement_; }
  PreciseReal& relative_displacement () { return relative_displacement_; }
  const PreciseReal& relative_velocity () const { return relative_velocity_; }
  PreciseReal& relative_velocity () { return relative_velocity_; }
  const std::ofstream& ofs_log () const { return ofs_log_; }
  std::ofstream& ofs_log () { return ofs_log_; }
  const std::ofstream& ofs_bench () const { return ofs_bench_; }
  std::ofstream& ofs_bench () { return ofs_bench_; }
  const size_t& computational_size () const { return computational_size_; }
  const int& computational_size_x () const { return computational_size_x_; }
  const int& margin_x () const { return margin_x_; }
  const Real& computational_corner_x () const { return computational_corner_x_; }
  const Real& computational_extent_x () const { return computational_extent_x_; }
  const Real& dX () const { return dX_; }
  const int& computational_size_y () const { return computational_size_y_; }
  const int& margin_y () const { return margin_y_; }
  const Real& computational_corner_y () const { return computational_corner_y_; }
  const Real& computational_extent_y () const { return computational_extent_y_; }
  const Real& dY () const { return dY_; }
  const int& computational_size_z () const { return computational_size_z_; }
  const int& margin_z () const { return margin_z_; }
  const Real& computational_corner_z () const { return computational_corner_z_; }
  const Real& computational_extent_z () const { return computational_extent_z_; }
  const Real& dZ () const { return dZ_; }
  

  
  LocalSetup () :
mpi_rank_(0),output_directory_(),mpi_prefix_(),computational_extent_(Triplet<Real>(kExtentX, kExtentY, kExtentZ)),relative_displacement_(0),relative_velocity_(0),ofs_log_(NULL),ofs_bench_(NULL),computational_size_(kSize),computational_size_x_(kSizeX),margin_x_(kMarginSizeX),computational_corner_x_(0),computational_extent_x_(kExtentX),dX_(kDX),computational_size_y_(kSizeY),margin_y_(kMarginSizeY),computational_corner_y_(0),computational_extent_y_(kExtentY),dY_(kDY),computational_size_z_(kSizeZ),margin_z_(kMarginSizeZ),computational_corner_z_(0),computational_extent_z_(kExtentZ),dZ_(kDZ)  {}

LocalSetup (int mpi_rank0, std::string output_dir0) :  mpi_rank_(0),output_directory_(),mpi_prefix_(),computational_extent_(Triplet<Real>(kExtentX, kExtentY, kExtentZ)),relative_displacement_(0),relative_velocity_(0),ofs_log_(NULL),ofs_bench_(NULL),computational_size_(kSize),computational_size_x_(kSizeX),margin_x_(kMarginSizeX),computational_corner_x_(0),computational_extent_x_(kExtentX),dX_(kDX),computational_size_y_(kSizeY),margin_y_(kMarginSizeY),computational_corner_y_(0),computational_extent_y_(kExtentY),dY_(kDY),computational_size_z_(kSizeZ),margin_z_(kMarginSizeZ),computational_corner_z_(0),computational_extent_z_(kExtentZ),dZ_(kDZ) {
  mpi_rank_ = mpi_rank0;
  int sx, sy, sz, a;
  sx = mpi_rank0 % kParallelSizeX; a = mpi_rank0 / kParallelSizeX;
  sy = a % kParallelSizeY;
  sz = a / kParallelSizeY;
  
  mpi_rank_3d() = Triplet<int>(sx, sy, sz);
  computational_corner_x_ = computational_corner()[X] 
			      = kGlobalCornerX + (kLocalSizeX * mpi_rank_3d()[X] - kMarginSizeX) * kDX;
  computational_corner_y_ = computational_corner()[Y] 
			      = kGlobalCornerY + (kLocalSizeY * mpi_rank_3d()[Y] - kMarginSizeY) * kDY;
  computational_corner_z_ = computational_corner()[Z] 
			      = kGlobalCornerZ + (kLocalSizeZ * mpi_rank_3d()[Z] - kMarginSizeZ) * kDZ;

  {
    char work[256]; std::sprintf(work, "%04d",  mpi_rank0);
    output_directory_ = output_dir0;
    mpi_prefix_ = work;
  }
  ofs_log().open(log_filename().c_str(), std::ios_base::out | std::ios_base::app);
  ofs_bench().open(bench_filename().c_str(), std::ios_base::out | std::ios_base::app);
}
  
  bool operator==(const LocalSetup &other) const {
    if (mpi_rank() != other.mpi_rank()) return false;
    if (output_directory() != other.output_directory()) return false;
    if (mpi_prefix() != other.mpi_prefix()) return false;
    if (mpi_rank_3d() != other.mpi_rank_3d()) return false;
    if (computational_corner() != other.computational_corner()) return false;
    if (computational_extent() != other.computational_extent()) return false;
    if (relative_displacement() != other.relative_displacement()) return false;
    if (relative_velocity() != other.relative_velocity()) return false;
    if (computational_size() != other.computational_size()) return false;
    if (computational_size_x() != other.computational_size_x()) return false;
    if (margin_x() != other.margin_x()) return false;
    if (computational_corner_x() != other.computational_corner_x()) return false;
    if (computational_extent_x() != other.computational_extent_x()) return false;
    if (dX() != other.dX()) return false;
    if (computational_size_y() != other.computational_size_y()) return false;
    if (margin_y() != other.margin_y()) return false;
    if (computational_corner_y() != other.computational_corner_y()) return false;
    if (computational_extent_y() != other.computational_extent_y()) return false;
    if (dY() != other.dY()) return false;
    if (computational_size_z() != other.computational_size_z()) return false;
    if (margin_z() != other.margin_z()) return false;
    if (computational_corner_z() != other.computational_corner_z()) return false;
    if (computational_extent_z() != other.computational_extent_z()) return false;
    if (dZ() != other.dZ()) return false;
    return true;
  }


  std::string encode_member (const int &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, int &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const double &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, double &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const float &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, float &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  std::string encode_member (const size_t &x) const {
    std::ostringstream ostr; ostr << x; return ostr.str();
  }
  void decode_member (const std::string &str, size_t &ret) const {
    std::istringstream istr(str); istr >> ret;
  }
  
  std::string encode_member (const std::string &x) const {
    std::ostringstream ret;
    char buf[3];
    for (int i = 0; i < x.size(); ++i) {
      unsigned char c = x[i];
      if (33<=c && c <= 126 && c!='%') {
	ret << c;
      } else {
	sprintf(buf, "%02x", c);
	ret << '%' << buf; 
      }
    }
    return ret.str();
  }
  void decode_member (const std::string &x, std::string &ret0) const {
    std::ostringstream ret;
    char buf[3]; buf[2] = 0;

    for (int i = 0; i < x.size(); ++i) {
      char c = x[i];
      if (c!='%') {
	ret << c;
      } else {
	int ic;
	buf[0] = x[++i]; buf[1] = x[++i];
	sscanf(buf, "%02x", &ic);
	ret << char(ic);
      }
    }
    ret0 = ret.str();
  }

  
PositioningSetup positioning_setup () const {
  return PositioningSetup(computational_corner(), mpi_rank_3d()[Y]);
}

std::string mkdir_prefix (std::string tag) const {
  std::string dir_name = output_directory() + "/" + tag;
  std::system(("mkdir -p " + dir_name).c_str());
  return dir_name + "/" + mpi_prefix();
}

std::string log_filename () const {
  return mkdir_prefix("log") + ".txt";
}
std::string bench_filename () const {
  return mkdir_prefix("bench") + ".txt";
}



};


namespace {
std::ostream& operator<<(std::ostream& ostr, const LocalSetup & x) {

  ostr << 21 << std::endl;
  ostr << "mpi_rank" << " " << x.encode_member(x.mpi_rank()) << std::endl;
  ostr << "output_directory" << " " << x.encode_member(x.output_directory()) << std::endl;
  ostr << "mpi_prefix" << " " << x.encode_member(x.mpi_prefix()) << std::endl;
  ostr << "relative_displacement" << " " << x.encode_member(x.relative_displacement()) << std::endl;
  ostr << "relative_velocity" << " " << x.encode_member(x.relative_velocity()) << std::endl;
  ostr << "computational_size" << " " << x.encode_member(x.computational_size()) << std::endl;
  ostr << "computational_size_x" << " " << x.encode_member(x.computational_size_x()) << std::endl;
  ostr << "margin_x" << " " << x.encode_member(x.margin_x()) << std::endl;
  ostr << "computational_corner_x" << " " << x.encode_member(x.computational_corner_x()) << std::endl;
  ostr << "computational_extent_x" << " " << x.encode_member(x.computational_extent_x()) << std::endl;
  ostr << "dX" << " " << x.encode_member(x.dX()) << std::endl;
  ostr << "computational_size_y" << " " << x.encode_member(x.computational_size_y()) << std::endl;
  ostr << "margin_y" << " " << x.encode_member(x.margin_y()) << std::endl;
  ostr << "computational_corner_y" << " " << x.encode_member(x.computational_corner_y()) << std::endl;
  ostr << "computational_extent_y" << " " << x.encode_member(x.computational_extent_y()) << std::endl;
  ostr << "dY" << " " << x.encode_member(x.dY()) << std::endl;
  ostr << "computational_size_z" << " " << x.encode_member(x.computational_size_z()) << std::endl;
  ostr << "margin_z" << " " << x.encode_member(x.margin_z()) << std::endl;
  ostr << "computational_corner_z" << " " << x.encode_member(x.computational_corner_z()) << std::endl;
  ostr << "computational_extent_z" << " " << x.encode_member(x.computational_extent_z()) << std::endl;
  ostr << "dZ" << " " << x.encode_member(x.dZ()) << std::endl;
  return ostr;
}

std::istream& operator>>(std::istream& istr, LocalSetup & x) {
  int n; istr >> n;
  std::map<std::string, std::string> buf;

  for (int i = 0; i < n; ++i) {
    std::string key, val;
    istr >> key >> val;
    buf[key] = val;
  }
    
  if (buf.count("mpi_rank")) {
    std::istringstream iss(buf["mpi_rank"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.mpi_rank());
  }
  
  if (buf.count("output_directory")) {
    std::istringstream iss(buf["output_directory"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.output_directory());
  }
  
  if (buf.count("mpi_prefix")) {
    std::istringstream iss(buf["mpi_prefix"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.mpi_prefix());
  }
  
  if (buf.count("relative_displacement")) {
    std::istringstream iss(buf["relative_displacement"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.relative_displacement());
  }
  
  if (buf.count("relative_velocity")) {
    std::istringstream iss(buf["relative_velocity"]);
    std::string tmp;
    iss >> tmp; x.decode_member(tmp, x.relative_velocity());
  }
  
  return istr;
}
}

namespace {
LocalSetup blancLocalSetup;
}


class MHDSolver {
protected:
  GlobalSetup global_setup_;
  LocalSetup &local_setup_;

  // device memory allocated for main calculation
  std::vector<thrust::device_vector<Real> > buffer_;
  std::map<Real*, int> buffer_dict_;

  // aliases for buffer
  BodyScalar density_;
  BodyVector velocity_;
  BodyScalar pressure_;
  FaceVector magnet_;

  BodyScalar energy_;
  BodyVector momentum_;

  BodyScalar cfl_time_;
  
  BodyScalar density_pace_;
  BodyVector momentum_pace_;
  BodyScalar energy_pace_;
  FaceVector magnet_pace_;

  BodyScalar density_flux_;
  BodyVector momentum_flux_;
  BodyScalar energy_flux_;
  BodyScalar magnet_flux_ey_;
  BodyScalar magnet_flux_ez_;

  BodyScalar density_predictor_;
  BodyVector velocity_predictor_;
  BodyScalar pressure_predictor_;
  FaceVector magnet_predictor_;

  void initial_condition ();
  void boundary_condition ();
#ifdef USE_MPI
  void communicate ();
#endif
  
  // buffer for mpi communication
  bool mpi_buffer_initialized;
  thrust::host_vector<Real> mpi_buffer_y_inc_send;
  thrust::host_vector<Real> mpi_buffer_y_inc_recv;
  thrust::host_vector<Real> mpi_buffer_y_dec_send;
  thrust::host_vector<Real> mpi_buffer_y_dec_recv;

  thrust::device_vector<Real> mpi_buffer_y_inc_device;
  thrust::device_vector<Real> mpi_buffer_y_dec_device;
  
  thrust::host_vector<Real> mpi_buffer_z_inc_send;
  thrust::host_vector<Real> mpi_buffer_z_inc_recv;
  thrust::host_vector<Real> mpi_buffer_z_dec_send;
  thrust::host_vector<Real> mpi_buffer_z_dec_recv;
  
  template <int Bit>
  thrust::device_vector<Real> &buffer_dict(Mesh<Bit, Real> mesh) {
    return buffer_[buffer_dict_[mesh.ptr_]];
  }
  
  void add_flux_kernel_caller
  (Real dt, BodyScalar density, BodyVector velocity, BodyScalar pressure, FaceVector magnet,
   BodyScalar density_pace, BodyVector momentum_pace, BodyScalar energy_pace, FaceVector magnet_pace) ;

  thrust::host_vector<Real> host_buffer_;
  void fread_device(FILE *fp, thrust::device_vector<Real> &xs) {
    std::fread(&host_buffer_[0], sizeof(Real), kSize, fp);  
    xs = host_buffer_;
  }
  void fwrite_device(FILE *fp, thrust::device_vector<Real> &xs) {
    host_buffer_ = xs;
    std::fwrite(&host_buffer_[0], sizeof(Real), kSize, fp);  
  }
  template<class T>
  bool fread_singlet(FILE *fp, const T x) {
    T y;
    std::fread(&y, sizeof(T), 1, fp);
    return x==y;
  }
  template<class T>
  bool fwrite_singlet(FILE *fp, const T x) {
    std::fwrite(&x, sizeof(T), 1, fp);
    return true;
  }

public:
  MHDSolver(LocalSetup &local_setup = blancLocalSetup);

  GlobalSetup &global_setup()
  { return global_setup_; }
  const GlobalSetup &global_setup() const
  { return global_setup_; }
  LocalSetup &local_setup()
  { return local_setup_; }
  const LocalSetup &local_setup() const
  { return local_setup_; }
  template<class T>
  Triplet<Real> positioning (const T& positioner) const {
    return positioner.position(local_setup().positioning_setup());
  }
  
  PreciseReal current_time () { return global_setup().current_time(); }
  int generation () { return global_setup().generation(); }

  // try to integrate until the goal time,
  // returns true if goal_time is reached,
  // exits and returns false if the wallclock time is reached earlier.
  // ignores wallclock_time_limit if it's negative.
  bool proceed (const Real goal_time, const PreciseReal wallclock_time_limit = -1);
  
  void dump (std::string);
  bool read (std::string); // returns true if successful
  bool write(std::string);

  void initialize_state ();
  void load_state (std::string);
  void save_state (std::string);
  
  int monitor (); // called every timestep. return non-zero to break
  int initial_event (); // called just before calculation starts 
  int final_event ();  // called after all calculation ends

  void relative_displacement (PreciseReal t, int index, Real &normalized_displacement, int &int_displacement, Real &minor_displacement);
  void update_displacement (); 
};
