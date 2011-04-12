#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

typedef float Real;
typedef double PreciseReal;
typedef map<string,string> Hash;
Hash global_hash, local_hash;

template <class T>
T read (const string &str) {
  T ret;
  istringstream iss(str);
  iss >> ret;
  return ret;
}

template<class T>
T fread_singlet(FILE *fp) {
  T y;
  std::fread(&y, sizeof(T), 1, fp);
  return y;
}

template<class T>
void fwrite_singlet(FILE *fp, const T x) {
  std::fwrite(&x, sizeof(T), 1, fp);
}

template <class T>
T read_hash (const string &key) {
  if (local_hash.count(key) > 0) {
    return read<T>(local_hash[key]);
  }
  if (global_hash.count(key) > 0) {
    return read<T>(global_hash[key]);
  }
  cerr << "hash key not found: " << key << endl;
  exit(-1);
}

istream& operator>> (istream& istr, Hash &h) {
  int n; istr >> n;
  for (int i = 0; i < n; ++i) {
    string k,v; istr >> k >> v;
    h[k] = v;
  }
  return istr;
}

ostream& operator<< (ostream& ostr, Hash &h) {
  for (Hash::iterator it = h.begin(); it != h.end(); ++it) {
    ostr << it->first << " " << it->second << endl;
  }
  return ostr;
}

const int max_resolution = 120;

struct Accessor {
  int size_x, size_y, size_z;
  Real corner_x, corner_y, corner_z;
  Real extent_x, extent_y, extent_z;
  size_t size;

  
  size_t encode (int x0, int y0, int z0, int &carry) {
    int x = ((x0 % size_x) + size_x) % size_x;
    int y = ((y0 % size_y) + size_y) % size_y;
    int z = ((z0 % size_z) + size_z) % size_z;
    carry = (y-y0)/size_y;
    return (size_t(z) * size_y + y) * size_x + x;
  }
  size_t encode_real (Real rx, Real ry, Real rz, int &carry) {
    int x = int(floor((rx-corner_x)/extent_x*size_x));
    int y = int(floor((ry-corner_y)/extent_y*size_y));
    int z = int(floor((rz-corner_z)/extent_z*size_z));
    return encode(x,y,z,carry);
  }
  
  size_t decode (size_t addr, int &x, int &y, int &z) {
    x = addr % size_x;
    size_t addr2 = addr/ size_x;
    y = addr2 % size_y;
    z = addr2 / size_y;
  }
  size_t decode_real (size_t addr, Real &rx, Real &ry, Real &rz) {
    int x,y,z;
    decode(addr, x, y, z);
    rx = corner_x + (Real(x)+Real(0.5)) / size_x * extent_x;
    ry = corner_y + (Real(y)+Real(0.5)) / size_y * extent_y;
    rz = corner_z + (Real(z)+Real(0.5)) / size_z * extent_z;
  }
};

template <class T>
struct Accumulator {
  vector<T> numerators, denominators;
  Accumulator (size_t sz) :
    numerators(sz,0), denominators(sz,0) {}
  void add (size_t addr, const T &val) {
    numerators[addr] += val;
    denominators[addr] += 1;
  }
  T get (size_t addr) {
    return numerators[addr] / denominators[addr];
  }
};

const int item_size = 8;

int main (int argc, char **argv) {
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " folder/of/snapshot output/folder/filename_prefix" << endl;
    return -1;
  }
  const string ifn_prefix = string(argv[1]) + "/";
  const string ofn_prefix = string(argv[2]);
  {
    string setup0_fn = ifn_prefix + "0000_setup.txt";
    
    ifstream ifs(setup0_fn.c_str());
    if(!ifs) {
      cerr << "cannot open: " << setup0_fn << endl;
      return -1;
    }
    ifs >> global_hash;
    string spacer; ifs >> spacer;
    ifs >> local_hash;

    assert(sizeof(Real) == read_hash<int>("real_size"));
    
  }

  Accessor output;
  output.size_x = min(max_resolution, read_hash<int>("global_size_x"));
  output.size_y = min(max_resolution, read_hash<int>("global_size_y"));
  output.size_z = min(max_resolution, read_hash<int>("global_size_z"));
  output.size = size_t(output.size_x) * size_t(output.size_y) * size_t(output.size_z);

  output.corner_x = read_hash<Real>("global_corner_x");
  output.corner_y = read_hash<Real>("global_corner_y");
  output.corner_z = read_hash<Real>("global_corner_z");
  output.extent_x = read_hash<Real>("global_extent_x");
  output.extent_y = read_hash<Real>("global_extent_y");
  output.extent_z = read_hash<Real>("global_extent_z");


  vector<Accumulator<PreciseReal> > accums(item_size,
					   Accumulator<PreciseReal>(output.size));

  const int parallel_size_x = read_hash<int>("parallel_size_x");
  const int parallel_size_y = read_hash<int>("parallel_size_y");
  const int parallel_size_z = read_hash<int>("parallel_size_z");
  vector<PreciseReal> integrated_displacement(parallel_size_y);
  vector<PreciseReal> integrated_velocity(parallel_size_y);
  
  const int parallel_size = read_hash<int>("parallel_size");
  for (int mpi_rank = 0; mpi_rank < parallel_size; ++mpi_rank) {
    { // read local setting for the mpi_rank
      char tmp[256]; sprintf(tmp, "%04d", mpi_rank);
      string setup_fn = ifn_prefix + tmp + "_setup.txt";
      ifstream ifs(setup_fn.c_str());
      if(!ifs) {
	cerr << "cannot open: " << setup_fn << endl;
	return -1;
      }
      Hash global_hash1;
      ifs >> global_hash1;
      if (global_hash != global_hash1) {
	cerr << "global setup differes: " << setup_fn << endl;
      return -1;
      }
      string spacer; ifs >> spacer;
      ifs >> local_hash;
    }
    
    Accessor input;
    input.size_x = read_hash<int>("computational_size_x");
    input.size_y = read_hash<int>("computational_size_y");
    input.size_z = read_hash<int>("computational_size_z");
    input.size = size_t(input.size_x) * size_t(input.size_y) * size_t(input.size_z);
    assert(input.size == read_hash<int>("computational_size"));
 
    
    input.corner_x = read_hash<Real>("computational_corner_x");
    input.corner_y = read_hash<Real>("computational_corner_y");
    input.corner_z = read_hash<Real>("computational_corner_z");
    input.extent_x = read_hash<Real>("computational_extent_x");
    input.extent_y = read_hash<Real>("computational_extent_y");
    input.extent_z = read_hash<Real>("computational_extent_z");

    const int mpi_rank = read_hash<int>("mpi_rank");
    const int mpi_rank_y = mpi_rank % parallel_size_y;

    if (mpi_rank_y < parallel_size_y-1) {
      integrated_displacement[mpi_rank_y+1] =
	integrated_displacement[mpi_rank_y]
	+ read_hash<Real>("relative_displacement");
      integrated_velocity[mpi_rank_y+1] =
	integrated_velocity[mpi_rank_y]
	+ read_hash<Real>("relative_velocity");
    }
    
    { // read binary dump for the mpi_rank
      char tmp[256]; sprintf(tmp, "%04d", mpi_rank);
      string setup_fn = ifn_prefix + tmp + ".bin";
      FILE *fp = fopen(setup_fn.c_str(), "r");
      int sx = fread_singlet<int>(fp);
      int sy = fread_singlet<int>(fp);
      int sz = fread_singlet<int>(fp);
      assert(sx == read_hash<int>("computational_size_x"));
      assert(sy == read_hash<int>("computational_size_y"));
      assert(sz == read_hash<int>("computational_size_z"));

      vector<Real> buf(input.size);
      for (int item = 0; item < item_size; ++item) {
	fread(&buf[0], sizeof(Real), input.size, fp);

	for (int input_addr = 0; input_addr < input.size; ++input_addr) {
	  Real x,y,z; size_t output_addr;
	  input.decode_real(input_addr, x, y, z);
	  x += integrated_displacement[mpi_rank_y];
	  Real boost = 0;
	  if(item==1) { //vx
	    boost = integrated_velocity[mpi_rank_y];
	  }
	  int carry = 0;
	  output_addr = output.encode_real(x, y, z, carry);
	  if(item==1) { //vx
	    boost = integrated_velocity[mpi_rank_y]
	      + parallel_size_y * read_hash<Real>("relative_velocity") * carry;
	  }
	  if (carry == 0) // hiyotta
	    accums[item].add(output_addr, buf[input_addr] + boost);
	}
      }
      
      fclose(fp);
    }
  } // end of an mpi_rank

  // output cross sections of the accumulated data
  for (int axis = 0; axis < 3; ++axis) {
    ostringstream oss_fn; 
    
    int x0 = 0, x1 = output.size_x;
    int y0 = 0, y1 = output.size_y;
    int z0 = 0, z1 = output.size_z;
    
    string axis_str;
    switch (axis) {
    case 0:
      x0 = x1 = output.size_x/2; ++x1; axis_str = "yz";
      break;
    case 1:
      y0 = y1 = output.size_y/2; ++y1; axis_str = "xz";
      break;
    case 2:
      z0 = z1 = output.size_z/2; ++z1; axis_str = "xy";
      break;
    }
    oss_fn << ofn_prefix << "_" << axis_str << ".txt";
    
    ofstream ofs(oss_fn.str().c_str());
    
    for (int z = z0; z < z1; ++z) {
      for (int y = y0; y < y1; ++y) {
	for (int x = x0; x < x1; ++x) {
	  Real rx, ry, rz; int carry;
	  int addr = output.encode(x,y,z,carry);
	  output.decode_real(addr, rx, ry, rz);
	  ofs << rx << " " << ry << " " << rz;
	  for (int i = 0; i < item_size; ++i) {
	    ofs << " " << accums[i].get(addr);
	  }
	  if (x1-x0 > 1) ofs << endl;
	}
	if (y1-y0 > 1) ofs << endl;
      }
      if (z1-z0 > 1) ofs << endl;
    }
  }
}
