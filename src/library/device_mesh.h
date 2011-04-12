#pragma once
#include "library/crystal.h"
#include "library/thrust_vector.h"

// a pair of device vector and reduced-coordinat flavored pointer to it
template <int Bit, class T>
class DeviceMesh {
public:
  typedef typename thrust::device_vector<T>::size_type size_type;
  typedef typename thrust::device_vector<T> vector_type;
  typedef Mesh<Bit, T> mesh_type;
private:
  vector_type device_;
  vector_type *p_device_;
  mesh_type mesh_;
public:
  DeviceMesh (size_type size) : // reserve your own device memory
    device_(size), p_device_(&device_), mesh_(raw(device_)) {}
  DeviceMesh (vector_type &vec) : // borrow the device memory
    device_(0), p_device_(&vec), mesh_(raw(vec)) {}


  vector_type &device() { return *p_device_; }
  const mesh_type &mesh() const { return mesh_; }
};
