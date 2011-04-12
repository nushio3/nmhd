#pragma once
#include <iostream>
struct XAxis;
struct XHalf;
struct YAxis;
struct YHalf;
struct ZAxis;
struct ZHalf;
struct XAxis
{
  enum { mask = 1 };
  typedef XHalf Half;
  typedef ZAxis Prev;
  typedef YAxis Next;
};
struct XHalf
{
  enum { mask = 1 };
  typedef XAxis Whole;
  typedef ZHalf Prev;
  typedef YHalf Next;
};
__device__ XAxis X;
__device__ XHalf HX;
struct YAxis
{
  enum { mask = 2 };
  typedef YHalf Half;
  typedef XAxis Prev;
  typedef ZAxis Next;
};
struct YHalf
{
  enum { mask = 2 };
  typedef YAxis Whole;
  typedef XHalf Prev;
  typedef ZHalf Next;
};
__device__ YAxis Y;
__device__ YHalf HY;
struct ZAxis
{
  enum { mask = 4 };
  typedef ZHalf Half;
  typedef YAxis Prev;
  typedef XAxis Next;
};
struct ZHalf
{
  enum { mask = 4 };
  typedef ZAxis Whole;
  typedef YHalf Prev;
  typedef XHalf Next;
};
__device__ ZAxis Z;
__device__ ZHalf HZ;

namespace direction{
template<class T>
__device__ __host__ typename T::Half half(T t) { typename T::Half ret; return ret; }
template<class T>
__device__ __host__ typename T::Whole whole(T t) { typename T::Whole ret; return ret; }
template<class T>
__device__ __host__ typename T::Prev prev(T t) { typename T::Prev ret; return ret; }
template<class T>
__device__ __host__ typename T::Next next(T t) { typename T::Next ret; return ret; }
}

/*
  Triplet : tuple of three object of same type
  (essentially, a three dimensional vector)  */

template <class T>
class Triplet {
public:
  T x_, y_, z_;
  __device__ __host__ Triplet () : x_(), y_(), z_() {}
  __device__ __host__ Triplet (const T& x,const T& y,const T& z) :
    x_(x), y_(y), z_(z) {}

  __device__ __host__
  T &operator[](XAxis) { return x_; }
  __device__ __host__
  const T &operator[](XAxis) const { return x_; }
  __device__ __host__
  T &operator[](YAxis) { return y_; }
  __device__ __host__
  const T &operator[](YAxis) const { return y_; }
  __device__ __host__
  T &operator[](ZAxis) { return z_; }
  __device__ __host__
  const T &operator[](ZAxis) const { return z_; }

  // vector addition between triplets

  __device__ __host__
  Triplet<T> operator+ (const Triplet<T> &other) const {
    return Triplet<T>(x_ + other.x_, y_ + other.y_, z_ + other.z_);
  }

  __device__ __host__
  Triplet<T> operator+= (const Triplet<T> &other) {
    x_ += other.x_; y_ += other.y_; z_ += other.z_;
    return *this;
  }

  __device__ __host__
  Triplet<T> operator- (const Triplet<T> &other) const {
    return Triplet<T>(x_ - other.x_, y_ - other.y_, z_ - other.z_);
  }

  __device__ __host__
  Triplet<T> operator-= (const Triplet<T> &other) {
    x_ -= other.x_; y_ -= other.y_; z_ -= other.z_;
    return *this;
  }

  // comparison
  __device__ __host__
  bool operator== (const Triplet<T> &other) const {
    return x_ == other.x_ && y_ == other.y_ && z_ == other.z_;
  }

  __device__ __host__
  bool operator!= (const Triplet<T> &other) const {
    return x_ != other.x_ || y_ != other.y_ || z_ != other.z_;
  }
};

// inner and outer product
template<class T>
__device__ __host__ T inner_prod (const Triplet<T> &a, const Triplet<T> &b) {
  return a.x_ * b.x_ +  a.y_ * b.y_ +  a.z_ * b.z_;
}
template<class T>
__device__ __host__ Triplet<T> outer_prod (const Triplet<T> &a, const Triplet<T> &b) {
  return Triplet<T>
    (
     a[Y] * b[Z] - a[Z] * b[Y]
,     a[Z] * b[X] - a[X] * b[Z]
,     a[X] * b[Y] - a[Y] * b[X]
);
}

// absolute and norm
template<class T>
__device__ __host__ T norm (const Triplet<T> &a) {
  return inner_prod(a,a);
}

namespace { // triplet function and operators
template<class T>
__device__ __host__ T abs (const Triplet<T> &a) {
  return sqrt(norm(a));
}
template<>
__device__ __host__ float abs (const Triplet<float> &a) {
  return sqrtf(norm(a));
}
}


template <class T>
std::ostream& operator<<(std::ostream& cout, Triplet<T> t) {
  return cout << t.x_ << " " << t.y_ << " " << t.z_;
}
template <class T>
std::istream& operator>>(std::istream& cin, Triplet<T> t) {
  return cin >> t.x_ >> t.y_ >> t.z_;
}




/*
  Triad : tuple of three object of different type */
 

template <class TX, class TY, class TZ>
class Triad {
public:
  TX x_; TY y_; TZ z_;
  __device__ __host__ Triad () : x_(), y_(), z_() {}
  __device__ __host__ Triad (const TX& x,const TY& y,const TZ& z) :
    x_(x), y_(y), z_(z) {}
  
  __device__ __host__
  TX &operator[](XAxis) { return x_; }
  __device__ __host__
  const TX &operator[](XAxis) const { return x_; }
  __device__ __host__
  TY &operator[](YAxis) { return y_; }
  __device__ __host__
  const TY &operator[](YAxis) const { return y_; }
  __device__ __host__
  TZ &operator[](ZAxis) { return z_; }
  __device__ __host__
  const TZ &operator[](ZAxis) const { return z_; }

  // vector addition between triplets

  __device__ __host__
  Triad<TX,TY,TZ> operator+ (const Triad<TX,TY,TZ> &other) const {
    return Triad<TX,TY,TZ>(x_ + other.x_, y_ + other.y_, z_ + other.z_);
  }

  __device__ __host__
  Triad<TX,TY,TZ> operator+= (const Triad<TX,TY,TZ> &other) {
    x_ += other.x_; y_ += other.y_; z_ += other.z_;
    return *this;
  }

  __device__ __host__
  Triad<TX,TY,TZ> operator- (const Triad<TX,TY,TZ> &other) const {
    return Triad<TX,TY,TZ>(x_ - other.x_, y_ - other.y_, z_ - other.z_);
  }

  __device__ __host__
  Triad<TX,TY,TZ> operator-= (const Triad<TX,TY,TZ> &other) {
    x_ -= other.x_; y_ -= other.y_; z_ -= other.z_;
    return *this;
  }

  // comparison
  __device__ __host__
  bool operator== (const Triad<TX,TY,TZ> &other) const {
    return x_ == other.x_ && y_ == other.y_ && z_ == other.z_;
  }

  __device__ __host__
  bool operator!= (const Triad<TX,TY,TZ> &other) const {
    return x_ != other.x_ || y_ != other.y_ || z_ != other.z_;
  }
};


template <class TX, class TY, class TZ>
std::ostream& operator<<(std::ostream& cout, Triad<TX,TY,TZ> t) {
  return cout << t.x_ << " " << t.y_ << " " << t.z_;
}
template <class TX, class TY, class TZ>
std::istream& operator>>(std::istream& cin, Triad<TX,TY,TZ> t) {
  return cin >> t.x_ >> t.y_ >> t.z_;
}

