#include <thrust/device_vector.h>
#include <iostream>
using namespace std;

#include "conster.h"

#include "conster.inl"
int main () {
  cout << kX << endl;
  thrust::device_vector<Real> ret(1);
  kernel<<< 1, 1 >>> (thrust::raw_pointer_cast(&*ret.begin()));
  cout << ret[0] << endl;  
}
