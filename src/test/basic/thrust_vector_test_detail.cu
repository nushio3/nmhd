#include "library/thrust_vector.h"
#include <vector>
using namespace std;

const int kMaxT = 1000;

__global__ void in_gpu_kernel (int *xs, int *ys) {
  int tid = threadIdx.x;
  int x = xs[tid];
  for (int t = 0; t < kMaxT; ++t) {
    x = x%2==0 ? x/2 : 3*x+1;
  }
  ys[tid] = x;
}

void in_gpu (vector<int> &xs, vector<int> &ys) {
  thrust::thrust_vector<int> txs(xs.size());
  thrust::thrust_vector<int> tys(xs.size());
  for (int i = 0; i < xs.size(); ++i) {
    txs[i] = xs[i];
  }
  in_gpu_kernel<<<1, xs.size()>>>(raw(txs), raw(tys));
  for (int i = 0; i < xs.size(); ++i) {
    ys[i] = tys[i];
  }
}

int in_cpu_kernel (int x) {
  for (int t = 0; t < kMaxT; ++t) {
    x = x%2==0 ? x/2 : 3*x+1;
  }  
  return x;
}

void in_cpu (std::vector<int> &xs, std::vector<int> &ys) {
  for (int i = 0; i < xs.size(); ++i) {
    ys[i] = in_cpu_kernel(xs[i]);
  }
}
