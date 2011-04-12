#include "library/mersenne.h"

#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
using namespace std;

TEST(mersenne, uniformness) {
  unsigned long seed = time(NULL);
  for (int t = 0; t < 10; ++t) {
    MersenneTwister gen(seed + t);
    double lo = gen.drand(0, 100); 
    double hi = gen.drand(lo + 1, 1000);
    vector<int> hist(100, 0);
    const int n = 10000;
    for (int i = 0; i < n * hist.size(); ++i) {
      double r = gen.drand(lo, hi);
      EXPECT_TRUE(lo <= r);
      EXPECT_TRUE(r < hi);
      int box = int((r - lo)/(hi - lo) * hist.size());
      EXPECT_TRUE(box >=0 && box < hist.size());
      hist[box]++;
    }
    
    for (int i = 0; i < hist.size(); ++i) {
      EXPECT_TRUE(hist[i] > n - 5.0 * sqrt(n));
      EXPECT_TRUE(hist[i] < n + 5.0 * sqrt(n));
    }
  }
}
