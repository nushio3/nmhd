#include <gtest/gtest.h>
#include <vector>
using namespace std;

const int kN = 100;

void in_cpu (std::vector<int> &xs, std::vector<int> &ys);
void in_gpu (std::vector<int> &xs, std::vector<int> &ys);

TEST(thrust_vector, consistency) {
  std::vector<int> xs(kN);
  std::vector<int> ys1(kN);
  std::vector<int> ys2(kN);
  
  for (int i = 0; i < kN; ++i) {
    xs[i] = i*i*i+1;
  }
  
  in_cpu(xs, ys1);
  in_gpu(xs, ys2);
  for (int i = 0; i < kN; ++i) {
    EXPECT_EQ(ys1[i], ys2[i]);
  }
}
