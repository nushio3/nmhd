#include <iostream>
using namespace std;

template<class T> struct Triplet {
  T x,y,z;
  Triplet(T i) : x(i), y(i), z(i) {}
};

Triplet<int> a = Triplet<int>(42);

int main () {
  cout << "hi " << a.z << endl;
  return 0;
}
