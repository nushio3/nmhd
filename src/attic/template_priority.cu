#include <iostream>
using namespace std;

template<class T> void f(T x) {
  cout << "generic : " << x << endl;
}
template<> void f(int x){
  cout << "int : " << x << endl;
}

int main () {
  f(1); f(2.3);
}
