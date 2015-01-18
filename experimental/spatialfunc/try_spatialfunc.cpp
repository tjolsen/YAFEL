#include <iostream>
#include "yafel.hpp"

using namespace yafel;

template<typename T>
class SpatialFunc {

private:
  T (*func)(const Vector &x);

public:
  SpatialFunc(T (*fp)(const Vector &x)) : func(fp) {}
  T operator()(const Vector &x) {return func(x);}
};


int main() {

  SpatialFunc<double> sf( [](const Vector &x){return x.norm();} );    
  Vector x(3, 4.0);
  
  std::cout << sf(x) << std::endl;
  
  return 0;
}
