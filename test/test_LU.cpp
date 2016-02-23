#include "yafel_globals.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/operators.hpp"
#include "lin_alg/solver/direct/LUDecomposition.hpp"
#include <iostream>
#include <random>

using namespace yafel;

bool test_1() {
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> dist(0,10);
  
  std::size_t N = 100;

  Matrix<double> A(N,N);
  Vector<double> x(N);

  for(std::size_t i=0; i<N; ++i) {
    x(i) = dist(generator);
    for(std::size_t j=0; j<N; ++j) {
      A(i,j) = dist(generator);
    }
  }

  Vector<double> b = A*x;

  LUDecomposition<double> LU(A);

  Vector<double> y = LU.linsolve(b);

  return (x-y).dot(x-y) < 1.0e-6;
}



int main() {

  int retval = 0;
  
  if(!test_1()) {
    std::cout << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  return retval;
}
