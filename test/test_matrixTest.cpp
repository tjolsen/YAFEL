#include "lin_alg/Matrix.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>

using namespace yafel;

int main(int argc, char **argv) {
  
  std::size_t N=10, M=20;
  double val = 1.0;
  
  Matrix<double> a(N);
  Matrix<double> b(M,N);
  Matrix<double> c(M,N,val);

  // Test proper size
  assert(a.rows()==N && a.cols()==N && "TEST:Matrix<double> constructed to correct square size");
  assert(b.rows()==M && b.cols()==N && "TEST:Matrix<double> constructed to correct rectangular size");

  // Test proper value construction
  double sum=0;
  for(std::size_t i=0; i<c.rows(); ++i) {
    for(std::size_t j=0; j<c.cols(); ++j) {
      sum += c(i,j);
    }
  }
  assert(sum == M*N*val && "TEST:Matrix<double> constructed with proper values");
  
  return 0;
}
