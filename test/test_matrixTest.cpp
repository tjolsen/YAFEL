#include "old_handmade_linalg/Matrix.hpp"
#include "old_handmade_linalg/operators.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>

using namespace yafel;

int main() {
  
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
  

  // Test copy constructor
  Matrix<> d(c);
  assert(d==c && "TEST: Matrix equality using copy ctor");

  // Test addition and scalar multiplication
  Matrix<> e = c+d;
  Matrix<> f = 2*c;
  assert(e==f && "TEST: Matrix doubling using addition and scalar multiplication");

  //Test subtraction
  Matrix<> g(M,N,0);
  Matrix<> h = c-d;
  assert(g==h && "TEST: Matrix subtraction");

  return 0;
}
