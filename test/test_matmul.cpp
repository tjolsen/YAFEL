#include "yafel_globals.hpp"
#include "lin_alg/matmul.hpp"
#include "lin_alg/Matrix.hpp"

#include <iostream>
#include <cstdlib>
#include <typeinfo>

using namespace yafel;


void test_1() {

  
  Matrix<double> A(10, 10, 1); 
  Matrix<double> B(10, 9, 1);
  
  // Will streamline with operator* and take advantage of default copy ctor
  Matrix<double> C = matmul(A,B);
  
  // ensure that dimensions are correct
  assert(C.rows() == A.rows() &&
  	 C.cols() == B.cols() &&
  	 "TEST: matmul dimensions");


  // ensure correct values
  // Each element of C should be A.cols()
  for(std::size_t i=0; i<C.rows(); ++i) {
    for(std::size_t j=0; j<C.cols(); ++j) {
        assert(C(i,j) == (double)(A.cols()) &&
  	 "TEST: Matmul value of C(i,j)");
    }
  }

}


void test_2() {
  
  std::size_t N = 10;

  Matrix<std::size_t> A(N,N,0);
  Matrix<std::size_t> B(N,N,0);
  
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      A(i,j) = i+j+1;
      B(i,j) = (i+1)*(j+1);
    }
  }

  // compute with matmul
  auto C = matmul(A,B);

  //compute using naive Aik Bkj (easy to be correct
  Matrix<std::size_t> D(N,N,0);
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      for(std::size_t k=0; k<N; ++k) {
	D(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  
  // assert correctness
  assert(C == D &&
	 "TEST: non-trivial matmul correctness");
}


void test_3() {
  
  std::size_t N = 10;
  Matrix<int> A(N,N,-3);
  Matrix<int> B(N,N,1);
  
  auto MEa = 4*A + 3*B;
  auto MEb = 3*A + 4*B - B;
  
  // compute with matmul
  auto C = matmul(4*A+3*B, 3*A + 4*B - B);

  //compute using naive Aik Bkj (easy to be correct
  Matrix<int> D(N,N,0);
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      for(std::size_t k=0; k<N; ++k) {
	D(i,j) += (4*A + 3*B)(i,k) * (3*A + 4*B - B)(k,j);
      }
    }
  }
  
  // assert correctness
  assert(C == D &&
	 "TEST: non-trivial MatrixExpression matmul correctness");
  
  
}

int main(int argc, char **argv) {
  
  // simple construction and uniform matrix multiplication correctness
  test_1();
  
  // more complicated (non-uniform) matmul correctness, using size_t matrices to
  // circumvent non-associativity of floating-point multiplication
  test_2();

  // Test matmul with composed MatrixExpressions
  test_3();
  
  return 0;
}
