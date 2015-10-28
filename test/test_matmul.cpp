#include "yafel_globals.hpp"
#include "lin_alg/matmul.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/operators.hpp"

#include <iostream>
#include <cstdlib>
#include <typeinfo>

using namespace yafel;


void test_1() {
  
  printf("test_matmul: test_1\n");
  std::size_t N=1000, N2=999;
  Matrix<> A(N,N, 1); 
  Matrix<> B(N, N2, 1);
  
  // Will streamline with operator* and take advantage of default copy ctor
  Matrix<> C = matmul(A,B);
  
  // ensure that dimensions are correct
  assert(C.rows() == A.rows() &&
  	 C.cols() == B.cols() &&
  	 "TEST: matmul dimensions");


  // ensure correct values
  // Each element of C should be A.cols()
  for(std::size_t i=0; i<C.rows(); ++i) {
    for(std::size_t j=0; j<C.cols(); ++j) {
      if(C(i,j) != (double)(A.cols()))
	std::cout << "C(" << i <<"," << j << ")=" << C(i,j) << std::endl;
      assert(C(i,j) == (double)(A.cols()) &&
	     "TEST: Matmul value of C(i,j)");
    }
  }

}


void test_2() {
  printf("test_matmul: test_2\n");
  std::size_t N = 100;

  Matrix<std::size_t> A(N,N,1);
  Matrix<std::size_t> B(N,N,1);
  
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      //A(i,j) = i+j+1;
      //B(i,j) = (i+1)*(j+1);
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
      if(C(i,j) != D(i,j)) {
	std::cout << "C(" << i <<"," << j << ")=" << C(i,j) <<
	  ", "<< "D(" << i <<"," << j << ")=" << D(i,j) << std::endl;
      }
    }
  }


  // assert correctness
  assert(C == D &&
	 "TEST: non-trivial matmul correctness");
}


void test_3() {
  printf("test_matmul: test_3\n");
  std::size_t N = 100;
  Matrix<double> A(N,N,-3);
  Matrix<double> B(N,N,1);
  
  // compute with matmul
  auto C = matmul(4*A+3*B, 3*A + 4*B - B);

  //compute using naive Aik Bkj (easy to be correct
  Matrix<double> D(N,N,0);
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

int main() {
  
  // simple construction and uniform matrix multiplication correctness
  test_1();
  
  // more complicated (non-uniform) matmul correctness, using size_t matrices to
  // circumvent non-associativity of floating-point multiplication
  test_2();

  // Test matmul with composed MatrixExpressions
  test_3();
  
  return 0;
}
