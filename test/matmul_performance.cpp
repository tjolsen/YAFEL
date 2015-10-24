#include "yafel_globals.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/MatrixExpression.hpp"
#include "lin_alg/matmul.hpp"
#include "output/MatrixVisualization.hpp"

#include <chrono>

//horrible horrible hack. need to fix library compilation
#include "../src/output/MatrixVisualization.cpp"


using namespace yafel;


std::chrono::steady_clock::duration time_matmul(std::size_t N, Matrix<double> &A, Matrix<double> &B, Matrix<double> &C) {

  //Matrix<double> A(N);
  //Matrix<double> B(N);
  //Matrix<double> C(N);

  auto t1 = std::chrono::steady_clock::now();
  
  divconq_matmul(C,A,B, 0,0, 0,0, N, N, N);
  //auto C = matmul(A,B);


  auto t2 = std::chrono::steady_clock::now();

  return t2-t1;
}


std::chrono::steady_clock::duration time_naive_matmul(std::size_t N, Matrix<double> &A, Matrix<double> &B, Matrix<double> &C) {

  auto t1 = std::chrono::steady_clock::now();
  
  for(size_t i=0; i<N; ++i) {
    for(size_t j=0; j<N; ++j) {
      C(i,j) = 0;
      for(size_t k=0; k<N; ++k) {
	C(i,j) += A(i,k)*B(k,j);
      }
    }
  }

  auto t2 = std::chrono::steady_clock::now();

  return t2-t1;
}




int main() {

  std::vector<double> Nvec;
  std::vector<double> Gflops1;
  std::vector<double> Gflops2;

  std::size_t Nmax = 500;
  
  Matrix<double> A(Nmax);
  Matrix<double> B(Nmax);
  Matrix<double> C(Nmax);
  
  for(std::size_t n=10; n<Nmax; n += 25) {
    std::cout << n << ": ";
    Nvec.push_back(n);
    
    //2n^3 add/multiplies for non-strassen-like matrix mult
    double work = 2*n*n*n;
    
    auto deltaTime1 = time_matmul(n,A,B,C);
    auto deltaTime2 = time_naive_matmul(n,A,B,C);
    
    //time in seconds
    double dt1 = ((double)deltaTime1.count()*std::chrono::steady_clock::period::num)/std::chrono::steady_clock::period::den;
    double dt2 = ((double)deltaTime2.count()*std::chrono::steady_clock::period::num)/std::chrono::steady_clock::period::den;
    
    std::cout << dt1 << ", " << dt2 << std::endl;
    Gflops1.push_back(work/dt1);
    Gflops2.push_back(work/dt2);
  }
  
  MatrixVisualization MV;
  
  MV.scatter_xy(Nvec, Gflops1);
  MV.scatter_xy(Nvec, Gflops2);
  return 0;
}

