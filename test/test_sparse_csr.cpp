#include "yafel_globals.hpp"
#include "old_handmade_linalg/sparse_csr.hpp"
#include "old_handmade_linalg/sparse_coo.hpp"

#include <iostream>
#include <tuple>
#include <random>
#include <chrono>

using namespace yafel;


bool test_1() {
  
  std::size_t N = 10;
  // construct sparse_coo diagonal matrix of 1's
  // call sparse_csr ctor, test rows(), cols(), nnz()
  sparse_coo<double> coo;
  for(std::size_t i=0; i<N; ++i) {
    coo.add(i,i,1);
  }
  sparse_csr<double> csr(coo);

  return csr.rows()==N && csr.cols()==N && csr.nnz()==N;
}


// Test operator()
bool test_2() {

 std::size_t N = 10;
  // construct sparse_coo diagonal matrix of 1's
  // call sparse_csr ctor, test rows(), cols(), nnz()
  sparse_coo<double> coo;
  for(std::size_t i=0; i<N; ++i) {
    coo.add(i,i,1);
  }
  sparse_csr<double> csr(coo);
 
  double sum = 0;
  for(std::size_t i=0; i<N; ++i) {
    sum += csr(i,i);
  }
  
  return sum==(double)N;
}


// Test accuulation of values into 2 nonzeros
bool test_3() {

  std::size_t N = 10;
  std::vector<typename sparse_csr<>::triplet> ts;

  for(std::size_t i=0; i<N; ++i) {
    ts.push_back(std::make_tuple(N-1,N-1,1.0));
    ts.push_back(std::make_tuple(0,0,1.0));
  }

  sparse_csr<> csr(ts);
  
  return csr.rows()==N && csr.cols()==N && csr.nnz()==2;
}

// Test proper matrix sum
bool test_4() {
  
  std::size_t nvals = 500;
  std::size_t N = 20;
  sparse_coo<int> coo;
  int maxval = 10;
  
  unsigned seed0 = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine valgenerator(seed0);
  auto valdist = std::uniform_int_distribution<int>(-maxval, maxval);

  unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine ijgenerator(seed1);
  auto ijdist = std::uniform_int_distribution<std::size_t>(0, N-1);

  std::vector<std::tuple<std::size_t, std::size_t>> ijs;
  
  
  int sum = 0, truesum=0;
  for(std::size_t k=0; k<nvals; ++k) {

    int v = valdist(valgenerator);
    
    std::size_t i = ijdist(ijgenerator);
    std::size_t j = ijdist(ijgenerator);
    
    ijs.push_back(std::tuple<std::size_t,std::size_t>(i,j));
    coo.add(i,j,v);
    truesum += v;
  }

  sparse_csr<int> csr(coo);

  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      sum += csr(i,j);
    }
  }
  /*
  for(std::size_t i=0; i<ijs.size(); ++i) {
    int v = csr(std::get<0>(ijs[i]), std::get<1>(ijs[i]));
    sum += v;
  }
  */  
  return sum==truesum;
}

int main() {
  
  int retval = 0;



  if(!test_1()) {
    std::cout << "Failed test_1" << "\n";
    retval |= 1<<0;
  }
  if(!test_2()) {
    std::cout << "Failed test_2" << "\n";
    retval |= 1<<1;
  }
  if(!test_3()) {
    std::cout << "Failed test_3" << "\n";
    retval |= 1<<2;
  }

  if(!test_4()) {
    std::cout << "Failed test_4" << "\n";
    retval |= 1<<3;
  }


 return retval;
}
