#include "yafel_globals.hpp"
#include "lin_alg/sparse_bcsr.hpp"
#include "lin_alg/sparse_coo.hpp"

#include <iostream>
#include <tuple>
#include <random>
#include <chrono>

using namespace yafel;



/*
 * Test operator() for correctness.
 * Test fills a matrix with integers and then reads
 * them back out to ensure correct retrieval
 */
bool test_1() {
  
  std::size_t N = 60;
  
  sparse_coo<int> coo;
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      coo.add(i,j,int(i*N+j));
    }
  }
  
  sparse_bcsr<3,int> bcsr(coo);

  // test for equality with known values 
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      if( bcsr(i,j) != int(i*N+j) )
        return false;
    }
  }

  return true;
}

/*
 * Test that bcsr.rows() and bcsr.cols() work properly.
 * Recall that bcsr rounds up to the nearest multiple of
 * BLOCK. This functionality may or may not change, depending
 * on if matrices of non-BLOCK-multiple sizes become useful.
 */
bool test_2() {
  std::size_t N = 60;

  sparse_coo<int> coo;
  for(std::size_t i=0; i<N; ++i) {
    coo.add(i,i,int(i));
  }
  
  sparse_bcsr<3,int> bcsr(coo);

  return bcsr.rows()==3*(N/3) && bcsr.cols()==3*(N/3);
}

/*
 * Construct directly from triplets, test accumulation.
 * test for rows() and cols() correctness.
 */
bool test_3() {
  
  std::size_t N = 60;
  std::vector<typename sparse_bcsr<3>::triplet> ts;

  for(std::size_t i=0; i<N; ++i) {
    ts.push_back(std::make_tuple(N-1, N-1, 1.0));
    ts.push_back(std::make_tuple(0, 0, 1.0));
  }

  sparse_bcsr<3> bcsr(ts);

  return bcsr.rows()==N && bcsr.cols()==N;

}


/*
 * Construct directly from triplets, test accumulation.
 * test for operator() correctness.
 */
bool test_4() {
  
  std::size_t N = 60;
  std::vector<typename sparse_bcsr<3>::triplet> ts;

  for(std::size_t i=0; i<N; ++i) {
    ts.push_back(std::make_tuple(N-1, N-1, 1.0));
    ts.push_back(std::make_tuple(0, 0, 1.0));
  }

  sparse_bcsr<3> bcsr(ts);

  bool retval = true;
  retval = retval && bcsr(0,0)==double(N);
  retval = retval && bcsr(N-1,N-1)==double(N);
  for(std::size_t i=1; i<N-1; ++i) {
    retval = retval && bcsr(i,i) == 0;
  }

  return retval;

}



int main() {

  std::size_t retval=0;
  
  if(!test_1()) {
    std::cout << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  if(!test_2()) {
    std::cout << "Failed test_2()" << std::endl;
    retval |= 1<<1;
  }

  if(!test_3()) {
    std::cout << "Failed test_3()" << std::endl;
    retval |= 1<<2;
  }

  if(!test_4()) {
    std::cout << "Failed test_4()" << std::endl;
    retval |= 1<<3;
  }
  
  
  
  return retval;
}
