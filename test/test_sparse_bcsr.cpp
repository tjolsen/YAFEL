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
  
  
  
  return retval;
}
