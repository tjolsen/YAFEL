#include "yafel_globals.hpp"
#include "lin_alg/sparse_bcsr.hpp"
#include "lin_alg/sparse_coo.hpp"

#include <iostream>
#include <tuple>
#include <random>
#include <chrono>

using namespace yafel;


bool test_1() {
  
  std::size_t N = 6;
  
  sparse_coo<int> coo;
  for(std::size_t i=0; i<N; ++i) {
    for(std::size_t j=0; j<N; ++j) {
      coo.add(i,j,int(i*N+j));
    }
  }
  
  sparse_bcsr<3,int> bcsr(coo);

  for(auto biptr : bcsr.brow_ptr)
    std::cout << biptr << " ";
  std::cout << std::endl;

  for(auto biptr : bcsr.bcol_index)
    std::cout << biptr << " ";
  std::cout << std::endl;

  for(auto biptr : bcsr.data)
    std::cout << biptr << " ";
  std::cout << std::endl;
  
  return true;
}




int main() {

  std::size_t retval=0;
  
  if(!test_1()) {
    std::cout << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }
  
  
  
  return retval;
}
