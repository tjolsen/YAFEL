#include "yafel_globals.hpp"
#include "lin_alg/ILUPreconditioner.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/sparse_coo.hpp"

#include <iostream>
#include <vector>

using namespace yafel;

bool test_1() {

  std::size_t N = 5;
  std::size_t NN = N*N;
  
  sparse_coo<double> coo;

  for(std::size_t i=0; i<NN; ++i) {
    coo.add(i,i,-4);
    if(i >= 1) {
      coo.add(i, i-1, 1);
    }
    if(i >= N) {
      coo.add(i, i-N, 1);
    }
    if(i < NN-1) {
      coo.add(i, i+1, 1);
    }
    if(i < NN-N) {
      coo.add(i, i+N, 1);
    }
  }
  

  sparse_csr<double> csr(coo);
  
  ILUPreconditioner<double> ilu(csr);

  std::vector<typename sparse_csr<double>::triplet> ts = ilu.getILU().copy_triplets();
  
  //for(std::tuple<std::size_t, std::size_t, double> t : ts) {
  //std::cout << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t)<<std::endl;
  ///}
  
  return true;
}



int main() {

  int retval = 0;
  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }

  return retval;
}
