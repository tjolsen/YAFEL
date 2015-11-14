#include "yafel_globals.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/sparse_coo.hpp"

#include <iostream>

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


int main() {
  
  int retval = 0;
  
  if(!test_1()) {
    std::cout << "Failed test_1" << "\n";
    retval |= 1<<0;
  }

  return retval;
}
