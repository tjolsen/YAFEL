#include "yafel_globals.hpp"
#include "old_handmade_linalg/sparse_csr.hpp"
#include "old_handmade_linalg/sparse_coo.hpp"
#include "old_handmade_linalg/Vector.hpp"
#include "old_handmade_linalg/csr_spmv.hpp"
#include "old_handmade_linalg/operators.hpp"

#include <iostream>

using namespace yafel;

// Construct sparse_csr identity matrix, 
bool test_1() {
  
  std::size_t N=10;
  sparse_coo<int> coo;
  
  Vector<int> u(N);

  for(unsigned i=0; i<N; ++i) {
    coo.add(i,i,1);
    u(i) = i+1;
  }
  
  sparse_csr<int> csr(coo);
  
  Vector<int> v = csr_spmv(csr, u);
  
  return u==v;
}

int main() {

  int retval = 0;

  if(!test_1()) {
    std::cout << "Failed test_1\n";
    retval |= (1<<0);
  }

  return retval;
}
