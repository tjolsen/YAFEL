#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/solver/iterative/cg_solve.hpp"
#include <cstdlib>
#include <iostream>

using namespace yafel;

// solve an identity matrix system with random rhs
bool test_1() {
  
  // seed rand w/ 0 for deterministic results
  // maybe switch to time(0) ?
  srand(0);
  
  std::size_t N = 10;
  sparse_coo<> coo;
  Vector<> rhs(N,0);
  for(std::size_t i=0; i<N; ++i) {
    coo.add(i,i,1.0);
    rhs(i) = double(rand()%100);
  }
  
  sparse_csr<> csr(coo);

  auto res = cg_solve(csr, rhs);

  return res == rhs;
}






int main() {
  
  int retval = 0;

  if(!test_1()) {
    std::cerr << "Failed test_1()" << std::endl;
    retval |= 1<<0;
  }
  
  return retval;
}
