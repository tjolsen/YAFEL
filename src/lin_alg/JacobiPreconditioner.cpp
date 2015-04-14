#include "lin_alg/JacobiPreconditioner.hpp"
#include <iostream>

YAFEL_NAMESPACE_OPEN

JacobiPreconditioner::JacobiPreconditioner(const sparse_csr &A) :
  AdiagInv(A.getRows())
{
  
  bool flag;
  //at this stage, it's not going to do anything about zeros on diag
  //other than issue a warning
  for(unsigned i=0; i<A.getRows(); ++i) {
    double aii = A(i,i,flag);
    if(!flag) {
      std::cerr << "Warning: JacobiPreconditioner found 0.0 on diag of A." << std::endl;
      std::cerr << "Warning: Results using this JacobiPreconditioner are bad." << std::endl;
    }
    AdiagInv(i) = 1.0/aii;
  }
  
}

Vector JacobiPreconditioner::MinvV(const Vector &rhs) const {

  unsigned N = rhs.getLength();
  Vector retvec(N,0.0);
  
  for(unsigned i=0; i<N; ++i) {
    retvec(i) = rhs(i)*AdiagInv(i);
  }
  
  return retvec;
}



YAFEL_NAMESPACE_CLOSE
