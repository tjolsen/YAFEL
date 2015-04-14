#ifndef _YAFEL_ILUPRECONDITIONER_HPP
#define _YAFEL_ILUPRECONDITIONER_HPP

#include "yafel_globals.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/Preconditioner.hpp"

YAFEL_NAMESPACE_OPEN

class ILUPreconditioner : public Preconditioner {
private:
  sparse_csr ILU;

public:
  // construct an incomplete LU factorization with no extra fill-in.
  // filling options could be added later, for S&G.
  // At first, no pivoting will be used, as this tends to increase the
  // bandwidth of resulting matrices arising from PDE numerical methods.
  ILUPreconditioner(const sparse_csr &A);
  Vector MinvV(const Vector &rhs) const;
  sparse_csr & getILU() {return ILU;}
};


YAFEL_NAMESPACE_CLOSE

#endif
