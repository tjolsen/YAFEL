#ifndef _YAFEL_JACOBIPRECONDITIONER_HPP
#define _YAFEL_JACOBIPRECONDITIONER_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Preconditioner.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/sparse_csr.hpp"

YAFEL_NAMESPACE_OPEN

class JacobiPreconditioner : public Preconditioner {

private:
  Vector AdiagInv;

public:
  JacobiPreconditioner(const sparse_csr &A);
  
  Vector MinvV(const Vector &rhs) const;
};

YAFEL_NAMESPACE_CLOSE

#endif
