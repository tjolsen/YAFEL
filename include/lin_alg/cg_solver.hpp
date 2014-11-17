#ifndef _YAFEL_CG_SOLVER_HPP
#define _YAFEL_CG_SOLVER_HPP
/*
This follows the algorithm laid out on the Wikipedia page for Conjugate Gradient,
but some steps have been broken apart to minimize memory reallocations.
 */

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/sparse_csr.hpp"

#define CG_SOLVER_TOL 1.0e-14

YAFEL_NAMESPACE_OPEN

Vector cg_solve(const sparse_csr & A, const Vector & rhs);
Vector cg_solve(const sparse_csr & A, const Vector & rhs, const Vector & x);
Vector pcg_solve(const sparse_csr & A, const Vector & rhs, const Vector & x);
Vector cg_solve(const FullMatrix & A, const Vector & rhs);  
  
YAFEL_NAMESPACE_CLOSE

#endif
