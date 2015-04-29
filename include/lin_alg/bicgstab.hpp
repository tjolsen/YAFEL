#ifndef _YAFEL_BICGSTAB_HPP
#define _YAFEL_BICGSTAB_HPP

/*
Follows the algorithm on wikipedia for the stabilized biconjugate gradient method
 */

#include "yafel_globals.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/sparse_csr.hpp"

#define BICGSTAB_TOL 1.0e-10

YAFEL_NAMESPACE_OPEN

Vector bicgstab(const sparse_csr &A, const Vector &rhs);
Vector bicgstab(const FullMatrix &A, const Vector &rhs);

YAFEL_NAMESPACE_CLOSE


#endif
