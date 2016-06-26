#ifndef _YAFEL_IDENTITYPRECONDITIONER_HPP
#define _YAFEL_IDENTITYPRECONDITIONER_HPP

#include "yafel_globals.hpp"
#include "lin_alg/solver/iterative/Preconditioner.hpp"
#include "lin_alg/Vector.hpp"
#include "lin_alg/access_sparse_matrix.hpp"

YAFEL_NAMESPACE_OPEN


/*
 * A do-nothing preconditioner to test whether the preconditioned iterative
 * solvers effectively fall back on the non-preconditioned case in the absense
 * of any work being done in the solve() method.
 *
 * Ctor is given this signature to maintain a consistent API with the other preconditioners
 */

template<typename MT, typename dataType>
class IdentityPreconditioner : public Preconditioner<IdentityPreconditioner<MT,dataType>, dataType>
{

public:
    IdentityPreconditioner(const access_sparse_matrix<MT,dataType>&) {}

    void solve(Vector<dataType>&) const {}

};


YAFEL_NAMESPACE_CLOSE


#endif
