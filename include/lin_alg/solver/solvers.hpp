#ifndef __YAFEL_SOLVERS_HPP

/*
 * Convenience header for including all linear system solvers
 */ 

// Iterative Solvers and preconditioners
#include "lin_alg/iterative/cg_solve.hpp"
#include "lin_alg/iterative/bicgstab.hpp"
#include "lin_alg/iterative/Preconditioner.hpp"
#include "lin_alg/iterative/JacobiPreconditioner.hpp"
#include "lin_alg/iterative/ILUPreconditioner.hpp"

// Direct Solvers
#include "lin_alg/direct/LUDecomposition.hpp"


#endif
