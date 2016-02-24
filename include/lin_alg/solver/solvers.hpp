#ifndef __YAFEL_SOLVERS_HPP

/*
 * Convenience header for including all linear system solvers
 */ 

// Iterative Solvers and preconditioners
#include "lin_alg/solver/iterative/cg_solve.hpp"
#include "lin_alg/solver/iterative/bicgstab_solve.hpp"
#include "lin_alg/solver/iterative/Preconditioner.hpp"
#include "lin_alg/solver/iterative/JacobiPreconditioner.hpp"
#include "lin_alg/solver/iterative/ILUPreconditioner.hpp"

// Direct Solvers
#include "lin_alg/solver/direct/LUDecomposition.hpp"


#endif
