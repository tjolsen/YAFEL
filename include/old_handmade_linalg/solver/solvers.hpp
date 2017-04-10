#ifndef __YAFEL_SOLVERS_HPP

/*
 * Convenience header for including all linear system solvers
 */ 

// Iterative Solvers and preconditioners
#include "old_handmade_linalg/solver/iterative/cg_solve.hpp"
#include "old_handmade_linalg/solver/iterative/pcg_solve.hpp"
#include "old_handmade_linalg/solver/iterative/bicgstab_solve.hpp"
#include "old_handmade_linalg/solver/iterative/Preconditioner.hpp"
#include "old_handmade_linalg/solver/iterative/IdentityPreconditioner.hpp"
#include "old_handmade_linalg/solver/iterative/JacobiPreconditioner.hpp"
#include "old_handmade_linalg/solver/iterative/ILUPreconditioner.hpp"

// Direct Solvers
#include "old_handmade_linalg/solver/direct/LUDecomposition.hpp"


#endif
