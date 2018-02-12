//
// Created by tyler on 2/11/18.
//

#ifndef YAFEL_EIGENBICGSTAB_HPP
#define YAFEL_EIGENBICGSTAB_HPP


#include "yafel_globals.hpp"
#include <Eigen/IterativeLinearSolvers>

YAFEL_NAMESPACE_OPEN

namespace LinearSolve {

/**
 * Tag type to dispatch the Eigen Conjugate Gradient solver
 */
struct EigenBICGSTABTag
{
};

namespace detail {

template<typename MatrixType, typename VectorType>
void solve_impl(VectorType &result, MatrixType const &A, VectorType const &b, EigenBICGSTABTag&)
{
    Eigen::BiCGSTAB<MatrixType> solver;
    solver.compute(A);
    result = solver.solveWithGuess(b,result);
};

}//end namespace detail


}//end namespace LinearSolve

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_EIGENBICGSTAB_HPP
