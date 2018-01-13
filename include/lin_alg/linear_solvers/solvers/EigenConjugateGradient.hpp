//
// Created by tyler on 1/12/18.
//

#ifndef YAFEL_EIGENCONJUGATEGRADIENT_HPP
#define YAFEL_EIGENCONJUGATEGRADIENT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/linear_solvers/LinearSolve.hpp"

#include <Eigen/IterativeLinearSolvers>

YAFEL_NAMESPACE_OPEN

namespace LinearSolve {

/**
 * Tag type to dispatch the Eigen Conjugate Gradient solver
 */
struct EigenConjugateGradientTag
{
};

namespace detail {

template<typename MatrixType, typename VectorType>
void solve_impl(VectorType &result, MatrixType const &A, VectorType const &b, EigenConjugateGradientTag)
{
    Eigen::ConjugateGradient<MatrixType, Eigen::Upper | Eigen::Lower> solver;
    solver.compute(A);
    result = solver.solveWithGuess(b,result);
};

}//end namespace detail


}//end namespace LinearSolve

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_EIGENCONJUGATEGRADIENT_HPP
