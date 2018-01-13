//
// Created by tyler on 1/12/18.
//

#ifndef YAFEL_LINEARSOLVER_HPP
#define YAFEL_LINEARSOLVER_HPP

#include "yafel_globals.hpp"

#include "lin_alg/linear_solvers/solvers/EigenConjugateGradient.hpp"
#include "lin_alg/linear_solvers/solvers/EigenCholesky.hpp"
#include "lin_alg/linear_solvers/solvers/VCLConjugateGradient.hpp"

#include <cassert>

YAFEL_NAMESPACE_OPEN

namespace LinearSolve {

/**
 * Solve a linear system using a pre-allocated result buffer.
 * Throws std::runtime_error if result is incorrectly-sized.
 * (i.e., it assumes you know what you're doing. Use solve() if you don't.)
 *
 * If relevant, the passed-in value of result is used as an initial guess.
 *
 * @tparam MatrixType
 * @tparam VectorType
 * @tparam TagType
 * @param result
 * @param A
 * @param b
 * @param tag
 */
template<typename MatrixType, typename VectorType, typename TagType>
void solve_inplace(VectorType &result, MatrixType const &A, VectorType const &b, TagType &tag)
{
    assert(A.rows() == b.rows() && "Incorrect linear system dimensions");
    assert(A.rows() == A.cols() && "A must be a square matrix");

    detail::solve_impl(result, A, b, tag);
};


/**
 * Solve a linear system, allocating a new vector for the result
 * @tparam MatrixType
 * @tparam VectorType
 * @tparam TagType
 * @param A
 * @param b
 * @param tag
 * @return
 */
template<typename MatrixType, typename VectorType, typename TagType>
VectorType solve(MatrixType const &A, VectorType const &b, TagType &tag)
{
    VectorType result(b.rows());
    for (int i = 0; i < b.rows(); ++i) {
        result(i) = 0;
    }

    solve_inplace(result, A, b, tag);

    return result;
}

}//end namespace LinearSolve

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_LINEARSOLVER_HPP
