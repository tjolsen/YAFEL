//
// Created by tyler on 1/12/18.
//

#ifndef YAFEL_EIGENCHOLESKY_HPP
#define YAFEL_EIGENCHOLESKY_HPP

#include "yafel_globals.hpp"

#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <stdexcept>

YAFEL_NAMESPACE_OPEN

namespace LinearSolve {

/**
 * Tag type to dispatch the Eigen Conjugate Gradient solver
 */
struct EigenCholeskyLLT
{};
struct EigenCholeskyLDLT
{};

namespace detail {

template<typename S, int O>
constexpr auto getUpLo(const Eigen::SparseMatrix<S,O>&) {
    if constexpr(O == Eigen::RowMajor) {
        return Eigen::Upper;
    }
    else {
        return Eigen::Lower;
    }
};

template<typename MatrixType, typename VectorType>
void solve_impl(VectorType &result, MatrixType const &A, VectorType const &b, EigenCholeskyLLT)
{
    constexpr int UpLo = getUpLo(A);

    Eigen::SimplicialLLT<MatrixType, UpLo> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        if(solver.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Eigen::SimplicialLLT failure: NumericalIssue");
        }
        else {
            throw std::runtime_error("Eigen::SimplicialLLT failure: Unknown");
        }

    }

    result = solver.solve(b);
};

template<typename MatrixType, typename VectorType>
void solve_impl(VectorType &result, MatrixType const &A, VectorType const &b, EigenCholeskyLDLT)
{
    constexpr auto UpLo = getUpLo(A);

    Eigen::SimplicialLDLT<MatrixType, UpLo> solver;
    solver.compute(A);
    if(solver.info() != Eigen::Success) {
        if(solver.info() == Eigen::NumericalIssue) {
            throw std::runtime_error("Eigen::SimplicialLLT failure: NumericalIssue");
        }
        else {
            throw std::runtime_error("Eigen::SimplicialLLT failure: Unknown");
        }

    }

    result = solver.solve(b);
};

}//end namespace detail


}//end namespace LinearSolve

YAFEL_NAMESPACE_CLOSE



#endif //YAFEL_EIGENCHOLESKY_HPP
