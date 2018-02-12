//
// Created by tyler on 1/12/18.
//

#ifndef YAFEL_VCLCONJUGATEGRADIENT_HPP
#define YAFEL_VCLCONJUGATEGRADIENT_HPP

#include "yafel_globals.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <viennacl/linalg/cg.hpp>

YAFEL_NAMESPACE_OPEN

namespace LinearSolve {

/**
 * Tag type to dispatch the Eigen Conjugate Gradient solver
 */
struct VCLConjugateGradientTag
{
    int verbose{1};
};

namespace detail {

template<typename MatrixType, typename VectorType>
void solve_impl(VectorType &result, MatrixType const &A, VectorType const &b, VCLConjugateGradientTag &vclcgtag)
{

    viennacl::vector<double> vcl_rhs(b.rows());
    viennacl::compressed_matrix<double> vcl_A(A.rows(), A.cols());
    viennacl::copy(b, vcl_rhs);
    viennacl::copy(A, vcl_A);

    //solve
    viennacl::linalg::cg_tag tag(1.0e-14, 2 * b.rows());
    viennacl::vector<double> vcl_result = viennacl::linalg::solve(vcl_A, vcl_rhs, tag);

    if(vclcgtag.verbose) {
        std::cout << ": System solved.\n\t Iterations: " << tag.iters() << "\n\tError: " << tag.error() << std::endl;
    }
    //Copy back
    viennacl::copy(vcl_result, result);

};

}//end namespace detail


}//end namespace LinearSolve

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_VCLCONJUGATEGRADIENT_HPP
