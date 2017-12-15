//
// Created by tyler on 3/23/17.
//

#ifndef YAFEL_SHAPEFUNCTIONUTILS_HPP
#define YAFEL_SHAPEFUNCTIONUTILS_HPP

/**
 * \file
 * This file contains functions for computing shape functions
 * on tensor product and simplex elements.
 *
 * - Jacobi polynomials (Legendre polynomials for tensor product elements)
 * - Koornwinder polynomials (for 2d/3d simplices)
 */


#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "utils/Range.hpp"
#include <Eigen/Core>
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Evaluate a polynomial using Horner's Method
 *
 * Polynomial coeffs should be in decreasing order.
 * i.e.: polynomials of degree "n" are evaluated as:
 *
 * \begin{equation}
 * P(x) = coeffs[0]*x^{n} + coeffs[1]^x^{n-1} + ... + coeffs[n]
 * \end{equation}
 *
 * @tparam T datatype of "x"
 * @param coeffs coefficients of polynomial:
 * @param x
 * @return
 */
template<typename T>
T poly_eval(const std::vector<double> &coeffs, T x)
{
    T retval{0};
    for (int i = 0; i < static_cast<int>(coeffs.size()); ++i) {
        retval *= x;
        retval += coeffs[i];
    }
    return retval;
}

/**
 * \brief Create coefficients of n-th degree Jacobi polynomial
 *
 * Recurrance relation given in http://lsec.cc.ac.cn/~hyu/teaching/shonm2013/STWchap3.2p.pdf
 *
 * @param n Jacobi polynomial degree
 * @param alpha Jacobi alhpa
 * @param beta Jacobi beta
 * @return polynomial coefficients suitable for use in poly_eval
 */
std::vector<double> jacobi(int n, double alpha, double beta);


void tensor_product_shape_functions(const std::vector<coordinate<>> &localPoints,
                                    const std::vector<coordinate<>> &quadPoints,
                                    int topoDim,
                                    std::vector<Eigen::VectorXd> &shapeValue,
                                    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &shapeGrad);


void triangle_shape_functions(const std::vector<coordinate<>> &localPoints,
                              const std::vector<coordinate<>> &quadPoints,
                              int polyOrder,
                              std::vector<Eigen::VectorXd> &shapeValue,
                              std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &shapeGrad);

void tetrahedron_shape_functions(const std::vector<coordinate<>> &localPoints,
                                 const std::vector<coordinate<>> &quadPoints,
                                 int polyOrder,
                                 std::vector<Eigen::VectorXd> &shapeValue,
                                 std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> &shapeGrad);


template<typename FUNC>
Eigen::MatrixXd make_vandermonde(const std::vector<coordinate<>> &localPoints, FUNC &func)
{

    int N = localPoints.size();
    Eigen::MatrixXd V(N, N);

    for (auto i : IRange(0, N)) {
        for (auto j : IRange(0, N)) {
            auto vij = func(j, localPoints[i]);
            V(i, j) = vij;
        }
    }

    return V;
}

YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_SHAPEFUNCTIONUTILS_HPP
