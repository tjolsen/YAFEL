//
// Created by tyler on 3/24/17.
//


#include "element/ShapeFunctionUtils.hpp"
#include "utils/DualNumber.hpp"
#include "utils/Range.hpp"

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cmath>

YAFEL_NAMESPACE_OPEN


std::vector<double> jacobi(int n, double alpha, double beta)
{
    std::vector<double> coeffs_1{(alpha + beta + 2) / 2, (alpha - beta) / 2};
    std::vector<double> coeffs_2{1.0};
    std::vector<double> coeffs;

    if (n == 0) {
        return coeffs_2;
    } else if (n == 1) {
        return coeffs_1;
    }

    coeffs.reserve(n + 1);
    coeffs.assign(coeffs_1.begin(), coeffs_1.end());
    coeffs_1.reserve(n + 1);
    coeffs_2.reserve(n + 1);

    for (auto i : IRange(1, n)) {
        auto tmp = i + alpha + beta + 1;
        double a = (2 * i + alpha + beta + 1) * (2 * i + alpha + beta + 2) / (2 * (i + 1) * tmp);

        double b = (beta * beta - alpha * alpha) * (2 * i + alpha + beta + 1) /
                   (2 * (i + 1) * tmp * (2 * i + alpha + beta));

        double c = (i + alpha) * (i + beta) * (2 * i + alpha + beta + 2) / ((i + 1) * tmp * (2 * i + alpha + beta));

        coeffs.push_back(0);
        for (auto &c : coeffs)
            c *= a;

        auto rit = coeffs.rbegin();

        for (auto rit_1 = coeffs_1.rbegin(); rit_1 != coeffs_1.rend(); ++rit, ++rit_1) {
            *rit -= b * (*rit_1);
        }

        rit = coeffs.rbegin();
        for (auto rit_2 = coeffs_2.rbegin(); rit_2 != coeffs_2.rend(); ++rit, ++rit_2) {
            *rit -= c * (*rit_2);
        }

        coeffs_2.swap(coeffs_1);
        coeffs_1.assign(coeffs.begin(), coeffs.end());
    }

    return coeffs;
}

void tensor_product_shape_functions(const std::vector<coordinate<>> &localPoints,
                                    const std::vector<coordinate<>> &quadPoints,
                                    int topoDim,
                                    std::vector<Eigen::VectorXd> &shapeValue,
                                    std::vector<Eigen::MatrixXd> &shapeGrad)
{
    shapeValue.resize(quadPoints.size());
    shapeGrad.resize(quadPoints.size());
    std::vector<std::vector<double>> jacobi_coeffs;
    int p = static_cast<int>(std::round(pow(localPoints.size(), 1.0 / topoDim)));
    for (auto n : IRange(0, p)) {
        jacobi_coeffs.push_back(jacobi(n, 0, 0));
    }

    Eigen::MatrixXd V;
    Eigen::VectorXd rhs;
    Eigen::MatrixXd gradRhs;

    std::function<double(int, coordinate<>)> shape_func;
    std::function<DualNumber<double>(int, coordinate<DualNumber<double>>)> shape_grad_func;

    if (topoDim == 1) {
        auto shape_func_lambda = [p, &jacobi_coeffs](int n, auto xi) {
            return poly_eval(jacobi_coeffs[n], xi(0));
        };
        shape_func = shape_func_lambda;
        shape_grad_func = shape_func_lambda;

    } else if (topoDim == 2) {
        auto shape_func_lambda = [p, &jacobi_coeffs](int n, auto xi) {
            int i = n / p;
            int j = n % p;
            return poly_eval(jacobi_coeffs[j], xi(0)) * poly_eval(jacobi_coeffs[i], xi(1));
        };
        shape_func = shape_func_lambda;
        shape_grad_func = shape_func_lambda;

    } else if (topoDim == 3) {
        auto shape_func_lambda = [p, &jacobi_coeffs](int n, auto xi) {
            int i = n / (p * p);
            int j = (n % (p * p)) / p;
            int k = n % p;
            return poly_eval(jacobi_coeffs[k], xi(0)) * poly_eval(jacobi_coeffs[j], xi(1)) *
                   poly_eval(jacobi_coeffs[i], xi(2));
        };
        shape_func = shape_func_lambda;
        shape_grad_func = shape_func_lambda;
    }


    V = make_vandermonde(localPoints, shape_func);
    rhs.resize(V.rows());
    gradRhs.resize(V.rows(), topoDim);
    auto VTLU = V.transpose().lu();

    for (auto i : IRange(0, static_cast<int>(quadPoints.size()))) {

        for (auto j : IRange(0, (int) V.rows())) {
            rhs(j) = shape_func(j, quadPoints[i]);
            for (auto d : IRange(0, topoDim)) {
                coordinate<DualNumber<double>> xi;
                for (auto di : IRange(0, xi.dim())) {
                    xi(di) = make_dual(quadPoints[i](di));
                }
                xi(d).second = 1;
                gradRhs(j, d) = shape_grad_func(j, xi).second;
            }
        }

        shapeValue[i] = VTLU.solve(rhs);
        shapeGrad[i] = VTLU.solve(gradRhs);
    }

}

YAFEL_NAMESPACE_CLOSE