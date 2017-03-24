//
// Created by tyler on 3/24/17.
//


#include "element/ShapeFunctionUtils.hpp"
#include "utils/Range.hpp"


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
            *rit -= c*(*rit_2);
        }


        coeffs_2 = coeffs_1;
        coeffs_1 = coeffs;
        //coeffs_2.swap(coeffs_1);
        //coeffs_1.assign(coeffs.begin(), coeffs.end());
    }

    return coeffs;
}


YAFEL_NAMESPACE_CLOSE