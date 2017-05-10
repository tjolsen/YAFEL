//
// Created by tyler on 5/9/17.
//

#ifndef YAFEL_SIMPLYSUPPORTEDSQUARE_HPP
#define YAFEL_SIMPLYSUPPORTEDSQUARE_HPP


#include "boundary_conditions/DirichletBC.hpp"

YAFEL_NAMESPACE_OPEN

auto SimplySupportedSquare(const DoFManager &dofm, double L = 1.0)
{

    double TOL = 1.0e-6;
    std::vector<DirichletBC> bcs;

    std::vector<int> comps{0, 1};
    std::vector<double> xvals{0, L};

    for (auto comp : comps) {
        for (auto xval : xvals) {

            auto func = [comp, xval, TOL](const coordinate<> &x) {
                return std::abs(x(comp) - xval) < TOL && std::abs(x(2)) < TOL;
            };

            for (int i = 0; i < 3; ++i) {
                bcs.emplace_back(dofm, 0, i);
                bcs.back().selectByFunction(func);
            }
        }
    }


    return bcs;
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SIMPLYSUPPORTEDSQUARE_HPP
