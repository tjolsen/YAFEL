//
// Created by tyler on 6/21/17.
//

#ifndef YAFEL_UNIAXIALTENSION_HPP
#define YAFEL_UNIAXIALTENSION_HPP

#include "yafel.hpp"

YAFEL_NAMESPACE_OPEN

auto UniaxialTension(DoFManager &dofm, double L, double strainValue)
{
    double TOL = 1.0e-3;
    DirichletBC x0(dofm,0,0);
    x0.selectByFunction([L,TOL](auto x){return abs(x)(0) < TOL;});
    DirichletBC y0(dofm,0,1);
    y0.selectByFunction([L,TOL](auto x){return abs(x)(1) < TOL;});
    DirichletBC z0(dofm,0,2);
    z0.selectByFunction([L,TOL](auto x){return abs(x)(2) < TOL;});

    DirichletBC xL(dofm, L*strainValue, 0);
    xL.selectByFunction([L,TOL](auto x){return std::abs(x(0)-L) < TOL;});

    return std::vector<DirichletBC>{x0,y0,z0,xL};
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_UNIAXIALTENSION_HPP
