//
// Created by tyler on 3/14/17.
//


#include "quadrature/LobattoPoints1D.hpp"
#include <cmath>
#include <cassert>

YAFEL_NAMESPACE_OPEN

std::vector<double> make_LobattoPoints_1D(int polyOrder)
{
    assert(polyOrder >= 1 && "make_LobattoPoints_1D: polyOrder must be >= 1");


    int Npoints = polyOrder + 1;
    std::vector<double> nodes(Npoints,0);

    double PI = std::atan2(1., 1.) * 4;
    double TOL = 1.0e-15;
    //newton-raphson solve for quadrature points
    for (unsigned xi = 0; xi < Npoints; ++xi) {
        double x = -std::cos((PI * xi) / polyOrder);
        double xold = x;
        double Pold, Pn, Pnew;

        do {
            xold = x;
            Pold = 1;
            Pn = x;

            //compute legendre polynomial value at current x
            for (unsigned n = 2; n < Npoints; ++n) {
                Pnew = ((2 * n - 1) * x * Pn - (n - 1) * Pold) / n;
                Pold = Pn;
                Pn = Pnew;
            }

            x = x - (x - Pold / Pn) / Npoints;

        } while (std::abs(x - xold) > TOL);

        nodes[xi] = x;
    }

    return nodes;
}

YAFEL_NAMESPACE_CLOSE