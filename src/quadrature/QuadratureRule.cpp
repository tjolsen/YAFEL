//
// Created by tyler on 3/25/17.
//

#include "quadrature/QuadratureRule.hpp"
#include "utils/Range.hpp"
#include <cmath>

YAFEL_NAMESPACE_OPEN

QuadratureRule::QuadratureRule(int polyOrder, QuadratureType quadratureType)
        : nodes(), weights()
{

    if (quadratureType == QuadratureType::GAUSS_LOBATTO) {
        make_gauss_lobatto_1D(polyOrder);
    } else if (quadratureType == QuadratureType::GAUSS_LEGENDRE) {
        make_gauss_legendre_1D(polyOrder);
    }


}


void QuadratureRule::make_gauss_lobatto_1D(int Npoints)
{
    //int Npoints = polyOrder + 1;
    nodes.resize(Npoints);
    weights.resize(Npoints);

    double PI = std::atan2(1., 1.) * 4;
    double TOL = 1.0e-15;
    //newton-raphson solve for quadrature points
    for (int xi = 0; xi < Npoints; ++xi) {
        double x = -std::cos((PI * xi) / (Npoints - 1));
        double xold = x;
        double Pold, Pn, Pnew;

        do {
            xold = x;
            Pold = 1;
            Pn = x;

            //compute legendre polynomial value at current x
            for (int n = 2; n < Npoints; ++n) {
                Pnew = ((2 * n - 1) * x * Pn - (n - 1) * Pold) / n;

                Pold = Pn;
                Pn = Pnew;
            }

            x = x - (x - Pold / Pn) / Npoints;

        } while (std::abs(x - xold) > TOL);

        nodes[xi](0) = x;
        //compute weight from formula
        double w = 2 / (Npoints * (Npoints - 1) * Pn * Pn);
        weights[xi] = w;
    }


}


void QuadratureRule::make_gauss_legendre_1D(int Npoints)
{

    int polyOrder = Npoints;
    double PI = std::atan2(1., 1.) * 4;
    double TOL = 1.0e-14;

    nodes.resize(polyOrder);
    weights.resize(polyOrder);

    // Newton-Raphson solve for xi-th root of Legendre polynomial
    for (int xi = 0; xi < polyOrder; ++xi) {

        double x = -std::cos(PI * (1 + xi - 0.25) / (polyOrder + 0.5));
        double Pold, Pn, Pnew, dP;

        do {
            //compute legendre polynomial
            Pold = 1;
            Pn = x;
            for (int n = 2; n <= polyOrder; ++n) {
                Pnew = ((2 * n - 1) * x * Pn - (n - 1) * Pold) / n;

                Pold = Pn;
                Pn = Pnew;
            }

            //legendre polynomial derivative
            dP = (polyOrder / (x * x - 1)) * (x * Pn - Pold);

            x -= Pn / dP;
        } while (std::abs(Pn) > TOL);

        nodes[xi](0) = x;
        //compute weight
        double w = 2.0 / ((1 - x * x) * (dP * dP));
        weights[xi] = w;
    }


}


QuadratureRule QuadratureRule::make_tensor_product(QuadratureType qt, int topoDim, int polyOrder)
{
    QuadratureRule Q1d;
    if (qt == QuadratureType::GAUSS_LEGENDRE) {
        int npts_1d = (polyOrder + 1) / 2 + (polyOrder + 1) % 2;
        Q1d.make_gauss_legendre_1D(npts_1d);
    } else if (qt == QuadratureType::GAUSS_LOBATTO) {
        int npts_1d = (polyOrder + 3) / 2 + (polyOrder + 3) % 2;
        Q1d.make_gauss_lobatto_1D(npts_1d);
    }


    if(topoDim == 1) {
        return Q1d;
    }

    if(topoDim == 2) {
        int N = Q1d.nodes.size();
        QuadratureRule Q2d;
        Q2d.nodes.reserve(N*N);
        Q2d.weights.reserve(N*N);

        for(auto i : IRange(0,N)) {
            for(auto j : IRange(0,N)) {
                coordinate<> x;
                x(0) = Q1d.nodes[j](0);
                x(1) = Q1d.nodes[i](0);
                Q2d.nodes.push_back(x);
                Q2d.weights.push_back(Q1d.weights[i]*Q1d.weights[j]);
            }
        }

        return Q2d;
    }

    if(topoDim == 3) {
        int N = Q1d.nodes.size();
        QuadratureRule Q3d;
        Q3d.nodes.reserve(N*N);
        Q3d.weights.reserve(N*N);

        for(auto i : IRange(0,N)) {
            for(auto j : IRange(0,N)) {
                for(auto k : IRange(0,N)) {
                    coordinate<> x;
                    x(0) = Q1d.nodes[k](0);
                    x(1) = Q1d.nodes[j](0);
                    x(2) = Q1d.nodes[i](0);
                    Q3d.nodes.push_back(x);
                    Q3d.weights.push_back(Q1d.weights[i] * Q1d.weights[j] * Q1d.weights[k]);
                }
            }
        }

        return Q3d;
    }

    return {};
}


YAFEL_NAMESPACE_CLOSE