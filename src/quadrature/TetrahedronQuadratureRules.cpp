//
// Created by tyler on 4/13/17.
//

#include "quadrature/QuadratureRule.hpp"

YAFEL_NAMESPACE_OPEN

void QuadratureRule::get_tetrahedron_quadrature(int porder)
{

    switch(porder) {
        default:
        case 1:
            nodes = {{1./4., 1./4., 1./4.}};
            weights = {1./6.};
            return;
    }


}


YAFEL_NAMESPACE_CLOSE