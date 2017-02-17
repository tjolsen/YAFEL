//
// Created by tyler on 2/17/17.
//

#ifndef YAFEL_ELEMENT_HPP
#define YAFEL_ELEMENT_HPP

#include "yafel_globals.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class Element
 */
class Element {

public:
    int poly_order; ///< spatial interpolation order
    int spatialDims; ///< number of spatial dimensions of element/simulation
    int topoDims; ///< topological dimension of element (line=1, quad=2, hex=3)

    //Vector of coordinates of local points
    std::vector<coordinate<>> localPoints_xi;

    //Quadrature rule (ie, points and weights)
    std::vector<double> quadrature_weights;
    std::vector<coordinate<>> quadrature_points;

    std::vector<Vector<double>> shape_values;
    std::vector<Matrix<double>> shape_gradients_xi;
};



YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_HPP
