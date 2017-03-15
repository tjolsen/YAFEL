//
// Created by tyler on 2/17/17.
//

#ifndef YAFEL_ELEMENT_HPP
#define YAFEL_ELEMENT_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include "ElementType.hpp"
#include "new_mesh/Mesh.hpp"
#include "lin_alg/Matrix.hpp"
#include "lin_alg/Vector.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

/**
 * \class Element
 */
class Element {

public:
    Element(ElementType={ElementClass::None, 0, 0});

    ElementType elementType;

    // Mesh of the local "Master" element
    Mesh localMesh;

    //Quadrature rule (ie, points and weights)
    //std::vector<double> quadrature_weights;
    //std::vector<coordinate<>> quadrature_points;

    //std::vector<Vector<double>> shape_values;
    //std::vector<Matrix<double>> shape_gradients_xi;


private:
    void make_simplex();
    void make_tensorProduct();
};



YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_HPP
