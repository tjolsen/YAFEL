//
// Created by tyler on 3/14/17.
//

#include "element/ShapeFunctionUtils.hpp"
#include "utils/Range.hpp"
#include <iostream>
#include <element/Element.hpp>

using namespace yafel;


int main()
{

    Element e({ElementClass::TensorProduct, 1, 3});

    std::vector<std::vector<double>> shapeVals;
    std::vector<std::vector<coordinate<>>> shapeGrads;

    tensor_product_shape_functions(e.localMesh.getGeometryNodes(),
                                   e.quadratureRule.nodes,
                                   e.elementType.topoDim,
                                   shapeVals,
                                   shapeGrads);


    return 0;
}


