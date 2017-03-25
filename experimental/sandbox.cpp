//
// Created by tyler on 3/14/17.
//

#include "element/ShapeFunctionUtils.hpp"
#include "utils/Range.hpp"
#include <iostream>

using namespace yafel;


int main() {

    auto coeffs = jacobi(8,0,0);

    for(auto c : coeffs)
        std::cout << c << std::endl;

    //for(auto x : Range<double>(-1.0,1.01,.01))
    //std::cout << x << "  " << poly_eval(coeffs, x) << std::endl;


    return 0;
}


