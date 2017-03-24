//
// Created by tyler on 3/14/17.
//

#include "element/ShapeFunctionUtils.hpp"
#include <iostream>

using namespace yafel;


int main() {

    auto coeffs = jacobi(5,1,0);

    for(auto c : coeffs)
        std::cout << c << std::endl;

    return 0;
}


