//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorAddition.hpp"
#include "lin_alg/new_tensor/mp_utils/sequences.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"

#include <iostream>

using namespace yafel;


int main()
{

    Tensor<3,2,int> x,y;

    for(auto &xi : x) {
        xi = 1;
    }

    for(auto &yi : y) {
        yi = 2;
    }

    Tensor<3,2,int> z = x+y;

    int s{0};
    for(auto zi : z) {
        s += zi;
    }
    std::cout << s << std::endl;

    return 0;
}