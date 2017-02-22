//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorAddition.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSubtraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorScaled.hpp"

using namespace yafel;


int main()
{
    constexpr int N = 3;


    Tensor<N,2,int> x,y;

    int count{1};
    for(auto &xi : x) {
        xi = count++;
    }

    count = 1;
    for(auto &yi : y) {
        yi = count++;
    }

    Tensor<N,2,int> z = 3*(x+y);

    int s{0};
    for(auto zi : z) {
        s += zi;
    }
    //std::cout << s << std::endl;

    return s;
}