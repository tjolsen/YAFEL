//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorAddition.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSubtraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>

using namespace yafel;
using namespace std;

int main(int argc, char **argv)
{
    int val;
    if(argc > 1)
        val = atoi(argv[1]);
    else
        val = 1;

    Tensor<3,3,int> x;

    //int idx=0;
    for(auto& xi : x)
        xi = val;

    //zero out the x(0,:,:) slice
    x(0,colon(),colon()) = Tensor<3,2,int>();

    int s{0};
    for(auto xit : x)
        s += xit;

    return s;
}