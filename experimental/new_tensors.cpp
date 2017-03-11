//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/tensor_expression_types/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_functions/operators.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorContraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorFunctor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorMap.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorDyadicProduct.hpp"
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/tensor_functions/UnaryFunctions.hpp"
#include "lin_alg/new_tensor/tensor_functions/update_assignment.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>

using namespace yafel;

int main()
{

    Tensor<3,2,int> x(1),y(0);

    y(0,1) = 1;

    x += y.perm<1,0>();
    x += y;

    for(int i=0; i<x.dim(); ++i) {
        for (auto xi : x(i, colon()))
            std::cout << xi << " ";
        std::cout << std::endl;
    }


    //for(auto xi : x(0,colon()))
    //std::cout << xi << "\n";

    return dot(x,x);
}
