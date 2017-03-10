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
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"


#include <iostream>

using namespace yafel;
//using namespace std;


int main()
{

    Tensor<10,1,double> x,y;
    int count = 1;
    int sgn=1;
    auto yit = y.begin();
    for(auto &xi : x) {
        *yit = xi = sgn*count++;
        ++yit;
        sgn *= -1;
    }

    auto z = sign(x).eval();


    for(auto zi : z) {
        std::cout << zi << std::endl;
    }

    auto yy = dot(y,z);

    return yy>0;
}
