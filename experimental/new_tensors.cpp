//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_functions/operators.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorContraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorFunctor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorMap.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorDyadicProduct.hpp"
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>
#include <cmath>
#include <vector>

using namespace yafel;
using namespace std;



int main()
{

    Tensor<4,1,int> x,y;
    int count = 1;
    auto yit = y.begin();
    for(auto &xi : x) {
        *yit = xi = count++;
        ++yit;
    }

    auto A = otimes(x,y).eval();

    auto C=otimes(A,x).eval();
    auto D = otimes(y,A).eval();

    auto E = (C - D.perm<1,2,0>().eval()).eval();

    return dot(E,E);
}
