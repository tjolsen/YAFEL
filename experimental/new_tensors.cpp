//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorAddition.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSubtraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>

using namespace yafel;
using namespace std;

int main()//(int argc, char **argv)
{
    constexpr int N = 3;
    Tensor<N,3,int> x(1);

    int count = 0;
    for(auto &xi : x)
        xi = count++;


    auto xT = permute(x, sequence<1,0,2>());

    for(int i=0; i<N; ++i) {
        for(int j=0; j<N; ++j){
            for(int k=0; k<N; ++k) {
                cout << xT(i, j, k) << "  ";
            }
            cout << endl;
        }
        cout << endl << endl;
    }

    return 0;
}