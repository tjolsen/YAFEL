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

int main()
{

    Tensor<3,3> x;

    int idx=0;
    for(auto& xi : x)
        xi = idx++;


    auto xslice = x(slice_sentinel(), 0,slice_sentinel());

    for(int i=0; i<3; ++i) {
        cout << "i = " << i <<":"<<std::endl;
        for(int j=0; j<3; ++j) {
            for(int k=0; k<3; ++k) {
                cout << x(i, j, k) << "  ";
            }
            cout << endl;
        }
        cout << endl;
    }


    cout << "Slice:"<<endl;
    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            cout << xslice(i,j) << "  ";
        }
        cout << endl;
    }
    return 0;
}