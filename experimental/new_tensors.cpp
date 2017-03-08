//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/Tensor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorAddition.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSubtraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorScaled.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorSlice.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorPermutation.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorContraction.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorFunctor.hpp"
#include "lin_alg/new_tensor/tensor_expression_types/TensorMap.hpp"
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>

using namespace yafel;
using namespace std;

/*
template<typename T, bool b>
Tensor<3, 2, double> Hookean_stress(const TensorExpression<T, 3, 2, double, b> &strain)
{


};


double do_contract2(double a)
{
    Tensor<3, 1> lhs;
    Tensor<3, 3> rhs;

    for (auto &li : lhs) {
        li = a;
        a *= 1.1;
    }
    for (auto &ri : rhs) {
        ri = a;
        a /= 1.2;
    }

    Tensor<3, 2, double> res;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            res(i, j) = dot(lhs(colon()), rhs(colon(), i, j));
        }
    }

    return dot(res, res);
}

double do_contract(double a)
{

    Tensor<3, 1> lhs;
    Tensor<3, 3> rhs;

    for (auto &li : lhs) {
        li = a;
        a *= 1.1;
    }
    for (auto &ri : rhs) {
        ri = a;
        a /= 1.2;
    }

    Tensor<3, 2, double> res = contract(lhs, rhs, sequence<1>());

    return contract(res, res, sequence<2>());
}

*/

template<typename T>
void const_fancy_func(T const* ptr)
{
    auto xmap = make_TensorMap<3,2>(ptr);
    //xmap(1,1) *= 3;
    cout << dot(xmap,xmap) << endl;
}

template<typename T>
void fancy_func(T* ptr)
{
    auto xmap = make_TensorMap<3,2>(ptr);

    xmap(0,0) *= 2;

    const_fancy_func(ptr);
}

int main()//(int argc, char **argv)
{

    array<int,81> buffer;
    buffer.fill(1);

    //const_fancy_func(buffer.data());
    fancy_func(buffer.data());

    return 0;
}
