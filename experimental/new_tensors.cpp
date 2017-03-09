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
#include "lin_alg/new_tensor/tensor_functions/tensor_dot.hpp"
#include "lin_alg/new_tensor/mp_utils/sequence_functions.hpp"
#include "lin_alg/new_tensor/mp_utils/slice_mp_utils.hpp"

#include <iostream>
#include <cmath>

using namespace yafel;
using namespace std;



int main()//(int argc, char **argv)
{
    double lambda=1.0e9;
    double mu=1.0e8;

    auto Cijkl = make_TensorFunctor<3,4,double>([lambda,mu](int i, int j, int k, int l)
                                                {return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));});


    Tensor<3,2> gradU(0); gradU(0,1) = .01;
    auto strain = ((gradU + gradU.perm<1,0>())/2);

    auto stress = contract<2>(Cijkl,strain).eval();

    int s = dot(stress,stress)>0;
    for(auto s : stress)
        cout << s << endl;
    cout << endl;
    cout << sqrt(dot(stress,stress)) << endl;
    return s;
}
