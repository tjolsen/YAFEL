//
// Created by tyler on 2/21/17.
//

#include "lin_alg/new_tensor/tensors.hpp"

#include <iostream>

using namespace yafel;

int main()
{

    Tensor<3,2> Q = TensorEye<3>();
    Q(1,0) = .1;
    Tensor<3,2> A;

    double lam[3] = {1,2,3};

    for(auto i : IRange(0,3))
        A += otimes(lam[i]*Q(colon(),i), Q(colon(),i));

    Tensor<3,2> A0(A);
    spectral_decomposition(A,Q); // warning, this modifies A, Q

    // Compute ||A0 - QAQ^T|| --> Check accuracy of spectral decomposition
    double err = norm((A0 - ((Q*A).eval()*Q.perm<1,0>().eval()).eval()));

    return err < 1.0e-15;
}
