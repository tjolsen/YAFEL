//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensors.hpp"
#include "utils/DualNumber.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;

int main()
{

    auto C = [](auto F){return (F.template perm<1,0>()*F).eval();};

    auto dCdF_ij = [](auto F, int k, int l) {
        auto d_il = TensorEye<3,double>()(colon(),l);
        auto d_jl = TensorEye<3,double>()(colon(),l);
        return (otimes(d_il, F(k,colon())) + otimes(F(k,colon()),d_jl)).eval();
    };

    Tensor<3,2,DualNumber<double>> F = TensorEye<3,DualNumber<double>>();
    F(0,1) = 0.1;


    int k=1, l=1;
    F(k,l).second = 1;

    for(auto f : F) {
        cout << f << endl;
    }
    cout << endl;

    for(auto c : C(F)) {
        cout << c.second << endl;
    }

    cout << endl;

    for(auto dcdf : dCdF_ij(F,k,l)) {
        cout << dcdf.first << endl;
    }



    return 0;
}


