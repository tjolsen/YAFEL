//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/SmallVector.hpp"
#include "utils/parallel/yafel_parallel.hpp"
#include "utils/parallel/pardot.hpp"
#include "lin_alg/tensor/tensors.hpp"
#include <iostream>
#include <Eigen/Dense>


using namespace yafel;
using std::cout;
using std::endl;

int main()
{
    int N = 100;
    Eigen::VectorXd A, B;
    A = Eigen::VectorXd::Constant(N, 1.0);
    B = Eigen::VectorXd::Constant(N, 1.0);

    auto val = pardot(A,B);


    std::cout << val << std::endl;
    return 0;
}


