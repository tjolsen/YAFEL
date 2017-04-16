//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"


#include <iostream>

using namespace yafel;

struct PoissonEquation
{
    static void LocalResidual(const Element&, double, Eigen::Map<double>&) {}
    static void LocalTangent(const Element&, double, Eigen::Map<double>&) {}
};

int main()
{

    CGAssembly<PoissonEquation>();
    std::cout << "Success!" << std::endl;

    return 0;
}


