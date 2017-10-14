//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/SmallVector.hpp"
#include "utils/parallel/yafel_parallel.hpp"
#include "lin_alg/tensor/tensors.hpp"
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;

int main()
{
    cout << "Num cores: " << config::num_cores << endl
         << "VERSION " << config::VERSION << endl
         << "BUILD TYPE: " << config::BUILD_TYPE << endl
         << "Git Revision: " << config::GIT_REVISION << endl;

    return 0;
}


