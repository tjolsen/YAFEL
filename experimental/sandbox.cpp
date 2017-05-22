//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensors.hpp"
#include "utils/DualNumber.hpp"
#include <element/ElementFactory.hpp>

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>



using namespace yafel;
using std::cout;
using std::endl;

int main()
{

    ElementFactory EF(1);
    ElementType et(ElementTopology::Simplex,3,1);

    auto &E = EF.getElement(et);

    for(auto &FN : E.face_nodes) {
        for(auto n : FN) {
            cout << n << "  ";
        }
        cout << endl;
    }
    return 0;
}


