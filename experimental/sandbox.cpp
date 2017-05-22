//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "lin_alg/tensor/tensors.hpp"
#include "utils/DualNumber.hpp"
#include <element/ElementFactory.hpp>
#include "utils/DoFManager.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>



using namespace yafel;
using std::cout;
using std::endl;

int main()
{

    Mesh M("twoQuads.msh");
    M.buildInternalFaces();
    DoFManager dofm(M,DoFManager::ManagerType::DG,1);


    for(auto &F : dofm.interior_faces) {



    }


    return 0;
}


