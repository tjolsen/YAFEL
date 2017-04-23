//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "assembly/CGAssembly.hpp"
#include "element/element_boundary_nodes.hpp"
#include "boundary_conditions/DirichletBC.hpp"

#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"
#include "output/VTUBackend.hpp"

#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>


using namespace yafel;

int main()
{
    /*
    std::vector<coordinate<>> M_nodes{{0,0,0},
                                      {1,0,0},
                                      {2,0,0},
                                      {0,1,0},
                                      {1,1,0},
                                      {2,1,0}};
    std::vector<int> cells {0,1,4,3,1,2,5,4};
    std::vector<int> cell_offsets{0,4,8};
    std::vector<CellType> cell_types(2,CellType::Quad4);
    Mesh M(Mesh::DefinitionScheme::Explicit,M_nodes,cells,cell_offsets,cell_types);
    M.buildInternalFaces();*/


    Mesh M("twoQuads.msh");
    M.buildInternalFaces();

    int p = 5;

    DoFManager dofm(M, DoFManager::ManagerType::CG, p);

    OutputMesh om(dofm);
    OutputFrame frame(om);
    VTUBackend vtu;
    vtu.initialize("test");
    vtu.write_frame(frame);
    vtu.finalize();

    ElementFactory EF(1);

    for (auto f : M.getInternalFaces()) {
        if (!(f.left < 0 || f.right < 0)) {
            std::cout << f << std::endl;

            int L = f.left;
            int R = f.right;
            int fL = f.left_flocal;
            int fR = f.right_flocal;
            int rL = f.left_rot;
            int rR = f.right_rot;

            auto etL = dofm.element_types[L];
            auto etR = dofm.element_types[R];

            auto &EL = EF.getElement(etL);
            auto &ER = EF.getElement(etR);

            std::vector<int> container;
            dofm.getGlobalNodes(L, container);
            for (auto n : EL.face_perm[fL][rL][0]) {
                auto x = dofm.dof_nodes[container[n]];
                std::cout << container[n] << ":  " << x(0) << "  " << x(1) << "  " << x(2) << std::endl;
            }
            std::cout << std::endl << std::endl;
            dofm.getGlobalNodes(R, container);
            for (auto n : ER.face_perm[fR][rR][1]) {
                auto x = dofm.dof_nodes[container[n]];
                std::cout << container[n] << ":  " << x(0) << "  " << x(1) << "  " << x(2) << std::endl;
            }
            std::cout << std::endl << std::endl;
        }
    }


    std::cout << std::endl;
    return 0;
}


