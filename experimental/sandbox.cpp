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


    Mesh M("twoTets.msh");
    M.buildInternalFaces();

    int p = 4;

    DoFManager dofm(M, DoFManager::ManagerType::CG, p);

    OutputMesh om(dofm);
    OutputFrame frame(om);
    VTUBackend vtu;
    vtu.initialize("test");
    vtu.write_frame(frame);
    vtu.finalize();

    ElementFactory EF(1);

    std::cout << dofm.dof_nodes[36](0) << "  " << dofm.dof_nodes[36](1) << "  " << dofm.dof_nodes[36](2) << std::endl;
    std::cout << dofm.dof_nodes[39](0) << "  " << dofm.dof_nodes[39](1) << "  " << dofm.dof_nodes[39](2) << std::endl;
    std::cout << dofm.dof_nodes[50](0) << "  " << dofm.dof_nodes[50](1) << "  " << dofm.dof_nodes[50](2) << std::endl;
    std::cout << dofm.dof_nodes[32](0) << "  " << dofm.dof_nodes[32](1) << "  " << dofm.dof_nodes[32](2) << std::endl;
    std::cout << dofm.dof_nodes[37](0) << "  " << dofm.dof_nodes[37](1) << "  " << dofm.dof_nodes[37](2) << std::endl;

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
                std::cout << container[n] << "  ";
            }
            std::cout << std::endl;
            dofm.getGlobalNodes(R, container);
            for (auto n : ER.face_perm[fR][rR][1]) {
                std::cout << container[n] << "  ";
            }
            std::cout << std::endl << std::endl;
        }
    }

    /*
    ElementType et{ElementTopology::TensorProduct, 2, p};
    Element E(et);

    auto idxs_l = make_element_boundary_nodes(E.localMesh.getGeometryNodes(), {1,0,0}, {0,1,0},
                                              [](auto xi) { return std::abs(xi(0) - 1) < 1.0e-6; });


    auto idxs_r = make_element_boundary_nodes(E.localMesh.getGeometryNodes(), {-1,0,0}, {0,1,0},
            [](auto xi) { return std::abs(xi(0) + 1) < 1.0e-6; });


    std::vector<int> left_nodes, right_nodes;
    dofm.getGlobalNodes(0,left_nodes);
    dofm.getGlobalNodes(1, right_nodes);

    for(auto il : idxs_l) {
        std::cout << left_nodes[il] << "  ";
    }
    std::cout << "\n";
    for(auto ir : idxs_r) {
        std::cout << right_nodes[ir] << "  ";
    }
    std::cout << std::endl;
    for(auto l : left_nodes) {
        std::cout << l << "  ";
    }
    std::cout << std::endl;
    for(auto r : right_nodes) {
        std::cout << r << "  ";
    }*/

    std::cout << std::endl;
    return 0;
}


