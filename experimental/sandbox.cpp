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

    int p = 4;

    DoFManager dofm(M,DoFManager::ManagerType::CG, p);



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
    }

    std::cout << std::endl;
    return 0;
}


