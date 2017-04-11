//
// Created by tyler on 3/14/17.
//

#include "mesh/Mesh.hpp"
#include "utils/Range.hpp"
#include <iostream>
#include <element/Element.hpp>
#include <element/ElementFactory.hpp>
#include "utils/DoFManager.hpp"

#include "output/VTUBackend.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"

using namespace yafel;

int main()
{
    /*
    //---------------------------------------------------
    std::vector<coordinate<>> X;
    int N = 3;
    double dx = 1.0;
    for(auto i : IRange(0,N)) {
        for(auto j : IRange(0,N)) {
            coordinate<> x;
            x(0) = j*dx;
            x(1) = i*dx;
            X.push_back(x);
        }
    }

    std::vector<int> cells, cell_offsets;
    std::vector<CellType> celltypes;
    int offset = 0;
    for(auto i : IRange(0,N-1)) {
        for(auto j : IRange(0,N-1)) {
            int idx = i*N + j;
            cells.push_back(idx);
            cells.push_back(idx+1);
            cells.push_back(idx+1+N);
            cells.push_back(idx+N);
            cell_offsets.push_back(offset);
            offset += 4;
            celltypes.push_back(CellType::Quad4);
        }
    }
    cell_offsets.push_back(offset);

    //---------------------------------------------------

    Mesh M(Mesh::DefinitionScheme::Explicit,
           X,cells,cell_offsets,celltypes);
    */
    Mesh M("minsquare.msh");
    int p = 2;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, 1);


    VTUBackend VTU;
    OutputMesh outputMesh(dofm);
    OutputFrame outputFrame(outputMesh);

    Eigen::VectorXd x_coord(dofm.dof_nodes.size());
    Eigen::VectorXd y_coord(dofm.dof_nodes.size());
    Eigen::VectorXd xy_coord(2 * dofm.dof_nodes.size());
    for (auto i : IRange(0, static_cast<int>(dofm.dof_nodes.size()))) {
        x_coord(i) = dofm.dof_nodes[i](0);
        y_coord(i) = dofm.dof_nodes[i](1);
        xy_coord(2 * i) = dofm.dof_nodes[i](0);
        xy_coord(2 * i + 1) = dofm.dof_nodes[i](1);
    }

    Eigen::VectorXd cell_id(dofm.element_types.size());
    double id{0};
    for(auto i : IRange(0,static_cast<int>(cell_id.rows()))) {
        cell_id(i) = id;
        id += 1;
    }

    OutputData Xdata(x_coord, "X");
    OutputData Ydata(y_coord, "Y");
    OutputData XYdata(xy_coord, "XY", OutputData::DataLocation::Point, OutputData::DataType::Vector, {1, 1});
    OutputData cellid_data(cell_id,"ID",OutputData::DataLocation::Cell,OutputData::DataType::Scalar, {1});
    outputFrame.point_data.push_back(&Xdata);
    outputFrame.point_data.push_back(&Ydata);
    outputFrame.point_data.push_back(&XYdata);
    outputFrame.cell_data.push_back(&cellid_data);

    VTU.initialize("test_output", 0);
    VTU.write_frame(outputFrame);
    VTU.finalize();

    return 0;
}


