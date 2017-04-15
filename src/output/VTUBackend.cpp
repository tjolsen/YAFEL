//
// Created by tyler on 4/11/17.
//

#include "output/OutputFrame.hpp"
#include "output/VTUBackend.hpp"
#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "element/ElementFactory.hpp"
#include "utils/Range.hpp"

#include <string>


YAFEL_NAMESPACE_OPEN
//forward declare
static int ElementType_to_VTKType(const ElementType &et);


void VTUBackend::initialize(const std::string &fname_base, double time)
{
    if (is_initialized) {
        return;
    }

    std::string fname = fname_base + "_" + std::to_string(time) + ".vtu";
    outfile.open(fname);

    if (outfile.good()) {
        is_initialized = true;
    }

    outfile << "<?xml version=\"1.0\"?>\n";
    outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";

    cleanup_stack.push_back([this]() {
        outfile << "</VTKFile>";
        outfile << std::endl;
        outfile.close();
    });
}


void VTUBackend::finalize()
{
    while (cleanup_stack.size() > 0) {
        auto &&f = cleanup_stack.back();
        f();
        cleanup_stack.pop_back();
    }
}

void VTUBackend::finalize_frame()
{
    while (frame_stack.size() > 0) {
        auto &&f = frame_stack.back();
        f();
        frame_stack.pop_back();
    }
}


void VTUBackend::write_frame(OutputFrame &frame)
{
    //write mesh
    write_mesh(frame.outputMesh);

    //write point data
    if (frame.point_data.size() > 0) {
        outfile << "<PointData>\n";
        for (auto &pd : frame.point_data) {
            write_data(*pd);
        }
        outfile << "</PointData>\n";
    }

    //write cell data
    int n_local_cells{0};
    for(auto n : frame.outputMesh->local_cells_per_cell) {
        n_local_cells += n;
    }

    if (frame.cell_data.size() > 0) {
        outfile << "<CellData>\n";
        for (auto &pd : frame.cell_data) {
            Eigen::VectorXd tmp(n_local_cells);
            int idx{0};
            for(int c : IRange(0,static_cast<int>(frame.outputMesh->local_cells_per_cell.size()))) {
                for(int i : IRange(0,frame.outputMesh->local_cells_per_cell[c])) {
                    tmp(idx++) = (*(pd->data))(c);
                }
            }
            OutputData tmp_data(*pd);
            tmp_data.data = &tmp;
            write_data(tmp_data);
        }
        outfile << "</CellData>\n";
    }

    finalize_frame();
}

void VTUBackend::write_data(const OutputData &data)
{

    std::vector<int> comps;
    for (auto i : IRange(0, static_cast<int>(data.dof_mask.size()))) {
        if (data.dof_mask[i] != 0)
            comps.push_back(i);
    }

    if (comps.size() > 1 && data.dataType == OutputData::DataType::Scalar) {
        throw (std::runtime_error(
                "VTUBackend::write_data: Incompatible OutputData::dof_mask and OutputData::DataType: data = " +
                data.name));
    }
    if (comps.size() > 3 && data.dataType == OutputData::DataType::Vector) {
        throw (std::runtime_error(
                "VTUBackend::write_data: Incompatible OutputData::dof_mask and OutputData::DataType: data = " +
                data.name));
    }
    if (comps.size() > 9) {
        throw (std::runtime_error(
                "VTUBackend::write_data: Too many components in OutputData::dof_mask: data = " + data.name));
    }

    if (data.data->rows() % data.dof_mask.size() != 0) {
        throw (std::runtime_error(
                "VTUBackend::write_data: Incompatible dof_mask length. "
                        "Data length must be a multiple of dof_mask length: data = " +
                data.name));
    }

    //number of data locations (points/cells)
    int n_locs = data.data->rows() / data.dof_mask.size();
    int dof_per_loc = data.dof_mask.size();
    auto &V = *data.data;

    //characters for separating vector/tensor components
    std::vector<char> v_sep{' ', ' ', '\n'};
    std::vector<char> tensor_sep(9, ' ');
    tensor_sep.back() = '\n';


    //Write DataArray header
    int ncomps = 1 * (data.dataType == OutputData::DataType::Scalar)
                 + 3 * (data.dataType == OutputData::DataType::Vector)
                 + 9 * (data.dataType == OutputData::DataType::Tensor);

    outfile << "<DataArray type=\"Float32\" format=\"ascii\" Name=\""
               + data.name
               + "\" NumberOfComponents=\""
               + std::to_string(ncomps) + "\">\n";

    //loop over locations writing data components to file
    for (auto n : IRange(0, n_locs)) {

        switch (data.dataType) {
            case OutputData::DataType::Scalar: {
                outfile << V(n * dof_per_loc + comps[0]) << "\n";
            }
                break;
            case OutputData::DataType::Vector: {
                int vcomp{0};
                for (auto c : comps) {
                    outfile << V(n * dof_per_loc + c) << v_sep[vcomp++];
                }
                for (auto i : IRange(vcomp, 3)) {
                    outfile << 0 << v_sep[i];
                }
            }
                break;
            case OutputData::DataType::Tensor: {
                int tcomp{0};
                for (auto c : comps) {
                    outfile << V(n * dof_per_loc + c) << v_sep[tcomp++];
                }
                for (auto i : IRange(tcomp, 3)) {
                    outfile << 0 << v_sep[i];
                }
            }
                break;
        }

    } // end n_locs loop

    outfile << "</DataArray>\n";

}


void VTUBackend::write_mesh(OutputMesh *outputMesh)
{
    cleanup_stack.push_back([this]() {
        outfile << "</Piece>\n";
        outfile << "</UnstructuredGrid>\n";
    });

    auto &dofm = *outputMesh->dofm;
    int n_points = outputMesh->dofm->dof_nodes.size();
    int n_parent_cells = outputMesh->dofm->element_offsets.size() - 1;

    //Need an element factory to build visualization topology
    ElementFactory EF(1);

    std::vector<int> expanded_cells;
    std::vector<int> expanded_cell_offsets;
    std::vector<int> expanded_cell_vtk_type;
    outputMesh->local_cells_per_cell.resize(n_parent_cells);

    int offset{0};
    std::vector<int> local_cell;
    std::vector<int> global_nodes;
    for (auto e : IRange(0, n_parent_cells)) {
        auto et = dofm.element_types[e];
        if (et.elementTopology == ElementTopology::None) {
            continue;
        }
        auto &E = EF.getElement(et);
        dofm.getGlobalNodes(e, global_nodes);

        outputMesh->local_cells_per_cell[e] = E.localMesh.nCells();



        for (auto lc : IRange(0, E.localMesh.nCells())) {
            int vtk_type = ElementType_to_VTKType(dofm.CellType_to_ElementType(E.localMesh.getCellType(lc), 1));
            E.localMesh.getCellNodes(lc, local_cell);
            expanded_cell_offsets.push_back(offset);
            offset += local_cell.size();

            expanded_cell_vtk_type.push_back(vtk_type);

            for (auto n : local_cell) {
                expanded_cells.push_back(global_nodes[n]);
            }
        }
    }
    expanded_cell_offsets.push_back(offset);
    int n_cells = expanded_cell_offsets.size() - 1;


    //Begin writing
    outfile << "<UnstructuredGrid>\n"
            << "<Piece NumberOfPoints=\""
               + std::to_string(n_points)
               + "\" NumberOfCells=\""
               + std::to_string(n_cells)
               + "\">\n";

    //write points
    outfile << "<Points>\n"
            << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto &x : dofm.dof_nodes) {
        outfile << x(0) << ' ' << x(1) << ' ' << x(2) << '\n';
    }
    outfile << "</DataArray>\n"
            << "</Points>\n";


    outfile << "<Cells>\n"
            << "<DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">\n";

    //write cell connectivity. each cell written on a separate line only for human-readability. not important for format.
    for (auto c : IRange(0, n_cells)) {
        auto start = expanded_cell_offsets[c];
        auto end = expanded_cell_offsets[c + 1];
        for (auto idx : IRange(start, end - 1)) {
            outfile << expanded_cells[idx] << ' ';
        }
        outfile << expanded_cells[end - 1] << '\n';
    }
    outfile << "</DataArray>\n";

    //write cell offsets
    outfile << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    for (auto c : IRange(1, n_cells + 1)) {
        outfile << expanded_cell_offsets[c] << '\n';
    }
    outfile << "</DataArray>\n";

    //write cell types
    outfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (auto vtktype : expanded_cell_vtk_type) {
        outfile << vtktype << '\n';
    }
    outfile << "</DataArray>\n";
    outfile << "</Cells>\n";
}

static int ElementType_to_VTKType(const ElementType &et)
{
    switch (et.elementTopology) {
        case ElementTopology::Simplex:
            return 5 * (et.topoDim == 2) + 10 * (et.topoDim == 3); //tri or tet
        case ElementTopology::TensorProduct:
            return 3 * (et.topoDim == 1) + 9 * (et.topoDim == 2) + 12 * (et.topoDim == 3); //quad or hex
        case ElementTopology::None:
            return 2; //poly-vertex, draw it as dots and see how it goes
    }
}


YAFEL_NAMESPACE_CLOSE