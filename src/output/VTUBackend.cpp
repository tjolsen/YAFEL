//
// Created by tyler on 4/11/17.
//

#include "output/VTUBackend.hpp"
#include "output/OutputFrame.hpp"
#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "element/ElementFactory.hpp"
#include "utils/Range.hpp"

#include <string>


YAFEL_NAMESPACE_OPEN
//forward declare
static int ElementType_to_VTKType(const ElementType &et);


void VTUBackend::initialize(const std::string &fname_base)
{
    if (is_initialized) {
        return;
    }

    std::string fname = fname_base + ".vtu";
    outfile.open(fname);

    if (outfile.good()) {
        is_initialized = true;
    }

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
        auto &f = frame_stack.back();
        f();
        frame_stack.pop_back();
    }
}


void VTUBackend::write_frame(OutputFrame &frame)
{
    outfile << "<?xml version=\"1.0\"?>\n";
    outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";

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
    for (auto n : frame.outputMesh->local_cells_per_cell) {
        n_local_cells += n;
    }

    if (frame.cell_data.size() > 0) {
        outfile << "<CellData>\n";
        for (auto &cd : frame.cell_data) {
            Eigen::VectorXd tmp(n_local_cells);
            int idx{0};
            for (int c : IRange(0, static_cast<int>(frame.outputMesh->local_cells_per_cell.size()))) {
                for (int i : IRange(0, frame.outputMesh->local_cells_per_cell[c])) {
                    tmp(idx++) = (cd->data)(c);
                }
            }
            OutputData tmp_data(*cd);
            tmp_data.data = Eigen::Map<Eigen::VectorXd>(tmp.data(), tmp.size());
            write_data(tmp_data);
        }
        outfile << "</CellData>\n";
    }

    finalize_frame();
}

void VTUBackend::write_data(const OutputData &data,double)
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

    if (data.data.rows() % data.dof_mask.size() != 0) {
        throw (std::runtime_error(
                "VTUBackend::write_data: Incompatible dof_mask length. "
                        "Data length must be a multiple of dof_mask length: data = " +
                data.name));
    }

    //number of data locations (points/cells)
    int n_locs = data.data.rows() / data.dof_mask.size();
    int dof_per_loc = data.dof_mask.size();
    auto &V = data.data;

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

    int n_cells = outputMesh->expanded_cell_offsets.size() - 1;
    int n_points = outputMesh->dofm->dof_nodes.size();

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
        auto start = outputMesh->expanded_cell_offsets[c];
        auto end = outputMesh->expanded_cell_offsets[c + 1];
        for (auto idx : IRange(start, end - 1)) {
            outfile << outputMesh->expanded_cells[idx] << ' ';
        }
        outfile << outputMesh->expanded_cells[end - 1] << '\n';
    }
    outfile << "</DataArray>\n";

    //write cell offsets
    outfile << "<DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">\n";
    for (auto c : IRange(1, n_cells + 1)) {
        outfile << outputMesh->expanded_cell_offsets[c] << '\n';
    }
    outfile << "</DataArray>\n";

    //write cell types
    outfile << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (auto et : outputMesh->expanded_cell_element_type) {
        int vtktype = ElementType_to_VTKType(et);
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

    // gcc can't tell that this is unreachable code, for some reason,
    // so it issues a warning without this line.
    return -1;
}


YAFEL_NAMESPACE_CLOSE