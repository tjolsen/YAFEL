//
// Created by tyler on 4/30/17.
//

#include <output/OutputFrame.hpp>
#include <output/OutputMesh.hpp>
#include "output/HDFBackend.hpp"

#include <hdf5/serial/hdf5.h>
#include <hdf5/serial/hdf5_hl.h>


YAFEL_NAMESPACE_OPEN

//forward declaration
static int ElementType_to_XDMFType(ElementType et);


void HDFBackend::initialize(const std::string &fname_base)
{
    if (is_initialized) {
        return;
    }

    // Open the xmf file
    xdmf_fname = fname_base + ".xmf";
    xdmf_outfile.open(xdmf_fname);
    if (!xdmf_outfile.is_open()) {
        throw std::runtime_error("Could not open XDMF file " + xdmf_fname);
    }


    // Open the HDF5 file
    h5_fname = fname_base + ".h5";
    h5_file = H5Fcreate(h5_fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (h5_file < 0) {
        throw std::runtime_error("Could not open h5 file " + h5_fname);
    }


    xdmf_outfile << "<?xml version=\"1.0\" ?>\n"
                 << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
                 << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
                 << "<Domain>\n"
                 << "<Grid Type=\"Collection\" CollectionType=\"Temporal\">\n";


    cleanup_stack.push_back(
            [this]() {
                xdmf_outfile << "</Grid>\n"
                             << "</Domain>\n"
                             << "</Xdmf>\n";
                xdmf_outfile.close();
                H5Fclose(h5_file);
            }
    );

    is_initialized = true;
}


void HDFBackend::finalize()
{
    while (cleanup_stack.size() > 0) {
        auto &&f = cleanup_stack.back();
        f();
        cleanup_stack.pop_back();
    }
}


void HDFBackend::write_frame(OutputFrame &frame)
{
    auto time = frame.time;


    xdmf_outfile << "<Grid Name=\"grid\" Type=\"Uniform\">\n"
                 << "<Time Type=\"Single\" Value=\"" + std::to_string(time) + "\"/>\n";

    write_mesh(frame.outputMesh);


    if (frame.point_data.size() > 0) {
        for (auto &pd : frame.point_data) {
            write_data(*pd, time);
        }
    }

    /*
    // Cell data not supported right now

    //write cell data
    int n_local_cells{0};
    for (auto n : frame.outputMesh->local_cells_per_cell) {
        n_local_cells += n;
    }

    if (frame.cell_data.size() > 0) {
        for (auto &cd : frame.cell_data) {
            Eigen::VectorXd tmp(n_local_cells);
            int idx{0};
            for (int c : IRange(0, static_cast<int>(frame.outputMesh->local_cells_per_cell.size()))) {
                for (int i : IRange(0, frame.outputMesh->local_cells_per_cell[c])) {
                    tmp(idx++) = (*(cd->data))(c);
                }
            }
            OutputData tmp_data(*cd);
            tmp_data.data = &tmp;
            write_data(tmp_data);
        }
    }*/


    xdmf_outfile << "</Grid>\n";
};

void HDFBackend::write_data(const OutputData &data, double time)
{

    std::vector<int> comps;
    int dof_per_loc = data.dof_mask.size();
    for (auto i : IRange(0, static_cast<int>(data.dof_mask.size()))) {
        if (data.dof_mask[i] != 0)
            comps.push_back(i);
    }

    if (comps.size() > 1 && data.dataType == OutputData::DataType::Scalar) {
        throw (std::runtime_error(
                "HDFBackend::write_data: Incompatible OutputData::dof_mask and OutputData::DataType: data = " +
                data.name));
    }
    if (comps.size() > 3 && data.dataType == OutputData::DataType::Vector) {
        throw (std::runtime_error(
                "HDFBackend::write_data: Incompatible OutputData::dof_mask and OutputData::DataType: data = " +
                data.name));
    }
    if (comps.size() > 9) {
        throw (std::runtime_error(
                "HDFBackend::write_data: Too many components in OutputData::dof_mask: data = " + data.name));
    }

    if (data.data->rows() % data.dof_mask.size() != 0) {
        throw (std::runtime_error(
                "HDFBackend::write_data: Incompatible dof_mask length. "
                        "Data length must be a multiple of dof_mask length: data = " +
                data.name));
    }


    int ncomps{0};
    int n_locs = data.data->rows() / data.dof_mask.size();
    std::string AttrType = "";
    std::string DataLoc = "";
    if (data.dataType == OutputData::DataType::Scalar) {
        ncomps = 1;
        AttrType = "Scalar";
    } else if (data.dataType == OutputData::DataType::Vector) {
        ncomps = 3;
        AttrType = "Vector";
    } else if (data.dataType == OutputData::DataType::Tensor) {
        ncomps = 9;
        AttrType = "Tensor";
    }

    if (data.dataLocation == OutputData::DataLocation::Point) {
        DataLoc = "Node";
    } else {
        DataLoc = "Cell";
    }

    std::string dataName = data.name + "_" + std::to_string(time);

    xdmf_outfile << "<Attribute Name=\""
                    + data.name
                    + "\" AttributeType=\""
                    + AttrType
                    + "\" Center=\""
                    + DataLoc
                    + "\">\n"
                 << "<DataItem Dimensions=\"" + std::to_string(n_locs * ncomps)
                    + "\" Format=\"HDF\">\n"
                 << h5_fname + ":/" + dataName + "\n"
                 << "</DataItem>\n"
                 << "</Attribute>\n";


    std::vector<double> tmp_data(n_locs * ncomps, 0);

    int idx{0};
    for (int loc = 0; loc < n_locs; ++loc) {
        for(auto c : comps) {
            tmp_data[idx++] = (*data.data)(loc*dof_per_loc + c);
        }
    }

    hsize_t tmp_size = tmp_data.size();
    H5LTmake_dataset(h5_file, ("/"+dataName).c_str(), 1, &tmp_size, H5T_NATIVE_DOUBLE, tmp_data.data());
}

void HDFBackend::write_mesh(OutputMesh *outputMesh)
{
    int nNodes = static_cast<int>(outputMesh->dofm->dof_nodes.size());
    int ncells = static_cast<int>(outputMesh->expanded_cell_element_type.size());
    int elem_len = static_cast<int>(outputMesh->expanded_cells.size());

    // length of elements vector plus an integer for the element type.
    // A triangle would be "4 n0 n1 n2"
    int cell_vec_len = elem_len + ncells;
    xdmf_outfile << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" + std::to_string(ncells) + "\">\n"
                 << "<DataItem Dimensions=\"" + std::to_string(cell_vec_len) + "\" NumberType=\"Int\" Format=\"HDF\">\n"
                 << h5_fname << ":/cells\n"
                 << "</DataItem>\n"
                 << "</Topology>"
                 << "<Geometry GeometryType=\"XYZ\">\n"
                 << "<DataItem Dimensions=\"" + std::to_string(nNodes * 3) + "\" Format=\"HDF\">\n"
                 << h5_fname << ":/XYZ\n"
                 << "</DataItem>\n"
                 << "</Geometry>\n";


    std::vector<int> cell_data(cell_vec_len, 0);
    int idx{0};
    for (int c = 0; c < ncells; ++c) {
        cell_data[idx++] = ElementType_to_XDMFType(outputMesh->expanded_cell_element_type[c]);
        for (int i = outputMesh->expanded_cell_offsets[c]; i < outputMesh->expanded_cell_offsets[c + 1]; ++i) {
            cell_data[idx++] = outputMesh->expanded_cells[i];
        }
    }

    hsize_t tmp_size = cell_vec_len;
    H5LTmake_dataset(h5_file, "/cells", 1, &tmp_size, H5T_NATIVE_INT, cell_data.data());


    std::vector<double> xyz_data(3 * nNodes);
    idx = 0;
    for (int i = 0; i < nNodes; ++i) {
        xyz_data[idx + 0] = outputMesh->dofm->dof_nodes[i](0);
        xyz_data[idx + 1] = outputMesh->dofm->dof_nodes[i](1);
        xyz_data[idx + 2] = outputMesh->dofm->dof_nodes[i](2);
        idx += 3;
    }

    tmp_size = 3 * nNodes;
    H5LTmake_dataset(h5_file, "/XYZ", 1, &tmp_size, H5T_NATIVE_DOUBLE, xyz_data.data());
}


// For element type numbering, see documentation
// http://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology
int ElementType_to_XDMFType(ElementType et)
{

    switch (et.elementTopology) {

        case ElementTopology::Simplex:
            if (et.topoDim == 2)
                return 4;
            if (et.topoDim == 3)
                return 6;
        case ElementTopology::TensorProduct:
            if (et.topoDim == 2)
                return 5;
            if (et.topoDim == 3)
                return 9;
        case ElementTopology::None:
            return 1;
    }

}

YAFEL_NAMESPACE_CLOSE