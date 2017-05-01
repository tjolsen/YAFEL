//
// Created by tyler on 4/30/17.
//

#include <output/OutputFrame.hpp>
#include <output/OutputMesh.hpp>
#include "output/HDFBackend.hpp"


YAFEL_NAMESPACE_OPEN

//forward declaration
static int ElementType_to_XDMFType(const ElementType &et);


void HDFBackend::initialize(const std::string &fname_base)
{
    if (is_initialized) {
        return;
    }


    xdmf_outfile << "<?xml version=\"1.0\" ?>\n"
                 << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
                 << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
                 << "<Domain>"
                 << "<Grid Type=\"Collection\" CollectionType=\"Temporal\">\n";


    cleanup_stack.push_back(
            [this]() {
                xdmf_outfile << "</Grid>\n"
                             << "</Domain>\n"
                             << "</Xdmf>\n";
                xdmf_outfile.close();
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
    int ncells = frame.outputMesh->dofm->nCells();
    int elem_len = static_cast<int>(frame.outputMesh->dofm->elements.size());

    // length of elements vector plus an integer for the element type.
    // A triangle would be "4 n0 n1 n2"
    int cell_vec_len = elem_len + ncells;


    /*
     * Write the XDMF Parts:
     *
     * - Topology
     * - Geometry
     * - Data
     */
    xdmf_outfile << "<Grid Name=\"grid\" Type=\"Uniform\">\n"
                 << "<Time Type=\"Single\" Value=\"" + std::to_string(time) +">\n";

    xdmf_outfile << "<Topology TopologyType=\"Mixed\" NumberOfElements=\"" + std::to_string(ncells) + "\">\n"
                 << "DataItem Dimensions=\"" + std::to_string(cell_vec_len) + "\" NumberType=\"Int\" Format=\"HDF\">\n"
                 << h5_fname << ":/cells";



    /*
     * Write the HDF5 Parts:
     * - Mesh (Topo/Geom) if not already written
     * - Data
     */
};


// For element type numbering, see documentation
// http://www.xdmf.org/index.php/XDMF_Model_and_Format#Topology
int ElementType_to_XDMFType(ElementType et)
{

    switch(et.elementTopology) {

        case ElementTopology::Simplex:
            if(et.topoDim == 2)
                return 4;
            if(et.topoDim == 3)
                return 6;
        case ElementTopology::TensorProduct:
            if(et.topoDim == 2)
                return 5;
            if(et.topoDim == 3)
                return 9;
        case ElementTopology::None:
            return 1;
    }

}

YAFEL_NAMESPACE_CLOSE