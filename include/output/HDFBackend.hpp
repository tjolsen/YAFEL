//
// Created by tyler on 4/30/17.
//

#ifndef YAFEL_HDFBACKEND_HPP
#define YAFEL_HDFBACKEND_HPP

#include "yafel_globals.hpp"
#include "output/OutputBackend.hpp"

#ifdef USE_HDF5
#include <hdf5/serial/hdf5.h>
#else
#define hid_t int
#endif
#include <fstream>

YAFEL_NAMESPACE_OPEN

/**
 * \class HDFBackend
 * Output backend for the XDMF + HDF5 supported by Paraview.
 * The backend will be specialized for a non-changing mesh (modulo
 * displacements accounted for in a data field).
 *
 */
class HDFBackend : public OutputBackend
{
public:
    virtual void initialize(const std::string &fname_base) override;
    virtual void finalize() override;
    virtual void write_frame(OutputFrame &frame) override;
    virtual void write_data(const OutputData &data, double time) override;
    virtual void write_mesh(OutputMesh *outputMesh) override;

private:
    std::string xdmf_fname;
    std::string h5_fname;

    std::ofstream xdmf_outfile;
    hid_t h5_file;

    std::vector<std::function<void()>> cleanup_stack;
    bool is_initialized{false};
    bool mesh_is_written{false};
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_HDFBACKEND_HPP
