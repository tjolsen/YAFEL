//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_VTUBACKEND_HPP
#define YAFEL_VTUBACKEND_HPP

#include "yafel_globals.hpp"
#include "output/OutputBackend.hpp"
#include <fstream>
#include <functional>

YAFEL_NAMESPACE_OPEN


/**
 * \class VTUBackend
 * \brief XML-based ascii backend for data output
 *
 * This backend implements the XML-based VTU file format,
 * and is intended to be compatible with Paraview (or other
 * compatible VTK readers)
 *
 * The format is designed to be one frame per file,
 * so the mesh is re-written at each time step.
 * For time-dependent outputs with many frames, consider
 * using an alternative format.
 */
class VTUBackend : public OutputBackend
{
public:
    //Public Interface from OutputBackend
    virtual void initialize(const std::string &fname_base) override;
    virtual void finalize() override;
    virtual void write_frame(OutputFrame &frame) override;
    virtual void write_data(const OutputData &data, double time=0) override;
    virtual void write_mesh(OutputMesh *outputMesh) override;

private:
    std::ofstream outfile;
    std::vector<std::function<void()>> cleanup_stack;
    std::vector<std::function<void()>> frame_stack;

    void finalize_frame();
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_VTUBACKEND_HPP
