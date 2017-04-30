//
// Created by tyler on 4/10/17.
//

#ifndef YAFEL_OUTPUTBACKEND_HPP
#define YAFEL_OUTPUTBACKEND_HPP

#include "yafel_globals.hpp"

#include <eigen3/Eigen/Core>
#include <vector>


YAFEL_NAMESPACE_OPEN
// Forward Declarations
class OutputData;
class OutputFrame;
class OutputMesh;


/**
 * \class OutputBackend
 * \brief Base class defining interface for output backends
 */
class OutputBackend
{
public:
    virtual void initialize(const std::string & fname_base) = 0;
    virtual void finalize() = 0;
    //virtual void finalize_frame() = 0;

    /// Write an entire frame
    virtual void write_frame(OutputFrame &frame) = 0;

    /// Write a data vector associated with a dof manager.
    virtual void write_data(const OutputData &data) = 0;

    /// Write a mesh via a DoFManager (constructs an ElementFactory under the hood for topology)
    virtual void write_mesh(OutputMesh *outputMesh) = 0;

protected:
    bool is_initialized{false};
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTBACKEND_HPP
