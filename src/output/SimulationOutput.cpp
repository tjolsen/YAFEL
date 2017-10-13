//
// Created by tyler on 4/30/17.
//

#include "output/SimulationOutput.hpp"
#include "output/VTUBackend.hpp"
#include "output/HDFBackend.hpp"

YAFEL_NAMESPACE_OPEN

SimulationOutput::SimulationOutput(std::string fname_base, BackendType backendType)
        : fname_base(fname_base)
{

    switch (backendType) {
        case BackendType::VTU:
            backend = std::make_unique<VTUBackend>();
            break;
#ifdef USE_HDF5
        case BackendType::HDF5:
            backend = std::make_unique<HDFBackend>();
            break;
#endif
    }


    backend->initialize(fname_base);
}

SimulationOutput::~SimulationOutput()
{
    backend->finalize();
}

YAFEL_NAMESPACE_CLOSE