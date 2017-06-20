//
// Created by tyler on 4/29/17.
//

#ifndef YAFEL_SIMULATIONOUTPUT_HPP
#define YAFEL_SIMULATIONOUTPUT_HPP

#include "yafel_globals.hpp"
#include "output/OutputBackend.hpp"
#include "output/OutputData.hpp"
#include "output/OutputMesh.hpp"
#include "output/OutputFrame.hpp"

#include <string>
#include <memory>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Enum to pass into SimulationOutput constructor to select
 * an output backend
 */
enum class BackendType
{
    VTU,
    HDF5
};

/**
 * \class SimulationOutput
 * \brief Handles writing simulation output data using the specified backend.
 *
 */
class SimulationOutput
{
public:
    SimulationOutput(
            std::string fname_base = "output",
            BackendType backendType = BackendType::VTU
    );

    ~SimulationOutput();

    /**
     * \brief Capture an output frame and write it to file
     */
    template<typename T>
    void captureFrame(T &fesystem_maybe, std::function<void(T&,OutputFrame&)> &captureFunction = default_capture<T>());


    /**
     * Default frame capture function. Copy each component into scalar data
     * and add separately. AtNSD * tempts to detect cell/point data automatically,
     * defaulting to point data.
     * @tparam T FESystem type. Leaving the door open for non-virtual polymorphism
     * @return return a function that will capture an output frame
     */
    template<typename T>
    static auto default_capture() //std::function<void(T&,OutputFrame&)>
    {
        return [](T &feSys, OutputFrame &frame) -> void {
            frame.time = feSys.currentTime();
            auto &dofm = feSys.getDoFManager();
            int dofpn = dofm.dof_per_node;
            auto &solution = feSys.getSolution();
            std::vector<int> dof_mask(dofpn, 0);

            //Default data location is point data
            OutputData::DataLocation dataLocation = OutputData::DataLocation::Point;
            if (dofpn * dofm.nCells() == solution.rows()) {
                //Odds are, if you've triggered this, you've got Cell data, not Point data
                // If not, then don't be lazy and write your own frame capture function...
                dataLocation = OutputData::DataLocation::Cell;
            }

            std::string name_base("Solution_");
            OutputData::DataType dt = OutputData::DataType::Scalar;
            for (auto i : IRange(0, dofpn)) {
                dof_mask[i] = 1;

                frame.addData(
                        std::make_shared<OutputData>(
                                solution,
                                name_base + std::to_string(i),
                                dataLocation, dt,
                                dof_mask));

                dof_mask[i] = 0;
            }


            //frame

        };
    }


private:
    std::string fname_base;

    std::unique_ptr<OutputBackend> backend;
};


//--------------------------------------------------------
// Implementation of captureFrame
//--------------------------------------------------------
template<typename T>
void SimulationOutput::captureFrame(T &feSystem, std::function<void(T&,OutputFrame&)> &captureFunction)
{

    OutputMesh outputMesh(feSystem.getDoFManager());
    OutputFrame frame(outputMesh);

    captureFunction(feSystem, frame);

    backend->write_frame(frame);
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_SIMULATIONOUTPUT_HPP
