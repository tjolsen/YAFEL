//
// Created by tyler on 9/11/17.
//
#include "yafel.hpp"
#include "cell_centroids.hpp"
#include "FVDofm.hpp"
#include "CellToNodeMap.hpp"

using namespace yafel;


int main()
{

    BasicTimer timer;
    Mesh M("channel.msh");


    M.buildInternalFaces();

    constexpr int NSD = 2;




    // p=1 DoFManager, since p=0 doesn't really work. I will
    // manage my own DoFs so that a "node" lives at the centroid
    // of a cell.
    int p = 1;
    int dofpn = 1;
    DoFManager dofm(M, DoFManager::ManagerType::CG, p, dofpn);
    FESystem feSystem(dofm);
    auto[centroids, volumes] = CellCentroidsVolumes<NSD>(dofm);

    FVDofm fvDofm(dofm, NSD);

    Eigen::VectorXd phi_cells = Eigen::VectorXd::Constant(fvDofm.centroids.size(), 0.0);
    CellToNodeMap<NSD> cellToNodeMap(fvDofm);

    auto phi_func = [](coordinate<> &x) { return x(0)*x(1); };
    int idx{0};
    for (auto &x : fvDofm.centroids) {
        phi_cells(idx++) = phi_func(x);
    }

    Eigen::VectorXd exactSolution = Eigen::VectorXd::Constant(dofm.nNodes(), 0.0);
    Eigen::VectorXd error;
    idx = 0;
    for (auto &x : dofm.dof_nodes) {
        exactSolution(idx++) = phi_func(x);
    }

    cellToNodeMap.interpolateToNodes(phi_cells, feSystem.getSolution());

    SimulationOutput simulationOutput("output", BackendType::VTU);
    simulationOutput.captureFrame(feSystem,
                                  std::function<void(FESystem&, OutputFrame&)>(
                                  [&exactSolution, &error](auto &sys, auto &frame) {

                                      OutputData::DataType dt = OutputData::DataType::Scalar;
                                      OutputData::DataLocation dL = OutputData::DataLocation::Point;

                                      error = sys.getSolution() - exactSolution;

                                      frame.addData(std::make_shared<OutputData>(sys.getSolution(), "interpolated"));
                                      frame.addData(std::make_shared<OutputData>(exactSolution, "exact"));
                                      frame.addData(std::make_shared<OutputData>(error, "error"));

                                  }));

}