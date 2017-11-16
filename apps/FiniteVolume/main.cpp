//
// Created by tyler on 9/11/17.
//
#include "yafel.hpp"
#include "cell_centroids.hpp"

using namespace yafel;



int main() {

    BasicTimer timer;
    timer.tic();
    Mesh M("channel.msh");
    timer.toc();

    std::cout << "Time to read mesh: " << timer.duration<>() << " ms\n";

    timer.tic();
    M.buildInternalFaces();
    timer.toc();
    std::cout << "Internal face time: " << timer.duration<>() << "ms\n";

    constexpr int NSD = 2;


    // p=1 DoFManager, since p=0 doesn't really work. I will
    // manage my own DoFs so that a "node" lives at the centroid
    // of a cell.
    DoFManager dofm(M, DoFManager::ManagerType::CG, 1, NSD+1);
    FESystem feSystem(dofm);


    timer.tic();
    auto [centroids, volumes] = CellCentroidsVolumes<NSD>(dofm);
    timer.toc();


    std::cout << timer.duration<std::chrono::milliseconds>() << std::endl;
}