//
// Created by tyler on 9/11/17.
//
#include "yafel.hpp"
#include "cell_centroids.hpp"

using namespace yafel;



int main() {

    Mesh M("channel.msh");
    constexpr int NSD = 2;


    // p=1 DoFManager, since p=0 doesn't really work. I will
    // manage my own DoFs so that a "node" lives at the centroid
    // of a cell.
    DoFManager dofm(M, DoFManager::ManagerType::CG, 1, NSD+1);
    FESystem feSystem(dofm);

    auto [centroids, volumes] = CellCentroidsVolumes<NSD>(dofm);

    int i{0};
    for(auto c : centroids) {

        std::cout << c(0) << ' ' << c(1) << ' ' << volumes[i++] << '\n';

    }
}