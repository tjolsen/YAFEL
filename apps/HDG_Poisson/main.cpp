//
// Created by tyler on 5/22/17.
//

#include "yafel.hpp"
#include "Poisson.hpp"

using namespace yafel;

int main()
{
    constexpr int NSD = 2;
    Mesh M("mesh.msh");

    int polyOrder = 12;
    int dof_per_volume_node = NSD + 1;

    DoFManager bulk_dofm(M, DoFManager::ManagerType::DG, polyOrder, dof_per_volume_node);

}