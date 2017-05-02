//
// Created by tyler on 5/2/17.
//

#ifndef YAFEL_DGASSEMBLY_HPP
#define YAFEL_DGASSEMBLY_HPP

#include "yafel_globals.hpp"
#include "element/ElementFactory.hpp"
#include "fe_system/FESystem.hpp"
#include "assembly/AssemblyRequirement.hpp"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>


YAFEL_NAMESPACE_OPEN


/**
 * Assembly procedure for Discontinuous Galerkin FEM.
 * The procedure involves looping over mesh faces and
 * computing fluxes between elements.
 *
 * @tparam Physics
 * @param feSystem
 * @param requirements
 */
template<typename Physics>
void DGAssembly(FESystem &feSystem,
                std::vector<AssemblyRequirement> requirements = {AssemblyRequirement::Residual})
{

    constexpr int simulation_dimension = Physics::nsd();

    //Unpack the FESystem
    auto &GlobalResidual = feSystem.getGlobalResidual();
    auto &dofm = feSystem.getDoFManager();
    auto dof_per_node = dofm.dof_per_node;
    double time = feSystem.currentTime();

    bool assemble_tangent{false};
    bool assemble_residual{false};
    bool assemble_dt_mass{false};
    bool assemble_dtdt_mass{false};
    for (auto req : requirements) {
        switch (req) {
            case AssemblyRequirement::Residual:
                assemble_residual = true;
                break;
            case AssemblyRequirement::Tangent:
                assemble_tangent = true;
                break;
            case AssemblyRequirement::DtMass:
                assemble_dt_mass = true;
                break;
            case AssemblyRequirement::DtDtMass:
                assemble_dtdt_mass = true;
                break;
        }
    }

#pragma omp parallel
    {
        //Define thread-local variables

        //storage buffers
        std::vector<double> local_tangent_buffer;
        std::vector<double> local_dtmass_buffer;
        std::vector<double> local_residual_buffer;
        std::vector<int> global_dof_buffer;



        // Boundary Fluxes

        for()


        // Element-level fluxes


    }


}




YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_DGASSEMBLY_HPP
