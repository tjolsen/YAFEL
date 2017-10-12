//
// Created by tyler on 5/8/17.
//

#ifndef YAFEL_DGRK4_HPP
#define YAFEL_DGRK4_HPP

#include "yafel_globals.hpp"
#include "fe_system/FESystem.hpp"
#include "assembly/DGAssembly.hpp"
#include <eigen3/Eigen/Dense>

YAFEL_NAMESPACE_OPEN

/**
 * \class DGRK4
 * \brief Class to perform RK4 integration,
 * designed for use with explicit DG methods (which provide M^{-1}*R quickly)
 *
 * Consider creating a base timestepping class that provides an interface like:
 * - initialize()
 * - step()
 */
class DGRK4
{

public:
    DGRK4() = delete; //need to supply a timestep
    DGRK4(double dt);


    template<typename Physics>
    void step(FESystem &feSystem, Physics &P)
    {
        auto U0 = feSystem.getSolution();
        double time0 = feSystem.currentTime();
        DGAssembly(feSystem, P,
                   {AssemblyRequirement::Residual,
                    AssemblyRequirement::DtMass});
        auto k1 = feSystem.getGlobalResidual();


        feSystem.getSolution() = U0 + (dt / 2) * k1;
        feSystem.currentTime() = time0 + dt / 2;
        DGAssembly(feSystem, P,
                   {AssemblyRequirement::Residual,
                    AssemblyRequirement::DtMass});
        auto k2 = feSystem.getGlobalResidual();

        feSystem.getSolution() = U0 + (dt / 2) * k2;
        feSystem.currentTime() = time0 + dt / 2;
        DGAssembly(feSystem, P,
                   {AssemblyRequirement::Residual,
                    AssemblyRequirement::DtMass});

        auto k3 = feSystem.getGlobalResidual();


        feSystem.getSolution() = U0 + dt * k3;
        feSystem.currentTime() = time0 + dt;
        DGAssembly(feSystem, P,
                   {AssemblyRequirement::Residual,
                    AssemblyRequirement::DtMass});
        auto k4 = feSystem.getGlobalResidual();

        feSystem.getSolution() = U0 + (dt / 6) * k1 + (dt / 3) * k2 + (dt / 3) * k3 + (dt / 6) * k4;

        feSystem.currentTime() = time0 + dt;
    }


private:
    double dt;

};


YAFEL_NAMESPACE_CLOSE


#endif //YAFEL_DGRK4_HPP
