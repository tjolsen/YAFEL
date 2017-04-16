//
// Created by tyler on 4/16/17.
//

#ifndef YAFEL_CGASSEMBLY_HPP
#define YAFEL_CGASSEMBLY_HPP

#include "yafel_globals.hpp"
#include "element/Element.hpp"
#include "fe_system/FESystem.hpp"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>

YAFEL_NAMESPACE_OPEN

enum class AssemblyRequirement
{
    Residual,
    Tangent,
    DtMass,
    DtDtMass
};

template<typename Physics>
constexpr bool hasLocalResidual()
{
    return std::is_same<void, typename std::result_of<decltype(&Physics::LocalResidual)(Element const &, double, Eigen::Map<double>&)>::type>::value;
}

template<typename Physics>
constexpr bool hasLocalTangent()
{
    return std::is_same<void, typename std::result_of<decltype(&Physics::LocalTangent)(Element const &, double, Eigen::Map<double>&)>::type>::value;
}






/**
 *
 * \brief General-purpose Continuous-Galerkin finite element assembly
 *
 * @tparam Physics Class that defines the local element matrix/vector construction in static void methods
 */
template<typename Physics>
void CGAssembly(FESystem &feSystem)
{
    static_assert(hasLocalResidual<Physics>(), "Must Implement static void method 'LocalResidual' in Physics class");
    static_assert(hasLocalTangent<Physics>(), "Must Implement static void method 'LocalTangent' in Physics class");


    // Unpack teh FESystem
    auto &GlobalTangent = feSystem.getGlobalTangent();
    auto &GlobalResidual = feSystem.getGlobalResidual();
    auto &dofm = feSystem.getDoFManager();


    //storage for local element matrix/vectors
    std::vector<double> local_tangent_buffer;
    std::vector<double> local_residual_buffer;


    for(auto elnum : IRange(0,dofm.nCells())) {

    }

}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_CGASSEMBLY_HPP
