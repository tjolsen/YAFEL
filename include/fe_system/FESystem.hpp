//
// Created by tyler on 4/16/17.
//

#ifndef YAFEL_FESYSTEM_HPP
#define YAFEL_FESYSTEM_HPP

#include "yafel_globals.hpp"
#include "utils/DoFManager.hpp"
#include <eigen3/Eigen/Sparse>

YAFEL_NAMESPACE_OPEN

/**
 * \class FESystem
 * \brief Class to hold all of the data structures that comprise a FEM simulation
 */
class FESystem
{

public:
    inline FESystem(DoFManager &dofm, int dim = 0)
            : dofm(dofm), simulation_dimension(dim)
    {
        int ndofs = dofm.dof_nodes.size() * dofm.dof_per_node;
        global_residual.resize(ndofs);
        global_tangent.resize(ndofs, ndofs);
    }

    inline auto &getGlobalTangent() { return global_tangent; }

    inline auto &getGlobalResidual() { return global_residual; }

    inline auto &getDoFManager() { return dofm; }

    inline auto &getDimension() { return simulation_dimension; }

    inline auto &getTime() {return time; }
protected:
    DoFManager &dofm;
    Eigen::SparseMatrix<double, Eigen::ColMajor> global_tangent;
    Eigen::VectorXd global_residual;

    int simulation_dimension;
    double time;
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_FESYSTEM_HPP
