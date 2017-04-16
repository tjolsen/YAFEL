//
// Created by tyler on 4/16/17.
//

#ifndef YAFEL_FESYSTEM_HPP
#define YAFEL_FESYSTEM_HPP

#include "yafel_globals.hpp"
#include <eigen3/Eigen/Sparse>

/**
 * \class FESystem
 * \brief Class to hold all of the data structures that comprise a FEM simulation
 */
class FESystem
{

public:
    inline auto &getGlobalTangent() { return global_tangent; }

    inline auto &getGlobalResidual() { return global_residual; }

protected:
    Eigen::SparseMatrix<double, Eigen::ColMajor> global_tangent;
    Eigen::VectorXd global_residual;

};

#endif //YAFEL_FESYSTEM_HPP
