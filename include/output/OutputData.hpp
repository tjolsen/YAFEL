//
// Created by tyler on 4/10/17.
//

#ifndef YAFEL_OUTPUTDATA_HPP
#define YAFEL_OUTPUTDATA_HPP

#include "yafel_globals.hpp"

#include <eigen3/Eigen/Core>
#include <vector>
#include <string>

YAFEL_NAMESPACE_OPEN

/**
 * /class OutputData
 * /brief Wrapper around an Eigen::VectorXd containing metadata relevant for output
 *
 * Very little data integrity checking is (or can be) done at this level,
 * so the user is expected to provide a well-formed data vector.
 * The data is considered well-formed if it has \f$dof_mask.size()*N\f$
 * components, where \f$N\f$ is the number of nodes (or cells) in the mesh.
 *
 * The dof_mask
 */
class OutputData
{
public:
    enum class DataType : int {
        Scalar,
        Vector,
        Tensor
    };

    enum class DataLocation {
        Point,
        Cell
    };

    OutputData(const Eigen::VectorXd &d, const std::string &name);

    OutputData(const Eigen::VectorXd &d, const std::string &name, const std::vector<int> &mask);


    const Eigen::VectorXd *data;
    std::vector<int> dof_mask;
    std::string name;
    DataType dataType;
    DataLocation dataLocation;
};

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_OUTPUTDATA_HPP
