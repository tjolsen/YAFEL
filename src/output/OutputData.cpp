//
// Created by tyler on 4/11/17.
//


#include "output/OutputData.hpp"

YAFEL_NAMESPACE_OPEN


OutputData::OutputData(const Eigen::VectorXd &d,
                       const std::string &name,
                       DataLocation dL,
                       DataType dT,
                       const std::vector<int> &mask)
        : data(&d), name(name), dataLocation(dL), dataType(dT),dof_mask(mask) {}

OutputData::OutputData(const Eigen::VectorXd &d, const std::string &name)
        : OutputData(d, name, DataLocation::Point, DataType::Scalar, {1}) {}


YAFEL_NAMESPACE_CLOSE