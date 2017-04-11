//
// Created by tyler on 4/11/17.
//


#include "output/OutputData.hpp"

YAFEL_NAMESPACE_OPEN


OutputData::OutputData(const Eigen::VectorXd &d, const std::string &name, const std::vector<int> &mask)
        : data(&d), name(name), dof_mask(mask) {}

OutputData::OutputData(const Eigen::VectorXd &d, const std::string &name)
        : OutputData(d, name, {1}) {}


YAFEL_NAMESPACE_CLOSE