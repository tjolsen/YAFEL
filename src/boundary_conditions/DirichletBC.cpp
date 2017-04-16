//
// Created by tyler on 4/15/17.
//

#include "boundary_conditions/DirichletBC.hpp"
#include <algorithm>

YAFEL_NAMESPACE_OPEN

void DirichletBC::selectByRegionID(int region_id)
{
    int N = dofm.nCells();
    std::vector<int> cellNodes;
    for (auto c : IRange(0, N)) {
        if (dofm.cell_region_idx[c] == region_id) {
            dofm.getGlobalNodes(c, cellNodes);
            for (auto n : cellNodes) {
                bc_nodes.push_back(n);
            }
        }
    }
    auto it = std::unique(bc_nodes.begin(), bc_nodes.end());
    bc_nodes.resize(std::distance(bc_nodes.begin(), it));
}


template<>
void DirichletBC::apply(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &rhs, double time)
{

    int *row_ptr = A.outerIndexPtr();
    int *col_ptr = A.innerIndexPtr();
    double *value_ptr = A.valuePtr();
    Eigen::VectorXd bc_values = Eigen::VectorXd::Constant(dofm.dof_per_node * dofm.dof_nodes.size(), 0.0);
    std::vector<bool> bc_mask(dofm.dof_nodes.size() * dofm.dof_per_node, false);
    for (auto n : bc_nodes) {
        bc_mask[n * dofm.dof_per_node + component] = true;
        bc_values(n * dofm.dof_per_node + component) = value_func(dofm.dof_nodes[n], time);
    }

    rhs -= A * bc_values;

    for (auto r : IRange(0, static_cast<int>(A.rows()))) {
        for (auto idx : IRange(row_ptr[r], row_ptr[r + 1])) {
            int c = col_ptr[idx];

            if (bc_mask[r] || bc_mask[c]) {
                value_ptr[idx] = 0;
                if (r == c) {
                    rhs(r) = bc_values(c);
                    value_ptr[idx] = 1;
                }
            }
        }
    }
}

template<>
void DirichletBC::apply(Eigen::SparseMatrix<double, Eigen::ColMajor> &A, Eigen::VectorXd &rhs, double time)
{

    int *col_ptr = A.outerIndexPtr();
    int *row_ptr = A.innerIndexPtr();
    double *value_ptr = A.valuePtr();
    Eigen::VectorXd bc_values = Eigen::VectorXd::Constant(dofm.dof_per_node * dofm.dof_nodes.size(), 0.0);
    std::vector<bool> bc_mask(dofm.dof_nodes.size() * dofm.dof_per_node, false);
    for (auto n : bc_nodes) {
        bc_mask[n * dofm.dof_per_node + component] = true;
        bc_values(n * dofm.dof_per_node + component) = value_func(dofm.dof_nodes[n], time);
    }

    rhs -= A * bc_values;

    for (auto c : IRange(0, static_cast<int>(A.cols()))) {
        for (auto idx : IRange(col_ptr[c], col_ptr[c + 1])) {
            int r = row_ptr[idx];

            if (bc_mask[c] || bc_mask[r]) {
                value_ptr[idx] = 0;
                if (r == c) {
                    rhs(r) = bc_values(c);
                    value_ptr[idx] = 1;
                }
            }
        }
    }


}

YAFEL_NAMESPACE_CLOSE