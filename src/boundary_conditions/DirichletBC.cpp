//
// Created by tyler on 4/15/17.
//

#include "boundary_conditions/DirichletBC.hpp"
#include <algorithm>

YAFEL_NAMESPACE_OPEN

void DirichletBC::selectByRegionID(int region_id, int component)
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
void DirichletBC::apply(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::VectorXd &rhs)
{

    int *row_ptr = A.outerIndexPtr();
    int *col_ptr = A.innerIndexPtr();
    double *value_ptr = A.valuePtr();

    std::vector<bool> bc_mask(dofm.dof_nodes.size()*dofm.dof_per_node,false);
    for(auto n : bc_nodes) {
        
    }


    for(auto r : IRange(0,A.rows())) {
        for(auto idx : IRange(row_ptr[r], row_ptr[r+1])) {
            int c = col_ptr[idx];
            if(bc_mask[r]) {

            }
        }
    }
}

template<>
void DirichletBC::apply(Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::VectorXd &rhs)
{


}

YAFEL_NAMESPACE_CLOSE