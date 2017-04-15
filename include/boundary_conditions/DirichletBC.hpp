//
// Created by tyler on 4/11/17.
//

#ifndef YAFEL_DIRICHLETBC_HPP
#define YAFEL_DIRICHLETBC_HPP

#include "yafel_globals.hpp"
#include "utils/DoFManager.hpp"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <functional>
#include <type_traits>

YAFEL_NAMESPACE_OPEN

/**
 * \brief Class to represent (and apply) a dirichlet boundary condition.
 *
 */
class DirichletBC
{
public:

    template<typename T,
            typename=typename std::enable_if<std::is_fundamental<T>::value>::type>
    DirichletBC(const DoFManager &dofm, T bcval);

    template<typename T,
            typename=typename std::enable_if<std::is_constructible<std::function<double(const coordinate<> &, double)>, T>::value>::type>

    DirichletBC(const DoFManager &dofm, T &&bcval);

    template<int RCMajor,
            typename=typename std::enable_if<RCMajor == Eigen::ColMajor || RCMajor == Eigen::RowMajor>::type>
    void apply(Eigen::SparseMatrix<double, RCMajor>, Eigen::VectorXd &rhs);

    void selectByRegionID(int region_id, int component);

    template<typename Lambda>

    void selectByFunction(Lambda &&func, int component);

    // for debugging purposes
    inline double test_value(const coordinate<> &x = coordinate<>(), double t = 0) const { return value_func(x, t); }

    inline auto &getMask() { return bc_mask; }

private:
    const DoFManager &dofm;
    std::vector<bool> bc_mask;
    std::function<double(const coordinate<> &, double)>
    value_func;
};


//-------------------------------------------------
// Constructor implementation for fundamental types
//-------------------------------------------------
template<typename T, typename>
DirichletBC::DirichletBC(const DoFManager &dofm, T bcval)
        : dofm(dofm),
          bc_mask(dofm.dof_nodes.size() * dofm.dof_per_node, false)
{
    value_func = [val = double(bcval)](const coordinate<> &, double) { return val; };
}

//----------------------------------------------
// Constructor implementation for function types
//----------------------------------------------
template<typename T, typename>
DirichletBC::DirichletBC(const DoFManager &dofm, T &&bcval)
        :dofm(dofm), bc_mask(dofm.dof_nodes.size() * dofm.dof_per_node, false), value_func(bcval) {}


template<typename Lambda>
void DirichletBC::selectByFunction(Lambda &&func, int component)
{
    int N = dofm.dof_nodes.size();
    for (auto i : IRange(0, N)) {
        if (func(dofm.dof_nodes[i])) {
            bc_mask[dofm.dof_per_node * i + component] = true;
        }
    }
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_DIRICHLETBC_HPP
