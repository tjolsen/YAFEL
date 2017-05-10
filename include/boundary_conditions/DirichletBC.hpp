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
            typename = typename std::enable_if<std::is_fundamental<T>::value>::type>
    DirichletBC(const DoFManager &dofm, T bcval, int comp = 0);

    template<typename T,
            typename = typename std::enable_if<std::is_constructible<std::function<double(coordinate<>, double)>, T>::value>::type>
    DirichletBC(const DoFManager &dofm, T &&bcval, int comp = 0);

    template<int RCMajor,
            typename=typename std::enable_if<RCMajor == Eigen::ColMajor || RCMajor == Eigen::RowMajor>::type>
    void apply(Eigen::SparseMatrix<double, RCMajor> &A, Eigen::VectorXd &rhs, double time = 0.0);

    void selectByRegionID(int region_id);

    template<typename Lambda>
    void selectByFunction(Lambda &&func);

    // for debugging purposes
    inline double test_value(const coordinate<> &x = coordinate<>(), double t = 0) const { return value_func(x, t); }

    inline auto &getNodes() { return bc_nodes; }

private:
    const DoFManager &dofm;
    int component;
    std::vector<int> bc_nodes;
    std::function<double(const coordinate<> &, double)> value_func;
};


//-------------------------------------------------
// Constructor implementation for fundamental types
//-------------------------------------------------
template<typename T, typename>
DirichletBC::DirichletBC(const DoFManager &dofm, T bcval, int comp)
        : dofm(dofm),
          component(comp),
          bc_nodes()
{
    value_func = [val = double(bcval)](const coordinate<> &, double) { return val; };
}

//----------------------------------------------
// Constructor implementation for function types
//----------------------------------------------
template<typename T, typename>
DirichletBC::DirichletBC(const DoFManager &dofm, T &&bcval, int comp)
        :dofm(dofm),
         component(comp),
         bc_nodes(),
         value_func(bcval) {}


template<typename Lambda>
void DirichletBC::selectByFunction(Lambda &&func)
{
    int N = dofm.dof_nodes.size();
    for (auto i : IRange(0, N)) {
        if (func(dofm.dof_nodes[i])) {
            bc_nodes.push_back(i);
        }
    }
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_DIRICHLETBC_HPP
