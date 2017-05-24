//
// Created by tyler on 5/23/17.
//

#ifndef YAFEL_POISSON_HPP
#define YAFEL_POISSON_HPP

#include "yafel_globals.hpp"
#include <eigen3/Eigen/Dense>

/**
 * Element matrix and local solve computation for
 * an HDG implementation of the Poisson equation.
 */
template<int NSD>
struct Poisson : PDEBase<NSD>
{

    template<typename T>
    static void LocalAssembly(
            const Element &E,
            Eigen::DenseBase<T> &VolVol,
            Eigen::DenseBase<T> &VolEdge,
            Eigen::DenseBase<T> &EdgeVol,
            Eigen::DenseBase<T> &EdgeEdge,
            Eigen::DenseBase<T> &ElemMatrix,
            Eigen::DenseBase<T> &R1,
            Eigen::DenseBase<T> &R2,
            Eigen::DenseBase<T> &R3);

    static void LocalSolve();

protected:
    static void SchurComplement();
};


//---------------------------------------------
// Impelementation
//---------------------------------------------

template<typename T>
static void LocalAssembly(
        const Element &E,
        Eigen::DenseBase<T> &VolVol,
        Eigen::DenseBase<T> &VolEdge,
        Eigen::DenseBase<T> &EdgeVol,
        Eigen::DenseBase<T> &EdgeEdge,
        Eigen::DenseBase<T> &ElemMatrix,
        Eigen::DenseBase<T> &R1,
        Eigen::DenseBase<T> &R2,
        Eigen::DenseBase<T> &R3)
{


}


#endif //YAFEL_POISSON_HPP
