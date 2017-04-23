//
// Created by tyler on 4/17/17.
//

#ifndef YAFEL_ELEMENT_BOUNDARY_NODES_HPP
#define YAFEL_ELEMENT_BOUNDARY_NODES_HPP

#include "yafel_globals.hpp"
#include "yafel_typedefs.hpp"
#include <vector>
#include <iostream>

YAFEL_NAMESPACE_OPEN

template<typename Lambda>
std::vector<int> make_element_boundary_nodes(
        const std::vector<coordinate<>> &xi_all,
        Tensor<3, 1, double> boundary_normal,
        Tensor<3, 1, double> edge_direction_1,
        Tensor<3, 1, double> edge_direction_2,
        Lambda &&boundary_mark)
{
    int npts = static_cast<int>(xi_all.size());
    std::vector<int> boundary_idxs;
    std::vector<int> idx_perm;
    std::vector<coordinate<>> xi_boundary;
    int idx{0};
    for (auto i : IRange(0, npts)) {
        if (boundary_mark(xi_all[i])) {
            boundary_idxs.push_back(i);
            idx_perm.push_back(idx++);
            xi_boundary.push_back(xi_all[i]);
        }
    }


    //generate new basis
    Tensor<3, 2> basis;
    basis(colon(), 0) = edge_direction_1;
    basis(colon(), 1) = edge_direction_2;
    basis(colon(), 2) = boundary_normal;


    Tensor<3, 2> Ainv = inverse(basis);

    for (auto &x : xi_boundary) {
        x = (Ainv * x);
        x(2) = 0;
    }

    std::sort(idx_perm.begin(), idx_perm.end(),
              [&xi_boundary](int L, int R) -> bool {
                  return std::lexicographical_compare(xi_boundary[L].data.rbegin(), xi_boundary[L].data.rend(),
                                                      xi_boundary[R].data.rbegin(), xi_boundary[R].data.rend());
              });

    for (auto &i : idx_perm) {
        i = boundary_idxs[i];
    }
    return idx_perm;
}


YAFEL_NAMESPACE_CLOSE
#endif //YAFEL_ELEMENT_BOUNDARY_NODES_HPP
