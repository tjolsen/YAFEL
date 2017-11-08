//
// Created by tyler on 11/8/17.
//

#ifndef YAFEL_PARASSIGN_HPP
#define YAFEL_PARASSIGN_HPP


#include "yafel_globals.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/parfor.hpp"
#include <eigen3/Eigen/Core>
#include <type_traits>

YAFEL_NAMESPACE_OPEN

/**
 * \brief parallel assignment function for eigen vectorXd's
 * using the parfor infrastructure.
 * Warning: Assumes no aliasing of A and B.
 */
template<typename V1, typename V2, int BLK = 4>
auto parassign(Eigen::MatrixBase<V1> & A, Eigen::MatrixBase<V2> const& B, int chunkSize = 1024) {


    if(A.rows() < 10*chunkSize) {
        A.noalias() = B;
        return A;
    }

    parfor(0, A.rows(), [&](auto i) {
        A(i) = B(i);
    },getGlobalScheduler(), chunkSize);

    return A;
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_PARASSIGN_HPP
