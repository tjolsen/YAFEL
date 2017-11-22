//
// Created by tyler on 11/8/17.
//

#ifndef YAFEL_PARDOT_HPP
#define YAFEL_PARDOT_HPP

#include "yafel_globals.hpp"
#include "utils/parallel/ReductionVariable.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/parfor.hpp"
#include "utils/SmallVector.hpp"
#include <Eigen/Core>
#include <type_traits>

YAFEL_NAMESPACE_OPEN

/**
 * \brief parallel dot() function for eigen vectorXd's
 * using the parfor infrastructure
 */
template<typename V1, typename V2, int BLK = 4>
auto pardot(Eigen::MatrixBase<V1> const& A, Eigen::MatrixBase<V2> const& B, int chunkSize = 1024) {

    if(A.rows() < 10*chunkSize) {
        return A.dot(B);
    }

    using T = std::decay_t<decltype(A(0))>;
    struct r_struct{T val{0};};

    SmallVector<ReductionVariable<r_struct>, 8> S(yafel::config::num_cores, r_struct{0});

    parfor(0, A.rows(), [&](auto i) {
        auto id = worker_global::worker_id;
        S[id].val += A(i)*B(i);
    },getGlobalScheduler(), chunkSize);

    T total{0};
    for(auto &s : S){
        total += s.val;
    }
    return total;
}


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_PARDOT_HPP
