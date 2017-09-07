//
// Created by tyler on 9/6/17.
//

#ifndef YAFEL_BARRIER_HPP_HPP
#define YAFEL_BARRIER_HPP_HPP

#include "yafel_globals.hpp"
#include <future>
#include <vector>


/**
 * \brief Wait on each of a vector of futures.
 * Since all are being waited on, the order is irrelevant.
 * Hence, scanning linearly through each is an acceptable order.
 *
 * Some utilities (eg: parfor) return a std::vector of futures.
 * this provides a synchronization primitive to wait until all
 * futures in the vector return.
 *
 * WARNING: This function is guaranteed to block until all
 * futures are complete. Do not use if there is a possibility
 * of a non-terminating task. (Also, don't abuse the TaskScheduler
 * with that. Spin up your own thread...)
 *
 * @tparam T future type
 * @param futures vector of future<T> objects
 */
template<typename T>
void wait_all(std::vector<std::future<T>> &futures) {
    for(auto &f : futures) {
        f.wait();
    }
}


#endif //YAFEL_BARRIER_HPP_HPP
