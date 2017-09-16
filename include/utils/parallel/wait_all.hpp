//
// Created by tyler on 9/6/17.
//

#ifndef YAFEL_WAIT_ALL_HPP
#define YAFEL_WAIT_ALL_HPP

#include "yafel_globals.hpp"
#include <future>
#include <initializer_list>
#include <vector>

YAFEL_NAMESPACE_OPEN

namespace detail{

template<typename Item>
inline void wait(Item &&item) {
    item.wait();
}

}

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
    for(auto &f: futures) {
        f.wait();
    }
}


/**
 * \brief An implementation of wait all that allows you
 * to pass in an arbitrary number of (potentially heterogeneous)
 * std::futures.
 */
inline void wait_all () {}
template<typename Item, typename ...Items>
inline void wait_all(Item &&item, Items &&...items) {
    detail::wait(item);
    wait_all(items...);
}

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_WAIT_ALL_HPP
