//
// Created by tyler on 9/6/17.
//

#ifndef YAFEL_WAIT_ANY_HPP
#define YAFEL_WAIT_ANY_HPP

#include "yafel_globals.hpp"
#include <future>
#include <vector>
#include <chrono>

/**
 * \brief Wait for any of a vector of futures.
 * Returns a reference to the first (encountered) completed future.
 *
 * Some utilities (eg: parfor) return a std::vector of futures.
 * this provides a synchronization primitive to wait until all
 * futures in the vector return.
 *
 * @tparam T future type
 * @param futures vector of future<T> objects
 */
template<typename T>
auto& wait_any(std::vector<std::future<T>> &futures) {
    while(true) {
        for (auto &f : futures) {
            if(f.wait_for(std::chrono::seconds(0)) == std::future_status::ready) {
                return f;
            }
        }
    }
}


#endif //YAFEL_WAIT_ANY_HPP
