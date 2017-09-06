//
// Created by tyler on 9/6/17.
//

#ifndef YAFEL_BARRIER_HPP_HPP
#define YAFEL_BARRIER_HPP_HPP

#include "yafel_globals.hpp"
#include <future>
#include <vector>

template<typename T>
void barrier(std::vector<std::future<T>> &futures) {
    for(auto &f : futures) {
        f.wait();
    }
}


#endif //YAFEL_BARRIER_HPP_HPP
