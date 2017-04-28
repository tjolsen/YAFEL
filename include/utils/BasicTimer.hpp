//
// Created by tyler on 4/28/17.
//

#ifndef YAFEL_BASICTIMER_HPP
#define YAFEL_BASICTIMER_HPP

#include "yafel_globals.hpp"
#include <chrono>

YAFEL_NAMESPACE_OPEN


class BasicTimer
{

public:

    inline void tic() { tic_point = std::chrono::high_resolution_clock::now(); }

    inline void toc() { toc_point = std::chrono::high_resolution_clock::now(); }

    template<typename T = std::chrono::milliseconds>
    auto duration() const { return std::chrono::duration_cast<T>(toc_point - tic_point).count(); }

private:

    std::chrono::high_resolution_clock::time_point tic_point;
    std::chrono::high_resolution_clock::time_point toc_point;
};


YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_BASICTIMER_HPP
