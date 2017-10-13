//
// Created by tyler on 9/29/17.
//

#ifndef YAFEL_YAFEL_CONFIG_HPP
#define YAFEL_YAFEL_CONFIG_HPP

#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

namespace config {

/**
 * Number of bytes per cache line. 64 is a reasonable guess,
 * but configure to taste.
 */
constexpr int cacheline_bytes = 64;

/**
 * Number of cores on the machine. Used in initializing the
 * global TaskScheduler.
 */
#ifndef NUM_CPU_CORES
#define NUM_CPU_CORES 8
constexpr int num_cores = 8;
#else
constexpr int num_cores = NUM_CPU_CORES;
#endif

} //end namespace config

YAFEL_NAMESPACE_CLOSE

#endif //YAFEL_YAFEL_CONFIG_HPP
