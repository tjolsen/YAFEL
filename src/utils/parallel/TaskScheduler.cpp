//
// Created by tyler on 10/9/17.
//

#include "utils/parallel/TaskScheduler.hpp"

YAFEL_NAMESPACE_OPEN

namespace worker_global {
thread_local int worker_id = 0;
}

namespace {

TaskScheduler globalTS(config::num_cores);

} //end anon namespace

TaskScheduler& getGlobalScheduler() {
    return globalTS;
}

YAFEL_NAMESPACE_CLOSE