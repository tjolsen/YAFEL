//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/parfor.hpp"
#include "utils/parallel/wait_all.hpp"
#include "utils/BasicTimer.hpp"
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;


int main()
{

    TaskScheduler TS(8);
    BasicTimer timer;
    int N = 1000000;
    std::size_t BS = 32;
    std::vector<std::size_t> partial_sums(TS.workers.size(), 0);

    /*
    std::vector<std::future<size_t>> futs;
    timer.tic();
    for(size_t i=0; i<N; ++i) {

        futs.push_back(TS.enqueue([](size_t i1, size_t i2) {
            size_t sum{0};
            for(size_t i=i1; i<i2; ++i) {
                sum += i;
            }

            return sum;
        }, i, i+N));

    }
    long unsigned int bigsum{0};
    for(auto &fut : futs) {
        bigsum += fut.get();
    }
    timer.toc();*/


    timer.tic();
    auto futs = parfor(0, N,
                       [&partial_sums,N](std::size_t idx) {
                           auto me = worker_global::worker_id;
                           std::size_t sum{0};
                           for (std::size_t i = idx; i < idx+N; ++i) {
                               sum += i;
                           }
                           partial_sums[me] += sum;
                       },
                       TS, BS);

    wait_all(futs);

    std::size_t bigsum{0};
    for(auto ps : partial_sums) {
        bigsum += ps;
    }
    timer.toc();

    std::cout << "Sum = " << bigsum << "  timer = "
              << timer.duration() << " ms" << std::endl;

    return 0;
}


