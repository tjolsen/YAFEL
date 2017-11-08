//
// Created by tyler on 7/31/17.
//

#ifndef YAFEL_PARFOR_HPP
#define YAFEL_PARFOR_HPP

#include "yafel_globals.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/wait_all.hpp"
#include <vector>
#include <iostream>

YAFEL_NAMESPACE_OPEN


template<typename Lambda>
struct parfor_body
{
    Lambda loopBody;
};

/**
 * Break a "for" loop with known extent into a statically-scheduled
 * set of tasks. Loop blocks must be able to be run simultaneously
 * and/or out of order, while loop iterations within a block are
 * guaranteed to run sequentially.
 *
 * The loop body function must take a single std::size_t parameter,
 * which is the loop index.
 *
 * A std::vector of loop block futures is returned to the caller.
 *
 * @tparam Lambda
 * @param scheduler
 * @param blockSize
 */
template<typename Lambda>
void parfor(std::size_t idx_start,
            std::size_t idx_end,
            Lambda &&loopBody,
            TaskScheduler &scheduler,
            std::size_t blockSize = 32)
{
    static_assert(std::is_same<void, typename std::result_of<Lambda(std::size_t)>::type>::value,
                  "Loop body must return void");

    auto loopLen = (idx_end - idx_start);

    //Execute short loops in serial
    if(loopLen < 10*blockSize) {
        for(std::size_t i=0; i<loopLen; ++i) {
            loopBody(i);
        }
        return;
    }

    auto full_blocks = loopLen / blockSize;
    auto cleanup_start = full_blocks * blockSize + idx_start;
    auto Nblocks = full_blocks + ((cleanup_start < idx_end) ? 1 : 0);


    std::vector<std::future<void>> futs;
    futs.reserve(Nblocks);

    for (std::size_t iblock = 0; iblock < Nblocks; ++iblock) {

        std::size_t i_start = idx_start + iblock * blockSize;
        std::size_t i_end = i_start + blockSize;
        i_end = (i_end < idx_end) ? i_end : idx_end;

        auto [task, fut] = scheduler.createTask([&loopBody, i_start, i_end]() {
            for (auto i = i_start; i < i_end; ++i) {
                loopBody(i);
            }
        });
        scheduler.enqueue(task);
        futs.push_back(std::move(fut));

    }

    wait_all(futs);
    //return futs;
}


YAFEL_NAMESPACE_CLOSE

#endif //PARTICLEDYN_PARFOR_HPP
