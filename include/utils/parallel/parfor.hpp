//
// Created by tyler on 7/31/17.
//

#ifndef YAFEL_PARFOR_HPP
#define YAFEL_PARFOR_HPP

#include "yafel_globals.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include <vector>

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
auto parfor(std::size_t idx_start,
            std::size_t idx_end,
            Lambda &&loopBody,
            TaskScheduler &scheduler,
            std::size_t blockSize = 32)
{
    static_assert(std::is_same<void, typename std::result_of<Lambda(std::size_t)>::type>::value,
                  "Loop body must return void");

    auto loopLen = (idx_end - idx_start);
    std::size_t Nblocks = loopLen / blockSize + ((loopLen % blockSize == 0) ? 0 : 1);

    //lambda to execute a chunk of the loop
    auto loopChunk = [&loopBody](auto istart, auto iend) {
        for (auto i = istart; i < iend; ++i) {
            loopBody(i);
        }
    };


    std::vector<std::future<void>> futs;
    futs.reserve(Nblocks);

    for (std::size_t iblock = 0; iblock < Nblocks; ++iblock) {
        futs.push_back(
                scheduler.enqueue(loopChunk,
                                  iblock * blockSize,
                                  std::min((iblock + 1) * blockSize, idx_end))
        );

    }


    return futs;
}


YAFEL_NAMESPACE_CLOSE

#endif //PARTICLEDYN_PARFOR_HPP
