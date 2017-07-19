//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/TaskScheduler.hpp"
#include "utils/BasicTimer.hpp"
#include <iostream>



using namespace yafel;
using std::cout;
using std::endl;


int main()
{

    TaskScheduler TS;
    BasicTimer timer;
    int N = 1000000;
    std::vector<std::future<int>> futs;
    timer.tic();
    for(int i=0; i<N; ++i) {

        futs.push_back(TS.enqueue([](int i1, int i2) {
            int sum{0};
            for(int i=i1; i<i2; ++i) {
                sum += i;
            }
            return sum;
        }, i, i+N));

    }
    long unsigned int bigsum{0};
    for(auto &fut : futs) {
        bigsum += fut.get();
    }
    timer.toc();

    std::cout << "Sum = " << bigsum << "  timer = "
              << timer.duration() << " ms" << std::endl;

    return 0;
}


