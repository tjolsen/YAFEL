//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/SmallVector.hpp"
#include "utils/parallel/yafel_parallel.hpp"
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;


int main()
{

    auto& TS = getGlobalScheduler();


    auto [task,fut] = TS.createTask([](){std::cout << "In Task" << std::endl;});

    TS.enqueue(task);

    fut.get();
    return 0;
}


