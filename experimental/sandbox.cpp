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

    auto [task, fut] = TS.createTask([](){cout << "Parent!\n";});

    auto [kid, kidfut] = task->addChild([](){cout << "Child!\n"; return 1;});

    auto [kid2, kidfut2] = task->addChild([](){cout << "Grandchild!\n"; return 1;});

    TS.enqueue(task);

    wait_all(fut,kidfut,kidfut2);

    return 0;
}


