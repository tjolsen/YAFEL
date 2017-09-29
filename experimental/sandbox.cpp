//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/parallel/TaskScheduler.hpp"
#include "utils/parallel/parfor.hpp"
#include "utils/parallel/ReductionVariable.hpp"
#include "utils/BasicTimer.hpp"
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;


int main()
{
    struct alignas(32) my_reducers {
        double x[8];
        //int y;
    };

    ReductionVariable<my_reducers> RV;

    cout << sizeof(RV) << endl << alignof(RV) << endl;

    return 0;
}


