//
// Created by tyler on 3/14/17.
//

#include "yafel_globals.hpp"
#include "utils/SmallVector.hpp"
#include <iostream>


using namespace yafel;
using std::cout;
using std::endl;


int main()
{
    struct alignas(16) Tp
    {
        double a, b, c;

        Tp(double x) : a(x), b(x), c(x)
        {}
    };
    SmallVector<Tp, 2> V;

    V.push_back(1.0);
    V.push_back(2.0);

    int a = 0;

    V.push_back(3.0);
    V.push_back(4.0);

    cout << a << endl;

    for (auto x : V) {
        cout << x.a << endl;
    }

    return 0;
}


