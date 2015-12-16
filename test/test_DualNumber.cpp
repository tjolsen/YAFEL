#include "utils/DualNumber.hpp"

#include <iostream>
#include <cmath>

using namespace yafel;
template class yafel::DualNumber<double>; // help lcov generate useful results

template <typename T>
std::ostream &operator<<(std::ostream &os, const yafel::DualNumber<T> &A) {
    os << "DualNumber(" << A.first << ", " << A.second << ")";
    return os;
}

bool (*tests[])() = {
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        yafel::DualNumber<double> y; 
        auto f = x + y;
        if (f.first != 5 || f.second != 1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = x + 5;
        if (y.first != 10 || y.second != 1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = x - 5;
        if (y.first != 0 || y.second != 1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = x * 2;
        if (y.first != 10 || y.second != 2) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = x / 2;
        if (y.first != 2 || y.second != 0.5) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        yafel::DualNumber<double> y(5, 0.5); 
        return (x < y);
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        yafel::DualNumber<double> y(5, 0.5); 
        return (y > x);
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = -x;
        if (y.first != -4 || y.second != -1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = 5 + x;
        if (y.first != 10 || y.second != 1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = -(5 - x);
        if (y.first != 0 || y.second != 1) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(5, 1); 
        auto y = 2 * x;
        if (y.first != 10 || y.second != 2) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = 1.0 / (2.0 / x);
        if (y.first != 2 || y.second != 0.5) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = cos(x);
        if (y.first != cos(4.0) || y.second != -sin(4.0)) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = sin(x);
        if (y.first != sin(4.) || y.second != cos(4.0)) {
            return false;
        }
        return true;
    },
    []() -> bool {
        yafel::DualNumber<double> x(4, 1); 
        auto y = exp(x);
        if (y.first != exp(4.0) || y.second != exp(4.0)) {
            return false;
        }
        return true;
    },
};

int run_test(int which_test)
{
    bool pass = (*tests[which_test])();
    return pass?(0):(1 << which_test);
}

int main()
{
    int retval = 0;
    const size_t num_tests = sizeof(tests) / sizeof(tests[0]);
    for (size_t i = 0; i < num_tests; i++) {
        retval |= run_test(i);
    }
    return retval;

    return 0;
}
