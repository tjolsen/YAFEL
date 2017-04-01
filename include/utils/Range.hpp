#ifndef YAFEL_RANGE_HPP
#define YAFEL_RANGE_HPP

#include "yafel_globals.hpp"
#include <type_traits>

YAFEL_NAMESPACE_OPEN

template<typename T>
class Range
{
public:
    template<typename U, typename = typename std::enable_if<std::is_convertible<U, T>::value>::type>
    Range(U first, U last, U inc = 1, bool inclusive = false)
            : first(first), last(last+inclusive*inc), inc(inc)
    {}

    inline T getFirst() const { return first; }

    inline T getLast() const { return last; }

    inline T getInc() const { return inc; }

private:
    T first;
    T last;
    T inc;
};


template<typename T>//, typename=typename std::enable_if<std::is_integral<T>::value>::type>
class RangeIterator
{
public:
    RangeIterator() = delete;

    RangeIterator(T val, T inc) : val(val), inc(inc) {}

    RangeIterator<T> operator++()
    {
        val += inc;
        return *this;
    }

    T operator*()
    {
        return val;
    }

    bool operator<(RangeIterator<T> rhs) const
    {
        return (val - *rhs) / inc < 0;
    }

    bool operator==(RangeIterator<T> rhs) const
    {
        return val == *rhs;
    }

    bool operator!=(RangeIterator<T> rhs) const
    {
        return *this < rhs;
    }


private:
    T val;
    T inc;
};


template<typename T>
RangeIterator<T> begin(Range<T> R)
{
    return RangeIterator<T>(R.getFirst(), R.getInc());
}

template<typename T>
RangeIterator<T> end(Range<T> R)
{
    return RangeIterator<T>(R.getLast(), R.getInc());
}


//Some Typedefs
using IRange = Range<int>;


YAFEL_NAMESPACE_CLOSE

#endif
