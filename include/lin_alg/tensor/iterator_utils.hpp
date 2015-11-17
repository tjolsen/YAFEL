#ifndef _ITERATOR_UTILS
#define _ITERATOR_UTILS

#include <tuple>

// Sequence generator
template <int ...S>
struct seq{};

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> {};

template<int ...S>
struct gens<0, S...> {
  typedef seq<S...> type;
};


// Tuple-type generator
template<int N, typename ...Args>
struct type_gen : type_gen<N-1, unsigned, Args...> {};

template<typename ...Args>
struct type_gen<0, Args...> {
  typedef std::tuple<Args...> type;
};

#endif
