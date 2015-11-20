#ifndef __YAFEL_GENERIC_INDEX_ITERATOR_HPP
#define __YAFEL_GENERIC_INDEX_ITERATOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/iterator_utils.hpp"

YAFEL_NAMESPACE_OPEN

template<unsigned DIM, unsigned RANK, int... IDXS>
class generic_index_iterator {

protected:
  typename gens<RANK>::type sequence;
  typename type_gen<RANK>::type indices;

private:
  bool done;

  template <int I>
  void inc(const seq<I> &) {
    if(++std::get<I>(indices) == DIM) {
      done = true;
    }
  }
  
  template <int I, int... II>
  void inc(const seq<I, II...> &) {
    if(++std::get<I>(indices) == DIM) {
      std::get<I>(indices) = 0;
      inc(seq<II...>());
   } 
  }

  template <int I>
  void zero_level(const seq<I> &) {
    std::get<I>(indices) = 0;
  }
  
  template <int I, int ...II>
  void zero_level(const seq<I,II...> &) {
    std::get<I>(indices) = 0;
    zero_level(seq<II...>());
  }
  
  template<int ...S>
  std::string make_string(const seq<S...> &) {
    return std::string("(") + stringify(std::get<S>(indices)...) + std::string(")");
  }
  template<typename TT, typename ...Args>
  std::string stringify(TT t, Args ...args) {
    return std::to_string(t) + std::string(", ") + stringify(args...);
  }
  template<typename TT>
  std::string stringify(TT t) {
    return std::to_string(t);
  }

public:
  template<typename ...Args>
  generic_index_iterator(Args... _indices) : indices(std::make_tuple(_indices...)), done(false)
  {
    static_assert(RANK>0, "Must have at least rank 1 tensor");
    static_assert(sizeof...(Args) == RANK, "RANK/nargs mismatch");
    reset();
  }
  
  void next() {
    inc(seq<IDXS...>());
  }
  
  void reset() {
    done = false;
    zero_level(seq<IDXS...>());
  }
  
  bool end() {
    return done;
  }

  std::string to_string() {
    return make_string(sequence);
  }

}; //end class


YAFEL_NAMESPACE_CLOSE

#endif
