#ifndef _GENERIC_TENSOR_ITERATOR_HPP
#define _GENERIC_TENSOR_ITERATOR_HPP

#include "iterator_utils.hpp"

template<unsigned DIM, unsigned RANK>
class generic_tensor_iterator {

protected:
  typename gens<RANK>::type sequence;
  typename type_gen<RANK>::type indices;

private:
  bool done;

  void inc(const seq<0> &) {
    if(++std::get<0>(indices) == DIM) {
      done = true;
    }
  }
  
  template <int I>
  void inc(const seq<I> &) {
    if(++std::get<I>(indices) == DIM) {
      std::get<I>(indices) = 0;
      inc(seq<I-1>());
   } 
  }

  template <int I>
  void zero_level(const seq<I> &) {
    std::get<I>(indices) = 0;
    zero_level(seq<I-1>());
  }
  void zero_level(const seq<0>&) {
    std::get<0>(indices) = 0;
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


  generic_tensor_iterator() : done(false)
  {
    static_assert(RANK>0, "Must have at least rank 1 tensor");
    reset();
  }
  
  void next() {
    inc(seq<RANK-1>());
  }
  
  void reset() {
    done = false;
    zero_level(seq<RANK-1>());
  }
  
  bool end() {
    return done;
  }

  std::string to_string() {
    return make_string(sequence);
  }

}; //end class

#endif
