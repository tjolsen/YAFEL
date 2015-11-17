#ifndef _TENSOR_ITERATOR_HPP
#define _TENSOR_ITERATOR_HPP

#include "iterator_utils.hpp"
#include "et_tensor.hpp"
#include <tuple>
#include <string>


//forward declaration
template<unsigned DIM, unsigned RANK, typename dataType>
class Tensor;


template<typename TENSTYPE, unsigned DIM, unsigned RANK, typename dataType=double>
class tensor_iterator {
private:
  bool all_done;
  const TENSTYPE &_t;

  typename gens<RANK>::type sequence; // <-- empty struct w/ type seq<0,1,...RANK-1>
  typename type_gen<RANK>::type indices; // <-- a tuple with RANK 'unsigned int'
  
  template<int I>
  void inc_level(const seq<I> &) {
    if(++std::get<I>(indices)==DIM) {
      std::get<I>(indices)=0;
      inc_level(seq<I-1>());
    }
  }
  void inc_level(const seq<0> &) {
    if(++std::get<0>(indices) == DIM) {
      all_done = true;
    }
  }
  template<int I>
  void zero_level(const seq<I> &) {
    std::get<I>(indices)=0;
    zero_level(seq<I-1>());
  }
  void zero_level(const seq<0>&) {
    std::get<0>(indices)=0;
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
  
  template<int ...S, typename ...Args>
  dataType apply_indices(const seq<S...>&,
			 const std::tuple<Args...> &indices) const {
    return _t(std::get<S>(indices)...);
  }
  template<int ...S, typename ...Args>
  dataType& apply_indices(const seq<S...>&,
			 const std::tuple<Args...> &indices) {
    return _t(std::get<S>(indices)...);
  }

public:

  template<int ...S>
  tensor_iterator(const TENSTYPE &t) : all_done(false), _t(t)
  {
    static_assert(RANK > 0, 
		  "Must have at least rank 1 tensor_iterator."
		  "Something else also probably blew up.");
  }


  dataType operator*() const {
    return apply_indices(sequence,indices);
  }

  dataType& operator*() {
    return apply_indices(sequence,indices);
  }

  
  //void operator++()
  
  void next() {
    inc_level( seq<RANK-1>() );
  }
  
  void reset() {
    all_done = false;
    zero_level(seq<RANK-1>());
  }
  
  bool end() {
    return all_done;
  }

  std::string to_string() {
    return make_string(sequence);
  }
};


#endif
