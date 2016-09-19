#ifndef _YAFEL_GENERIC_TENSOR_ITERATOR_HPP
#define _YAFEL_GENERIC_TENSOR_ITERATOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/tensor/iterator_utils.hpp"

YAFEL_NAMESPACE_OPEN

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
	    set_level(seq<RANK-1>(),DIM); //special "end" value where all slots are set to DIM
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
    void set_level(const seq<I> &, unsigned val) {
	std::get<I>(indices) = val;
	set_level(seq<I-1>(), val);
    }
    void set_level(const seq<0>&, unsigned val) {
	std::get<0>(indices) = val;
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
    generic_tensor_iterator(unsigned val=0) : done(false)
	{
	    static_assert(RANK>0, "Must have at least rank 1 tensor");
	    reset(val);
	}
  
    void next() {
	inc(seq<RANK-1>());
    }
  
    void reset(int val) {
	done = false;
	set_level(seq<RANK-1>(),val);
    }
  
    bool end() {
	return done;
    }

    const typename type_gen<RANK>::type & get_indices() const {
	return indices;
    }

    typename type_gen<RANK>::type & get_indices() {
	return indices;
    }

    bool operator!=(const generic_tensor_iterator<DIM,RANK>&rhs) {
	static_assert(RANK>0 && DIM>0, "Must have nonzero dimension, rank");

	return indices != rhs.get_indices();
    }

    generic_tensor_iterator<DIM,RANK>& operator++() {
	next();
	return *this;
    }
    
    std::string to_string() {
	return make_string(sequence);
    }
    
}; //end class


YAFEL_NAMESPACE_CLOSE

#endif
