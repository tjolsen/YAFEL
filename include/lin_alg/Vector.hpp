#ifndef _YAFEL_VECTOR_HPP
#define _YAFEL_VECTOR_HPP

#include "yafel_globals.hpp"
#include "lin_alg/VectorExpression.hpp"

YAFEL_NAMESPACE_OPEN

template<typename dataType=double>
class Vector: public VectorExpression<Vector<dataType>, dataType> {
  
public:
    typedef typename VectorExpression<Vector<dataType>,dataType>::container_type container_type;
    typedef typename VectorExpression<Vector<dataType>,dataType>::value_type value_type;
    typedef typename VectorExpression<Vector<dataType>,dataType>::size_type size_type;
    typedef typename VectorExpression<Vector<dataType>,dataType>::reference reference;


    reference operator()(size_type i) {return _data[i];}
    value_type operator()(size_type i) const {return _data[i];}
    size_type size() const {return _data.size();}
  
    /*
     * Constructors
     */
    // Construct Vector of length N
    Vector(size_type N) : _data(N) {}

    // Construct Vector of length N and fill with value
    Vector(size_type N, value_type val) : _data(N, val) {}
  
    // Construct from VectorExpression<T,dataType>
    template<typename T>
    Vector(const VectorExpression<T,dataType> & v) {
        _data.resize(v.size());
        for(size_type i=0; i<v.size(); ++i) {
            _data[i] = v(i);
        }
    }


    /*
     * Vector-specific operators. Cannot be used with VectorExpressions as
     * left-hand argument.
     */
    template<typename T>
    Vector<dataType>& operator+=(const VectorExpression<T,dataType> &rhs) {
#ifndef _OPTIMIZED
        assert(size() == rhs.size() && "Vector::operator+= dimension mismatch");
#endif

        constexpr size_type UNROLL = 8;

        auto N = size()/UNROLL;
        auto NN = UNROLL*N;
        for(size_type i=0; i<NN; i+=UNROLL) {
            _data[i] += rhs(i);
            _data[i+1] += rhs(i+1);
            _data[i+2] += rhs(i+2);
            _data[i+3] += rhs(i+3);
            _data[i+4] += rhs(i+4);
            _data[i+5] += rhs(i+5);
            _data[i+6] += rhs(i+6);
            _data[i+7] += rhs(i+7);
        }
        for(size_type i=NN; i<size(); ++i) {
            _data[i] += rhs(i);
        }
    
        return *this;
    }

    template<typename T>
    Vector<dataType>& operator-=(const VectorExpression<T,dataType> &rhs) {
#ifndef _OPTIMIZED
        assert(size() == rhs.size() && "Vector::operator-= dimension mismatch");
#endif
        for(size_type i=0; i<size(); ++i) {
            _data[i] -= rhs(i);
        }
    
        return *this;
    }

    template<typename T>
    Vector<dataType>& operator*=(const VectorExpression<T,dataType> &rhs) {
#ifndef _OPTIMIZED
        assert(size() == rhs.size() && "Vector::operator*= dimension mismatch");
#endif
        for(size_type i=0; i<size(); ++i) {
            _data[i] *= rhs(i);
        }
    
        return *this;
    }

    Vector<dataType>& operator*=(dataType rhs) {


        constexpr size_type UNROLL = 8;
      
        auto N = size()/UNROLL;
        auto NN = UNROLL*N;
        for(size_type i=0; i<NN; i+=UNROLL) {
            _data[i] *= rhs;
            _data[i+1] *= rhs;
            _data[i+2] *= rhs;
            _data[i+3] *= rhs;
            _data[i+4] *= rhs;
            _data[i+5] *= rhs;
            _data[i+6] *= rhs;
            _data[i+7] *= rhs;
        }
        for(size_type i=NN; i<size(); ++i) {
            _data[i] *= rhs;
        }

    
        return *this;
    }




private:
    container_type _data;

};






YAFEL_NAMESPACE_CLOSE

#endif
