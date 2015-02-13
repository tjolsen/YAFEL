#ifndef _MATRIXVISUALIZATION_HPP
#define _MATRIXVISUALIZATION_HPP

#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/FullMatrix.hpp"

YAFEL_NAMESPACE_OPEN

namespace MatrixVisualization {

  void spy(sparse_coo &coo);
  void spy(const sparse_csr &csr);

  void contour(const FullMatrix &Z);
  void contour(const Vector &X, const Vector &Y, const FullMatrix &Z);
  
}
YAFEL_NAMESPACE_CLOSE

#endif
