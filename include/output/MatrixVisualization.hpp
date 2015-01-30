#ifndef _MATRIXVISUALIZATION_HPP
#define _MATRIXVISUALIZATION_HPP

#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

namespace MatrixVisualization {

  void spy(sparse_coo &coo);
  void spy(const sparse_csr &csr);
  
}
YAFEL_NAMESPACE_CLOSE

#endif
