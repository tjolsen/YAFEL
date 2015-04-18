#ifndef _MATRIXVISUALIZATION_HPP
#define _MATRIXVISUALIZATION_HPP

#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/FullMatrix.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

class MatrixVisualization {
private:
  FILE *gp;

public:

  MatrixVisualization();
  ~MatrixVisualization();
  void spy(sparse_coo &coo);
  void spy(const sparse_csr &csr);

  void contour(const FullMatrix &Z);
  void contour(const Vector &X, const Vector &Y, const FullMatrix &Z);
  
  void scatter_xy(const std::vector<double> &x, 
		  const std::vector<double> &y);
  void contour_xyz(const std::vector<double> &x, 
		   const std::vector<double> &y,
		   const std::vector<double> &z);
};

YAFEL_NAMESPACE_CLOSE

#endif
