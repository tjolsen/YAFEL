#ifndef _MATRIXVISUALIZATION_HPP
#define _MATRIXVISUALIZATION_HPP


/*
 * MatrixVisualization:
 *
 * This file provides functions for runtime visualization of matrices.
 * The functions are self-contained, and they are contained with the
 * MatrixVisualization namespace.
 *
 */ 

#include "yafel_globals.hpp"
#include "lin_alg/sparse_matrix.hpp"
#include <vector>
#include <sstream>
#include <tuple>
#include <exception>


YAFEL_NAMESPACE_OPEN

namespace MatrixVisualization {
  
  template<typename T, typename dataType>
  void spy(sparse_matrix<T,dataType> &A);
  
  
//  //  void contour(const FullMatrix &Z);
//  //  void contour(const Vector &X, const Vector &Y, const FullMatrix &Z);
//  
//  void scatter_xy(const std::vector<double> &x, 
//		  const std::vector<double> &y);
//  //void scatter_xy(const Vector &x, 
//  //		  const Vector &y);
//  void contour_xyz(const std::vector<double> &x, 
//		   const std::vector<double> &y,
//		   const std::vector<double> &z);
//


} // end namespace MatrixVisualization



/*
 * Template function definition
 */
template<typename T, typename dataType>
void MatrixVisualization::spy(sparse_matrix<T,dataType> &A) {
  auto triplets = A.copy_triplets();

  FILE *gp = popen("gnuplot --persist", "w");
  if(gp == nullptr) {
    throw(std::runtime_error("Could not open gnuplot pipe"));
  }
  
  std::stringstream ss;
  ss << "reset\n" 
     << "set terminal wxt\n"
     << "set style data lines\n"
     << "set xrange [0:"<<A.cols()-1<<"]\n"
     << "set yrange [0:"<<A.rows()-1<<"] reverse\n"
     << "set size square\n"
     << "set key off\n"
     << "plot '-' using 1:2 with points pt 7 ps 0.5\n";

  for(auto t : triplets) {
    ss << std::get<1>(t) << " " << std::get<0>(t) << "\n";
  }

  fprintf(gp, "%s", ss.str().c_str());

  pclose(gp);
}

YAFEL_NAMESPACE_CLOSE

#endif
