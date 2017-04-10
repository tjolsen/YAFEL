#ifndef _YAFEL_MATRIXVISUALIZATION_HPP
#define _YAFEL_MATRIXVISUALIZATION_HPP


/*
 * MatrixVisualization:
 *
 * This file provides functions for runtime visualization of matrices.
 * The functions are self-contained, and they are contained with the
 * MatrixVisualization namespace.
 *
 */ 

#include "yafel_globals.hpp"
#include "old_handmade_linalg/sparse_matrix.hpp"
#include <vector>
#include <sstream>
#include <tuple>
#include <exception>


YAFEL_NAMESPACE_OPEN

namespace MatrixVisualization {
  
    template<typename T, typename dataType>
    void spy(sparse_matrix<T,dataType> &A);
  
  
    template<typename dataType>
    void scatter_xy(const std::vector<dataType> &x, 
                    const std::vector<dataType> &y);

    //void scatter_xy(const Vector &x, 
    //const Vector &y);
    //  void contour_xyz(const std::vector<double> &x, 
    //		   const std::vector<double> &y,
    //		   const std::vector<double> &z);
    //

    //  //  void contour(const FullMatrix &Z);
    //  //  void contour(const Vector &X, const Vector &Y, const FullMatrix &Z);
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
        if(std::get<2>(t) != 0) { // <-- semantics of "spy" func are to look at nonzeros
            ss << std::get<1>(t) << " " << std::get<0>(t) << "\n";
        }
    }

    fprintf(gp, "%s", ss.str().c_str());

    pclose(gp);
}



template<typename dataType>
void MatrixVisualization::scatter_xy(const std::vector<dataType> &x,
                                     const std::vector<dataType> &y)
{
    FILE *gp = popen("gnuplot --persist", "w");
    if(gp == nullptr) {
        throw(std::runtime_error("Could not open gnuplot pipe"));
    }
    
    fprintf(gp, "reset\n");
    fprintf(gp, "set terminal wxt\n");
    fprintf(gp, "set style data lines\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "plot '-' using 1:2 with points pt 7 ps 0.5\n");
    
    auto x_it = x.begin();
    auto y_it = y.begin();
    for(x_it = x.begin(); x_it<x.end() && y_it<y.end(); ++x_it, ++y_it) {
        fprintf(gp, "%f %f\n", *x_it, *y_it);
    }
    fprintf(gp, "e\n");
    
}




YAFEL_NAMESPACE_CLOSE

#endif
