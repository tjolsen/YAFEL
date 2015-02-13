#include "output/MatrixVisualization.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

YAFEL_NAMESPACE_OPEN
namespace MatrixVisualization {
  //========================================================================================
  void scatter_xy(const std::vector<double> &x, const std::vector<double> &y) {
    
    FILE *gp = popen("gnuplot --persist", "w");
    if(gp == NULL) {
      fprintf(stderr, "Popen failure\n");
      exit(1);
    }

    fprintf(gp, "set terminal wxt\n");
    fprintf(gp, "set style data lines\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "plot '-' using 1:2 with points pt 7 ps 0.5\n");
    
    auto x_it = x.begin();
    auto y_it = y.begin();
    for(x_it = x.begin(); x_it<x.end() && y_it<y.end(); ++x_it, ++y_it) {
      fprintf(gp, "%f %f\n", *x_it, *y_it);
    }
    
    pclose(gp);
  }

  //========================================================================================
  void contour_xyz(const std::vector<double> &x, 
		   const std::vector<double> &y,
		   const std::vector<double> &z) {

    FILE *gp = popen("gnuplot --persist", "w");
    if(gp == NULL) {
      fprintf(stderr, "Popen failure\n");
      exit(1);
    }  
    fprintf(gp, "set terminal wxt\n");
    fprintf(gp, "set pm3d map\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "splot '-' with pm3d\n");
    
    //write contour data....
    for(unsigned yi=0; yi<y.size(); ++yi) {
      for(unsigned xi=0; xi<x.size(); ++xi) {
	fprintf(gp, "%f %f %f\n", x[xi], y[yi], z[yi*y.size() + xi]);
      }
      fprintf(gp, "\n");
    }
    
    pclose(gp);
  }
  
  //========================================================================================
  void spy(sparse_coo &coo) {
    
    std::vector<double> x(coo.nnz(), 0);
    std::vector<double> y(coo.nnz(), 0);
    
    unsigned *rp = coo.get_row_ptr();
    unsigned *cp = coo.get_col_ptr();
    
    
    
    for(unsigned i=0; i<coo.nnz(); ++i) {
      x[i] = (double)cp[i];
      y[i] = (double)(coo.getRows() - rp[i]);
    }
    
    scatter_xy(x,y);
  }

  //========================================================================================
  void spy(const sparse_csr &csr) {
    sparse_coo coo(csr);
    spy(coo);
  }
  
  //========================================================================================
  void contour(const FullMatrix &Z) {
    
    std::vector<double> x(Z.getCols(),0.0);
    std::vector<double> y(Z.getRows(),0.0);
    std::vector<double> z(Z.getRows()*Z.getCols(),0.0);

    for(unsigned i=0; i<x.size(); ++i) {
      x[i] = i;
    }
    for(unsigned i=0; i<y.size(); ++i) {
      y[i] = i;
    }
    
    for(unsigned i=0; i<Z.getRows(); ++i) {
      for(unsigned j=0; j<Z.getCols(); ++j) {
	z[i*Z.getRows() + j] = Z(i,j);
      }
    }
    
    contour_xyz(x,y,z);
  }
  
} //end namespace MatrixVisualization
YAFEL_NAMESPACE_CLOSE
