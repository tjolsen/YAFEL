#include "output/MatrixVisualization.hpp"
#include <cstdio>
#include <cstdlib>

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
}
YAFEL_NAMESPACE_CLOSE
