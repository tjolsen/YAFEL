#ifndef __YAFEL_LINEARADVECTION_HPP
#define __YAFEL_LINEARADVECTION_HPP

/*
 * Linear Advection:
 *
 *
 * This header sets up the program structure for a linear advection solver.
 * It also defines a utility class AdvectionParameters that holds 
 * the problem parameters. This is convenient for purposes of instantiating
 * a "LinearAdvection" class.
 */

#include "yafel.hpp"
#include <vector>
#include <string>

YAFEL_NAMESPACE_OPEN

class AdvectionParameters {
public:
  // Mesh parameters
  std::vector<double> mesh_dims;
  std::vector<std::size_t> dir_elems;
  
  // Advection Velocity
  Tensor<2,1,double> v_advection;

  // DG parameters
  std::size_t polyOrder;
  QuadratureRule<2> & Q2D; //volume quadrature
  QuadratureRule<1> & Q1D; //face quadrature

  // Time parameters
  double dt;
  double Tfinal;

  //initial conditions
  SpatialFunction<2,double> u0_func;
  
  // Output parameters
  std::string output_file_base;
};



//======================================================================
//======================================================================
//======================================================================
class LinearAdvection {

public:
  
  LinearAdvection()=delete;
  LinearAdvection(const AdvectionParameters &AP);


  // Problem steps
  void run();
  void setup();
  void RK4();

  
  Vector<double> set_initial_condition();

  //return the solved residual M\R (block diagonal, so should be easy!)
  Vector<double> residual(const Vector<double> &u);

  //return the element residual, not solved by mass matrix.
  //this will be handled in calling function
  Vector<double> element_residual(const Vector<double> & u_elem);

  //function used to define numerical fluxes at element boundaries
  Tensor<2,1,double> flux_function(double u_in, double u_out, 
                                   Tensor<2,1,double> c, Tensor<2,1,double> n) const;

  //calculate element Mass, convection matrices
  void set_Me_Se();

  //write a timestep out to disk
  void write_output(const Vector<double> &u, double ti);



  // Member variables
  RectilinearMesh<2> M;
  AdvectionParameters AP;
  DG_DoFManager<RectilinearMesh<2>,2> dofm;
  DG_Quad<> DGQ;
  Matrix<double> Me;
  Matrix<double> Se;
  LUDecomposition<double> LU_Me;
};


YAFEL_NAMESPACE_CLOSE


#endif
