#include "LinearAdvection.hpp"

#include <cmath>

//==============================================================
LinearAdvection::LinearAdvection(const AdvectionParameters &ap)
  : M(ap.mesh_dims, ap.dir_elems),
    AP(ap),
    dofm(M, 1)
    DGQ(ap.polyOrder, 1, ap.Q2D, ap.Q1D)
    Me(DGQ.nodes_per_element, DGQ.nodes_per_element),
    Se(DGQ.nodes_per_element, DGQ.nodes_per_element),
    LU_Me(Me)
{}


//==============================================================
void LinearAdvection::run() {

  setup();

  //do time integration (handles output writing)
  RK4();
  
}

//==============================================================
Tensor<2,1,double> 
LinearAdvection::flux_function(double u_in, double u_out,
                               Tensor<2,1,double> c, Tensor<2,1,double> n,
                               bool is_boundary) const
{
  if(is_boundary) {
    //boundary conditions handled here
    if(contract<1>(c,n) < 0) {
      //inflow
      return Tensor<2,1,double>{0,0};
    }
    else {
      //outflow
      return u_in*c;
    }
  }
  else {
  //Upwind flux
  return c*((u_in + u_out)/2) + (1.0/2.0)*std::abs(contract<1>(c,n))*(u_out-u_in)*n;
  }
}

//==============================================================
void LinearAdvection::setup()
{
  
  M.build_faces();
  set_Me_Se();
  LU_Me.reinit(Me);

}

//==============================================================
void LinearAdvection::RK4()
{

  Vector<double> u = set_initial_condition();
  
  write_output(u, 0);
  
  //number of time steps
  std::size_t Nt = std::size_t(T/dt); //round down. no point being more precise now
  
  for(std::size_t ti=1; ti<=Nt; ++ti) {
    
    Vector<double> k1 = dt*residual(u);
    Vector<double> k2 = dt*residual(u + 0.5*k1);
    Vector<double> k3 = dt*residual(u + 0.5*k2);
    Vector<double> k4 = dt*residual(u + k3);

    u += (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    
    write_output(u, ti);
  }

}



//==============================================================
Vector<double> LinearAdvection::residual(const Vector<double> &u) {
  
  std::size_t nnodes = DGQ.nodes_per_element;
  
  Vector<double> u_elem(nnodes);
  
  for(std::size_type elnum=0; elnum < M.n_elements(); ++e) {
  
    

  }
}


//==============================================================
void LinearAdvection::set_Me_Se()
{
  
  //use first element to compute
  DGQ.update_element(M,0);

  Me *= 0;
  Se *= 0;

  for(std::size_t qpi=0; qpi<DGQ.Q2D.n_qp(); ++qpi) {
    double jxw = DGQ.detJ[qpi]*DGQ.Q2D.weight(qpi);

    for(std::size_t A=0; A<DGQ.nodes_per_element; ++A) {
      for(std::size_t B=0; B<DGQ.nodes_per_element; ++B) {
        Me(A,B) += DGQ.shape_values[qpi][A]*DGQ.shape_values[qpi][B]*jxw;

        for(std::size_t i=0; i<2; ++i) {
          Se(A,B) += DGQ.shape_gradients[qpi](A,i)*AP.v_advect(i)*DGQ.shape_values[qpi][B]*jxw;
        }
      }
    }
  }

}

//==============================================================
void LinearAdvection::write_output(const Vector<double> &u, double ti)
{

  


}
