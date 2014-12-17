#include "ChorinMethod.hpp"
#include "yafel.hpp"

using namespace yafel;

//=========== CONSTRUCTOR ================================
ChorinMethod::ChorinMethod(const char *meshFileName, const char *obase, unsigned nsd, double Tf) :
  meshFname(std::string(meshFileName)), outputBase(std::string(obase)),
  NSD(nsd), Tfinal(Tf) 
{
  //read mesh
  M = MeshReader::gmsh_read(meshFname.c_str());
  EFu = ElementFactory(M, NSD);
  EFp = ElementFactory(M, 1);
}

// ========= SETUP ========================================
void ChorinMethod::setup() {

  //set viscosity
  viscosity = 1.0;

  //set density
  density = 1000;
  
  //set velocity dirichlet boundary conditions
  for(unsigned i=0; i<M.el_tags.size(); ++i) {
    if(M.el_tags[i][0] == 100) {
      for(unsigned j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomp.push_back(0);
	bcvals.push_back(0.0);
      }
    }
    if(M.el_tags[i][0] == 100) {
      for(unsigned j=0; j<M.elements[i].size(); ++j) {
	bcnodes.push_back(M.elements[i][j]);
	bccomp.push_back(1);
	bcvals.push_back(0.0);
      }
    }
  }

  // set pressure boundary conditions (Dirichlet-like on corrector step)
  for(unsigned i=0; i<M.el_tags.size(); ++i) {
    
    
  }

  // set initial conditions for fields
  u_sol = Vector(EFu.get_n_dof(), 0.0);
  P_sol = Vector(EFp.get_n_dof(), 0.0);
}

void ChorinMethod::run() {
  
  
  double time = 0;
  double dt = Tfinal/100; // BS VALUE: WRITE ROUTINE TO ESTIMATE (Mesh will need to return a "characteristic length")
  
  //set up output structures
  VTKMesh vtkm(EFu);

  while(time < Tfinal) {
    if((Tfinal - time) < dt) {
      dt = Tfinal-time;
    }
    time += dt;
    
    bool goodstep = false;
    while(!goodstep) {

      goodstep = timestep(dt);
      
      //use return value from timestep() to adaptively change dt
      //if residuals are diverging... not sure how to do this yet,
      //but seems like a nice idea, so I'll leave room for it.
      if(!goodstep) {
	dt /= 2;
      }
      else {
	//good time step: copy solution vectors
	u_out.clear();
	for(unsigned i=0; i<M.get_n_nodes(); ++i) {
	  Vector v(3,0.0);
	  for(unsigned j=0; j<NSD; ++j) {
	    v(j) = u_sol(i*NSD + j);
	  }
	  u_out.push_back(v);

	}
      }            
    }

    VTKVectorData vtku(u_out, VTKObject::VTKPOINTDATA, std::string("Velocity"));
    VTKScalarData vtkp(P_sol, VTKObject::VTKPOINTDATA, std::string("Pressure"));
    
    //write output for this timestep
    VTKOutput VO;
    VO.addVTKObject(&vtkm);
    VO.addVTKObject(&vtku);
    VO.addVTKObject(&vtkp);
    char buf[256];
    sprintf(buf, "%s_t%04.4f.vtu", outputBase.c_str(), time);
    VO.write(buf);
  }
  
  
}

// ========= ASSEMBLE =======================================
void ChorinMethod::assemble_u(double dt) {
  

  for(unsigned elnum=0; elnum < M.get_n_elems(); ++elnum) {
    
    Element *e = EFu.getElement(elnum);
    if(e==NULL) continue;
    if(e->n_spaceDim != NSD) continue;
  
    e->update_element(M,elnum);
    
    unsigned Ndof = e->dof_per_el;
    FullMatrix K_el(Ndof, Ndof, 0.0);
    FullMatrix M_el(Ndof, Ndof, 0.0);
    Vector R_el(Ndof, 0.0);
    
    for(unsigned qpi = 0; qpi < e->n_quadPoints; ++qpi) {
      
      Vector qp = e->quad_points[qpi];
      
      //get local u values from global vector
      Vector u_loc(Ndof, 0.0);
      for(unsigned A=0; A<Ndof; ++A) {
	u_loc(A) = u_sol(e->global_dofs[A]);
      }
      
      //compute velocies at quadrature point
      Vector uqp(e->n_spaceDim, 0.0);
      for(unsigned A=0; A<Ndof; ++A) {
	int comp = e->getComp(A);
	uqp(comp) += e->vals[qpi](e->getBase(A)) * u_loc(A);
      }
      
      //compute velocity gradient at quadpoint
      FullMatrix ugradqp(e->n_spaceDim, e->n_spaceDim, 0.0);
      for(unsigned A=0; A<Ndof; ++A) {
	int abase = e->getBase(A);
	int i = e->getComp(A);
	for(unsigned j=0; j<e->n_spaceDim; ++j) {
	  ugradqp(i,j) += e->grads[qpi](abase, j)*u_loc(A);
	}
      }

      Vector gradu_u = ugradqp * uqp;
      
      // Element tangent
      for(unsigned A=0; A<Ndof; ++A) {
	unsigned i=e->getComp(A);
	unsigned abase = e->getBase(A);
	for(unsigned B=0; B<Ndof; ++B) {
	  unsigned j=e->getComp(B);
	  unsigned bbase = e->getBase(B);

	  K_el(A,B) += (dt*viscosity/density)*(e->grads[qpi](A,i))*(e->grads[qpi](B,j))*(e->JxW(qpi));

	  M_el(A,B) += (e->vals[qpi](abase))*(e->vals[qpi](bbase))*(e->JxW(qpi));
	}
      }
      
      // Element Residual
      for(unsigned A=0; A<Ndof; ++A) {
	unsigned n = e->getBase(A);
	unsigned i = e->getComp(A);
	
	R_el(A) += (uqp(i) - dt*gradu_u(i))*(e->vals[qpi](n))*(e->JxW(qpi));
      }
      
    } //end: quadpoints loop
    

    //assemble into global
    FullMatrix Kfull_el = M_el - K_el;
    
    
    
  }// end: elements loop
  
}

void ChorinMethod::assemble_p(double dt) {
  
  
  
  
  
  
}


// =========== TIMESTEP =====================================
bool ChorinMethod::timestep(double dt) {

  
  return true;
}
