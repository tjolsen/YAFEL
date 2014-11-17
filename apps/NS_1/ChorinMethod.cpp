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
  EF = ElementFactory(M, NSD);
}

// ========= SETUP ========================================
void ChorinMethod::setup() {
  
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
  u_sol = Vector(EF.get_n_dof(), 0.0);
  P_sol = Vector(M.get_n_nodes(), 0.0);
}

void ChorinMethod::run() {
  
  
  double time = 0;
  double dt = Tfinal/100; // BS VALUE: WRITE ROUTINE TO ESTIMATE (Mesh will need to return a "characteristic length")
  
  //set up output structures
  VTKMesh vtkm(EF);

  while(time < Tfinal) {
    if((Tfinal - time) < dt) {
      dt = Tfinal-time;
    }
    time += dt;
    
    bool goodstep = false;
    while(!goodstep) {

      goodstep = timestep();
      
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
void ChorinMethod::assemble() {






}


// =========== TIMESTEP =====================================
bool ChorinMethod::timestep() {

  
  return true;
}
