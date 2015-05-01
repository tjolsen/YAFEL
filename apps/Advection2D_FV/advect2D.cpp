#include "yafel.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <unistd.h>

using namespace yafel;

/*

2D Advection equation on periodic LxL grid

|---------------------------------------|(Nx,Ny)
|                                       |
|                                       |
|                                       |
(i,0)                 (i,j)             |
|                                       |
|                                       |
|                                       |
|                                       |
|                                       |
|(0,0)________________(0,j)_____________|

Index of (i,j) = (i*Nx) + j

 */

int main() {
  //parameters to set
  unsigned Nx = 500;
  unsigned Ny = 500;
  double L = 1.0;
  double T = 1.0;
  double cx = 1.0;  //ci = speed in ith direction
  double cy = 1.0;
  double omega = 4;

  //computed parameters
  double dx = L/Nx;
  double dy = L/Ny;
  double pi = 2*asin(1);

  FullMatrix U1(Ny, Nx, 0.0);
  
  //initial conditions
  for(unsigned i=0; i<Ny; ++i) {
    for(unsigned j=0; j<Nx; ++j) {
      U1(i,j) = sin(omega*pi*j*dx)*sin(omega*pi*i*dy);
      //1U(i,j) = (i>2*Ny/5) && (i<3*Ny/5) && (j>2*Nx/5) && (j<3*Nx/5);
    }
  }
  FullMatrix U2(U1);

  FullMatrix &U = U1;
  FullMatrix &Uold = U2;
  
  //set stable time step
  double dt = 0.5*std::min(dx/cx, dy/cy);
  std::cout << "dt = " << dt << std::endl;
  MatrixVisualization MV;

  MV.contour(U);

  double t=0;

  unsigned iter = 0;
  while(t < T) {
    ++iter;
#pragma omp parallel for schedule(static,50)
    for(unsigned i=0; i<Ny; ++i) {
	unsigned n = (i+1 + Ny)%Ny;
	unsigned nn = (i+2 + Ny)%Ny;
	unsigned s = (i-1 + Ny)%Ny;
	unsigned ss = (i-2 + Ny)%Ny;

      for(unsigned j=0; j<Nx; ++j) {
	double Fxe, Fxw, Fyn, Fys;
	
	unsigned e = (j+1 + Nx)%Nx;
	unsigned ee = (j+2 + Nx)%Nx;
	unsigned w = (j-1 + Nx)%Nx;
	unsigned ww = (j-2 + Nx)%Nx;
	
	/*
	//use upwind fluxes
	Fxe = (cx>0) ? cx*U(i,j) : cx*U(i,e);
	Fxw = (cx>0) ? cx*U(i,w) : cx*U(i,j);
	Fyn = (cy>0) ? cy*U(i,j) : cy*U(n,j);
	Fys = (cy>0) ? cy*U(s,j) : cy*U(i,j);
	*/
	
	//use QUICK fluxes
	Fxe = (cx>0) ? 
	  cx*(0.75*U(i,j) + 0.375*U(i,e) - 0.125*U(i,w)) :
	  cx*(0.75*U(i,e) + 0.375*U(i,j) - 0.125*U(i,ee));
	
	Fxw = (cx>0) ? 
	  cx*(0.75*U(i,w) + 0.375*U(i,j) - 0.125*U(i,ww)) :
	  cx*(0.75*U(i,j) + 0.375*U(i,w) - 0.125*U(i,e));
	
	Fyn = (cy>0) ? 
	  cy*(.75*U(i,j) + 0.375*U(n,j) - 0.125*U(s,j)) :
	  cy*(.75*U(n,j) + 0.375*U(i,j) - 0.125*U(nn,j));
	
	Fys = (cy>0) ?
	  cy*(0.75*U(s,j) + 0.375*U(i,j) - 0.125*U(ss,j)) :
	  cy*(0.75*U(i,j) + 0.375*U(s,j) - 0.125*U(n,j));

	U(i,j) = Uold(i,j) + dt*((Fxw-Fxe)/dx + (Fys-Fyn)/dy);
      }
    }
    
    if(iter % 10 == 0) {
      std::cout << "Calling contour" << std::endl;
      MV.contour(U);
    }
    

    FullMatrix &Utmp = Uold;
    Uold = U;
    U = Utmp;
    t += dt;
  }
  
  
  return 0;
}
