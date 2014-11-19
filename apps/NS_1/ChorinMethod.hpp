#ifndef _CHORINMETHOD_HPP
#define _CHORINMETHOD_HPP

#include <string>
#include <vector>
#include "yafel.hpp"

YAFEL_NAMESPACE_OPEN

class ChorinMethod {
  
private:
  std::string meshFname;
  std::string outputBase;
  unsigned NSD;
  double Tfinal;
  
  std::vector<Vector> u_out;
  std::vector<double> P_out;
  Vector u_sol;
  Vector P_sol;
  std::vector<int> bcnodes;
  std::vector<int> bccomp;
  std::vector<double> bcvals;

  Mesh M;
  ElementFactory EFu, EFp;
  
  //=======================
  bool timestep();
  void assemble_u();
  void assemble_p();
  void setup();

public:
  ChorinMethod(const char *meshFileName, const char *obase, unsigned nsd, double Tf);
  void run();
  
};

YAFEL_NAMESPACE_CLOSE

#endif
