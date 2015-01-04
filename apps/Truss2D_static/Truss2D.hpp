#include "yafel.hpp"
#include <string>
#include <vector>
#include <cmath>

using namespace yafel;

class Truss2D {

private:
  std::string MeshFilename;
  std::string OutputFilename;
  std::vector<int> bcnodes;
  std::vector<int> bccomps;
  std::vector<double> bcvals;
  double Eyoungs;
  double Axsection;
  double rho;
  Mesh M;
  sparse_coo Kcoo;
  Vector Fsys;
  Vector Usol;

public:
  Truss2D(const char *Mfname, const char *outFname);
  void setup();
  void assemble();
  void solve();
  void output();
  void run();
};
