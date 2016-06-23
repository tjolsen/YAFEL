#ifndef _POISSON_HPP
#define _POISSON_HPP

#include "yafel.hpp"
#include <string>
#include <vector>

YAFEL_NAMESPACE_OPEN

class Poisson {

private:
    GmshMesh M;
    DoFManager DOFM;
    ElementFactory EF;
    DirBC BC,BC2;
    Vector<double> Usol;
    Vector Fsys;
    sparse_coo Kcoo;
    std::vector<int> bcnodes;
    std::vector<double> bcvals;
    double fvol;
  

    void setup();//const std::string &inputFilename);
    void assemble();
    void solve();
    void output(const std::string &outputFilename);

public:
    Poisson(const char *fname);
    void run(const std::string & outputFilename);
};

YAFEL_NAMESPACE_CLOSE

#endif
