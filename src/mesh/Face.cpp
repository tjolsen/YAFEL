#include "mesh/Face.hpp"

YAFEL_NAMESPACE_OPEN


//constructors
Face::Face() : Face(std::vector<size_type>(2,0), 0, 0, false)
{}

Face::Face(const std::vector<size_type> &n, size_type i, size_type o, bool b)
  : nodes(n), inner(i), outer(o), boundary(b)
{}

//printer
std::ostream & operator<<(std::ostream& out, const Face &F) {
  
  out << "Face{ n={";
  for(std::size_t i=0; i<F.nodes.size(); ++i) {
    out << F.nodes[i];
    if(i<F.nodes.size()-1)
      out <<",";
  }
  out << "} i=" << F.inner 
      << " o=" << F.outer 
      << " b=" << (F.boundary?"true":"false")
      << " }";

  return out;
}


YAFEL_NAMESPACE_CLOSE
