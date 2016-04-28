#include "mesh/Face.hpp"

YAFEL_NAMESPACE_OPEN


//constructors
Face::Face() : Face(std::vector<size_type>(2,0), 0, 0, 0, 0, false)
{}

Face::Face(const std::vector<size_type> &n, size_type _i, size_type _if,
           size_type _o, size_type _of, bool _b)
    : nodes(n), 
      inner(_i),
      inner_face(_if),
      outer(_o),
      outer_face(_of),
      boundary(_b)
{}

//printer
std::ostream & operator<<(std::ostream& out, const Face &F) {
  
    out << "Face{ n={";
    for(std::size_t i=0; i<F.nodes.size(); ++i) {
        out << F.nodes[i];
        if(i<F.nodes.size()-1)
            out <<",";
    }
    out << "} i=" << F.inner << " if=" << F.inner_face
        << " o=" << F.outer << " of=" << F.outer_face
        << " b=" << (F.boundary?"true":"false")
        << " }";

    return out;
}


YAFEL_NAMESPACE_CLOSE
