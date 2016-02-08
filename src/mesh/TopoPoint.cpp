#include "mesh/TopoPoint.hpp"
#include "mesh/TopoLine.hpp"

YAFEL_NAMESPACE_OPEN

TopoPoint::TopoPoint() : TopoPoint(0) {}

TopoPoint::TopoPoint(std::size_t _id):
  id(_id), incoming(), outgoing() 
{
  incoming.clear();
  outgoing.clear();
}

YAFEL_NAMESPACE_CLOSE

std::ostream & operator<<(std::ostream &out, const yafel::TopoPoint &TP) {
  
  out << "Point: ID=" << TP.id;
  out << " IN={ ";
  for(auto l=TP.incoming.begin(); l!=TP.incoming.end(); ++l) {
    out << (*l)->id << " ";
  }
  out << "} OUT={ ";
  for(auto l=TP.outgoing.begin(); l!=TP.outgoing.end(); ++l) {
    out << (*l)->id << " ";
  }
  out << "}";
  return out;
}
