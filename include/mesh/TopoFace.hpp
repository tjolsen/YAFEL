#ifndef _YAFEL_TOPOFACE_HPP
#define _YAFEL_TOPOFACE_HPP

#include "yafel_globals.hpp"
#include <vector>

YAFEL_NAMESPACE_OPEN

//forward declaration
class TopoLine;
class TopoPoint;

class TopoFace {
  
public:
  std::size_t id;
  std::vector<TopoPoint*> vertices;
  std::vector<TopoLine*> boundary;

  TopoFace();
  TopoFace(std::size_t _id);
  
};

YAFEL_NAMESPACE_CLOSE

#endif
