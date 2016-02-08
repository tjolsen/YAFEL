#ifndef _YAFEL_TOPOPOINT_HPP
#define _YAFEL_TOPOPOINT_HPP

#include "yafel_globals.hpp"
#include <vector>
#include <iostream>

YAFEL_NAMESPACE_OPEN

//forward declaration of TopoLine
class TopoLine;


//definition of TopoPoint
class TopoPoint {

public:
  std::size_t id;
  std::vector<TopoLine*> incoming;
  std::vector<TopoLine*> outgoing;

  TopoPoint();
  TopoPoint(unsigned _id);
};

YAFEL_NAMESPACE_CLOSE

std::ostream &operator<<(std::ostream &out, const yafel::TopoPoint &TP);


#endif
