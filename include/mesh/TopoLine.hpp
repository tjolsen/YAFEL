#ifndef _YAFEL_TOPOLINE_HPP
#define _YAFEL_TOPOLINE_HPP

#include "yafel_globals.hpp"
#include <iostream>

YAFEL_NAMESPACE_OPEN

//necessary forward declarations
class TopoPoint;
class TopoFace;

class TopoLine {

public:
  unsigned id;
  TopoPoint *head;
  TopoPoint *tail;
  TopoFace *left;
  TopoFace *right;

  TopoLine();
  TopoLine(unsigned _id);
  TopoLine(unsigned _id, TopoPoint *h, TopoPoint *t);
  TopoLine(const TopoLine &TL);
};

YAFEL_NAMESPACE_CLOSE

std::ostream & operator<<(std::ostream &out, const yafel::TopoLine& TL);

#endif
