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
  std::size_t id;
  TopoPoint *head;
  TopoPoint *tail;
  TopoFace *left;
  TopoFace *right;

  TopoLine();
  TopoLine(std::size_t _id);
  TopoLine(std::size_t _id, TopoPoint *h, TopoPoint *t);
  TopoLine(const TopoLine &TL);
};

YAFEL_NAMESPACE_CLOSE

std::ostream & operator<<(std::ostream &out, const yafel::TopoLine& TL);

#endif
