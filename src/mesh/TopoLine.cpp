#include "mesh/TopoLine.hpp"
#include "mesh/TopoPoint.hpp"
#include "mesh/TopoFace.hpp"

YAFEL_NAMESPACE_OPEN

TopoLine::TopoLine(): TopoLine(0, nullptr, nullptr) {}

TopoLine::TopoLine(std::size_t _id) : TopoLine(_id, nullptr, nullptr) {}

TopoLine::TopoLine(std::size_t _id, TopoPoint *h, TopoPoint *t) :
  id(_id), head(h), tail(t), left(nullptr), right(nullptr)
{}

TopoLine::TopoLine(const TopoLine &TL) {
  id = TL.id;
  head = TL.head;
  tail = TL.tail;
  right = TL.right;
  left = TL.left;
}


YAFEL_NAMESPACE_CLOSE

std::ostream & operator<<(std::ostream &out, const yafel::TopoLine& TL) {
  out << "Line: ID=" << TL.id
      << " H=" << TL.head->id
      << " T=" << TL.tail->id
      << " L=" << (TL.left==nullptr ? -1 : (int)(TL.left->id))
      << " R=" << (TL.right==nullptr ? -1 : (int)(TL.right->id));
  return out;
}
