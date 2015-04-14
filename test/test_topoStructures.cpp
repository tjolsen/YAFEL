#include "yafel.hpp"
#include <iostream>
#include <vector>

using namespace yafel;

int main() {
  
  std::vector<TopoPoint> points;
  std::vector<TopoLine> lines;
  
  unsigned Npts = 4;
  
  for(unsigned i=0; i<Npts; ++i) {
    TopoPoint tp(i);
    points.push_back(tp);
  }
  
  for(unsigned i=0; i<Npts; ++i) {
    TopoLine tl;
    tl.id = i;
    tl.tail = &points[i];
    tl.head = &points[(i+1)%Npts];
    tl.left = nullptr;
    tl.right = nullptr;

    lines.push_back(tl);
    
    lines[i].head->incoming.push_back(&lines[i]);
    lines[i].tail->outgoing.push_back(&lines[i]);
  }


  for(auto p = points.begin(); p != points.end(); ++p) {
    std::cout << (*p) << std::endl;
  }

  for(auto l = lines.begin(); l != lines.end(); ++l) {
    std::cout << (*l) << std::endl;
  }

  return 0;
}
