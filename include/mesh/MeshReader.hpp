#ifndef _YAFEL_MESHREADER_HPP
#define _YAFEL_MESHREADER_HPP

#include <string>
#include "mesh/Mesh.hpp"
#include "yafel_globals.hpp"

YAFEL_NAMESPACE_OPEN

class MeshReader {
private:
  static void string_split(std::string s, char delim, 
			  std::vector<std::string> & elems);

public:
  static Mesh gmsh_read(std::string fname);

};

YAFEL_NAMESPACE_CLOSE

#endif // end #ifndef _YAFEL_MESHREADER_HPP
