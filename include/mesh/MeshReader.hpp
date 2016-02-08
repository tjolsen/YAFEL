#ifndef __YAFEL_MESHREADER_HPP
#define __YAFEL_MESHREADER_HPP

#include "yafel_globals.hpp"
#include "mesh/Mesh.hpp"

#include <string>

YAFEL_NAMESPACE_OPEN

class MeshReader {

public:
  using size_type = typename Mesh::size_type;
  using coordinate_type = typename Mesh::coordinate_type;

  static Mesh gmsh_read(std::string fname);

private:
  static void string_split(std::string s, char delim, 
			  std::vector<std::string> & elems);

};

YAFEL_NAMESPACE_CLOSE

#endif
