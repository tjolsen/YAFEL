#ifndef _YAFEL_HPP
#define _YAFEL_HPP

/*
  This file serves as the simplest interface to YAFEL by including everything.
  If desired, this file may be ignored in favor of more fine-grained control
  over the parts of the library that are included.
 */

// Library-wide global declarations
#include "yafel_globals.hpp"

// Linear Algebra
#include "lin_alg/Vector.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/cg_solver.hpp"
#include "lin_alg/bicgstab.hpp"
#include "lin_alg/LUDecomposition.hpp"

// Mesh
#include "mesh/Mesh.hpp"
#include "mesh/MeshReader.hpp"

// Elements
#include "Element.hpp"
#include "LinQuad.hpp"
#include "LinTri.hpp"
#include "LinTet.hpp"
#include "ElementFactory.hpp"

// Output
#include "VTKObject.hpp"
#include "VTKOutput.hpp"
#include "VTKScalarData.hpp"
#include "VTKVectorData.hpp"
#include "VTKTensorData.hpp"
#include "VTKMesh.hpp"

#endif
