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
#include "element/Element.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"
#include "element/LinTet.hpp"
#include "element/ElementFactory.hpp"

// Output
#include "output/VTKObject.hpp"
#include "output/VTKOutput.hpp"
#include "output/VTKScalarData.hpp"
#include "output/VTKVectorData.hpp"
#include "output/VTKTensorData.hpp"
#include "output/VTKMesh.hpp"
#include "output/VTKTimeOutput.hpp"

#endif
