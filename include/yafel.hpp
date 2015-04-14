#ifndef _YAFEL_HPP
#define _YAFEL_HPP

/*
  This file serves as the simplest interface to YAFEL by including everything.
  If desired, this file may be ignored in favor of more fine-grained control
  over the parts of the library that are included.
 */

// Library-wide global declarations
#include "yafel_globals.hpp"

// Utils
#include "utils/DirBC.hpp"
#include "utils/DoFManager.hpp"
#include "utils/SpatialFunction.hpp"

// Linear Algebra
#include "lin_alg/Vector.hpp"
#include "lin_alg/FullMatrix.hpp"
#include "lin_alg/sparse_coo.hpp"
#include "lin_alg/sparse_csr.hpp"
#include "lin_alg/cg_solver.hpp"
#include "lin_alg/bicgstab.hpp"
#include "lin_alg/LUDecomposition.hpp"
#include "lin_alg/Preconditioner.hpp"
#include "lin_alg/ILUPreconditioner.hpp"
#include "lin_alg/JacobiPreconditioner.hpp"

// Mesh
#include "mesh/Mesh.hpp"
#include "mesh/MeshReader.hpp"
#include "mesh/MeshGenerator.hpp"
//#include "mesh/DGMesh.hpp"
#include "mesh/MeshTopology.hpp"
#include "mesh/TopoPoint.hpp"
#include "mesh/TopoLine.hpp"
#include "mesh/TopoFace.hpp"
//#include "mesh/TopoCell.hpp"  // <--- planned addition

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
#include "output/MatrixVisualization.hpp"

#endif
