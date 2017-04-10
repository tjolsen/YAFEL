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
#include "utils/CG_DoFManager.hpp"
#include "utils/DG_DoFManager.hpp"
#include "utils/DualNumber.hpp"
#include "utils/ElementType.hpp"
#include "quadrature/QuadratureRule.hpp"
#include "utils/GaussLegendreQuadrature.hpp"
#include "utils/GaussLobattoQuadrature.hpp"
#include "utils/SpatialFunction.hpp"


// Linear Algebra
#include "old_handmade_linalg/Vector.hpp"
#include "old_handmade_linalg/Matrix.hpp"
#include "old_handmade_linalg/sparse_coo.hpp"
#include "old_handmade_linalg/sparse_csr.hpp"
#include "old_handmade_linalg/sparse_bcsr.hpp"
#include "lin_alg/solver/iterative/cg_solve.hpp"
#include "lin_alg/solver/iterative/bicgstab_solve.hpp"
#include "lin_alg/solver/iterative/Preconditioner.hpp"
#include "lin_alg/solver/iterative/ILUPreconditioner.hpp"
#include "lin_alg/solver/iterative/JacobiPreconditioner.hpp"
#include "lin_alg/solver/direct/LUDecomposition.hpp"

// Mesh
#include "mesh/GenericMesh.hpp"
#include "mesh/RectilinearMesh.hpp"
#include "mesh/GmshMesh.hpp"
#include "mesh/Face.hpp"

// Elements
#include "element/Element.hpp"
#include "element/LinQuad.hpp"
#include "element/LinTri.hpp"
#include "element/LinHex.hpp"
//#include "element/LinTet.hpp"
#include "element/ElementFactory.hpp"
#include "element/DG_Quad.hpp"

// Output
#include "output/VTKObject.hpp"
#include "output/VTKOutput.hpp"
#include "output/VTKScalarData.hpp"
#include "output/VTKVectorData.hpp"
#include "output/VTKTensorData.hpp"
#include "output/VTKMesh.hpp"
#include "output/VTKDGMesh.hpp"
#include "output/VTKTimeOutput.hpp"
#include "output/MatrixVisualization.hpp"

#endif
