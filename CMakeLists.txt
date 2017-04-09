cmake_minimum_required(VERSION 3.5)
project(YAFEL)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -ffast-math -march=native")

set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Og")

set(SOURCE_FILES
        include/element/Element.hpp
        include/element/ElementFactory.hpp
        include/element/ElementType.hpp
        include/element/ShapeFunctionUtils.hpp
        include/lin_alg/solver/direct/LUDecomposition.hpp
        include/lin_alg/solver/iterative/bicgstab_solve.hpp
        include/lin_alg/solver/iterative/cg_solve.hpp
        include/lin_alg/solver/iterative/IdentityPreconditioner.hpp
        include/lin_alg/solver/iterative/ILUPreconditioner.hpp
        include/lin_alg/solver/iterative/JacobiPreconditioner.hpp
        include/lin_alg/solver/iterative/pcg_solve.hpp
        include/lin_alg/solver/iterative/Preconditioner.hpp
        include/lin_alg/solver/solvers.hpp
        include/lin_alg/tensor/const_reference_index_iterator.hpp
        include/lin_alg/tensor/const_reference_tensor_iterator.hpp
        include/lin_alg/tensor/generic_index_iterator.hpp
        include/lin_alg/tensor/generic_tensor_iterator.hpp
        include/lin_alg/tensor/iterator_utils.hpp
        include/lin_alg/tensor/reference_index_iterator.hpp
        include/lin_alg/tensor/reference_tensor_iterator.hpp
        include/lin_alg/tensor/Tensor.hpp
        include/lin_alg/tensor/tensor_specializations.hpp
        include/lin_alg/tensor/TensorExpression.hpp
        include/lin_alg/access_sparse_matrix.hpp
        include/lin_alg/bcsr_spmv.hpp
        include/lin_alg/cg_solver.hpp
        include/lin_alg/construction_sparse_matrix.hpp
        include/lin_alg/csr_sparsity_pattern.hpp
        include/lin_alg/csr_spmv.hpp
        include/lin_alg/matmul.hpp
        include/lin_alg/Matrix.hpp
        include/lin_alg/MatrixExpression.hpp
        include/lin_alg/mvmul.hpp
        include/lin_alg/operators.hpp
        include/lin_alg/sparse_bcsr.hpp
        include/lin_alg/sparse_coo.hpp
        include/lin_alg/sparse_csr.hpp
        include/lin_alg/sparse_matrix.hpp
        include/lin_alg/sparse_utils.hpp
        include/lin_alg/sparsity_pattern.hpp
        include/lin_alg/Vector.hpp
        include/lin_alg/VectorExpression.hpp
        #include/mesh/Face.hpp
        #include/mesh/GenericMesh.hpp
        #include/mesh/GmshMesh.hpp
        #include/mesh/MeshIterator.hpp
        #include/mesh/MeshTopology.hpp
        #include/mesh/RectilinearMesh.hpp
        #include/mesh/TopoFace.hpp
        #include/mesh/TopoLine.hpp
        #include/mesh/TopoPoint.hpp

        include/new_mesh/ElementIterator.hpp
        include/quadrature/LobattoPoints1D.hpp
        include/new_mesh/Mesh.hpp
        include/new_mesh/mesh_typedefs.hpp

        include/output/MatrixVisualization.hpp
        include/output/VTKDGMesh.hpp
        include/output/VTKMesh.hpp
        include/output/VTKObject.hpp
        include/output/VTKOutput.hpp
        include/output/VTKScalarData.hpp
        include/output/VTKTensorData.hpp
        include/output/VTKTimeOutput.hpp
        include/output/VTKVectorData.hpp
        include/utils/CG_DoFManager.hpp
        include/utils/ConcurrentQueue.hpp
        include/utils/DG_DoFManager.hpp
        include/utils/DirBC.hpp
        include/utils/DoFManager.hpp
        include/utils/DualNumber.hpp
        include/utils/ElementType.hpp
        include/utils/ElementVtkType.hpp
        include/utils/GaussLegendreQuadrature.hpp
        include/utils/GaussLobattoQuadrature.hpp
        include/quadrature/QuadratureRule.hpp
        include/utils/Range.hpp
        include/utils/ScalarTraits.hpp
        include/utils/SpatialFunction.hpp
        include/utils/TaskScheduler.hpp
        include/utils/ThreadPool.hpp
        include/yafel.hpp
        include/yafel_globals.hpp
        include/yafel_typedefs.hpp



        src/element/Element.cpp
        src/element/ElementFactory.cpp
        src/element/ShapeFunctionUtils.cpp
        src/quadrature/LobattoPoints1D.cpp
        src/new_mesh/Mesh.cpp
        src/new_mesh/build_faces.cpp
        src/new_mesh/parse_gmsh.cpp
        #src/output/VTKObject.cpp
        #src/output/VTKOutput.cpp
        #src/output/VTKScalarData.cpp
        #src/output/VTKTimeOutput.cpp
        src/quadrature/QuadratureRule.cpp
        src/utils/DoFManager.cpp
        #src/utils/DirBC.cpp
        #src/utils/ElementType.cpp
        #src/utils/ElementVtkType.cpp
        )

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(yafel ${SOURCE_FILES})

add_executable(sandbox experimental/sandbox.cpp)
target_link_libraries(sandbox yafel)
