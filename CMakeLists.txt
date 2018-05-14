cmake_minimum_required(VERSION 3.8)
project(YAFEL)
#==========================================
set(YAFEL_VERSION_MAJOR 0)
set(YAFEL_VERSION_MINOR 1)
set(YAFEL_VERSION "${YAFEL_VERSION_MAJOR}.${YAFEL_VERSION_MINOR}")

set(CMAKE_CXX_STANDARD 17)


include(ProcessorCount)
ProcessorCount(NATIVE_CORE_COUNT)

option(ENABLE_HDF5 "Enable HDF5 Output" ON)
option(COMPILE_APPS "Compile application directory" ON)
set(VIENNACL_BACKEND OPENCL CACHE STRING "Valid ViennaCL backend options: OPENMP OPENCL")
set(CORE_COUNT ${NATIVE_CORE_COUNT} CACHE STRING "Number of cores to make YAFEL aware of")

#==========================================
# Get git revision version
execute_process(
        COMMAND "${CMAKE_SOURCE_DIR}/cmake/get_git_hash.sh"
        OUTPUT_VARIABLE YAFEL_GIT_REVISION
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
#strip trailing newline from echo command in get_git_hash.sh
string(REGEX REPLACE "\n$" "" YAFEL_GIT_REVISION "${YAFEL_GIT_REVISION}")

set(SOURCE_FILES
        include/assembly/AssemblyRequirement.hpp
        include/assembly/CGAssembly.hpp
        include/assembly/DGAssembly.hpp
        include/assembly/LocalSmoothingGradient.hpp
        include/assembly/ZZGradientRecovery.hpp

        include/boundary_conditions/DirichletBC.hpp

        include/element/Element.hpp
        include/element/ElementFactory.hpp
        include/element/ElementType.hpp
        include/element/ShapeFunctionUtils.hpp
        include/element/element_boundary_nodes.hpp

        include/fe_system/FESystem.hpp

        include/lin_alg/linear_solvers/LinearSolve.hpp
        include/lin_alg/linear_solvers/solvers/EigenConjugateGradient.hpp
        include/lin_alg/linear_solvers/solvers/EigenBICGSTAB.hpp
        include/lin_alg/linear_solvers/solvers/EigenCholesky.hpp
        include/lin_alg/linear_solvers/solvers/VCLConjugateGradient.hpp

        include/lin_alg/tensor/tensors.hpp
        include/lin_alg/tensor/TensorExpression.hpp
        include/lin_alg/tensor/mp_utils/contraction_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/map_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/mp_tuple_utils.hpp
        include/lin_alg/tensor/mp_utils/sequence_functions.hpp
        include/lin_alg/tensor/mp_utils/sequences.hpp
        include/lin_alg/tensor/mp_utils/slice_mp_utils.hpp
        include/lin_alg/tensor/mp_utils/TypeList.hpp
        include/lin_alg/tensor/tensor_expression_types/BinaryOperations.hpp
        include/lin_alg/tensor/tensor_expression_types/Tensor.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorConstant.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorContraction.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorCwiseBinaryOp.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorCwiseUnaryOp.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorDyadicProduct.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorFunctor.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorMap.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorPermutation.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorScaled.hpp
        include/lin_alg/tensor/tensor_expression_types/TensorSlice.hpp
        include/lin_alg/tensor/tensor_expression_types/UnaryOperations.hpp
        include/lin_alg/tensor/tensor_functions/BinaryFunctions.hpp
        include/lin_alg/tensor/tensor_functions/operators.hpp
        include/lin_alg/tensor/tensor_functions/rank2_specializations.hpp
        include/lin_alg/tensor/tensor_functions/rank1_specializations.hpp
        include/lin_alg/tensor/tensor_functions/Reductions.hpp
        include/lin_alg/tensor/tensor_functions/spectral_decomposition.hpp
        include/lin_alg/tensor/tensor_functions/tensor_dot.hpp
        include/lin_alg/tensor/tensor_functions/UnaryFunctions.hpp
        include/lin_alg/tensor/tensor_functions/update_assignment.hpp


        include/mesh/ElementIterator.hpp
        include/mesh/mesh_typedefs.hpp
        include/mesh/Mesh.hpp
        include/mesh/CellFace.hpp

        include/output/HDFBackend.hpp
        include/output/OutputBackend.hpp
        include/output/OutputData.hpp
        include/output/OutputFrame.hpp
        include/output/OutputMesh.hpp
        include/output/SimulationOutput.hpp
        include/output/VTUBackend.hpp

        include/PDE/PDEBase.hpp

        include/quadrature/LobattoPoints1D.hpp
        include/quadrature/QuadratureRule.hpp

        include/time_integration/DGRK4.hpp

        include/utils/BasicTimer.hpp
        include/utils/ConcurrentQueue.hpp
        include/utils/DoFManager.hpp
        include/utils/DualNumber.hpp
        include/utils/ElementVtkType.hpp
        include/utils/Printing.hpp
        include/utils/Range.hpp
        include/utils/ScalarTraits.hpp
        include/utils/SmallVector.hpp
        include/utils/SpatialFunction.hpp
        include/utils/parallel/TaskScheduler.hpp
        include/utils/parallel/Task.hpp
        include/utils/parallel/wait_all.hpp
        include/utils/parallel/wait_any.hpp
        include/utils/parallel/parfor.hpp
        include/utils/parallel/pardot.hpp
        include/utils/parallel/parassign.hpp
        include/utils/parallel/ReductionVariable.hpp
        include/utils/ThreadPool.hpp

        include/yafel.hpp
        include/yafel_globals.hpp
        include/yafel_typedefs.hpp
        include/yafel_config.hpp
        include/utils/parallel/yafel_parallel.hpp

        src/boundary_conditions/DirichletBC.cpp

        src/element/Element.cpp
        src/element/ElementFactory.cpp
        src/element/ShapeFunctionUtils.cpp
        src/element/element_boundary.cpp
        #src/element/build_tet_faces.cpp
        #src/element/build_topodim2_faces.cpp
        src/element/make_tensorproduct_element.cpp
        src/element/make_simplex_element.cpp

        src/mesh/CellFace.cpp
        src/mesh/Mesh.cpp
        src/mesh/build_faces.cpp
        src/mesh/parse_gmsh.cpp

        src/output/HDFBackend.cpp
        src/output/OutputBackend.cpp
        src/output/OutputData.cpp
        src/output/OutputMesh.cpp
        src/output/SimulationOutput.cpp
        src/output/VTUBackend.cpp

        src/quadrature/LobattoPoints1D.cpp
        src/quadrature/QuadratureRule.cpp
        src/quadrature/TriangleQuadratureRules.cpp
        src/quadrature/TetrahedronQuadratureRules.cpp

        src/time_integration/DGRK4.cpp

        src/utils/DoFManager.cpp
        src/utils/parallel/TaskScheduler.cpp
        )

#include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(yafel ${SOURCE_FILES})
target_include_directories(yafel PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/>
        $<INSTALL_INTERFACE:include/yafel>)

target_compile_options(yafel PUBLIC $<$<CONFIG:RELEASE>:-march=native; -ffast-math>)
target_compile_options(yafel PUBLIC $<$<CONFIG:DEBUG>:-fno-omit-frame-pointer>)

# Load OpenMP
find_package(OpenMP REQUIRED)
if(OPENMP_CXX_FOUND)
    target_compile_options(yafel PUBLIC ${OpenMP_CXX_FLAGS})
    target_link_libraries(yafel PUBLIC ${OpenMP_CXX_LIBS})
    target_link_libraries(yafel PUBLIC ${OpenMP_CXX_FLAGS})
endif()


target_link_libraries(yafel LINK_PUBLIC pthread)
set(CMAKE_VERBOSE_MAKEFILE 1)
set_target_properties(yafel PROPERTIES
        CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
#add_library(Eigen3 IMPORTED INTERFACE)
target_link_libraries(yafel PUBLIC Eigen3::Eigen)
target_compile_definitions(yafel PUBLIC "EIGEN_DEFAULT_INDEX_TYPE=int")
target_compile_definitions(yafel PUBLIC "VIENNACL_HAVE_EIGEN")

find_package(ViennaCL QUIET)
if(VIENNACL_FOUND)
    target_include_directories(yafel PUBLIC "${VIENNACL_INCLUDE_DIRS}")
    target_link_libraries(yafel PUBLIC "${VIENNACL_LIBRARIES}")
endif()


if (ENABLE_HDF5)
    find_package(HDF5 COMPONENTS C HL REQUIRED)
    add_library(HDF5 IMPORTED SHARED)
    #set_property(TARGET HDF5 PROPERTY
    #        INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
    target_include_directories(yafel PUBLIC "${HDF5_INCLUDE_DIRS}")
    target_compile_definitions(yafel PUBLIC USE_HDF5)
    target_link_libraries(yafel PUBLIC ${HDF5_LIBRARIES} ${HDF5_C_LIBRARY_hdf5_hl})
endif ()

if(VIENNACL_BACKEND STREQUAL "OPENCL")
    find_package(OpenCL)
    target_link_libraries(yafel INTERFACE OpenCL::OpenCL)
    target_compile_definitions(yafel INTERFACE VIENNACL_WITH_OPENCL)
elseif(VIENNACL_BACKEND STREQUAL "OPENMP")
    target_compile_definitions(yafel INTERFACE VIENNACL_WITH_OPENMP)
endif ()

configure_file(${CMAKE_SOURCE_DIR}/src/yafel_config.cpp.in src/yafel_config.cpp)
target_sources(yafel PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/src/yafel_config.cpp)


add_executable(sandbox experimental/sandbox.cpp)
target_link_libraries(sandbox PRIVATE yafel)


# Add any subdirectories to the project
if (COMPILE_APPS STREQUAL "ON")
    add_subdirectory(apps)
endif ()


install(TARGETS yafel EXPORT yafelConfig
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(DIRECTORY include/
        DESTINATION include/yafel
        FILES_MATCHING PATTERN "*.hpp")
install(EXPORT yafelConfig DESTINATION cmake)
