# Add Eigen3 as a required library
find_package(Eigen3 REQUIRED)
add_library(Eigen3 INTERFACE IMPORTED)
set_property(TARGET Eigen3 PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Eigen3_INCLUDE_DIR})

# Find ViennaCL. Technically optional, since Eigen solvers may be used.
find_package(ViennaCL OPTIONAL)
if(VIENNACL_FOUND)
    add_library(ViennaCL INTERFACE IMPORTED)
    set_property(TARGET ViennaCL PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES ${VIENNACL_INCLUDE_DIRS})
endif()

# Adjust compiler flags for OpenMP
find_package(OpenMP REUQIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()