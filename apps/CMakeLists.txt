include(CMakeDependentOption)


cmake_dependent_option(APP_MESH_AREA "Compile MeshArea" OFF
        "COMPILE_APPS" OFF)
if(APP_MESH_AREA)
    add_executable(meshArea MeshArea/main.cpp)
    target_link_libraries(meshArea yafel)
endif()


cmake_dependent_option(APP_BASIC_POISSON "Compile BasicPoisson" OFF
        "COMPILE_APPS" OFF)
if(APP_BASIC_POISSON)
    add_executable(basicPoisson BasicPoisson/main.cpp)
    target_link_libraries(basicPoisson yafel)
endif()

cmake_dependent_option(APP_BASIC_ELASTICITY "Compile BasicElasticity" OFF
        "COMPILE_APPS" OFF)
if(APP_BASIC_ELASTICITY)
    add_executable(basicElasticity BasicLinearElasticity/main.cpp)
    target_link_libraries(basicElasticity yafel)
endif()

cmake_dependent_option(
        APP_BASIC_ADVECTION "Compile basicAdvection" OFF
        "COMPILE_APPS" OFF)
if(APP_BASIC_ADVECTION)
    add_executable(basicAdvection BasicAdvection_DG/main.cpp)
    target_link_libraries(basicAdvection yafel)
endif()

cmake_dependent_option(APP_HDG_POISSON "Compile HDG_Poisson" OFF
        "COMPILE_APPS" OFF)
if(APP_HDG_POISSON)
    add_executable(hdg_poisson HDG_Poisson/main.cpp)
    target_link_libraries(hdg_poisson yafel)
endif()

cmake_dependent_option(APP_SHELL_PROJECT "Compile Shell Project" OFF
        "COMPILE_APPS" OFF)
if(APP_SHELL_PROJECT)
    add_executable(shellProject 2081_CG_Shell/main.cpp)
    target_link_libraries(shellProject yafel)
endif()

cmake_dependent_option(APP_FINITE_VOLUME "Compile FiniteVolume" ON
        "COMPILE_APPS" OFF)
if(APP_FINITE_VOLUME)
    add_executable(finiteVolume
            FiniteVolume/main.cpp
            FiniteVolume/FVDofm.cpp
            FiniteVolume/FVDofm.hpp
            FiniteVolume/CellToNodeMap.hpp
            FiniteVolume/cell_centroids.hpp
            )
    target_link_libraries(finiteVolume yafel)
    target_compile_options(finiteVolume PRIVATE -fopenmp)
endif()

cmake_dependent_option(APP_SUPG "Compile SUPG" OFF
        "COMPILE_APPS" OFF)
if(APP_SUPG)
    add_executable(supg supg_learning/main.cpp supg_learning/supg_advection.hpp)
    target_link_libraries(supg yafel)
endif()