add_executable(meshArea MeshArea/main.cpp)
target_link_libraries(meshArea yafel)

add_executable(basicPoisson BasicPoisson/main.cpp)
target_link_libraries(basicPoisson yafel)

add_executable(basicElasticity BasicLinearElasticity/main.cpp)
target_link_libraries(basicElasticity yafel)

add_executable(basicAdvection BasicAdvection_DG/main.cpp)
target_link_libraries(basicAdvection yafel)

add_executable(hdg_poisson HDG_Poisson/main.cpp)
target_link_libraries(hdg_poisson yafel)

add_executable(shellProject 2081_CG_Shell/main.cpp)
target_link_libraries(shellProject yafel)

add_executable(finiteVolume FiniteVolume/main.cpp)
target_link_libraries(finiteVolume yafel)