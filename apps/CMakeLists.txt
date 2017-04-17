add_executable(meshArea MeshArea/main.cpp)
target_link_libraries(meshArea yafel)

add_executable(basicPoisson BasicPoisson/main.cpp)
target_link_libraries(basicPoisson yafel)