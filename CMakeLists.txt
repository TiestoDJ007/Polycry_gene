cmake_minimum_required(VERSION 3.14)
project(Polycry_gene)

set(CMAKE_CXX_STANDARD 17)
include_directories(/home/kakaka/local_env/include)

add_executable(Atoms_Initialization Atoms_Initialization.cpp)
add_executable(Choosing_Atom Choosing_Atom.cpp)
add_executable(generator generator.cpp)
add_executable(test test.cpp)
add_executable(para_test para_test.cpp)
add_executable(generator_algo generator_algo.h generator_algo.cpp test_poly_vertex.cpp)