cmake_minimum_required(VERSION 3.14)
project(Polycry_gene)

set(CMAKE_CXX_STANDARD 17)
include_directories(/home/kakaka/local_env/include)

add_executable(Polycry_gene Read_Tess.cpp)
add_executable(Atoms_Initialization Atoms_Initialization.cpp)