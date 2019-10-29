cmake_minimum_required(VERSION 3.14)
project(Polycry_gene)

set(CMAKE_CXX_STANDARD 17)
include_directories(/home/kakaka/local_env/include)

add_executable( Ray_generator File_Func.cpp File_Func.h Vector_Generator.cpp test.cpp)