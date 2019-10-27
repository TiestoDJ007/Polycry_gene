cmake_minimum_required(VERSION 3.14)
project(Polycry_gene)

set(CMAKE_CXX_STANDARD 17)
include_directories(/home/kakaka/local_env/include)

add_executable(generator generator.cpp)
add_executable( Ray_generator Ray_File_Func.cpp Ray_File_Func.h Ray_Generator.cpp)