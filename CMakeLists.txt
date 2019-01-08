cmake_minimum_required(VERSION 3.13)
project(SparseSolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(SparseSolver main.cpp Utils.cpp Utils.h SparseSolver.cpp SparseSolver.h)