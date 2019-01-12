cmake_minimum_required(VERSION 3.13)
project(SparseSolver)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3")
#-ftree-vectorize -fopt-info-vec-missed
#-fprofile-generate
#-fprofile-use
#-ftree-vectorize -ftree-vectorizer-verbose=1
#set_target_properties(intrin-ompsimd PROPERTIES
#         -fopt-info-vec-missed"
#        LINKER_FLAGS")


add_executable(SparseSolver main.cpp Utils.cpp Utils.h SparseSolver.cpp SparseSolver.h TriangularSolver.cpp TriangularSolver.h NaiveTSolver.cpp NaiveTSolver.h DependencyGraph.cpp DependencyGraph.h DecoupledTSolver.cpp DecoupledTSolver.h SupernodalTSolver.cpp SupernodalTSolver.h)