//
// Created by Amir Masud on 1/8/19.
//

#include "TriangularSolver.h"

TriangularSolver::TriangularSolver(const utils::Matrix *L,
                                   const utils::Matrix *b) :
                                   m_L(L), m_b(b) {
}
