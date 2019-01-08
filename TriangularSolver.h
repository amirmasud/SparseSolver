//
// Created by Amir Masud on 1/8/19.
//

#ifndef SPARSE_SOLVER_TRIANGULAR_SOLVER_H
#define SPARSE_SOLVER_TRIANGULAR_SOLVER_H


#include "SparseSolver.h"
#include "Utils.h"

class TriangularSolver : public SparseSolver {
protected:
    TriangularSolver(const utils::Matrix *L, const utils::Matrix *b);

    const utils::Matrix *m_L, *m_b;
    utils::Matrix *m_result;
};


#endif //SPARSE_SOLVER_TRIANGULAR_SOLVER_H
