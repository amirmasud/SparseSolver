//
// Created by Amir Masud on 1/11/19.
//

#ifndef SPARSE_SOLVER_SUPERNODAL_TSOLVER_H
#define SPARSE_SOLVER_SUPERNODAL_TSOLVER_H


#include "TriangularSolver.h"

class SupernodalTSolver : TriangularSolver{
public:
    static SupernodalTSolver *create(const utils::Matrix *L, const utils::Matrix *b);
    void solve() override;
    bool evaluate() override;

private:
    SupernodalTSolver(const utils::Matrix *L, const utils::Matrix *b);
    void init();

    void solveWithCSCFormat();
};


#endif //SPARSE_SOLVER_SUPERNODAL_TSOLVER_H
