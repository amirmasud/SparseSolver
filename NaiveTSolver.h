//
// Created by Amir Masud on 1/8/19.
//

#ifndef SPARSE_SOLVER_NAIVE_TSOLVER_H
#define SPARSE_SOLVER_NAIVE_TSOLVER_H


#include "TriangularSolver.h"

class NaiveTSolver : public TriangularSolver {
public:
    static NaiveTSolver *create(const utils::Matrix *L, const utils::Matrix *b);
    void solve() override;
    bool evaluate() override;

private:
    NaiveTSolver(const utils::Matrix *L, const utils::Matrix *b);
    void init();

    void solveWithCSCFormat();
};


#endif //SPARSE_SOLVER_NAIVE_TSOLVER_H
