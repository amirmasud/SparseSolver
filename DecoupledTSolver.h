//
// Created by Amir Masud on 1/8/19.
//

#ifndef SPARSE_SOLVER_DECOUPLED_TSOLVER_H
#define SPARSE_SOLVER_DECOUPLED_TSOLVER_H


#include "TriangularSolver.h"
#include "DependencyGraph.h"

class DecoupledTSolver : TriangularSolver {
public:
    static DecoupledTSolver *create(const utils::Matrix *L, const utils::Matrix *b);
    void solve() override;
    bool evaluate() override;

private:
    DependencyGraph *m_dependencyGraph;
    DecoupledTSolver(const utils::Matrix *L, const utils::Matrix *b);
    void init();

    void solveWithCSCFormat();
};


#endif //SPARSE_SOLVER_DECOUPLED_TSOLVER_H
