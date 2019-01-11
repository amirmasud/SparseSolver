//
// Created by Amir Masud on 1/11/19.
//

#ifndef SPARSE_SOLVER_SUPERNODAL_TSOLVER_H
#define SPARSE_SOLVER_SUPERNODAL_TSOLVER_H


#include "TriangularSolver.h"

class SupernodalTSolver : TriangularSolver {
public:
    static SupernodalTSolver *create(const utils::Matrix *L, const utils::Matrix *b);

    void solve() override;

    bool evaluate() override;

private:
    SupernodalTSolver(const utils::Matrix *L, const utils::Matrix *b);

    void init();

    void solveWithCSCFormat();

    void makeContiniousX(double *x, double *XContinousMem, size_t compNo,
                         const vector<pair<size_t, size_t>> *comps);

    void copyBackContiniousX(double *x, double *XContinousMem, size_t compNo,
                             const vector<pair<size_t, size_t>> *comps);

    void solveRect(const size_t *Lp, const double *Lx, size_t currCol, size_t supWidth, size_t nSupR,
                    double *XContinousMem) const;
};


#endif //SPARSE_SOLVER_SUPERNODAL_TSOLVER_H
