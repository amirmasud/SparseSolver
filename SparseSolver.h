//
// Created by Amir Masud on 1/8/19.
//

#ifndef SPARSE_SOLVER_SPARSE_SOLVER_H
#define SPARSE_SOLVER_SPARSE_SOLVER_H


class SparseSolver {
public:
    virtual void solve() = 0;
    virtual bool evaluate() = 0;
};


#endif //SPARSE_SOLVER_SPARSE_SOLVER_H
