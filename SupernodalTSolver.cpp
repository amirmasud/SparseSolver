//
// Created by Amir Masud on 1/11/19.
//

#include <cassert>
#include <chrono>
#include <iostream>
#include "SupernodalTSolver.h"

SupernodalTSolver::SupernodalTSolver(const utils::Matrix *L, const utils::Matrix *b) :
        TriangularSolver(L, b) {
}
SupernodalTSolver* SupernodalTSolver::create(const utils::Matrix *L, const utils::Matrix *b) {
    auto solver = new SupernodalTSolver(L, b);
    solver->init();
    return solver;
}
void SupernodalTSolver::init() {
}

void SupernodalTSolver::solve() {
    solveWithCSCFormat();
}

void SupernodalTSolver::solveWithCSCFormat() {
    auto cscL = dynamic_cast<const utils::SupernodalCSCMatrix *>(m_L);
    auto cscB = dynamic_cast<const utils::CSCVector *>(m_b);

    m_result = cscB->clone();

    auto Li = cscL->Li;
    auto Lp = cscL->Lp;
    auto Lx = cscL->Lx;
    auto x = dynamic_cast<utils::CSCVector *>(m_result)->v;

    assert(Lp && Li && x);


    auto start = std::chrono::system_clock::now();
    for (int j = 0; j < cscL->rowCount; ++j) {
        if (x[j] != 0) {
            x[j] /= Lx[Lp[j]];
            for (size_t p = Lp[j] + 1; p < Lp[j + 1]; ++p)
                x[Li[p]] -= Lx[p] * x[j];
        }
    }
    auto end = std::chrono::system_clock::now();
    std::cout <<"Supernodal solve time:    "<< (end - start).count() << std::endl;
}

bool SupernodalTSolver::evaluate() {
    return evaluateWithCSCFormat();
}