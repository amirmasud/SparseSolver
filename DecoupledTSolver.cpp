//
// Created by Amir Masud on 1/8/19.
//

#include <iostream>
#include "DecoupledTSolver.h"

//TODO : test on graph with some unreached nodes

DecoupledTSolver::DecoupledTSolver(const utils::Matrix *L, const utils::Matrix *b) :
        TriangularSolver(L, b) {
}
DecoupledTSolver* DecoupledTSolver::create(const utils::Matrix *L,
                                           const utils::Matrix *b) {
    auto solver = new DecoupledTSolver(L, b);
    solver->init();
    return solver;
}
void DecoupledTSolver::init() {
    m_dependencyGraph = DependencyGraph::createWithCSCMatrix(
            dynamic_cast<const utils::CSCMatrix *>(m_L));
}

void DecoupledTSolver::solve() {
    solveWithCSCFormat();
}

void DecoupledTSolver::solveWithCSCFormat() {
    auto cscL = dynamic_cast<const utils::CSCMatrix *>(m_L);
    auto cscB = dynamic_cast<const utils::CSCVector *>(m_b);
    m_result = cscB->clone();
    auto reachSetPair = m_dependencyGraph->calculateReachSetFromCSCVector(cscB);
    auto reachSet = reachSetPair.first;
    auto reachSetSize = reachSetPair.second;

    auto Li = cscL->Li;
    auto Lp = cscL->Lp;
    auto Lx = cscL->Lx;
    auto x = dynamic_cast<utils::CSCVector *>(m_result)->v;

    for (size_t px = 0; px < reachSetSize; ++px) {
        size_t j = reachSet[px];
        x[j] /= Lx[Lp[j]];
        for (size_t p = Lp[j] + 1; p < Lp[j + 1]; ++p)
            x[Li[p]] -= Lx[p] * x[j];
    }
}

bool DecoupledTSolver::evaluate() {
    return evaluateWithCSCFormat();
}
