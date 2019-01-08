#include <iostream>

#include "Utils.h"
#include "NaiveTSolver.h"

int main() {
    auto *L = utils::CSCMatrix::createFromFile("af_0_k101.mtx");
    auto *b = utils::CSCVector::createFromFile("b.mtx");
    auto naiveTSolver = NaiveTSolver::create(L, b);
    naiveTSolver->solve();
    naiveTSolver->evaluate();
    return 0;
}
