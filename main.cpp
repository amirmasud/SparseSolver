#include <iostream>

#include "Utils.h"
#include "NaiveTSolver.h"
#include "DecoupledTSolver.h"

int main() {
    const auto *L = utils::CSCMatrix::createFromFile("af_0_k101.mtx");
    const auto *b = utils::CSCVector::createFromFile("b.mtx");
    /*auto naiveTSolver = NaiveTSolver::create(L, b);
    naiveTSolver->solve();
    naiveTSolver->evaluate();*/
    auto decoupledTSolver = DecoupledTSolver::create(L, b);
    decoupledTSolver->solve();
    decoupledTSolver->evaluate();
    return 0;
}
