#include <iostream>

#include "Utils.h"
#include "NaiveTSolver.h"
#include "DecoupledTSolver.h"
#include "SupernodalTSolver.h"


int main() {

//    const auto *L = utils::CSCMatrix::createFromFile("af_0_k101.mtx");
    const auto *b = utils::CSCVector::createFromFile("b_sparse_af_0_k101.mtx");
//    auto naiveTSolver = NaiveTSolver::create(L, b);
//    naiveTSolver->solve();
//    naiveTSolver->evaluate();
//    const auto *newL = utils::AugmentedCSCMatrix::createFromFile("af_0_k101.mtx");
//    newL->print();

    const auto *L = utils::SupernodalCSCMatrix::createFromFile("af_0_k101.mtx");
    auto supernodalTSovler = SupernodalTSolver::create(L,b);
    supernodalTSovler->solve();
    supernodalTSovler->evaluate();
//    auto decoupledTSolver = DecoupledTSolver::create(L, b);
//    decoupledTSolver->solve();
//    decoupledTSolver->evaluate();
    return 0;
}
