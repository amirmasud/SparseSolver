#include <iostream>

#include "Utils.h"
#include "NaiveTSolver.h"
#include "DecoupledTSolver.h"
#include "SupernodalTSolver.h"


int main(int argc, char *argv[]) {
    if(argc < 2){
        std::cout<<"The input matrices is missing!"<<std::endl;
        std::cout<<"First input is the path of  L matrix and the second is for b matrix!"<<std::endl;
        return -1;
    }
    std::string L_fileName = argv[1];
    std::string B_fileName = argv[2];

    const auto *L1 = utils::CSCMatrix::createFromFile(L_fileName);
    const auto *b = utils::CSCVector::createFromFile(B_fileName);
    auto naiveTSolver = NaiveTSolver::create(L1, b);
    naiveTSolver->solve();
    naiveTSolver->evaluate();
//    const auto *newL = utils::AugmentedCSCMatrix::createFromFile("af_0_k101.mtx");
//    newL->print();

    const auto *L2 = utils::SupernodalCSCMatrix::createFromFile(L_fileName);
    auto supernodalTSovler = SupernodalTSolver::create(L2,b);
    supernodalTSovler->solve();
    supernodalTSovler->evaluate();
//    auto decoupledTSolver = DecoupledTSolver::create(L, b);
//    decoupledTSolver->solve();
//    decoupledTSolver->evaluate();
    return 0;
}
