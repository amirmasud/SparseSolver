//
// Created by Amir Masud on 1/8/19.
//

//TODO: error handling

#include <cassert>
#include <iostream>
#include <chrono>
#include "NaiveTSolver.h"

NaiveTSolver::NaiveTSolver(utils::Matrix *L, utils::Matrix *b) :
        TriangularSolver(L, b) {
}
NaiveTSolver* NaiveTSolver::create(utils::Matrix *L, utils::Matrix *b) {
    auto solver = new NaiveTSolver(L, b);
    solver->init();
    return solver;
}
void NaiveTSolver::init() {
}

void NaiveTSolver::solve() {
    solveWithCSCFormat();

}
void NaiveTSolver::solveWithCSCFormat() {
    auto cscL = dynamic_cast<const utils::CSCMatrix *>(m_L);
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
    std::cout << (end - start).count() << std::endl;
}

bool NaiveTSolver::evaluate() {
    return evaluateWithCSCFormat();
}
bool NaiveTSolver::evaluateWithCSCFormat() {
    auto cscL = dynamic_cast<const utils::CSCMatrix *>(m_L);
    auto cscB = dynamic_cast<const utils::CSCVector *>(m_b);

    auto result = dynamic_cast<utils::CSCVector *>(m_result);
    auto ans = (*cscL) * (*result);

    if (ans.size() != result->rowCount) {
        std::cout << "oh! different size!" << std::endl;
        std::cout << "ans : " << ans.size() << "   ";
        std::cout << "res : " << result->rowCount << std::endl;
        return false;
    }

    for (size_t i = 0; i < ans.size(); ++i) {
        if (std::abs(ans[i] - cscB->v[i]) > 1e-13) {
            std::cout << "oh! not the same! with index " << i << " ";
            std::cout << "with values : ans=" << ans[i] << "  result=" <<
            cscB->v[i];
            return false;
        }
    }
    std::cout << "evaluation passed" << std::endl;
    return true;
}
