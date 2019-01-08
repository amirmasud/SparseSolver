//
// Created by Amir Masud on 1/8/19.
//

#include "TriangularSolver.h"

#include <iostream>

TriangularSolver::TriangularSolver(const utils::Matrix *L,
                                   const utils::Matrix *b) :
                                   m_L(L), m_b(b) {
}

bool TriangularSolver::evaluateWithCSCFormat() {
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
