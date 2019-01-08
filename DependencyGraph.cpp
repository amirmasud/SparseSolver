//
// Created by Amir Masud on 1/8/19.
//

#include <algorithm>
#include "DependencyGraph.h"

DependencyGraph* DependencyGraph::createWithCSCMatrix(
        const utils::CSCMatrix *cscMatrix) {
    auto dependencyGraph = new DependencyGraph();
    dependencyGraph->initWithCSCMatrix(cscMatrix);
    return dependencyGraph;
}
void DependencyGraph::initWithCSCMatrix(const utils::CSCMatrix *cscMatrix) {
    graph = new std::vector<size_t>[cscMatrix->rowCount];
    seen = new bool[cscMatrix->rowCount];
    reached = new size_t[cscMatrix->rowCount];
    for (int i = 0; i < cscMatrix->colCount - 1; ++i) {
        size_t dif = cscMatrix->Lp[i + 1] - cscMatrix->Lp[i] - 1;
        graph[i].resize(dif);
        int adjNumber = 0;
        for (size_t j = cscMatrix->Lp[i] + 1; j < cscMatrix->Lp[i + 1]; ++j) {
            graph[i][adjNumber] = cscMatrix->Li[j];
            ++adjNumber;
        }
    }
}

std::pair<size_t *, size_t> DependencyGraph::calculateReachSetFromCSCVector(
        const utils::CSCVector *cscVector) {
    std::fill_n(seen, cscVector->rowCount, false);
    size_t count = 0;
    for (int i = 0; i < cscVector->nnz; ++i)
        if (!seen[cscVector->Li[i]])
            dfs(cscVector->Li[i], count);
    std::sort(reached, reached + count);
    return {reached, count};
}
void DependencyGraph::dfs(size_t nodeIndex, size_t &count) {
    reached[count] = nodeIndex;
    ++count;
    seen[nodeIndex] = true;
    for (size_t i : graph[nodeIndex])
        if (!seen[i])
            dfs(i, count);
}
