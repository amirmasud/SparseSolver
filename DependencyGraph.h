//
// Created by Amir Masud on 1/8/19.
//

#ifndef SPARSE_SOLVER_DEPENDENCY_GRAPH_H
#define SPARSE_SOLVER_DEPENDENCY_GRAPH_H


#include "Utils.h"

class DependencyGraph {
public:
    static DependencyGraph *createWithCSCMatrix(const utils::CSCMatrix *cscMatrix);
    /**
     * @note For performance considerations we used pre-allocated array
     * (reached) instead of dynamic allocation.
     * @return A pointer to reached set array and size of reached set.
     */
    std::pair<size_t *, size_t> calculateReachSetFromCSCVector(
            const utils::CSCVector *cscVector);

private:
    std::vector<size_t> *graph;

    explicit DependencyGraph() = default;
    void initWithCSCMatrix(const utils::CSCMatrix *cscMatrix);

    void dfs(size_t nodeIndex, size_t &count);

    bool *seen;
    size_t *reached;
};


#endif //SPARSE_SOLVER_DEPENDENCY_GRAPH_H
