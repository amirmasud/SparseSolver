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
    auto Lcomp = cscL->Lcomp;
    auto sup2col = cscL->sup2col;
    auto supNo = cscL->supNo;

    assert(Lp && Li && x);
    size_t currCol, nxtCol, supWidth, nSupR;
    double *Ltrng;
    auto XContinousMem = new double[cscL->rowCount](); // magic number ! (can be reduced to max number of rows containng nonzero in each column)

    auto start = std::chrono::system_clock::now();

    for (size_t supIdx = 1; supIdx < supNo; ++supIdx) {

        currCol = supIdx==0? 0: sup2col[supIdx-1];
        nxtCol = sup2col[supIdx];
        supWidth = nxtCol-currCol;
        nSupR = Lp[currCol+1]-Lp[currCol]; // num of nonzero rows in the current column = num of nnz on that column
        Ltrng = &Lx[Lp[currCol]];//first nnz of current supernode
        switch (supWidth){
/*            case 1:
//                std::cout<<"salam"<<supWidth<<std::endl;
//                std::cout<<"salam"<<nSupR<<std::endl;
                if (x[currCol]!= 0) {
                    x[currCol] /= Lx[Lp[currCol]];
                    for (size_t p = Lp[currCol] + 1; p < Lp[currCol + 1]; ++p)
                        x[Li[p]] -= Lx[p] * x[currCol];
                }
                break;*/
            default:
                // solving diagonal values of the supernode
                matrix_hf::dlsolve_blas_nonUnit(nSupR,supWidth, Ltrng, &x[currCol]);

                Ltrng = &Lx[Lp[currCol]+supWidth];//first nnz of below diagonal

                size_t compNo = Lcomp[currCol].size();
                auto comps = &Lcomp[currCol]; // all component of current column
                // create a dense temp vector that includes corresponding indices of all components.
                makeContiniousX(x, XContinousMem, compNo, comps);
                solveRect(Lp, Lx, currCol, supWidth, nSupR, XContinousMem);
                copyBackContiniousX(x, XContinousMem, compNo, comps);



                

//                matrix_hf::dmatvec_blas(nSupR,nSupR-supWidth,supWidth,Ltrng,&x[currCol],XContinousMem);
//                if (supIdx ==1){
//                    for (int i = 840; i < 860; ++i) {
//                        std::cout<<i<<":"<<x[i]<<std::endl;
//                    }
//                }
//                for (size_t idx = Lp[currCol]+supWidth,k=0; idx < Lp[nxtCol]; ++idx, ++k) {
//                    x[Li[idx]]-=XContinousMem[k];
//                    XContinousMem[k]=0;
//                }
        }
    }
    delete XContinousMem;


    auto end = std::chrono::system_clock::now();
    std::cout <<"Supernodal solve time:    "<< (end - start).count() << std::endl;
}

void SupernodalTSolver::solveRect(const size_t *Lp, const double *Lx, size_t currCol, size_t supWidth, size_t nSupR,
                                   double *XContinousMem) const {
    for (size_t colIdx = 0; colIdx < supWidth; ++colIdx){
                    size_t triangleColumnHead = Lp[currCol + colIdx];
                    size_t rectColumnHead = triangleColumnHead + supWidth - colIdx;
                    size_t rectHeight = nSupR - supWidth;
                    for (size_t offset = 0; offset < rectHeight; ++offset) {
                        XContinousMem[offset + supWidth] -= Lx[rectColumnHead + offset] * XContinousMem[colIdx];
                    }
                }
}

void SupernodalTSolver::makeContiniousX(double *x, double *XContinousMem, size_t compNo,
                                                 const vector<pair<size_t, size_t>> *comps) {
    size_t tempVecSize = 0;
    for (size_t compIdx = 0; compIdx < compNo; ++compIdx) {
                    size_t startRow = (*comps)[compIdx].first;
                    size_t endRow = (*comps)[compIdx].second;
                    // component includes both first and last elements
                    std::copy(&x[startRow], &x[endRow + 1], &XContinousMem[tempVecSize]);
                    tempVecSize += endRow+1 - startRow;
                }
}

void SupernodalTSolver::copyBackContiniousX(double *x, double *XContinousMem, size_t compNo,
                                             const vector<pair<size_t, size_t>> *comps) {
    size_t tempVecSize = 0;
    for (size_t compIdx = 0; compIdx < compNo; ++compIdx) {
        size_t startRow = (*comps)[compIdx].first;
        size_t endRow = (*comps)[compIdx].second;
        size_t compHeight = endRow + 1 - startRow;
        // component includes both first and last elements
        std::copy(&XContinousMem[tempVecSize], &XContinousMem[tempVecSize+ compHeight], &x[startRow]);
        std::copy(&x[startRow], &x[endRow + 1], &XContinousMem[tempVecSize]);
        tempVecSize += compHeight;
    }

}

bool SupernodalTSolver::evaluate() {
    return evaluateWithCSCFormat();
}