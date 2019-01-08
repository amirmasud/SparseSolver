//
// Created by iman on 1/7/19.
//

#include "Utils.h"

#include <fstream>
#include <iostream>
#include <sstream>

namespace utils {

Matrix::Matrix(utils::Format format) : m_format(format) {
}
Format Matrix::getFormat() const {
    return m_format;
}

//template <class T>
CSCMatrix::CSCMatrix(size_t m_colCount, size_t m_nnz,
                     int *m_Lp, int *m_Li, double *m_Lx) :
        Matrix(Format::CSC), colCount(m_colCount), nnz(m_nnz),
        Lp(m_Lp), Li(m_Li), Lx(m_Lx) {
}
CSCMatrix::CSCMatrix() : Matrix(Format::CSC),
                         Lp(nullptr), Li(nullptr), Lx(nullptr) {
}
//template <class T>
CSCMatrix* CSCMatrix::create(size_t colCount, size_t nnz,
                             int *Lp, int *Li, double *Lx) {
    auto cscMatrix = new CSCMatrix(colCount, nnz, Lp, Li, Lx);
    cscMatrix->init();
    return cscMatrix;
}
//template <class T>
CSCMatrix* CSCMatrix::createFromFile(const std::string &name) {
    auto cscMatrix = new CSCMatrix();
    cscMatrix->initWithFile(name);
    return cscMatrix;
}
//template <class T>
void CSCMatrix::init() {
}
//template <class T>
void CSCMatrix::initWithFile(const std::string &name) {
    matrix_hf::readMatrix(name, rowCount, colCount, nnz, Lp, Li, Lx);
}

CSCMatrix::~CSCMatrix() {
    delete Lp;
    delete Li;
    delete Lx;
}

CSCVector::CSCVector() : CSCMatrix() {
}
CSCVector* CSCVector::createFromFile(const std::string &name) {
    auto cscVector = new CSCVector();
    cscVector->initWithFile(name);
    return cscVector;
}
void CSCVector::initWithFile(const std::string &name) {
    CSCMatrix::initWithFile(name);
    initVector();
}
void CSCVector::initVector() {
    v = new double[rowCount];
    std::fill_n(v, rowCount, 0);
    for (size_t i = 0; i < rowCount; ++i)
        v[Li[i]] = Lx[i];
}

CSCVector::~CSCVector() {
    delete v;
}

}


namespace matrix_hf {

bool readMatrix(const std::string &fName, size_t  &rowCount, size_t &colCount,
                size_t &NNZ, int* &col, int* &row, double* &val) {
    /*This function reads the input matrix from "fName" file and
     * allocate memory for matrix A, L and U.
     * - The input file is a coordinate version and e
     * ach row of the file shows (col, row, nnz)
     * - The matrices are zero-indexed
     */

    std::ifstream inFile;
    inFile.open(fName);
    std::string line,banner, mtx, crd, arith, sym;
    /*  File format:
     *    %%MatrixMarket matrix coordinate real general/symmetric/...
     *    % ...
     *    % (optional comments)
     *    % ...
     *    #rows    #non-zero
     *    Triplet in the rest of lines: row    col    value
     */
    std::getline(inFile,line);
    for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
    std::istringstream iss(line);
    if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
        std::cout<<"Invalid header (first line does not contain 5 tokens)\n";
        return false;
    }

    if(banner.compare("%%matrixmarket")) {
        std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")\n";
        return false;
    }
    if(mtx.compare("matrix")) {
        std::cout<<"Not a matrix; this driver cannot handle that.\"\n";
        return false;
    }
    if(crd.compare("coordinate")) {
        std::cout<<"Not in coordinate format; this driver cannot handle that.\"\n";
        return false;
    }
    if(arith.compare("real")) {
        if(!arith.compare("complex")) {
            std::cout<<"Complex matrix; use zreadMM instead!\n";
            return false;
        }
        else if(!arith.compare("pattern")) {
            std::cout<<"Pattern matrix; values are needed!\n";
            return false;
        }
        else {
            std::cout<<"Unknown arithmetic\n";
            return false;
        }
    }
    while (!line.compare(0,1,"%"))
    {
        std::getline(inFile, line);
    }
    std::istringstream issDim(line);
    if (!(issDim >> rowCount >> colCount >> NNZ)){
        std::cout<<"The matrix dimension is missing\n";
        return false;
    }
    if(colCount <= 0 || NNZ <= 0)
        return false;
    col = new int[colCount + 1]();
    // colL = new int[colCount + 1]; colU = new int[colCount + 1];
    row = new int[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    if(!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt=0, nnzCnt=0;
    double value;

    col[0]=0;
    for (int i = 0; nnzCnt<NNZ; ) {//Reading from file row by row
        inFile>>x;x--;
        inFile>>y;y--;//zero indexing
        inFile>>value;
        if(y > colCount)
            return false;
        if(y==i){
            val[nnzCnt]=value;
            row[nnzCnt]=x;
            colCnt++; nnzCnt++;
        }
        else{//New col
            col[i+1]=col[i]+colCnt;
            i++;//next iteration
            colCnt=1;
            val[nnzCnt]=value;
            row[nnzCnt]=x;
            nnzCnt++;
        }

    }
    col[colCount]= col[colCount - 1] + colCnt;//last col

    return true;
}

}