//
// Created by Amir Masud on 1/7/19.
//

#include "Utils.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <chrono>
#include <algorithm>
//TODO : think of better clone.



namespace utils {

    Matrix::Matrix(utils::Format format, size_t m_rowCount, size_t m_colCount) :
            m_format(format), rowCount(m_rowCount), colCount(m_colCount) {
    }

    Format Matrix::getFormat() const {
        return m_format;
    }

    Matrix *Matrix::clone() const {
        auto newMatrix = new Matrix(m_format, rowCount, colCount);
        newMatrix->cloneAllocations(this);
        return newMatrix;
    }

//template <class T>
    CSCMatrix::CSCMatrix(size_t m_rowCount, size_t m_colCount,
                         size_t m_nnz, size_t *m_Lp, size_t *m_Li, double *m_Lx) :
            Matrix(Format::CSC, m_rowCount, m_colCount),
            nnz(m_nnz), Lp(m_Lp), Li(m_Li), Lx(m_Lx) {
    }

    CSCMatrix::CSCMatrix() : Matrix(Format::CSC, 0, 0),
                             Lp(nullptr), Li(nullptr), Lx(nullptr) {
    }

    CSCMatrix::CSCMatrix(const utils::CSCMatrix &cscMatrix) : Matrix(cscMatrix) {
        Lp = new size_t[colCount];
        std::copy(cscMatrix.Lp, cscMatrix.Lp + colCount, Lp);
        Li = new size_t[nnz];
        std::copy(cscMatrix.Li, cscMatrix.Li + nnz, Li);
        Lx = new double[nnz];
        std::copy(cscMatrix.Lx, cscMatrix.Lx + nnz, Lx);
    }

//template <class T>
    CSCMatrix *CSCMatrix::create(size_t rowCount, size_t colCount, size_t nnz,
                                 size_t *Lp, size_t *Li, double *Lx) {
        auto cscMatrix = new CSCMatrix(rowCount, colCount, nnz, Lp, Li, Lx);
        cscMatrix->init();
        return cscMatrix;
    }

//template <class T>
    CSCMatrix *CSCMatrix::createFromFile(const std::string &name) {
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

    std::vector<double> CSCMatrix::operator*(const CSCVector &cscVector) const {
        auto x = cscVector.v;
        size_t n = rowCount;
        std::vector<double> ans(n, 0.0);
        assert(Lp && x);
        for (size_t j = 0; j < n; ++j)
            for (size_t p = Lp[j]; p < Lp[j + 1]; ++p)
                ans[Li[p]] += Lx[p] * x[j];
        return ans;
    }

    void CSCMatrix::cloneAllocations(const CSCMatrix *cscMatrix) {
        Matrix::cloneAllocations(cscMatrix);
        Lp = new size_t[cscMatrix->colCount];
        std::copy(cscMatrix->Lp, cscMatrix->Lp + colCount, Lp);
        Li = new size_t[cscMatrix->nnz];
        std::copy(cscMatrix->Li, cscMatrix->Li + nnz, Li);
        Lx = new double[cscMatrix->nnz];
        std::copy(cscMatrix->Lx, cscMatrix->Lx + nnz, Lx);
    }

    Matrix *CSCMatrix::clone() const {
        auto newCSCMatrix = CSCMatrix::create(rowCount, colCount, nnz,
                                              nullptr, nullptr, nullptr);
        newCSCMatrix->cloneAllocations(this);
        return newCSCMatrix;
    }

    CSCMatrix::~CSCMatrix() {
        delete Lp;
        delete Li;
        delete Lx;
    }

    CSCVector::CSCVector() : CSCMatrix() {
    }

    CSCVector::CSCVector(size_t rowCount, size_t colCount, size_t nnz, size_t *Lp,
                         size_t *Li, double *Lx, double *m_v) :
            CSCMatrix(rowCount, colCount, nnz, Lp, Li, Lx),
            v(m_v) {
    }

    CSCVector::CSCVector(const utils::CSCVector &cscVector) : CSCMatrix(cscVector) {
        v = new double[rowCount];
        std::copy(cscVector.v, cscVector.v + rowCount, v);
    }

    CSCVector *CSCVector::create(size_t rowCount, size_t colCount, size_t nnz,
                                 size_t *Lp, size_t *Li, double *Lx, double *m_v) {
        auto csvVector = new CSCVector(rowCount, colCount, nnz, Lp, Li, Lx, m_v);
        csvVector->init();
        return csvVector;
    }

    CSCVector *CSCVector::createFromFile(const std::string &name) {
        auto cscVector = new CSCVector();
        cscVector->initWithFile(name);
        return cscVector;
    }

    void CSCVector::init() {
        CSCMatrix::init();
    }

    void CSCVector::initWithFile(const std::string &name) {
        CSCMatrix::initWithFile(name);
        initVector();
    }

    void CSCVector::initVector() {
        v = new double[rowCount];
        std::fill_n(v, rowCount, 0);
        for (size_t i = 0; i < nnz; ++i)
            v[Li[i]] = Lx[i];
    }

    void CSCVector::cloneAllocations(const CSCVector *cscVector) {
        CSCMatrix::cloneAllocations(cscVector);
        v = new double[cscVector->rowCount];
        std::copy(cscVector->v, cscVector->v + rowCount, v);
    }

    Matrix *CSCVector::clone() const {
        auto newCSCVector = CSCVector::create(rowCount, colCount, nnz,
                                              nullptr, nullptr, nullptr, nullptr);
        newCSCVector->cloneAllocations(this);
        return newCSCVector;
    }

    void CSCVector::dump() {
        for (size_t i = 0; i < 100; ++i)
            std::cout << v[i] << std::endl;
    }

    CSCVector::~CSCVector() {
        delete v;
    }



    AugmentedCSCMatrix *AugmentedCSCMatrix::createFromFile(const std::string &name) {
        auto augmentedCSCMatrix = new AugmentedCSCMatrix();
        augmentedCSCMatrix->initWithFile(name);
        return augmentedCSCMatrix;
    }

    void AugmentedCSCMatrix::initWithFile(const std::string &name) {
        matrix_hf::readAugmentedMatrix(name, rowCount, colCount, nnz, Lp, Li, Lx, Lcomp);

    }

    void AugmentedCSCMatrix::print() const{
        for (int i = 0; i < 10; i++)
            for (int j = 0; j < Lcomp[i].size(); j++)
            std::cout<<i<<","<<j<<": "<<Lcomp[i][j].first<<" : "<<Lcomp[i][j].second<<std::endl;

    }


    SupernodalCSCMatrix* SupernodalCSCMatrix::createFromFile(const std::string &name) {
        auto supernodalCSCMatrix = new SupernodalCSCMatrix();
        supernodalCSCMatrix->initWithFile(name);
        return supernodalCSCMatrix;
    }

    void SupernodalCSCMatrix::initWithFile(const std::string &name) {
        AugmentedCSCMatrix::initWithFile(name);
        auto start = std::chrono::system_clock::now();
        // super node detection
//        initSupernodes();
        initSupernodesFaster();
        auto end = std::chrono::system_clock::now();
        std::cout<<"Supernode detection time: "<< (end - start).count() <<std::endl;
    }


    void SupernodalCSCMatrix::initSupernodes(){
        supNo = 0;
        sup2col.resize(colCount); // initially all column can be a separate node;
        sup2col[0] = 0;
        bool similar_col = true;
        size_t newColComponentsNo;
        for (size_t col = 1; col < colCount; ++col) {
            similar_col = true;
            newColComponentsNo = Lcomp[col].size();
            if (newColComponentsNo == Lcomp[col-1].size()){
                for (size_t comp = 1; comp < newColComponentsNo; ++comp) {
                    // in order to be a similar column, all components except first one should be equal
                    if (Lcomp[col][comp] !=  Lcomp[col-1][comp]){
                        similar_col = false;
                        break;
                    }
                }
                // also last row of first components should be equal
                if (Lcomp[col][0].second != Lcomp[col-1][0].second)
                    similar_col = false;
                if (!similar_col) {
                    supNo++;
                    sup2col[supNo] = col;
                }
            }
        }
        sup2col.resize(supNo); // free extra memory
        std::cout<<"num of supernodes: "<<supNo<<std::endl;
    }

    void SupernodalCSCMatrix::initSupernodesFaster(){
        supNo = 0;
        sup2col.resize(colCount); // initially all column can be a separate node;
        sup2col[0] = 0;
        bool similar_col;
        size_t newColComponentsNo;
        for (size_t col = 1; col < colCount; ++col) {
            similar_col = true;
            newColComponentsNo = Lcomp[col].size();

            if (newColComponentsNo == Lcomp[col-1].size()){

                // all components should be the same except start of the first component
                if (!std::equal(&(Lcomp[col-1][0])+1,&(Lcomp[col-1][newColComponentsNo]), &(Lcomp[col][0])+1))
                    similar_col = false;

                // also last row of first components should be equal
                if (Lcomp[col][0].second != Lcomp[col-1][0].second)
                    similar_col = false;
                if (!similar_col) {
                    supNo++;
                    sup2col[supNo] = col;
                }
            }
        }
        sup2col.resize(supNo); // free extra memory
        std::cout<<"num of supernodes: "<<supNo<<std::endl;
    }
}

namespace matrix_hf {


    bool readAugmentedMatrix(const std::string &fName, size_t &rowCount, size_t &colCount,
                             size_t &NNZ, size_t *&col, size_t *&row, double *&val,
                             vector<pair<size_t, size_t>> *&Lcomp) {
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

    col = new size_t[colCount + 1]();
    // colL = new int[colCount + 1]; colU = new int[colCount + 1];
    row = new size_t[NNZ];
    // rowL = new int[factorSize]; rowU = new int[factorSize];
    val = new double[NNZ];
    // valL = new double[factorSize]; valU = new double[factorSize];
    Lcomp = new vector<pair<size_t, size_t>>[NNZ];

    if(!val || !col || !row)
        return false;
    //Initializing the result vector
    int y, x, colCnt=0, nnzCnt=0;
    double value;

    Lcomp[0].push_back({0,0}); // start of the first component of the first column
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
            if (nnzCnt != 0 && row[nnzCnt-1] + 1 != row[nnzCnt]){ // for each space between the none zero positions

                    Lcomp[i].back().second = row[nnzCnt - 1]; // the end of the component
                    Lcomp[i].push_back({row[nnzCnt], row[nnzCnt]}); // start of a new component

            }
            colCnt++; nnzCnt++;
        }
        else{//New col
            col[i+1]=col[i]+colCnt;
            i++;//next iteration
            colCnt=1;
            val[nnzCnt]=value;
            row[nnzCnt]=x;
            Lcomp[i-1].back().second = row[nnzCnt - 1];  // end of the last component of the previous column
            Lcomp[i].push_back({row[nnzCnt], row[nnzCnt]});  // start of first component of the column is in this row
            nnzCnt++;
        }

    }
    Lcomp[colCount-1][0].second = Lcomp[colCount-1][0].first; // last column has just one value on its diagonal
    col[colCount]= col[colCount - 1] + colCnt;//last col

    return true;
    }

    bool readMatrix(const std::string &fName, size_t  &rowCount, size_t &colCount,
                size_t &NNZ, size_t * &col, size_t * &row, double* &val) {
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
    col = new size_t[colCount + 1]();
    // colL = new int[colCount + 1]; colU = new int[colCount + 1];
    row = new size_t[NNZ];
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
    // TODO: use revese memory layout
    void dlsolve_blas_nonUnit ( int ldm, int ncol, double *M, double *rhs ) {//general triangular solver
        int k;
        double x0, x1, x2, x3, x4, x5, x6, x7;
        double *M0;
        register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
        register int firstcol = 0;

        M0 = &M[0];

        while (firstcol < ncol - 7) { /* Do 8 columns */
            Mki0 = M0;
            Mki1 = Mki0 + ldm;
            Mki2 = Mki1 + ldm - 1;
            Mki3 = Mki2 + ldm - 2;
            Mki4 = Mki3 + ldm - 3;
            Mki5 = Mki4 + ldm - 4;
            Mki6 = Mki5 + ldm - 5;
            Mki7 = Mki6 + ldm - 6;

            x0 = rhs[firstcol] / *Mki0++;
            x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;
            x2 = (rhs[firstcol + 2] - x0 * *Mki0++ - x1 * *Mki1++) / *Mki2++;
            x3 = (rhs[firstcol + 3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++) / *Mki3++;
            x4 = (rhs[firstcol + 4] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
                  - x3 * *Mki3++) / *Mki4++;
            x5 = (rhs[firstcol + 5] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
                  - x3 * *Mki3++ - x4 * *Mki4++) / *Mki5++;
            x6 = (rhs[firstcol + 6] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
                  - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++) / *Mki6++;
            x7 = (rhs[firstcol + 7] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++
                  - x3 * *Mki3++ - x4 * *Mki4++ - x5 * *Mki5++
                  - x6 * *Mki6++) / *Mki7++;

            rhs[firstcol++] = x0;
            rhs[firstcol++] = x1;
            rhs[firstcol++] = x2;
            rhs[firstcol++] = x3;
            rhs[firstcol++] = x4;
            rhs[firstcol++] = x5;
            rhs[firstcol++] = x6;
            rhs[firstcol++] = x7;

            for (k = firstcol; k < ncol; k++)
                rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                         - x2 * *Mki2++ - x3 * *Mki3++
                         - x4 * *Mki4++ - x5 * *Mki5++
                         - x6 * *Mki6++ - x7 * *Mki7++;
            M0 += 8 * ldm - 28; // -1 -2 -3 -4 -5 -6 -7
            ldm = ldm - 8;

        }

        while (firstcol < ncol - 3) { /* Do 4 columns */
//            for (int shit = 0; shit< 10;shit++){
//                std::cout<<"shit: "<<*(M0+shit)<<std::endl;
//            }
            Mki0 = M0;
            Mki1 = Mki0 + ldm;
            Mki2 = Mki1 + ldm - 1;
            Mki3 = Mki2 + ldm - 2;
            Mki4 = Mki3 + ldm - 3; //next M0
            x0 = rhs[firstcol] / *Mki0++;
            x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;
            x2 = (rhs[firstcol + 2] - x0 * *Mki0++ - x1 * *Mki1++) / *Mki2++;
            x3 = (rhs[firstcol + 3] - x0 * *Mki0++ - x1 * *Mki1++ - x2 * *Mki2++) / *Mki3++;

            rhs[firstcol++] = x0;
            rhs[firstcol++] = x1;
            rhs[firstcol++] = x2;
            rhs[firstcol++] = x3;

            for (k = firstcol; k < ncol; k++)
                rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++
                         - x2 * *Mki2++ - x3 * *Mki3++;
            M0 = Mki4; //4 * ldm - 10; // -1 -2 -3 -4
            ldm -= 4;
        }

        if (firstcol < ncol - 1) { /* Do 2 columns */
            Mki0 = M0;
            Mki1 = Mki0 + ldm;
            Mki2 = Mki1 + ldm - 1; // next M0

            x0 = rhs[firstcol] / *Mki0++;
            x1 = (rhs[firstcol + 1] - x0 * *Mki0++) / *Mki1++;

            rhs[firstcol++] = x0;
            rhs[firstcol++] = x1;

            for (k = firstcol; k < ncol; k++)
                rhs[k] = rhs[k] - x0 * *Mki0++ - x1 * *Mki1++;
            M0 = Mki2;
            ldm -= 2;
        }

        if (firstcol == ncol - 1) { /* Do 1 columns */
            Mki0 = M0;
            x0 = rhs[firstcol] / *Mki0;
            rhs[firstcol] = x0;
        }
    }
//    void matVecMult(
//            size_t * Li,
//            size_t * Lp,
//            size_t * Lx,
//            size_t startCol,
//            size_t supWidth,
//            size_t blockLength,
//            double * x
//            ){
//        size_t curColFirstElementIdx;
//        size_t firstRow;
//        size_t lastCol = startCol + supWidth;
//        for (size_t col = startCol; col < lastCol; ++col) {
//            curColFirstElementIdx = Lp[col];
//            firstRow = Li[curColFirstElementIdx] + compOffset + col - startCol;
//            for (size_t offset = 0; offset < blockLength; ++offset) {
//                x[firstRow + offset] -= Lx[curColFirstElementIdx + offset] * x[col];
//            }
//        }
//    }

    void dmatvec_blas (
            int ldm,	/* in -- leading dimension of M */
            int nrow,	/* in */
            int ncol,	/* in */
            double *M,	/* in */
            double *vec,	/* in */
            double *Mxvec	/* in/out */
    ) {
        double vi0, vi1, vi2, vi3, vi4, vi5, vi6, vi7;
        double *M0;
        register double *Mki0, *Mki1, *Mki2, *Mki3, *Mki4, *Mki5, *Mki6, *Mki7;
        register int firstcol = 0;
        int k;

        M0 = &M[0];
        while ( firstcol < ncol - 7 ) {	/* Do 8 columns */

            Mki0 = M0;
            Mki1 = Mki0 + ldm;
            Mki2 = Mki1 + ldm;
            Mki3 = Mki2 + ldm;
            Mki4 = Mki3 + ldm;
            Mki5 = Mki4 + ldm;
            Mki6 = Mki5 + ldm;
            Mki7 = Mki6 + ldm;

            vi0 = vec[firstcol++];
            vi1 = vec[firstcol++];
            vi2 = vec[firstcol++];
            vi3 = vec[firstcol++];
            vi4 = vec[firstcol++];
            vi5 = vec[firstcol++];
            vi6 = vec[firstcol++];
            vi7 = vec[firstcol++];

            for (k = 0; k < nrow; k++)
                Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                            + vi2 * *Mki2++ + vi3 * *Mki3++
                            + vi4 * *Mki4++ + vi5 * *Mki5++
                            + vi6 * *Mki6++ + vi7 * *Mki7++;

            M0 += 8 * ldm;
        }

        while ( firstcol < ncol - 3 ) {	/* Do 4 columns */

            Mki0 = M0;
            Mki1 = Mki0 + ldm;
            Mki2 = Mki1 + ldm;
            Mki3 = Mki2 + ldm;

            vi0 = vec[firstcol++];
            vi1 = vec[firstcol++];
            vi2 = vec[firstcol++];
            vi3 = vec[firstcol++];
            for (k = 0; k < nrow; k++)
                Mxvec[k] += vi0 * *Mki0++ + vi1 * *Mki1++
                            + vi2 * *Mki2++ + vi3 * *Mki3++ ;

            M0 += 4 * ldm;
        }

        while ( firstcol < ncol ) {		/* Do 1 column */

            Mki0 = M0;
            vi0 = vec[firstcol++];
            for (k = 0; k < nrow; k++)
                Mxvec[k] += vi0 * *Mki0++;

            M0 += ldm;
        }

    }

}