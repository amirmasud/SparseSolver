//
// Created by Amir Masud on 2/7/19.
//

#ifndef SPARSE_SOLVER_UTILS_H
#define SPARSE_SOLVER_UTILS_H


#include <string>
#include <vector>


namespace utils {

enum class Format {
    CSC
};

class Matrix {
public:
    virtual Format getFormat() const;

    virtual void dump() {}

    virtual Matrix *clone() const;

    virtual ~Matrix() = default;

    size_t rowCount, colCount;

protected:
    void cloneAllocations(const Matrix *matrix) {}
    Matrix(const Matrix &matrix) = default;
    explicit Matrix(Format format, size_t m_rowCount, size_t m_colCount);

private:
    Format m_format;
};

class CSCVector;

//template <class T>
class CSCMatrix : public Matrix {
public:
    CSCMatrix(const CSCMatrix &cscMatrix);
    // Can be constructed from other formats too.
    static CSCMatrix *create(size_t rowCount, size_t colCount, size_t nnz,
                             size_t *Lp, size_t *Li, double *Lx);
    static CSCMatrix *createFromFile(const std::string &name);

    std::vector<double> operator*(const CSCVector &cscVector) const;

    Matrix *clone() const override;

    ~CSCMatrix() override;

    size_t nnz;
    size_t *Lp;
    size_t *Li;
    double *Lx;

protected:
    CSCMatrix();
    CSCMatrix(size_t m_rowCount, size_t m_colCount,
              size_t m_nnz,  size_t *m_Lp, size_t *m_Li, double *m_Lx);
    virtual void init();
    virtual void initWithFile(const std::string &name);
    void cloneAllocations(const CSCMatrix *cscMatrix);
};

class CSCVector : public CSCMatrix {
public:
    CSCVector(const CSCVector &cscVector);
    static CSCVector *create(size_t rowCount, size_t colCount, size_t nnz,
                             size_t *Lp, size_t *Li, double *Lx, double *m_v);
    static CSCVector *createFromFile(const std::string &name);

    Matrix *clone() const override;

    void dump() override;

    double *v;

    ~CSCVector() override;

private:
    CSCVector();
    CSCVector(size_t rowCount, size_t colCount, size_t nnz,
              size_t *Lp, size_t *Li, double *Lx, double *m_v);
    void init() override;
    void initWithFile(const std::string &name) override;
    void initVector();
    void cloneAllocations(const CSCVector *cscVector);
};

}


namespace matrix_hf {

bool readMatrix(const std::string &fName, size_t &rowCount, size_t &colCount,
                size_t &NNZ, size_t * &col, size_t * &row, double* &val);

}

#endif //SPARSE_SOLVER_UTILS_H
