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

    virtual ~Matrix() = default;

protected:
    explicit Matrix(Format format);

private:
    Format m_format;
};

//template <class T>
class CSCMatrix : public Matrix {
public:
    // Can be constructed from other formats too.
    static CSCMatrix *create(size_t colCount, size_t nnz,
                             int *Lp, int *Li, double *Lx);
    static CSCMatrix *createFromFile(const std::string &name);

    ~CSCMatrix() override;

    size_t rowCount, colCount;
    size_t nnz;
    int *Lp;
    int *Li;
    double *Lx;

protected:
    CSCMatrix();
    CSCMatrix(size_t m_colCount, size_t m_nnz, int *m_Lp, int *m_Li, double *m_Lx);
    virtual void init();
    virtual void initWithFile(const std::string &name);
};

class CSCVector : public CSCMatrix {
public:
    static CSCVector *createFromFile(const std::string &name);

    double *v;

    ~CSCVector() override;

private:
    CSCVector();
    void initWithFile(const std::string &name) override;
    void initVector();
};

}


namespace matrix_hf {

bool readMatrix(const std::string &fName, size_t &rowCount, size_t &colCount,
                size_t &NNZ, int* &col, int* &row, double* &val);

}

#endif //SPARSE_SOLVER_UTILS_H
