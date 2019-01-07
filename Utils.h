//
// Created by Amir Masud on 2/7/19.
//

#ifndef SPARSE_SOLVER_UTILS_H
#define SPARSE_SOLVER_UTILS_H


#include <string>


namespace utils {

enum class Format {
    CSC
};

class Matrix {
public:
    virtual Format getFormat() const;

protected:
    explicit Matrix(Format format);

private:
    Format m_format;
};

//template <class T>
class CSCMatrix : public Matrix {
public:
    // Can be constructed from other formats too.
    static CSCMatrix *create(unsigned int n, int *Lp, int *Li, double *Lx);
    static CSCMatrix *readFromFile(const std::string &name);
    void print();

    long unsigned int n;
    long unsigned int nnz;
    int *Lp;
    int *Li;
    double *Lx;

private:
    CSCMatrix();
    CSCMatrix(unsigned int m_n, int *m_Lp, int *m_Li, double *m_Lx);
    void init();
    void initWithFile(const std::string &name);
};

}

namespace matrix_hf {

bool readMatrix(const std::string &fName, size_t &n, size_t &NNZ, int* &col,
                int* &row, double* &val);

}

#endif //SPARSE_SOLVER_UTILS_H
