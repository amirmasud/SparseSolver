#include <iostream>
#include "Utils.h"

int main() {
    auto *cscMatrix = utils::CSCMatrix::readFromFile("b.mtx");
    cscMatrix->print();
    std::cout << "salam" << std::endl;
    return 0;
}
