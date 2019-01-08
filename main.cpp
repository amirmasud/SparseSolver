#include <iostream>
#include "Utils.h"

int main() {
    auto *cscMatrix = utils::CSCMatrix::createFromFile("b.mtx");
    return 0;
}
