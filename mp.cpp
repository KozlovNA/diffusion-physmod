#include <iostream>
#include <array>
#include "hdf5/serial/hdf5.h"
// #include "hdf5/serial/H5Cpp.h"
#include <fstream>
#include <string>
#include <omp.h>

int main()
{
    #pragma omp parallel
    std::cout << "hello" << '\n';
}