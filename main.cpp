#include <iostream>
#include <array>
#include <hdf5/serial/hdf5.h>

using namespace std;

#include "nlohmann/json.hpp"

using json = nlohmann::json;


using namespace std;

template<const int x_size, const int y_size>
array<array<double, y_size + 2>, x_size + 2> step(array<array<double, y_size + 2>, x_size + 2> phi, double dt, double dx, double dy,
 const double phi_0, const double phi_1)
//operator that counts field of phi configuration on next step in time
{
    array<array<double, y_size + 2>, x_size + 2> phi_new{0}; //new phi field (initialized with zeros)
    for (int i = 0; i < x_size; i++)
    {
        for (int j = 0; j < y_size; j++)
        {
            phi_new[i][j] = phi[i][j] + dt * ((phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(dx*dx) + (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(dy*dy));
        }
    }
    for (int j = 0; j < y_size; j++)
    {
        phi_new[0][j] = 2*phi_0 - phi[1][j];
    } 
    return phi_new;
}

int main(int argc, char* argv[])
{
    //initializing grid
    const int x_size = 10; //horizontal size of grid
    const int y_size = 10; //vertical size of grid
    const double phi_0 = 1;
    const double phi_1 = 0;
    array<array<double, y_size + 2>, x_size + 2> phi{0}; //+2 stands for buffer boundaries
    //initiate boundary conditions in buffers
    //<...>
    //
}