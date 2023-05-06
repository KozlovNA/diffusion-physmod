#include <iostream>
#include <array>
#include <hdf5/serial/hdf5.h>
#include <fstream>
#include <string>

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
        phi_new[0][j] = 2*phi_0 - phi[1][j];                    //vertical bounds
        phi_new[x_size+1][j] = - 2*phi_1 + phi[x_size][j];
    } 
    for (int i = 0; i < x_size; i++)
    {
        phi_new[i][1] = phi_new[i][0];              //horizontal bounds
        phi_new[i][y_size] = phi_new[i][y_size+1];
    }
    return phi_new;
}

int main(int argc, char* argv[])
{    
    const double phi_0 = 1;
    const double phi_1 = 0;
    const int x_size = 10; //horizontal size of grid
    const int y_size = 10;//vertical size of grid
    double dt = 0.01;
    double dx = 0.001;
    double dy = 0.001;
    array<array<double, y_size + 2>, x_size + 2> phi{0}; //+2 stands for buffer boundaries
    //initiate boundary conditions in buffers
    for (int j = 0; j < y_size; j++) 
    {
        phi[0][j] = 2*phi_0 - phi[1][j];                    //vertical bounds
        phi[x_size+1][j] = - 2*phi_1 + phi[x_size][j];
    } 
    for (int i = 0; i < x_size; i++)
    {
        phi[i][1] = phi[i][0];              //horizontal bounds
        phi[i][y_size] = phi[i][y_size+1];
    }
    phi = step<x_size, y_size>(phi, dt, dx, dy, phi_0, phi_1);
    printf("\nover\n");
}