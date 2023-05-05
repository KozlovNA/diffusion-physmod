#include<iostream>
#include<array>

using namespace std;

template<const int x_size, const int y_size>
array<array<double, y_size + 2>, x_size + 2> step(array<array<double, y_size + 2>, x_size + 2> phi, double dt, double dx, double dy)
//operator that counts field of phi configuration on next step in time
{
    array<array<double, y_size + 2>, x_size + 2> phi_new{0} //new phi field (initialized with zeros)
    for (int i = 0; i < x_size; i++)
    {
        for (int j = 0; j < y_size; j++)
        {
            phi_new[i][j] = phi[i][j] + dt * ((phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(dx*dx) + (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(dy*dy))
        }
    }
    
}

int main(int argc, char* argv[])
{
    //initializing grid
    const int x_size = 10; //horizontal size of grid
    const int y_size = 10; //vertical size of grid
    array<array<double, y_size + 2>, x_size + 2> phi{0}; //+2 stands for buffer boundaries
    //initiate boudary conditions in buffers
    //<...>
    //
}