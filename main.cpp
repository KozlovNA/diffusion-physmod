#include <iostream>
#include <array>
#include "hdf5/serial/hdf5.h"
// #include "hdf5/serial/H5Cpp.h"
#include <fstream>
#include <string>
#include <omp.h>
#include <sys/time.h>

using namespace std;

#include "nlohmann/json.hpp"

using json = nlohmann::json;


using namespace std;

template<const int x_size, const int y_size>
array<array<double, y_size + 2>, x_size + 2> step(array<array<double, y_size + 2>, x_size + 2>& phi, double& dt, double& dx, double& dy,
 const double& phi_0, const double& phi_1)
//operator that counts field of phi configuration on next step in time
{   
    array<array<double, y_size + 2>, x_size + 2> phi_new{0}; //new phi field (initialized with zeros)
    //#pragma omp parallel for num_threads(12)
        for (int i = 0; i < x_size+2; i++)
        {
            for (int j = 0; j < y_size+2; j++)
            {
                phi_new[i][j] = phi[i][j] + dt * ((phi[i+1][j] - 2*phi[i][j] + phi[i-1][j])/(dx*dx) + (phi[i][j+1] - 2*phi[i][j] + phi[i][j-1])/(dy*dy));
            }
        }
    //#pragma omp parallel for num_threads(12)       
    for (int j = 0; j < y_size+2; j++) 
    {
        phi_new[0][j] = 2*phi_0 - phi[1][j];                    //vertical bounds
        phi_new[x_size+1][j] = 2*phi_1 - phi[x_size][j];
    }
    //#pragma omp parallel for num_threads(12) 
    for (int i = 0; i < x_size+2; i++)
    {
        phi_new[i][0] = phi_new[i][1];              //horizontal bounds
        phi_new[i][y_size+1] = phi_new[i][y_size];
    }
    phi_new[x_size/2][y_size/2] = 2;
    return phi_new;
}

template<const int x_size, const int y_size>
void write_to_file(array<array<double, y_size + 2>, x_size + 2>& phi, string OUTPATH)
{
    ofstream out;
        out.open(OUTPATH);
        if (out.is_open())
        {
            for (int j = 1; j < y_size; j++){
                for (int i = 1; i < x_size; i++)
                out <<  phi[i][j] <<'\n';
            }
        }
        out.close();
} 

int main(int argc, char* argv[])
{    
    const double phi_0 = 2;
    const double phi_1 = 0;
    const int x_size = 100; //horizontal size of grid
    const int y_size = 100;//vertical size of grid
    double dt = 0.1;
    double dx = 1;
    double dy = 1;
    double num_of_iter = 100000;
    array<array<double, y_size + 2>, x_size + 2> phi{0}; //+2 stands for buffer boundaries
    //initiate boundary conditions in buffers
    #pragma omp parallel for num_threads(12)
    for (int j = 0; j < y_size+2; j++) 
    {
        phi[0][j] = 2*phi_0 - phi[1][j];                    //vertical bounds
        phi[x_size+1][j] = - 2*phi_1 + phi[x_size][j];
    }
    #pragma omp parallel for num_threads(12) 
    for (int i = 0; i < x_size+2; i++)
    {
        phi[i][0] = phi[i][1];              //horizontal bounds
        phi[i][y_size+1] = phi[i][y_size];
    }
    //count grid configuration
    double time = 0;
    struct timeval begin, end;
        gettimeofday(&begin, 0);
    for (int i = 0; i < num_of_iter; i++)
    {   
        phi = step<x_size, y_size>(phi, dt, dx, dy, phi_0, phi_1);
    }
    gettimeofday(&end, 0);
        long seconds = end.tv_sec - begin.tv_sec;
        long microseconds = end.tv_usec - begin.tv_usec;
        double elapsed = seconds + microseconds*1e-6;
        //time+=elapsed;
    cout << elapsed << '\n';

    write_to_file<x_size, y_size>(phi, "grid150000steps.txt");

    printf("\nover\n");
}