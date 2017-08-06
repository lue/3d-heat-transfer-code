//
// Created by maryhallow on 06.12.15.
//

#ifndef MAIN_CONTROL_PANEL_H
#define MAIN_CONTROL_PANEL_H
                           
                                                 
const unsigned int order    = 1;                    // order of basis functions
const unsigned int quadRule = 2;                    // order of quadrature rule
const int dimension = 3;                            // number of dimensions in the problem
double Euler_scheme = 1.0 ;                         // time integration method (out of Euler family)

#define CoolingFile "cooling_curve.dat"

const double t_max        = 4e6*yrtosec;            // time when computation stops [sec]
const double T_0          = 1.e10;                  // initial redshifted temperature of the star [K]
const double T_min        = 2.e5;
const double inner_radius = 5000;                   // lower boundary of the star in cm [for BSK21_1.40 NS model]
const double outer_radius = 1251100;                // upper boundary of the star in cm [for BSK21_1.40 NS model]

const unsigned int global_refinement = 3;           // refinement parameter
const unsigned int number_of_cells   = pow(2,global_refinement); // number of cells in the mesh
const unsigned int N_output = 100;                  // number of data points in output file containing cooling curves

// time step handler
//std::vector<double> time_steps =  {1,   10,     100,   1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9,  1.e10, 1.e10, 1.e10}; // array of time steps [sec]
std::vector<double> time_steps =  {10,   100,     1000,   10.e3, 10.e4, 10.e5, 10.e6, 10.e7, 10.e8, 10.e9,  10.e10, 10.e11, 10.e11}; // array of time steps [sec]
std::vector<double> time_points = {1.e-7, 1.e-6, 1e-4, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1, 3.e2, 1.e3, 2.e3, 1.e5, 1.e6};    // time values when we change time step according to the array above [yr]
// initial time step (at t=0.0 sec) is set to be 1 sec
                                                                                                          
#endif //MAIN_CONTROL_PANEL_H
