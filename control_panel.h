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

const double t_min        = 0.0;
const double t_max        = t_min + 1e5*yrtosec;    // time when computation stops [sec]
const double T_0          = 1.e8;                   // initial redshifted temperature of the star [K]
const double T_min        = 2.e5;
const double inner_radius = 5000;                   // lower boundary of the star in cm [for BSK21_1.40 NS model]
const double outer_radius = 1251100;                // upper boundary of the star in cm [for BSK21_1.40 NS model]
const double H = 1e18  ;

const unsigned int global_refinement = 4;           // refinement parameter
const unsigned int number_of_cells   = pow(2,global_refinement); // number of cells in the mesh
const unsigned int N_output = 200;                  // number of data points in output file containing cooling curves

// time step handler
//std::vector<double> time_steps =  {1,   10,     100,   1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9,  1.e10, 1.e10, 1.e10}; // array of time steps [sec]
//std::vector<double> time_steps =  {5,   50,     500,   5.e3, 5.e4, 5.e5, 5.e6, 5.e7, 5.e8, 5.e9,  5.e10, 5.e11, 10.e11}; // array of time steps [sec]
std::vector<double> time_steps =  {50,   500,     5000,   50.e3, 50.e4, 50.e5, 50.e6, 50.e7, 50.e8, 50.e9,  5.e10, 5.e11, 10.e11}; // array of time steps [sec]
std::vector<double> time_points = {1.e-7, 1.e-6, 1e-4, 1.e-3, 1.e-2, 1.e-1, 1.e0, 1.e1, 3.e2, 1.e3, 2.e3, 1.e5, 1.e6};    // time values when we change time step according to the array above [yr]
// initial time step (at t=0.0 sec) is set to be 1 sec
                                                                                                          
#endif //MAIN_CONTROL_PANEL_H
