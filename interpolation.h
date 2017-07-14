//
// Created by maryhallow on 02.12.15.
//

#ifndef MAIN_INTERPOLATION_H
#define MAIN_INTERPOLATION_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "constants.h"
#include "loaddata.cc"
#include <string>

double C(double T, double rho = 0.);
double kappa(double T, double rho = 0.);
double Q(double T, double rho = 0.);

double TiTe(double Ti);

double mass(double radius);
double phi(double radius);
double rho(double radius);

double rho_deform_grid(double radius);
double r_deform_grid(double rho);



#endif //MAIN_INTERPOLATION_H
