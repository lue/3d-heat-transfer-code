//
// Created by maryhallow on 02.12.15.
//
#ifndef MAIN_INTERPOLATION_CC
#define MAIN_INTERPOLATION_CC
#include "interpolation.h"

double Q(double T, double rho) // neutrino emissivity (redshifted temperature, density)
{
    return pow(10, gsl_interp2d_eval(Q_spline, grid1_interp, grid2_interp, Q_interp, log10(T), log10(rho), xacc, yacc));
}

double C(double T, double rho) // heat capacity (redshifted temperature, density)
{
    return pow(10, gsl_interp2d_eval(C_spline, grid1_interp, grid2_interp, C_interp, log10(T), log10(rho), xacc, yacc));
}

double kappa(double T, double rho) // thermal conductivity (redshifted temperature, density)
{
    return pow(10, gsl_interp2d_eval(k_spline, grid1_interp, grid2_interp, k_interp, log10(T), log10(rho), xacc, yacc));
}

double rho_deform_grid(double radius) // linear interpolation function "density (radius)". Only used at the start of the simulation for grid transformation.
{

    double f_2_max = model_rho[len_model-1];
    double f_1_max = model_rho[len_model-2];
    double x_2_max = model_r[len_model-1];
    double x_1_max = model_r[len_model-2];

    double f_2_min = model_rho[2];
    double f_1_min = model_rho[1];
    double x_2_min = model_r[2];
    double x_1_min = model_r[1];

    if(log10(radius)>x_2_max){
        return pow(10, f_1_max + (f_2_max - f_1_max) * (log10(radius) - x_1_max) / (x_2_max - x_1_max));
    }
    else if(log10(radius)<x_1_min){
        return pow(10, f_1_min + (f_2_min - f_1_min) * (log10(radius) - x_1_min) / (x_2_min - x_1_min));
    }
    else{
        return pow(10, gsl_spline_eval(rho_r_spline, log10(radius), racc));
    }
}

double r_deform_grid(double rho) // linear interpolation function "radius (density)". Only used at the start of the simulation for grid transformation.
{

    //double x_2 = pow(10,model_rho[len_model-1]);
    //double x_1 = pow(10,model_rho[1]);
    //double f_2 = pow(10,model_r[len_model-1]);
    //double f_1 = pow(10,model_r[1]);
    
    double x_2_max = model_rho[1];
    double x_1_max = model_rho[2];
    double f_2_max = model_r[1];
    double f_1_max = model_r[2];

    double x_2_min = model_rho[len_model-2];
    double x_1_min = model_rho[len_model-1];
    double f_2_min = model_r[len_model-2];
    double f_1_min = model_r[len_model-1];

    if(log10(rho)>x_2_max){
        return pow(10, f_1_max + (f_2_max - f_1_max) * (log10(rho) - x_1_max) / (x_2_max - x_1_max));
    }
    else if(log10(rho)<x_1_min){
        return pow(10, f_1_min + (f_2_min - f_1_min) * (log10(rho) - x_1_min) / (x_2_min - x_1_min));
    }
    else{
        return pow(10, gsl_spline_eval(r_rho_spline, log10(rho), rhoacc));
    }
}

double rho(double radius) // density (radius)
{
    return pow(10, gsl_spline_eval(rho_r_spline, log10(radius), racc));
}

double mass(double radius) // mass enclosed within a sphere with a radius r "mass (radius)"
{
    return pow(10, gsl_spline_eval(m_r_spline, log10(radius), racc));
}

double phi(double radius) // gravitational potential (radius)
{
    return gsl_spline_eval(phi_r_spline, log10(radius), racc);
}

double TiTe(double T) // external temperature (internal temperature)
{
    return pow(10, gsl_spline_eval(TiTe_spline, log10(T), TiTeacc));
}



#endif
