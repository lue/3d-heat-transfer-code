//
// Created by maryhallow on 29.11.15.
//
#ifndef MAIN_LOADDATA_CC
#define MAIN_LOADDATA_CC
#include "loaddata.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

std::vector<double> model_r;                                                  // This vector contains radius values (radius = distance form the center of the star)
std::vector<double> model_rho;                                                // This vector contains density values
std::vector<double> model_m;                                                  // This vector contains mass values that is enclosed within a sphere with a radius r if we write m(r)
std::vector<double> model_phi;                                                // This vector contains gravitational potential values

unsigned int len_model = 481;                                                 // The length of each vector out of four above 

// (same four vectors as above) created for gsl splines initialization
const size_t len_interp_model = 481;   
double *r_interp = (double *) malloc(len_interp_model  * sizeof(double));             
double *rho_interp = (double *) malloc(len_interp_model  * sizeof(double));          
double *m_interp = (double *) malloc(len_interp_model  * sizeof(double));   
double *phi_interp = (double *) malloc(len_interp_model  * sizeof(double));

double *r_interp_reversed = (double *) malloc(len_interp_model  * sizeof(double));             
double *rho_interp_reversed = (double *) malloc(len_interp_model  * sizeof(double)); 

// To create splines for NS model
const gsl_interp_type *interpolation_type_star_model = gsl_interp_linear;
gsl_interp_accel *racc = gsl_interp_accel_alloc();
gsl_spline *m_r_spline = gsl_spline_alloc(interpolation_type_star_model, len_interp_model);
gsl_spline *rho_r_spline = gsl_spline_alloc(interpolation_type_star_model, len_interp_model);
gsl_spline *phi_r_spline = gsl_spline_alloc(interpolation_type_star_model, len_interp_model);

gsl_interp_accel *rhoacc = gsl_interp_accel_alloc();
gsl_spline *r_rho_spline = gsl_spline_alloc(interpolation_type_star_model, len_interp_model);



// This is for TiTe spline initializaion
const size_t TiTe_interp_length = 59;                                         // The length of the Ti- and Te-vectors below [same as the length of TiTe file]
double *Ti_interp = (double *) malloc(TiTe_interp_length * sizeof(double));   // This vector contains internal temperature values
double *Te_interp = (double *) malloc(TiTe_interp_length * sizeof(double));   // This vector contains external temperature values    

// To create Ti-Te spline
const gsl_interp_type *interpolation_type_TiTe = gsl_interp_linear;
gsl_interp_accel *TiTeacc = gsl_interp_accel_alloc();
gsl_spline *TiTe_spline = gsl_spline_alloc(interpolation_type_TiTe, TiTe_interp_length);                                                  



// This is for splines Q, C and k (functions of density and temperature)
std::vector<double> grid1;                                                    // This vector contains redshifted temperature values
unsigned int len_x;                                                           // The length of grid1 vector

std::vector<double> grid2;                                                    // This vector contains density values
unsigned int len_y;                                                           // The length of grid2 vector

// All three matrices below have a size [len_x X len_y]
std::vector<double> Q_values;                                                 // This matrix contains redshifted neutrino emissivity values
std::vector<double> C_values;                                                 // This matrix contains specific heat capacity per unit volume values
std::vector<double> k_values;                                                 // This matrix contains thermal conductivity values

const gsl_interp2d_type *interpolation_type = gsl_interp2d_bilinear;

const size_t nx_interp = 150; /* x grid points along Temperature (length of file1.dat) */            
const size_t ny_interp = 481; /* y grid points along Density (length of file2.dat) */

double *grid1_interp = (double *) malloc(nx_interp * sizeof(double));
double *grid2_interp = (double *) malloc(ny_interp * sizeof(double));

double *Q_interp = (double *) malloc(nx_interp * ny_interp * sizeof(double));
double *C_interp = (double *) malloc(nx_interp * ny_interp * sizeof(double));
double *k_interp = (double *) malloc(nx_interp * ny_interp * sizeof(double));

gsl_interp2d *Q_spline = gsl_interp2d_alloc(interpolation_type, nx_interp, ny_interp);
gsl_interp2d *C_spline = gsl_interp2d_alloc(interpolation_type, nx_interp, ny_interp);
gsl_interp2d *k_spline = gsl_interp2d_alloc(interpolation_type, nx_interp, ny_interp);

gsl_interp_accel *xacc = gsl_interp_accel_alloc();
gsl_interp_accel *yacc = gsl_interp_accel_alloc();


void loadfiles(void)
{
    int i, j;

    std::ifstream inputFile1("./datafiles/file1.dat");
    if (inputFile1.is_open()) {
        double temp;
        i = 0;

        while (inputFile1 >> temp) {
            grid1.push_back(log10(temp));
            grid1_interp[i] = log10(temp);
            i++;
        }
        inputFile1.close();
    }
    else {
        printf("Temperature file is empty! Simulation is terminated.");
        exit(0);
    }
    len_x = grid1.size();


    std::ifstream inputFile2("./datafiles/file2.dat");
    if (inputFile2.is_open()) {
        double temp;
        i = 0;

        while (inputFile2 >> temp) {
            grid2.push_back(log10(temp));
            grid2_interp[i] = log10(temp);
            i++;
        }
        inputFile2.close();
    }
    else {
        printf("Density file is empty! Simulation is terminated.");
        exit(0);
    }
    len_y = grid2.size();


    std::ifstream inputFile3("./datafiles/file3.dat");
    if (inputFile3.is_open()) {
        double temp;
        i = 0;

        while(inputFile3 >> temp){
            Q_values.push_back(log10(temp));
            Q_interp[i] = log10(temp);
            i++;
        }
        inputFile3.close();
        gsl_interp2d_init(Q_spline, grid1_interp, grid2_interp, Q_interp, nx_interp, ny_interp);
    }
    else {
        printf("Qfile is empty! Simulation is terminated.");
        exit(0);
    }

    std::ifstream inputFile4("./datafiles/file4.dat");
    if (inputFile4.is_open()) {
        double temp;
        i = 0;

        while(inputFile4 >> temp){
            C_values.push_back(log10(temp));
            C_interp[i] = log10(temp);
            i++;
        }
        inputFile4.close();
        gsl_interp2d_init(C_spline, grid1_interp, grid2_interp, C_interp, nx_interp, ny_interp);
    }
    else {
        printf("Cfile is empty! Simulation is terminated.");
        exit(0);
    }


    std::ifstream inputFile5("./datafiles/file5.dat");
    if (inputFile5.is_open()) {
        double temp;
        i = 0;

        while (inputFile5 >> temp){
            k_values.push_back(log10(temp));
            k_interp[i] = log10(temp);
            i++;
        }
        inputFile5.close();
        gsl_interp2d_init(k_spline, grid1_interp, grid2_interp, k_interp, nx_interp, ny_interp);
    }
    else {
        printf("kfile is empty! Simulation is terminated.");
        exit(0);
    }


    std::ifstream inputFile6("./datafiles/tite.dat");
    if (inputFile6.is_open()) {
        double temp;
        i = 0;
        unsigned int swtch = 0;
        while (inputFile6 >> temp) {
            if(swtch){
                swtch = 0;
                Te_interp[i] = temp;
                i++;
            }
            else{
                swtch = 1;
                Ti_interp[i] = temp;
            }
        }
        inputFile6.close();
    }
    else {
        printf("TiTefile is empty! Simulation is terminated.");
        exit(0);
    }
    gsl_spline_init(TiTe_spline, Ti_interp, Te_interp, TiTe_interp_length);

    std::ifstream inputFile7("./datafiles/BSK21_1.40.dat");
    if (inputFile7.is_open()) {
        double temp;
        inputFile7 >> temp;
        inputFile7 >> temp;
        i = 0;
        unsigned int swtch = 0;
        while (inputFile7 >> temp) {
            if(swtch==0){
                model_m.push_back(log10(temp*MSun));
                m_interp[i] = log10(temp*MSun);
                swtch ++;
            }
            else if(swtch==1){
                model_r.push_back(log10(temp* 1.e5));
                r_interp[i] = log10(temp* 1.e5);
                r_interp_reversed[len_model-1-i] = log10(temp* 1.e5);
                swtch ++;
            }
            else if(swtch==3){
                model_rho.push_back(log10(temp));
                rho_interp[i] = log10(temp);
                rho_interp_reversed[len_model-1-i] = log10(temp);
                swtch ++;
            }
            else if(swtch==4){
                model_phi.push_back(temp); 
                phi_interp[i] = temp;
                swtch ++;
            }
            else if(swtch==7){
                swtch=0;
                i++;
            }
            else{
                swtch ++;
            }
        }
        inputFile7.close();
    }
    else {
        printf("Modelfile is empty! Simulation is terminated.");
        exit(0);
    }
    //len_model = model_rho.size();
    gsl_spline_init(m_r_spline, r_interp, m_interp, len_interp_model);
    gsl_spline_init(rho_r_spline, r_interp, rho_interp, len_interp_model);
    gsl_spline_init(phi_r_spline, r_interp, phi_interp, len_interp_model);
    gsl_spline_init(r_rho_spline, rho_interp_reversed, r_interp_reversed, len_interp_model);
}

#endif
