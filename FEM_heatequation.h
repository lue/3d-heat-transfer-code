//
// Created by maryhallow on 06.12.15.
//

#ifndef MAIN_FEM_HEATEQUATION_H
#define MAIN_FEM_HEATEQUATION_H

//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/manifold_lib.h>
#include "interpolation.cc"
#include "control_panel.h"

using namespace dealii;

#define r(i) sqrt(pow(nodeLocation[i][0],2)+pow(nodeLocation[i][1],2)+pow(nodeLocation[i][2],2))

std::vector<double> log_space(double start, double stop, unsigned int number_of_points){   // log space function

    std::vector<double> vector(number_of_points);
    double log_start = log10(start);
    double log_stop = log10(stop);
    double log_step = (log_stop - log_start)/(number_of_points-1);

    for (unsigned  i=0; i<number_of_points; i++) {
        vector[i] = pow(10., log_start + i * log_step);
    }

    return vector;
}

template <int dim>
class FEM
{
public:

    FEM (double Alpha);
    ~FEM();

    void generate_mesh();
    void setup_system();
    void assemble_system();
    void apply_initial_conditions();
    void solve_trans();
    void output_dat_file(unsigned int index);
    void output_trans_results(unsigned int index);

    std::vector<double> snapshot = log_space(2, t_max, N_output);

    Triangulation<dim> triangulation;                         // mesh
    FESystem<dim>      fe;                                    // FE element
    DoFHandler<dim>    dof_handler;                           // Connectivity matrices

    QGauss<dim>   quadrature_formula;                         // Quadrature
    QGauss<dim-1> face_quadrature_formula;                    // Face Quadrature

    //Data structures
    SparsityPattern      sparsity_pattern;                    // Sparse matrix pattern
    SparseMatrix<double> M, K, system_matrix;                 // Global stiffness matrix - Sparse matrix - used in the solver
    Vector<double>       D_trans, V_trans, F, RHS;            // Global vectors - Solution vector (D) and Global force vector (F)
    std::map<unsigned int,double> boundary_values_of_D;       // Map of dirichlet boundary conditions for the time derivative of temperature

    Table<2,double>	nodeLocation;	                      // Table of the coordinates of nodes by global dof number

    double alpha, eff_radius; 	                              // Specifies the Euler method, 0 <= alpha <= 1

    double redshift = sqrt(1 - (2*G*mass(outer_radius))/(outer_radius*c*c));

    //solution name array
    std::vector<std::string> nodal_solution_names;
    std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

struct grid_transform                                         // grid deform function ( quater hyper shell )
{
    Point<3> operator() (const Point<3> &in) const
    {
        double temp_r = sqrt(pow(in(0),2)+pow(in(1),2)+pow(in(2),2));
        double dr_lin = (outer_radius - inner_radius)/number_of_cells;
        double drho_log = (-(log10(rho_deform_grid(outer_radius)) - log10(rho_deform_grid(inner_radius))))/number_of_cells;
        double n = std::floor((temp_r-inner_radius)/dr_lin + 0.5);
        double r_new = r_deform_grid(pow(10,log10(rho_deform_grid(inner_radius))-n*drho_log));

        double r_ratio = r_new/temp_r;
        
        if (temp_r==inner_radius){
            return Point<3> (in(0), in(1),in(2));
        }
        else if (temp_r==outer_radius){
            return Point<3> (in(0), in(1),in(2));
        }
        else{
            return Point<3> (in(0)*r_ratio, in(1)*r_ratio,in(2)*r_ratio);
        }
    }
};


#endif //MAIN_FEM_HEATEQUATION_H
