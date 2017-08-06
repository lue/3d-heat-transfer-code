#ifndef MAIN_FEM_HEATEQUATION_H
#define MAIN_FEM_HEATEQUATION_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include "interpolation.cc"
#include "control_panel.h"

namespace LA
{
using namespace dealii::LinearAlgebraTrilinos;
}

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_solver.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <fstream>
#include <iostream>

using namespace dealii;

template <int dim>
class LaplaceProblem
{
    public:
    LaplaceProblem (double Alpha);
    ~LaplaceProblem ();
    
    void generate_mesh();
    void apply_initial_conditions();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    MPI_Comm                                  mpi_communicator;

    parallel::distributed::Triangulation<dim> triangulation;
 
    FESystem<dim>                             fe;
    DoFHandler<dim>                           dof_handler;

    QGauss<dim>   quadrature_formula;                       
    QGauss<dim-1> face_quadrature_formula;                

    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;

    ConstraintMatrix                          constraints;

    LA::MPI::SparseMatrix                     system_matrix, M, K;
    LA::MPI::Vector                           RHS, D_trans, V_trans, F, D_tilde;

    std::map<unsigned int,double> boundary_values_of_D;     
    Table<2,double>	nodeLocation;	                    
    double alpha, eff_radius; 	                            
    
    double redshift = sqrt(1 - (2*G*mass(outer_radius))/(outer_radius*c*c));

    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;
};


struct grid_transform                                     
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


std::vector<double> log_space(double start, double stop, unsigned int number_of_points){  

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
void LaplaceProblem<dim>::generate_mesh(){ 

    double  x_center = 0.,
            y_center = 0.,
            z_center = 0.;

    Point<dim,double> center(x_center,y_center,z_center);
    GridGenerator::quarter_hyper_shell (triangulation, center, inner_radius, outer_radius, 0, true);
    static const SphericalManifold<3> manifold_description(center);
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold (0, manifold_description);
    triangulation.refine_global(global_refinement);
    GridTools::transform(grid_transform(), triangulation);               
}


template <int dim>
void LaplaceProblem<dim>::apply_initial_conditions(){

    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    const unsigned int dofs_per_elem = fe.dofs_per_cell;                     
    std::vector<unsigned int> local_dof_indices (dofs_per_elem);             

    typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), endc = dof_handler.end();

    for (;elem!=endc; ++elem){
        if (elem->is_locally_owned()){
            elem->get_dof_indices (local_dof_indices);
            for(unsigned int C=0; C<dofs_per_elem; C++){
                completely_distributed_solution(local_dof_indices[C]) = T_0;
            }    
        }
    }
    constraints.distribute (completely_distributed_solution);

    D_trans = completely_distributed_solution;
}

template <int dim>
LaplaceProblem<dim>::LaplaceProblem (double Alpha)
    :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator),
    fe (FE_Q<dim>(order), 1),
    dof_handler (triangulation),
    quadrature_formula(quadRule),
    face_quadrature_formula(quadRule),
    pcout (std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator)
            == 0)),
    computing_timer (mpi_communicator,
                     pcout,
                     TimerOutput::summary,
                     TimerOutput::wall_times)
  {alpha = Alpha;}



template <int dim>
LaplaceProblem<dim>::~LaplaceProblem ()
{
    dof_handler.clear ();
}


template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
    TimerOutput::Scope t(computing_timer, "setup");

    dof_handler.distribute_dofs (fe);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);

    D_trans.reinit (locally_owned_dofs,locally_relevant_dofs, mpi_communicator) ;
    V_trans.reinit (locally_owned_dofs, mpi_communicator);
    D_tilde.reinit (locally_owned_dofs, mpi_communicator);
    RHS.reinit (locally_owned_dofs, mpi_communicator);
    F.reinit (locally_owned_dofs, mpi_communicator);

    pcout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    constraints.close ();

  
    DynamicSparsityPattern dsp (locally_relevant_dofs);

    DoFTools::make_sparsity_pattern (dof_handler, dsp,
                                     constraints, false);
    SparsityTools::distribute_sparsity_pattern (dsp,
                                                dof_handler.n_locally_owned_dofs_per_processor(),
                                                mpi_communicator,
                                                locally_relevant_dofs);

    system_matrix.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          dsp,
                          mpi_communicator);

    K.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          dsp,
                          mpi_communicator);
    M.reinit (locally_owned_dofs,
                          locally_owned_dofs,
                          dsp,
                          mpi_communicator);
}



template <int dim>
void LaplaceProblem<dim>::assemble_system(){

    TimerOutput::Scope t(computing_timer, "Assembly");

    M=0; K=0; F=0;

    FEValues<dim> fe_values(fe,                                               // we update fe values parameters
                            quadrature_formula,
                            update_values |
                            update_gradients |
                            update_quadrature_points |
                            update_JxW_values);

    FEFaceValues<dim> fe_face_values (fe,                                     // we update fe face values parameters
                                      face_quadrature_formula,
                                      update_values |
                                      update_quadrature_points |
                                      update_JxW_values);

    const unsigned int num_face_quad_pts = face_quadrature_formula.size();    // Total number of quad points in the face
    const unsigned int faces_per_elem = GeometryInfo<dim>::faces_per_cell;    // Faces per cell (for our geometry it is 6)

    const unsigned int dofs_per_elem = fe.dofs_per_cell;                      // This gives you dofs per element
    unsigned int       num_quad_pts = quadrature_formula.size();              // Total number of quad points in the element

    FullMatrix<double> Mlocal (dofs_per_elem, dofs_per_elem);                 // Mass matrix
    FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);                 // Stiffness matrix
    Vector<double>     Flocal (dofs_per_elem);                                // Force vector

    std::vector<unsigned int> local_dof_indices (dofs_per_elem);              // This relates local dof numbering to global dof numbering
    typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), endc = dof_handler.end();
    for (;elem!=endc; ++elem){
        if (elem->is_locally_owned()){

        fe_values.reinit(elem);
        elem->get_dof_indices (local_dof_indices);

        Mlocal = 0.;
        Klocal = 0.;
        Flocal = 0.;

        for(unsigned int q=0; q<num_quad_pts; q++){                           // We build our M and K matrices and F vector

            double D_trans_at_q = 0.;                                         // We interpolate temperature at quadrature point q using dealii basis functions

            for(unsigned int C=0; C<dofs_per_elem; C++){
                D_trans_at_q += D_trans(local_dof_indices[C])*fe_values.shape_value(C,q);
            }

            double r_q   = sqrt(pow(fe_values.quadrature_point(q)[0],2) + pow(fe_values.quadrature_point(q)[1],2) + pow(fe_values.quadrature_point(q)[2],2));
                                                                             
                                                                              // We get r values at quadrature points
            double rho_q = rho(r_q);                                          // We get density that corresponds to r at quadrature point q

            double phi_q = phi(r_q);                                          // Same with gravitational potential
            double m_q   = mass(r_q);                                         // Same with mass that is enclosed within a sphere with radius r_q
            double sqrt_q = sqrt(1 - (2*G*m_q)/(r_q*c*c));
            double exp_q = exp(phi_q);                                        // We use this factor when we, for example, relate the local temperature to the redshifted temperature

            double kappa_precalc = kappa(D_trans_at_q,rho_q);                 // Precalculating kappa outside of the loops
            double C_precalc = C(D_trans_at_q,rho_q);                         // Precalculating C outside of the loops

            for(unsigned int A=0; A<fe.dofs_per_cell; A++){                   // Loop over local DOFs to populate Mlocal, Klocal and Flocal at quadrature point q
                Flocal[A] -= fe_values.shape_value(A,q)*Q(D_trans_at_q,rho_q)/sqrt_q*fe_values.JxW(q);
                for(unsigned int B=0; B<fe.dofs_per_cell; B++){
                    Mlocal[A][B] += fe_values.shape_value(A,q)*fe_values.shape_value(B,q)*C_precalc/sqrt_q*fe_values.JxW(q);
                    for(unsigned int i=0; i<dim; i++){                        // Kappa is a tensor so we have to create additional loop over its own elements
                        for(unsigned int j=0; j<dim; j++){
                            if(i==j){                                         // We use this condition because all off diagonal elements of conductivity tensor kappa are zero. So
                                                                              // they do not contribute to Klocal
                                Klocal[A][B] += fe_values.shape_grad(A,q)[i]*kappa_precalc*sqrt_q*exp_q*fe_values.shape_grad(B,q)[j]*fe_values.JxW(q);
                            }
                        }
                    }
                }
            }
        }


        for (unsigned int f=0; f < faces_per_elem; f++){                      // This loop stands for Neumann condition at the boundary of the star
            fe_face_values.reinit (elem, f);
            if(elem->face(f)->at_boundary()){                                 // We check whether we are at the boundary or not. If we are not, nothing happens.
                if(elem->face(f)->boundary_id() == 1){
                    for (unsigned int q=0; q<num_face_quad_pts; ++q){
                        double T_face = 0.;                                   // We interpolate temperature at quadrature point q (we know temperature at face nodes).
                        for(unsigned int C=0; C<dofs_per_elem; C++){          // T at q = sum over C where C goes over all face nodes (C = 1 ... 6 ). In the sum we have dealii basis function "shape_value"
                                                                              // corresponding to C node at quadrature point q times T at C node)
                            T_face += D_trans[local_dof_indices[C]]*fe_face_values.shape_value(C,q);
                        }
                        for (unsigned int A=0; A<dofs_per_elem; A++){         // Photon radiation from the star surface
                            Flocal[A] -= sigma*pow(TiTe(T_face*exp(-phi(outer_radius))),4)*exp(2*phi(outer_radius))*fe_face_values.shape_value(A,q)*fe_face_values.JxW(q);
                        }
                    }
                }   
            }
        }


          constraints.distribute_local_to_global (Flocal,
                                                  local_dof_indices,
                                                  F);

          constraints.distribute_local_to_global (Mlocal,
                                                  local_dof_indices,
                                                  M);

          constraints.distribute_local_to_global (Klocal,
                                                  local_dof_indices,
                                                  K);
    }}

     M.compress (VectorOperation::add);
     K.compress (VectorOperation::add);
     F.compress (VectorOperation::add);
}


template <int dim>
void LaplaceProblem<dim>::solve ()
  { 

    double delta_t = 1;                                                       // initial time step
    double t_step = 0;       

    unsigned int trigger = 0;
    unsigned int N_proc_write = 0;

    unsigned int time_counter = 0;                                            // Just some parameters that help to handle the computation and output processes
    unsigned int snap_shot_counter = 0;
    
    if(Utilities::MPI::this_mpi_process(mpi_communicator) == N_proc_write){
        FILE* curve = std::fopen(CoolingFile, "w" );                                                                      
        std::fprintf(curve, "%s \n", "time [yr], Temperature at the boundary [K], Redshifted temperature at the surface [K]");
        std::fclose(curve);
    }

    LA::MPI::Vector completely_distributed_solution (locally_owned_dofs, mpi_communicator);
    std::vector<double> snapshot = log_space(2, t_max, N_output);
    const unsigned int dofs_per_elem = fe.dofs_per_face;                     
    std::vector<unsigned int> local_dof_indices (dofs_per_elem); 

    SolverControl solver_control (1e-12);

    dealii::TrilinosWrappers::SolverDirect solver(solver_control);

    apply_initial_conditions(); 

    while(t_step < t_max){
        
        assemble_system ();

        t_step += delta_t;     
        
        computing_timer.enter_subsection ("Solve"); 
                                            
        D_tilde = D_trans;
        D_tilde.add(delta_t*(1-alpha),V_trans);
        system_matrix.copy_from(M);
        system_matrix.add(alpha*delta_t,K);
        M.vmult(RHS,D_tilde);
        RHS.add(alpha*delta_t,F);
        
        solver.solve (system_matrix, completely_distributed_solution, RHS);
        
        constraints.distribute (completely_distributed_solution);

        D_trans = completely_distributed_solution;

        V_trans = D_trans;                                                    
        V_trans.add(-1.,D_tilde);
        V_trans /= alpha*delta_t;

        computing_timer.leave_subsection();

        if(t_step > snapshot[snap_shot_counter]){ 
            computing_timer.enter_subsection ("Output"); 
            output_results(snap_shot_counter);                           
            snap_shot_counter ++;
            pcout << "time = " << t_step << " sec" << std::endl;
            pcout << "time step = " << delta_t << " sec" << std::endl;
            pcout << "number of step = " << snap_shot_counter << " out of " << N_output << std::endl;
            
            typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), endc = dof_handler.end();

            if(Utilities::MPI::this_mpi_process(mpi_communicator) == N_proc_write){
                for (;elem!=endc; ++elem){
                    if (elem->is_locally_owned()){
                        for (unsigned int f=0; f < GeometryInfo<dim>::faces_per_cell; f++){
                            if(elem->face(f)->at_boundary()){ 
                                if(elem->face(f)->boundary_id()==1){
                                    elem->face(f)->get_dof_indices (local_dof_indices);
                                    FILE* curve = std::fopen(CoolingFile, "a+" );                                           
                                    std::fprintf(curve, "%14.7e %14.7e  %14.7e\n", t_step/yrtosec, D_trans[local_dof_indices[0]]*exp(-phi(outer_radius)), TiTe(D_trans[local_dof_indices[0]]*
                                    exp(-phi(outer_radius)))*redshift);
                                    std::fclose(curve);

                                    std::cout << "T_surface = " << TiTe(D_trans[local_dof_indices[0]]*exp(-phi(outer_radius)))*redshift << " K\n" << std::endl;
                               
                                    trigger = 1;
                                    break;
                                }
                            }
                        }
                    }
                
                    if (trigger == 1) {trigger=0; break;}
                }
            }
            computing_timer.leave_subsection(); 
        }

        if(time_counter<time_points.size()){                                                                                   
            if(t_step/yrtosec > time_points[time_counter]){
                delta_t = time_steps[time_counter];
                time_counter ++;
            }
        }
    }

    pcout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
    pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}


template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (D_trans, "u");

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches ();

    const std::string filename = ("output/solution-" +
                                  Utilities::int_to_string (cycle, 3) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));
    std::ofstream output ((filename + ".vtu").c_str());
    data_out.write_vtu (output);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          filenames.push_back ("solution-" +
                               Utilities::int_to_string (cycle, 3) +
                               "." +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output (("output/solution-" +
                                      Utilities::int_to_string (cycle, 3) +
                                      ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
    }
}


#endif //MAIN_FEM_HEATEQUATION_H


