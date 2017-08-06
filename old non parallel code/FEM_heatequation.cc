//
// Created by maryhallow on 06.12.15.
//

#include "FEM_heatequation.h"

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (double Alpha)
        :
        fe (FE_Q<dim>(order), 1),                                             // 1 means that we work with a scalar field (temperature)
        dof_handler (triangulation),
        quadrature_formula(quadRule),
        face_quadrature_formula(quadRule)
{
    alpha = Alpha;                                                            // parameter that defines the Euler method used for time discretization (e.g. if 1 then we deal with backwards Euler method)

    nodal_solution_names.push_back("D");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(){ // setting up the mesh

    double  x_center = 0.,
            y_center = 0.,
            z_center = 0.;

    Point<dim,double> center(x_center,y_center,z_center);
    GridGenerator::quarter_hyper_shell (triangulation, center, inner_radius, outer_radius, 0, true);
    static const SphericalManifold<3> manifold_description(center);
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold (0, manifold_description);
    triangulation.refine_global(global_refinement);
    GridTools::transform(grid_transform(), triangulation);                // we deform the quater hypershell manifold
   
}

// Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

    // Let deal.II organize degrees of freedom
    dof_handler.distribute_dofs (fe);

    MappingQ1<dim,dim> mapping;
    std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
    nodeLocation.reinit(dof_handler.n_dofs(),dim);
    DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
    for(unsigned int i=0; i<dof_coords.size(); i++){
        for(unsigned int j=0; j<dim; j++){
            nodeLocation[i][j] = dof_coords[i][j];
        }
    }

    // We Define the size of the global matrices and vectors
    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    sparsity_pattern.compress();
    K.reinit (sparsity_pattern);
    M.reinit (sparsity_pattern);
    system_matrix.reinit (sparsity_pattern);
    D_trans.reinit(dof_handler.n_dofs());
    V_trans.reinit(dof_handler.n_dofs());
    RHS.reinit(dof_handler.n_dofs());
    F.reinit(dof_handler.n_dofs());

    // Just some notes...
    std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

// In this function we form elemental vectors and matrices and assemble to the global vector (F) and matrices (K) and (M)
template <int dim>
void FEM<dim>::assemble_system(){

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
    unsigned int 	   num_quad_pts = quadrature_formula.size();              // Total number of quad points in the element

    FullMatrix<double> Mlocal (dofs_per_elem, dofs_per_elem);                 // Mass matrix
    FullMatrix<double> Klocal (dofs_per_elem, dofs_per_elem);                 // Stiffness matrix
    Vector<double>     Flocal (dofs_per_elem);                                // Force vector

    std::vector<unsigned int> local_dof_indices (dofs_per_elem);              // This relates local dof numbering to global dof numbering

    typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), endc = dof_handler.end();
    for (;elem!=endc; ++elem){

        fe_values.reinit(elem);
        elem->get_dof_indices (local_dof_indices);

        Mlocal = 0.;
        Klocal = 0.;
        Flocal = 0.;

        for(unsigned int q=0; q<num_quad_pts; q++){                           // We build our M and K matrices and F vector

            double D_trans_at_q = 0.;                                         // We interpolate temperature at quadrature point q using dealii basis functions

            for(unsigned int C=0; C<dofs_per_elem; C++){
                D_trans_at_q += D_trans[local_dof_indices[C]]*fe_values.shape_value(C,q);
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
            if(elem->face(f)->center()[0]*elem->face(f)->center()[1]*elem->face(f)->center()[2]>0 &&  
               sqrt(pow(elem->face(f)->center()[0],2)+pow(elem->face(f)->center()[1],2)+pow(elem->face(f)->center()[2],2)) > inner_radius &&
               elem->face(f)->at_boundary()){ // We check whether we are at the boundary or not. If we are not, nothing happens.

                for (unsigned int q=0; q<num_face_quad_pts; ++q){
                    double T_face = 0.;                                       // We interpolate temperature at quadrature point q (we know temperature at face nodes).
                    for(unsigned int C=0; C<dofs_per_elem; C++){              // T at q = sum over C where C goes over all face nodes (C = 1 ... 6 ). In the sum we have dealii basis function "shape_value"
                                                                              // corresponding to C node at quadrature point q times T at C node)
                        T_face += D_trans[local_dof_indices[C]]*fe_face_values.shape_value(C,q);
                    }
                    for (unsigned int A=0; A<dofs_per_elem; A++){             // Photon radiation from the star surface
                        Flocal[A] -= sigma*pow(TiTe(T_face*exp(-phi(outer_radius))),4)*exp(2*phi(outer_radius))*fe_face_values.shape_value(A,q)*fe_face_values.JxW(q);
                    }
                }
            }
        }


        for (unsigned int i=0; i<dofs_per_elem; ++i){                         // We already know local M,K and F. Every local M, K and F correspond to a particular element.
            F[local_dof_indices[i]] += Flocal[i];                             // From local M, K  matrices and F vector we construct global M, K matrices and global F vector.
            for (unsigned int j=0; j<dofs_per_elem; ++j){                     // This process has a name "assembly".
                K.add(local_dof_indices[i],local_dof_indices[j],Klocal[i][j]);
                M.add(local_dof_indices[i],local_dof_indices[j],Mlocal[i][j]);
            }
        }
    }
}


// Applying initial conditions
template <int dim>
void FEM<dim>::apply_initial_conditions(){

    const unsigned int totalNodes = dof_handler.n_dofs();
    for(unsigned int i=0; i<totalNodes; i++){
       D_trans[i] = T_0;                                                      // This initial condition is used for the real problem.
                                                                              // (Redshifted temperature is constant(r) at t = 0)
    }
}

template <int dim>
void FEM<dim>::solve_trans(){                                                 // The solver

    double delta_t = 0.1;                                                     // initial time step
    double t_step = 0;                                                        // initial time
    const unsigned int totalNodes = dof_handler.n_dofs();                     // Total number of nodes
    Vector<double>     D_tilde(totalNodes);

    apply_initial_conditions();                                               // applying initial conditions at t_step = 0

    unsigned int time_counter = 0;                                            // Just some parameters that help to handle the computation and output processes
    unsigned int snap_shot_counter = 0;

    FILE* curve = std::fopen(CoolingFile, "w" );                                                                      
    std::fprintf(curve, "%s \n", "time [yr], Temperature at the boundary [K], Redshifted temperature at the surface [K]");
    std::fclose(curve);


    while(t_step < t_max){                                                    // Loop over time

        t_step += delta_t;                                                    // updating time

        assemble_system();

        // Solving a matrix equation d_{n+1} = (M + alpha*delta_t*K)^{-1}*(alpha*delta_t*F_{n+1} + M*d_tilde_{n+1}),
        // where d_tilde_{n+1} = d_{n} + delta_t*(1-alpha)*v_{n}

        D_tilde = D_trans;
        D_tilde.add(delta_t*(1-alpha),V_trans);
        system_matrix.copy_from(M);
        system_matrix.add(alpha*delta_t,K);

        M.vmult(RHS,D_tilde);
        RHS.add(alpha*delta_t,F);

        SparseDirectUMFPACK  A;
        A.initialize(system_matrix);
        A.vmult (D_trans, RHS);

        V_trans = D_trans;                                                    // v_{n+1} = (d_{n+1} - d_tilde_{n+1})/(alpha*delta_t)
        V_trans.add(-1.,D_tilde);
        V_trans /= alpha*delta_t;

        if (D_trans[0]<T_min){
            std::cout << "Simulation has been stopped, because T star is less than T min.\n" << std::endl;
            break;
        }

        if(t_step > snapshot[snap_shot_counter]){                             // printing out some notes at snapshot time
            snap_shot_counter ++;
            std::cout << "time = " << t_step << " sec" << std::endl;
            std::cout << "time step = " << delta_t << " sec" << std::endl;
            std::cout << "number of step = " << snap_shot_counter << " out of " << N_output << std::endl;
                                                                              // writing out cooling data to the dat file

            for(unsigned int globalNode=0; globalNode<totalNodes; globalNode++) {
                if (sqrt(pow(nodeLocation[globalNode][0],2)+pow(nodeLocation[globalNode][1],2)+pow(nodeLocation[globalNode][2],2)) == outer_radius) {
                    FILE* curve = std::fopen(CoolingFile, "a+" );
                                                                              // writing out (time in years, (internal) redshifted temperature on the surface in K, (external) redshifted temperature on the surface in K)
                    std::fprintf(curve, "%14.7e %14.7e  %14.7e\n", t_step/yrtosec, D_trans[globalNode]*exp(-phi(outer_radius)), TiTe(D_trans[globalNode]*exp(-phi(outer_radius)))*redshift);
                    std::cout << "T_surface = " << TiTe(D_trans[globalNode]*exp(-phi(outer_radius)))*redshift << " K\n" << std::endl;
                    std::fclose(curve);
                    break;
                }
            }
        }

        if(time_counter<time_points.size()){                              // changing time step and writing out 3d data to the vtk file and
                                                                          // temperature profile data to the dat file
            if(t_step/yrtosec > time_points[time_counter]){
                delta_t = time_steps[time_counter];
                output_trans_results(time_counter);
                output_dat_file(time_counter);
                time_counter ++;
            }
        }
    }
}

template <int dim>
void FEM<dim>::output_trans_results (unsigned int index){

    //Write results to VTK file
    char filename[100];
    snprintf(filename, 100, "./output/solution_%d.vtk", index);
    std::ofstream output1 (filename);
    DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

    //Add nodal DOF data
    data_out.add_data_vector (D_trans, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
    data_out.build_patches (); data_out.write_vtk (output1); output1.close();
}


template <int dim>
void FEM<dim>::output_dat_file(unsigned int index){

    const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

    char buffer[100];

    snprintf(buffer, 100, "./output/temperature_profile_%i.dat", index);

    FILE* temperature_profile;
    temperature_profile = std::fopen(buffer, "w" );

    for(unsigned  int i=0; i<totalNodes; i++)
        std::fprintf(temperature_profile, "%14.7e  %14.7e\n", r(i), D_trans[i]);

    std::fclose(temperature_profile);
}
