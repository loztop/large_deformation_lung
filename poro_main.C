//iHat*s_u + jHat*s_v + kHat*s_w


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"


#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"

#include <petscksp.h>
#include <petsctime.h>

#include "poro.h"


#define SOLVER_NAME "mumps"
//#define SOLVER_NAME "superlu"
#define PC_TYPE PCLU
#define PETSC_MUMPS 1


//iHat*s_u + jHat*s_v + kHat*s_w

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Begin the main program.
int main (int argc, char** argv)
{
//
  
  
  static char help[] = "Test file for the PCFactorSetShiftType()\n";
   Mat                A,T;
   MPI_Comm           comm;
   PetscScalar        v;
    PetscInt           i,j;
   PetscErrorCode     ierr;

    PetscInitialize(&argc,&argv,0,help);
    comm = MPI_COMM_SELF;

   MatCreateSeqAIJ(comm,4,4,4,0,&A);
  i=3; j=3; v=2;
  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr); 
    i=2; j=3; v=12;
  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
    i=1; j=3; v=13;
  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
    i=2; j=2; v=23;
  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
	
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
  MatView(A,PETSC_VIEWER_STDOUT_WORLD);

   
  PetscInt       a_nrow,*a_rows,*a_cols;
  PetscBool      done;

  PetscFunctionBegin;
  ierr = MatGetRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&a_nrow,&a_rows,&a_cols,&done);CHKERRQ(ierr);

  PetscScalar *a_array;
  MatGetArray(A, &a_array);

	
	//Tree matrix
  MatCreateSeqAIJ(comm,3,3,3,0,&T);
  i=2; j=2; v=5;
  ierr = MatSetValues(T,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr); 
    i=1; j=2; v=55;
  ierr = MatSetValues(T,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
	
  ierr = MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
  MatView(T,PETSC_VIEWER_STDOUT_WORLD);

   
  PetscInt       t_nrow,*t_rows,*t_cols;
 
  PetscFunctionBegin;
  ierr = MatGetRowIJ(T,0,PETSC_FALSE,PETSC_FALSE,&t_nrow,&t_rows,&t_cols,&done);CHKERRQ(ierr);
  
  PetscScalar *t_array;
  MatGetArray(T, &t_array);
	

  
	//Build bigger matrix	
	//PetscInt      big_rows_a[a_nrow+1];
	PetscInt      big_rows_t[t_nrow+a_nrow+1];
	PetscInt      ext=a_nrow;

	PetscInt      big_nrows=t_nrow+ext;


		//MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows,big_cols,0,&Big);

	
	for ( i = 0; i < t_nrow+1 + a_nrow; i++ )
	{
		 if(i<a_nrow+1){
			big_rows_t[ i ] = a_rows[i];
		 }else{
			big_rows_t[ i ] = t_rows[i-ext]+big_rows_t[ a_nrow -1];			
		 }
	}
	
	
	PetscInt num_vals_a=big_rows_t[a_nrow];
	PetscInt num_vals_t=big_rows_t[t_nrow+a_nrow]+num_vals_a;

	PetscInt   big_cols_t[num_vals_t+num_vals_a];
	

	for ( i = 0; i < num_vals_a; i++ )
   {
      big_cols_t[ i ] = a_cols[ i ];
   }
   
	for ( i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
   {
      big_cols_t[ i ] = t_cols[ i ] + a_nrow;
   }
   
   
   //Create b_array
   	PetscScalar b_array[num_vals_t+num_vals_a];
   for ( i = 0; i < num_vals_a+1; i++ )
   {
      b_array[ i ] = a_array[ i ];
   }
   
	for ( i = num_vals_a+1; i < num_vals_a+num_vals_t+1; i++ )
   {
      b_array[ i ] = t_array[ i-num_vals_a-1 ];
   }
   
  
  /*
  //Print vectors
  
  for ( i = 0; i < num_vals_t+num_vals_t; i++ )
   {
			std::cout<< "b_array["<< i <<" ] " << b_array[ i ] << std::endl;
   }
  
  
	for ( i = 0; i < num_vals_a; i++ )
   {
			std::cout<< "big_cols_t["<<i<<"] " << big_cols_t[ i ] << std::endl;
   }
   for ( i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
   {
			std::cout<< "big_cols_t["<<i<<"] " << big_cols_t[ i ] << std::endl;
   }
   
   
  	for ( i = 0; i < a_nrow; i++ )
   {
			std::cout<< "big_rows_t["<<i<<"] " << big_rows_t[ i ] << std::endl;
   }
   for ( i = a_nrow; i < a_nrow+t_nrow; i++ )
   {
			std::cout<< "big_rows_t["<<i<<"] " << big_rows_t[ i ] << std::endl;
   }
  */
  
  //Insert the tree into the big matrix
  
 
	
	/*
	Mat   Big;

	//MatCreateSeqAIJ(comm,big_nrows,big_nrows,big_nrows,0,&Big);
  //MatSetValues(Big,big_nrows,big_rows_t,big_nrows,big_cols_t,t_array,INSERT_VALUES);	
 
  //MatSetValues(Big,t_nrow,big_rows_t,t_nrow,big_cols_t,t_array,INSERT_VALUES);	
 
 
 MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows_t,big_cols_t,b_array,&Big);
	
	
	std::cout << "Big  "<< std::endl;
	ierr = MatAssemblyBegin(Big,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Big,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
	MatView(Big,PETSC_VIEWER_STDOUT_WORLD);
	
	*/
	
	
	PetscMatrix<Number> Big_P(create_big_matrix(A,T)) ;

					

	
	
	
	

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //
LibMeshInit init (argc, argv);

Real time     = 0;
Real end_time     = 3;
const unsigned int n_nonlinear_steps = 15;
const Real nonlinear_tolerance       = 1.e-2;
const Real inita_rowsl_linear_solver_tol = 1.e-18;

Mesh mesh(3);
EquationSystems equation_systems (mesh);
read_parameters(equation_systems,argc,argv);
std::string result_file_name (equation_systems.parameters.get<std::string>("result_file_name"));

unsigned int n_timesteps = equation_systems.parameters.get<Real>("n_timesteps");
unsigned int N_eles=equation_systems.parameters.get<Real>("N_eles");


/*
std::string mesh_file_name ("MESHES/whole_lung_246.msh");
   
std::cout << mesh_file_name << std::endl;
GmshIO(mesh).read(mesh_file_name);
mesh.prepare_for_use();
mesh.print_info();
*/


MeshTools::Generation::build_cube (mesh,
                                   N_eles, N_eles, N_eles,
                                   0., 1.,
                                   0., 1.,
                                   0., 1.,
                                   HEX8);

mesh.all_first_order();
mesh.prepare_for_use();
//mesh.print_info();

setup_equationsystem(equation_systems);
equation_systems.init ();
//equation_systems.print_info();

Real dt = end_time/n_timesteps;
equation_systems.parameters.set<Real> ("dt")   = dt;

ExodusII_IO exo= ExodusII_IO(mesh);

TransientLinearImplicitSystem&  newton_update =   equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
TransientLinearImplicitSystem&  reference =   equation_systems.get_system<TransientLinearImplicitSystem>("Reference-Configuration");
TransientLinearImplicitSystem&  last_non_linear_soln = equation_systems.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");

// Loop over all nodes and copy the location from the current system to
// the auxila_rowsry system.

const MeshBase::const_node_iterator nd_end =
      equation_systems.get_mesh().local_nodes_end();
for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
  const Node *node = *nd; 

//Copy inita_rowsl mesh into reference.solution
  for (unsigned int d = 0; d < 3; ++d) {
    unsigned int dest_dof = node->dof_number(reference.number(), d, 0);
    Real value = (*node)(d);
    reference.current_local_solution->set(dest_dof, value);
    reference.solution->set(dest_dof, value);
  }
}

reference.solution->close();
reference.current_local_solution->close();
reference.update();   

AutoPtr<NumericVector<Number> > change_in_newton_update (newton_update.solution->clone());
AutoPtr<NumericVector<Number> > newton_solution (newton_update.solution->clone());
newton_solution->zero();

// Load in the reference mesh as an inita_rowsl guess (a neater way to do this ?)
last_non_linear_soln.solution->zero(); 
last_non_linear_soln.solution->add(1.,(*reference.solution));
last_non_linear_soln.solution->close();
last_non_linear_soln.current_local_solution->zero(); 
last_non_linear_soln.current_local_solution->add(1.,(*reference.current_local_solution));
last_non_linear_soln.current_local_solution->close();
last_non_linear_soln.old_local_solution->zero(); 
last_non_linear_soln.old_local_solution->add(1.,(*reference.current_local_solution));
last_non_linear_soln.old_local_solution->close();
last_non_linear_soln.update(); 



for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
{
  time += dt;
  equation_systems.parameters.set<Real> ("time") = time;
  double progress = (t_step+0.000000001) / (n_timesteps+0.000000001);
  equation_systems.parameters.set<Real>("progress") = progress;
  equation_systems.parameters.set<unsigned int>("step") = t_step; 
  std::cout << "\n\n*** Solving time step " << t_step << ", time = " << time <<  ", progress = " << progress << " ***" << std::endl;

  *last_non_linear_soln.old_local_solution = *last_non_linear_soln.current_local_solution;

  // Now we begin the nonlinear loop
  for (unsigned int l=0; l<n_nonlinear_steps; ++l)
  {
		
    clock_t begin_inside_nonlin=clock();

    equation_systems.parameters.set<Real> ("non_lin_step") = l;
    std::cout<<"\n Non-linear iteration " << l << std::endl;

    change_in_newton_update->zero();
    change_in_newton_update->add(*newton_update.solution);

    //Prepare the newton update system for it's linear solve
     *newton_update.old_local_solution = *newton_update.current_local_solution;
    
    newton_update.current_local_solution->zero();  
    newton_update.solution->zero();  
    newton_update.update();

	Tree tree;
	tree.read_tree();

	//Assemble the combined stiffness matrix
	PetscMatrix<Number> big_AP =assemble_coupled_stiffness (equation_systems,tree,mesh);

	//Assemble the combined rhs vector
	PetscVector<Number> big_rp =assemble_coupled_rhs (equation_systems,tree,mesh);
	
	//std::cout<< big_rp  <<std::endl;
	//std::cout<< *(newton_update.rhs) <<std::endl;

	//Create the solution vector
	PetscVector<Number> big_xp(big_rp) ;
	big_xp.zero();
	
	//Set up solver, proba can get rid of lower half
	PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
	PC pc = petsc_linear_solver->pc();
	int ierr = PCSetType(pc, PC_TYPE);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);
	ierr = PCFactorSetMatSolverPackage(pc,SOLVER_NAME);
	CHKERRABORT(libMesh::COMM_WORLD,ierr);
	
	
	//Finally solve the coupled poroelastic system	
	clock_t begin_solid_solve=clock();

	
	
	petsc_linear_solver->solve(big_AP, big_xp , big_rp, 1.e-15,25);

	
	equation_systems.get_system("Newton-update").solve();  

	
	clock_t end_solid_solve=clock();
	
	std::cout<<"Time to solve: "<<diffclock(end_solid_solve,begin_solid_solve) <<std::endl;
	
	//Output residual history  - Put this in convergence section
	std::vector< double > hist;
	petsc_linear_solver->get_residual_history(hist);
	std::cout<< "hist.size() " <<hist.size()<<std::endl;
	for (int i=0; i < hist.size(); i++) {
		std::cout<< "residual " << hist[i]<<std::endl;
	}
	
	/*
	//Copy FEM solution back to system solution vector
	for (int i=0; i<newton_update.rhs->size(); i++) {
		newton_update.solution->set(i,big_xp(i));
		newton_update.current_local_solution->set(i,big_xp(i));
	}
	newton_update.solution->close();
	newton_update.current_local_solution->close();
	newton_update.update();   
	*/
	
	
	equation_systems.reinit();

  //update the final solution
  // xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
  //Apply a full Newton-step
  Real K=1; //Newton step size

  last_non_linear_soln.solution->add(-1*K,*newton_update.solution);
  last_non_linear_soln.solution->close();
  last_non_linear_soln.current_local_solution->add(-1,*newton_update.current_local_solution);
   last_non_linear_soln.current_local_solution->close();
   last_non_linear_soln.update();

	MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 
    for ( ; el != end_el; ++el)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
    }
    
    equation_systems.update();
    equation_systems.allgather();
    
      /*****///Convergence and Norm Computing Stuff///***********/
      //Get the norm of the residual to check for convergence.
	  TransientLinearImplicitSystem & newton_update =equation_systems.get_system<TransientLinearImplicitSystem> ("Newton-update");
	  Real solid_residual=newton_update.rhs->l2_norm ();
      change_in_newton_update->add (-1., *newton_update.solution);
			change_in_newton_update->close();
			Real norm_delta = change_in_newton_update->l2_norm();
      change_in_newton_update->add (-1., *newton_update.solution);
      change_in_newton_update->close();
      norm_delta = change_in_newton_update->l2_norm();
      const unsigned int n_linear_iterations = newton_update.n_linear_iterations();
      const Real final_linear_residual = newton_update.final_linear_residual();   
      std::cout << "System: Linear conv at step: "
                    << n_linear_iterations
                    << ", resid: "
                    << final_linear_residual 
                    << ", time: " << double(diffclock(end_solid_solve,begin_solid_solve)) << " ms"
                    <<std::endl;             
      std::cout   << "System Nonlinear convergence: ||u - u_old|| = "
                    << norm_delta
                    << std::endl;
      if ((norm_delta < nonlinear_tolerance)&&(solid_residual < nonlinear_tolerance) ){
      std::cout << "Nonlinear solver converged after "<< l+1 <<" steps."<<std::endl;
      break;

			
  } 
 } // end nonlinear loop

 last_non_linear_soln.update();      
 newton_update.update(); 
 // Write out every nth timestep to file.
 const unsigned int write_interval = 1;

 
 #include "write_variable_results_and_mesh.cpp"
 
 
} // end timestep loop.
}
