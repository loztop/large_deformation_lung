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
#include <libmesh/tecplot_io.h>
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_matrix.h"


#include "poro.h"
#include "tree.h"


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

LibMeshInit init (argc, argv);



Mesh mesh(3);
EquationSystems equation_systems (mesh);
read_parameters(equation_systems,argc,argv);
std::string result_file_name (equation_systems.parameters.get<std::string>("result_file_name"));

unsigned int n_timesteps = equation_systems.parameters.get<Real>("n_timesteps");
unsigned int N_eles=equation_systems.parameters.get<Real>("N_eles");

 
Real time     = 0;
Real end_time     = equation_systems.parameters.get<Real>("end_time");
const unsigned int n_nonlinear_steps = 10;
const Real nonlinear_tolerance       = 1.e-3;
const Real initial_linear_solver_tol = 1.e-18;

if(!equation_systems.parameters.get<std::string>("problem").compare("cube")){

	 MeshTools::Generation::build_cube (mesh,
									N_eles, N_eles, N_eles,
									0., 1.,
									0., 1.,
									0., 1.,
									HEX8);
  mesh.all_first_order();

}	

if(!equation_systems.parameters.get<std::string>("problem").compare("lung")){
  
	//std::string mesh_file_name (equation_systems.parameters.get<std::string>("mesh_input"));
	//GmshIO(mesh).read(mesh_file_name);
	
	 MeshTools::Generation::build_cube (mesh,
                                       3, 3,4,
                                       0.0, 1.0,
                                       0.0, 1.0,
																				0.0, 2.0,
                                       TET4);
	
 }

mesh.prepare_for_use();
mesh.print_info();

setup_equationsystem(equation_systems);
equation_systems.init ();
equation_systems.print_info();

Real dt = end_time/n_timesteps;
equation_systems.parameters.set<Real> ("dt")   = dt;

ExodusII_IO exo= ExodusII_IO(mesh);

#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif

TransientLinearImplicitSystem&  newton_update =   equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
TransientLinearImplicitSystem&  reference =   equation_systems.get_system<TransientLinearImplicitSystem>("Reference-Configuration");
TransientLinearImplicitSystem&  last_non_linear_soln = equation_systems.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");

// Loop over all nodes and copy the location from the current system to
// the auxiliary system.

const MeshBase::const_node_iterator nd_end =
      equation_systems.get_mesh().local_nodes_end();
for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
  const Node *node = *nd; 

//Copy initial mesh into reference.solution
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

newton_update.assemble_before_solve=false;
newton_update.update();
	
// Load in the reference mesh as an initial guess (a neater way to do this ?)
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

//Create tree
Tree tree;
tree.read_tree(equation_systems);

//Write out 0th timestep
equation_systems.parameters.set<unsigned int>("step") = 0; 
tree.write_tree(equation_systems);
std::stringstream file_name;
file_name << equation_systems.parameters.get<std::string>("result_file_name");
file_name << std::setw(2) << std::setfill('0') << 0;
file_name << ".e-s.";
file_name << std::setw(3) << std::setfill('0') << 0;
exo.write_timestep(file_name.str(), equation_systems,1,0); 
std::cout<<"Write initial conditions "<< file_name.str() <<std::endl;
exo.write_element_data(equation_systems);

//Count dofs 
int size_fem=last_non_linear_soln.n_dofs();
int size_tree=tree.number_nodes+tree.number_edges;
std::cout<<"FEM dofs "<< size_fem <<std::endl;
std::cout<<"Tree dofs "<< size_tree <<std::endl;


size_tree=0;

//Create the big matrix to hold both systems
Mat big_A;            
MatCreate(PETSC_COMM_WORLD,&big_A);
MatSetSizes(big_A,PETSC_DECIDE,PETSC_DECIDE,size_fem+size_tree,size_fem+size_tree);
MatSetFromOptions(big_A);
MatSetUp(big_A);		
MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);


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
		
    equation_systems.parameters.set<Real> ("non_lin_step") = l;
    std::cout<<"\nNon-linear iteration " << l << std::endl;

    change_in_newton_update->zero();
    change_in_newton_update->add(*newton_update.solution);

    //Prepare the newton update system for it's linear solve
    *newton_update.old_local_solution = *newton_update.current_local_solution;
    
    newton_update.current_local_solution->zero();  
    newton_update.solution->zero();  
    newton_update.update();

    clock_t begin_big_assemble=clock();
		#include "update_big_matrix.cpp"
		PetscMatrix<Number> big_AP(big_A) ;
		big_AP.close();
		clock_t end_big_assemble=clock();
		std::cout<<"Assembly,"<< " total: " << double(diffclock(end_big_assemble,begin_big_assemble)) <<  " fem: " << double(diffclock(end_assemble_fem,begin_assemble_fem)) << " coupling: " << double(diffclock(end_assemble_coupling,begin_assemble_coupling)) <<  " ms"<<std::endl; 

	
		//Finally solve the poroelastic system
		clock_t begin_solid_solve=clock();
  
		//equation_systems.get_system("Newton-update").solve();
	
		//std::cout<< "Solving system "<<std::endl;
		PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
		//petsc_linear_solver->solve(big_AP, big_xp , big_rp, tolerance, m_its);
		petsc_linear_solver->solve( AP, big_xp , big_rp, 1.e-15,4);
		clock_t end_solid_solve=clock();
	
		
		//Copy FEM solution back to system solution vector, update mesh and tree positions
		#include "update_solution.cpp"

		//Get the norm of the residual to check for convergence.
		TransientLinearImplicitSystem & newton_update =equation_systems.get_system<TransientLinearImplicitSystem> ("Newton-update");
		
		//Real solid_residual=newton_update.rhs->l2_norm ();
    Real solid_residual=big_xp.l2_norm ();
		
    /*****///Convergence and Norm Computing Stuff///***********/
    change_in_newton_update->add (-1., *newton_update.solution);
		change_in_newton_update->close();
		Real norm_delta = change_in_newton_update->l2_norm();
    change_in_newton_update->add (-1., *newton_update.solution);
    change_in_newton_update->close();
    norm_delta = change_in_newton_update->l2_norm();
    const unsigned int n_linear_iterations = newton_update.n_linear_iterations();
    const Real final_linear_residual = newton_update.final_linear_residual();   
    std::cout << "Poro: Linear conv at step: "
                    << n_linear_iterations
                    << ", resid: "
                    << final_linear_residual 
                    << ", time: " << double(diffclock(end_solid_solve,begin_solid_solve)) << " ms"
                    <<std::endl;             
    std::cout   << "Poro Nonlinear convergence: ||u - u_old|| = "
                    << norm_delta << ", poro res:  " << solid_residual
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
