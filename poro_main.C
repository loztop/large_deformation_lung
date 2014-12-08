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
 
Real time = 0;
Real end_time = equation_systems.parameters.get<Real>("end_time");
const unsigned int n_nonlinear_steps = 6;
const Real nonlinear_tolerance = 1.e-3;
const Real initial_linear_solver_tol = 1.e-18;

if(!equation_systems.parameters.get<std::string>("problem").compare("cube")){

MeshTools::Generation::build_cube (
mesh,
                                       1,1,2,
                                       0.0, 1.0,
                                       0.0, 1.0,
																			 0.0, 1.0,
                                       HEX8);
                                      

}	


if((!equation_systems.parameters.get<std::string>("problem").compare("lung")) || (!equation_systems.parameters.get<std::string>("problem").compare("cylinder"))){
  
std::string mesh_file_name (equation_systems.parameters.get<std::string>("mesh_input"));
GmshIO(mesh).read(mesh_file_name);

}
 
//Create tree
Tree tree;
tree.read_tree(equation_systems);

mesh.prepare_for_use();

setup_equationsystem(equation_systems);
TransientLinearImplicitSystem& newton_update = equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
TransientLinearImplicitSystem& reference = equation_systems.get_system<TransientLinearImplicitSystem>("Reference-Configuration");
TransientLinearImplicitSystem& last_non_linear_soln = equation_systems.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
TransientLinearImplicitSystem& postvars = equation_systems.get_system<TransientLinearImplicitSystem>("postvars");




equation_systems.init ();
equation_systems.print_info();
 
/*
const MeshBase::const_node_iterator nd_end_tlc =
      equation_systems.get_mesh().local_nodes_end();
for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end_tlc; ++nd) {
  const Node *node = *nd;

//Copy initial mesh into reference.solution
  for (unsigned int d = 0; d < 3; ++d) {
		//Save the TLC reference
    unsigned int dest_dof = node->dof_number(reference.number(), d+4, 0);
    Real value = (*node)(d);
    reference.current_local_solution->set(dest_dof, value);
    reference.solution->set(dest_dof, value);

	}
}
 */

//Move mesh and tree from inspiration to registered expiration
  
 const Point Aexp = equation_systems.parameters.get<Point>("A");
const Point bexp = equation_systems.parameters.get<Point>("b");

Real facexp=1.4;
//Update the position of the airway tree
std::cout<<"Moving tree to expiratory FRC " <<std::endl;
    for (double j=0; j <tree.number_nodes ; j++) {
			for (unsigned int d = 0; d < 3; ++d) {
				tree.nodes(j)(d)=tree.nodes(j)(d)-facexp*( tree.nodes(j)(d)*Aexp(d)+bexp(d) );
			}
				
				
				//Move mesh a bit more to ecncompass tree on;y do this for N48
				//	if((!equation_systems.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048r_fine2881.msh")) || (!equation_systems.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048_node5036.msh"))){
					  	if( (!equation_systems.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048_node5036.msh"))){
						
						//tree.nodes(j)(0)= tree.nodes(j)(0)*0.98+2; //-paper don't know why ??
						//tree.nodes(j)(1)= tree.nodes(j)(1)*0.98+4; //-paper don't know why ??
						//tree.nodes(j)(2)= tree.nodes(j)(2)*0.98+2; //-paper don't know why ??
											
						tree.nodes_deformed(j)(0)=tree.nodes(j)(0);
						tree.nodes_deformed(j)(1)=tree.nodes(j)(1);
						tree.nodes_deformed(j)(2)=tree.nodes(j)(2);		
						
							
					}	
					
					
						tree.nodes_frc(j)(0)=tree.nodes(j)(0);
						tree.nodes_frc(j)(1)=tree.nodes(j)(1);
						tree.nodes_frc(j)(2)=tree.nodes(j)(2);	
					
     }


tree.fix_tree(equation_systems);
//tree.shrink_rad(equation_systems);

       
//Update the mesh position
std::cout<<"Moving mesh to expiratory FRC" <<std::endl;
Mesh::node_iterator it_node = mesh.nodes_begin();
const Mesh::node_iterator it_last_node = mesh.nodes_end();
for ( ; it_node != it_last_node ; ++it_node)
    {
      Node* node = *it_node;
      for (unsigned int d = 0; d < 3; ++d) {
           (*node)(d)= (*node)(d) - facexp*((*node)(d)*Aexp(d)+bexp(d)) ;
		   
      }
   }



//Calculate volume of frc
MeshBase::const_element_iterator el_jac_frc = mesh.active_local_elements_begin();
const MeshBase::const_element_iterator end_el_jac_frc = mesh.active_local_elements_end();
Real total_volume_frc=0;
for ( ; el_jac_frc != end_el_jac_frc; ++el_jac_frc)
{
const Elem* elem = *el_jac_frc;
Real elem_vol=elem->volume();
total_volume_frc=total_volume_frc+elem_vol;
}
std::cout<<"total_volume frc "<< total_volume_frc << std::endl;




const MeshBase::const_node_iterator nd_end_frc =
      equation_systems.get_mesh().local_nodes_end();
for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end_frc; ++nd) {
	const Node *node = *nd;

	  for (unsigned int d = 0; d < 3; ++d) {
		//Save the FRC reference mesh position
    unsigned int dest_dof = node->dof_number(reference.number(), d+4, 0);
    Real value = (*node)(d);
    reference.current_local_solution->set(dest_dof, value);
    reference.solution->set(dest_dof, value);

	}
}
//////////////////


std::cout<<"Moving tree to reference " <<std::endl;
    for (double j=0; j <tree.number_nodes ; j++) {
			for (unsigned int d = 0; d < 3; ++d) {
				tree.nodes(j)(d)=tree.nodes(j)(d)-( tree.nodes(j)(d)*REFSCALE );
				tree.nodes_deformed(j)(d)=tree.nodes(j)(d);
			}
				
					
     }

 
std::cout<<"Moving mesh to reference " <<std::endl;
Real ref_state=1;
Real NT_ref_final=REFNT;
Real NT_ref=0;
Mesh::node_iterator it_node_ref = mesh.nodes_begin();
const Mesh::node_iterator it_last_node_ref = mesh.nodes_end();
for ( ; it_node_ref != it_last_node_ref ; ++it_node_ref)
    {
      Node* node = *it_node_ref;
      for (unsigned int d = 0; d < 3; ++d) {
		 (*node)(d)= (*node)(d) - (*node)(d)*REFSCALE ;
      }
   }

   
//Copy REF to reference system
const MeshBase::const_node_iterator nd_end =
      equation_systems.get_mesh().local_nodes_end();
for (MeshBase::const_node_iterator nd = equation_systems.get_mesh().local_nodes_begin();
      nd != nd_end; ++nd) {
  const Node *node = *nd;

//Copy ref mesh into reference.solution REF MESH
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


Real dt = end_time/n_timesteps;
equation_systems.parameters.set<Real> ("dt") = dt;

ExodusII_IO exo= ExodusII_IO(mesh);

#if WRITE_TEC
TecplotIO tec= TecplotIO(equation_systems.get_mesh());
#endif


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


MeshBase::const_element_iterator el_jac = mesh.active_local_elements_begin();
const MeshBase::const_element_iterator end_el_jac = mesh.active_local_elements_end();
const DofMap & dof_map = reference .get_dof_map();
std::vector<unsigned int> dof_indices_p;
const unsigned int p_var = reference.variable_number ("vol_ref");
Real total_volume_ref=0;

int el_counter=0;
for ( ; el_jac != end_el_jac; ++el_jac)
{
const Elem* elem = *el_jac;
Real elem_vol=elem->volume();
total_volume_ref=total_volume_ref+elem_vol;
dof_map.dof_indices (elem, dof_indices_p, 3);
reference.current_local_solution->set(dof_indices_p[0], elem_vol);
reference.solution->set(dof_indices_p[0], elem_vol);
el_counter=el_counter+1;
}
std::cout<<"total_volume ref "<< total_volume_ref << std::endl;
std::cout<<"el_counter "<< el_counter << std::endl;


/*
//Save FRC to postvars /// CAN REMOVE ALL THIS ONLY NEEDED FOR AIRWAY RES
/////////////////////
equation_systems.parameters.set<Real>("ref_state")=0;
tree.add_constriction(equation_systems);
equation_systems.parameters.set<Real>("ref_state")=1;

tree.calculate_total_resistance(equation_systems);
//Do some post-processing (calculate stress etc)
///To get airway tree resistances
	postvars.assemble_before_solve=false;
	postvars.matrix->zero();
	postvars.rhs->zero();
	assemble_postvars(equation_systems,"postvars",tree);
	postvars.update();
	postvars.solve();
///
	///To get airway tree resistances at 0th time step !!
std::stringstream file_name_tec;
  file_name_tec <<result_file_name << "_"<< 0<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
  std::cout<<"Wrote "<< file_name_tec.str() <<std::endl;
	/////////
	*/
	


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

//Integrate the total outflow
Real total_outflow=0;

#if mats
#include "create_arrays_mats.cpp"
equation_systems.update();
equation_systems.allgather();
equation_systems.reinit();
#endif

PetscLinearSolver<Number>* petsc_linear_solver =dynamic_cast<PetscLinearSolver<Number>*>(newton_update.get_linear_solver());
	
for (unsigned int t_step=1; t_step<=n_timesteps; ++t_step)
{
  

  if(NT_ref>=NT_ref_final && ref_state>0){
	std::cout<<"Changing to normal breathing"<<std::endl;
	ref_state=0;
	time=0;
	equation_systems.parameters.set<Real> ("ref_state") = 0;
	
	tree.add_constriction(equation_systems);
	T=tree.make_tree_matrix( );

  }
	  
  if(ref_state>0 ){
		std::cout<<"Going from reference to FRC"<<std::endl;

    time += REFDT;

	NT_ref=NT_ref+1;
	t_step=t_step-1;
		  	
  }else{
  
  //dt = end_time/n_timesteps;
  double progress = (t_step+0.000000001) / (n_timesteps+0.000000001);
  equation_systems.parameters.set<Real>("progress") = progress;
  equation_systems.parameters.set<unsigned int>("step") = t_step;
  time = progress*end_time;

  std::cout << "\n\n*** Solving time step " << t_step << ", time = " << time << ", progress = " << progress << " ***" << std::endl;
  }
  
  equation_systems.parameters.set<Real> ("time") = time;

	
  *last_non_linear_soln.old_local_solution = *last_non_linear_soln.current_local_solution;

  // Now we begin the nonlinear loop
  for (unsigned int l=0; l<n_nonlinear_steps; ++l)
  {	
	clock_t begin_big_nonlin=clock();
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
	clock_t end_big_assemble=clock();

	std::cout<<"Assembly,"<< " total: " << double(diffclock(end_big_assemble,begin_big_assemble)) << " fem: " << double(diffclock(end_assemble_fem,begin_assemble_fem)) << " tree: " << double(diffclock(end_assemble_tree,begin_assemble_tree)) << " ms"<< " coupling: " << double(diffclock(end_assemble_coupling_fast,begin_assemble_coupling_fast)) << " ms"<<std::endl;

	//Finally solve the poroelastic system
	clock_t begin_solid_solve=clock();	
	//std::cout<< "Solving system "<<std::endl;
	
	//petsc_linear_solver->solve(big_AP, big_xp , big_rp, tolerance, m_its);
	petsc_linear_solver->solve( big_AP, big_xp , big_rp, 1.e-15,4);
	clock_t end_solid_solve=clock();

	//Copy FEM solution back to system solution vector, update mesh and tree positions
	clock_t begin_big_update=clock();
	#include "update_solution.cpp"
	clock_t end_big_update=clock();

	//std::cout<< "-------- Residual and convergence info ----------"<<std::endl;
	#include "residual_info.cpp"

	

	//end of nonlinear step computation
	clock_t end_big_nonlin=clock();

	//free memory
	#if mats
	//PetscFree(b_array); 
	//PetscFree(big_cols_t); 
	//PetscFree(vals_new); 
	//PetscFree(cols_new); 
	big_AP.clear();
	big_rp.clear();
	#endif
 
	std::cout<<"Total: " << double(diffclock(end_big_nonlin,begin_big_nonlin)) << " assembly: " << double(diffclock(end_big_assemble,begin_big_assemble)) << " update: " << double(diffclock(end_big_update,begin_big_update)) << " solve: " << double(diffclock(end_solid_solve,begin_solid_solve)) << " ms"<<std::endl;

	if ((norm_delta/l2_soln < nonlinear_tolerance)&&(solid_residual/l2_soln < nonlinear_tolerance) ){
	  std::cout << "Nonlinear solver converged after "<< l+1 <<" steps."<<std::endl;
	  break;
	}
	
  } // end nonlinear loop

 last_non_linear_soln.update();
 newton_update.update();
 // Write out every nth timestep to file.
 const unsigned int write_interval = 1;

 		if(!(ref_state>0)){
		total_outflow=total_outflow+tree.edges_flowrate(0)*dt;
		std::cout<<"total_outflow "<< total_outflow << std::endl;
		} 
		
 //update tree position
  clock_t begin_pos_update=clock();
    tree.update_positions(equation_systems);    
    clock_t end_pos_update=clock();
	std::cout<<"tree position_update: " << double(diffclock(end_pos_update,begin_pos_update)) <<  " ms"<<std::endl;
	
	//Do some post-processing (calculate stress etc)
	postvars.assemble_before_solve=false;
	postvars.matrix->zero();
	postvars.rhs->zero();
		std::cout<<"Assembling postvars !! " <<std::endl;

	assemble_postvars(equation_systems,"postvars",tree);
	postvars.update();
	postvars.solve();
	
 #include "write_variable_results_and_mesh.cpp"
 
} // end timestep loop.
} // end main function.

#include "assemble_postvars.cpp"


