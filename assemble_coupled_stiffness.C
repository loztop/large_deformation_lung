#include "poro.h"

PetscMatrix<Number> assemble_coupled_stiffness (EquationSystems& es,Tree& tree, Mesh& mesh)
{
 
	
  std::cout<< "Assembling coupled big matrix " <<std::endl;
	
	std::cout<<"Assembling FEM matrix "<<std::endl;
  assemble_solid(es,"Newton-update"); 
	TransientLinearImplicitSystem & system = 
    es.add_system<TransientLinearImplicitSystem> ("Newton-update");
 
	

	
	//Make a new and bigger matrix
  SparseMatrix< Number >  &matrix_in=*(system.matrix);
  NumericVector< Number > &solution_in= *(system.solution);
	NumericVector< Number > &rhs_in = *(system.rhs);
	Number size_mat = system.rhs->size();

  Mat            A;            
  PetscMPIInt    size;
  PetscBool      nonzeroguess = PETSC_FALSE;

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size_mat,size_mat);
  MatSetFromOptions(A);
  MatSetUp(A);
	
  PetscMatrix<Number> AP(A) ;
  AP.add(1,matrix_in);
  
	
	
	
	
	
  //Create the matrix for the tree
  Mat AT;
  AT=tree.make_tree_matrix( );

  PetscMatrix<Number> ATP(AT) ;
  std::cout<< "Tree matrix rows " << ATP.m()<<std::endl;
  Real size_t=ATP.m();
	
  //Create the big matrix to hold both systems
  Mat big_A;            
  MatCreate(PETSC_COMM_WORLD,&big_A);
  MatSetSizes(big_A,PETSC_DECIDE,PETSC_DECIDE,size_mat+size_t,size_mat+size_t);
	//MatSetSizes(big_A,PETSC_DECIDE,PETSC_DECIDE,size_mat,size_mat);
	
  MatSetFromOptions(big_A);
  MatSetUp(big_A);		
  MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);
  std::cout<< "FEM matrix rows " << AP.m()<<std::endl;
	
  std::cout<< "Copying element matrix into big matrix " <<std::endl;
  //Add the element matrix
  for (int i=0; i<AP.m(); i++) {
		//std::cout<< "FEMi " << i << std::endl;
		for (int j=0; j<AP.n(); j++) {
		//	if(abs(AP(i,j))>0){
			   MatSetValue(big_A, i, j  , AP(i,j) ,ADD_VALUES); 
		//	}
		}
	}
	
	


	//Add the network tree
	for (int i=AP.m(); i<AP.m()+size_t; i++) {
		for (int j=AP.n(); j<AP.n()+size_t; j++) {				
		//	if(abs(ATP(i-AP.m(),j-AP.n()))>0){
			   MatSetValue(big_A, i, j  , ATP(i-AP.m(),j-AP.n()) ,ADD_VALUES); 
		//	}
		}
	}

	 
	 //big_A=create_big_matrx(A,T)
	

	
	/*
	//Now couple the FEm system with the tree
	//This is p_nu not s_p !!
	std::vector<unsigned int> dof_indices_p;
	const unsigned int p_var = system.variable_number ("p_nu");
	const DofMap & dof_map = system.get_dof_map();

	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

	for ( ; el != end_el; ++el)
	{    	
	  const Elem* elem = *el;
		Point elem_pos=elem->centroid();
		
		//Find the closest distal tree_node to q_point[qp]
		int closest_j=0;
		Real dist_j=0;
		Real closest_dist=99999999;
		int closest_end_j=0;
		int closest_end_edge_j=0;
		int closest_end_node=0;

		for (int j=0; j < tree.number_nodes; j++) {	
		  
		  if(tree.nodes_type(j)==0){
				
				Real dist_j=pow(elem_pos(0)-  tree.nodes(j)(0),2)+pow(elem_pos(1)-  tree.nodes(j)(1),2)+pow(elem_pos(2)-  tree.nodes(j)(2),2);
											
				if(dist_j<closest_dist)
				{
					closest_dist=dist_j; 
					
					//Relies on assumption that lower node is +1 of its parent edge.
					closest_end_edge_j=j-1;
				}
			}
		}
				
		dof_map.dof_indices (elem, dof_indices_p, p_var);

		//Approx factor of distal airway resistance - 32
		Real distal_fac=1;
		
		//Need to check and xperiment with this
		Real dt  = es.parameters.get<Real>("dt");

		
		//Should also have some new contributions to the residual !!
		
		//Mass conservation coupling term (source)  //need to add 1/r
		MatSetValue(big_A, dof_indices_p[0], AP.m()+tree.number_nodes+closest_end_edge_j, - dt*1*(elem->volume()) ,ADD_VALUES); 
	
		//Distal flow rate coupling term (source)  //need to add 1/r
				
		//Q_j - (1/rd) P_{j} + ( (1/r_d)*(1/omega_j)*p_poro this in intro.C )
		MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , dof_indices_p[0]  , distal_fac*(1.0/tree.edges_resistance(closest_end_edge_j)) ,ADD_VALUES); 
		
		//We might get away with ignoring the residual contribution of the tree since this should always be zero.
		
		
	}

	*/
	
  
   //  Update the matrix big_A to reflect the changes made above
 // MatAssemblyBegin(big_A, MAT_FINAL_ASSEMBLY);
 // MatAssemblyEnd(big_A, MAT_FINAL_ASSEMBLY);

	//Print matrix
//	MatView(big_A,PETSC_VIEWER_STDOUT_WORLD);
  
  
  PetscMatrix<Number> big_AP(big_A) ;
  
  
  //	const DofMap & dof_map = system .get_dof_map();

	//std::cout<< dof_map.get_send_list () << std::endl;
	
   // big_AP.add_matrix(matrix_in,dof_map.get_send_list () );

	
  big_AP.close();
  

			
  return big_AP;
  
}