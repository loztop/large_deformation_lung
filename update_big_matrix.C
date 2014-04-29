#include "poro.h"

void update_big_matrix(Mat& big_A, EquationSystems& es, Tree& tree) {

  std::cout<< "Updating big coupled stiffness matrix "  <<std::endl;

  //Setup up some variables
  const MeshBase& mesh = es.get_mesh();
  const Real dt    = es.parameters.get<Real>("dt"); 

  TransientLinearImplicitSystem&  system = es.get_system<TransientLinearImplicitSystem>("Newton-update");
  //system.assemble_before_solve=false;
  system.update();
  
  SparseMatrix< Number > &matrix_in=*(system.matrix);
  NumericVector< Number > &solution_in= *(system.solution);
  NumericVector< Number > &rhs_in = *(system.rhs);
  Number size_mat = system.rhs->size();
	
  //Create the matrix for fem
  Mat            A;            
  PetscMPIInt    size;
  PetscBool      nonzeroguess = PETSC_FALSE;

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size_mat,size_mat);
  MatSetFromOptions(A);
  MatSetUp(A);
	
  PetscMatrix<Number>& AP = *libmesh_cast_ptr<PetscMatrix<Number>*>( system.matrix);

  //Create the matrix for the tree
  Mat T;
  MatCreate(PETSC_COMM_WORLD,&T);
  MatSetSizes(T,PETSC_DECIDE,PETSC_DECIDE,tree.number_nodes+tree.number_edges,tree.number_nodes+tree.number_edges);
  MatSetFromOptions(T);
  MatSetUp(T);		
  MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY);
  //Need to allocate size properly
  T=tree.make_tree_matrix( );
  PetscMatrix<Number> ATP(T) ;
  Real size_t=ATP.m();
  
  
  //Add the coupling between the poromedium and the tree
  //std::cout<< "Create big matrix " <<std::endl;
  #include "create_include.cpp"

  //std::cout<< "Adding coupling terms "<<std::endl;
  #include "add_couplingv4.cpp"
  
  
}


