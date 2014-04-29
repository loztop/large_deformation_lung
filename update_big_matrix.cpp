 // std::cout<< "Updating big coupled stiffness matrix "  <<std::endl;

 		clock_t begin_assemble_fem=clock();
		assemble_solid(equation_systems,"Newton-update");
		clock_t end_assemble_fem=clock();
	
  //Setup up some variables
  //const MeshBase& mesh = es.get_mesh();
  //const Real dt    = es.parameters.get<Real>("dt"); 

  TransientLinearImplicitSystem&  system = equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
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
	AP.close();
	
	 
  //Create the matrix for the tree
  Mat T;
  MatCreate(PETSC_COMM_WORLD,&T);
  MatSetSizes(T,PETSC_DECIDE,PETSC_DECIDE,tree.number_nodes+tree.number_edges,tree.number_nodes+tree.number_edges);
  MatSetFromOptions(T);
  MatSetUp(T);		
  MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY);
  //Need to allocate size properly
	clock_t begin_assemble_tree=clock();
  T=tree.make_tree_matrix( );
	clock_t end_assemble_tree=clock();

	PetscMatrix<Number> ATP(T) ;
  Real size_t=ATP.m();
  
  
  //Add the coupling between the poromedium and the tree
 // std::cout<< "Create big matrix " <<std::endl;
  #include "create_include.cpp"

 // std::cout<< "Adding coupling terms "<<std::endl;
	clock_t begin_assemble_coupling=clock();
  //#include "add_couplingv4.cpp"
	clock_t end_assemble_coupling=clock();
 

	
	
	//Sort out rhs	- put into big vector
	Vec x;
	VecCreate(PETSC_COMM_WORLD,&x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x,PETSC_DECIDE,size_fem);
  VecSetFromOptions(x);
	PetscVector<Number> xp(x) ;
  //NumericVector< Number > &solution_in= *(newton_update.solution);
	xp.add(1,solution_in);

	
	//NumericVector< Number > &rhs_in = *(newton_update.rhs);
	Vec r;
	VecCreate(PETSC_COMM_WORLD,&r);
  PetscObjectSetName((PetscObject) r, "rhs");
  VecSetSizes(r,PETSC_DECIDE,size_fem);
  VecSetFromOptions(r);
	PetscVector<Number> rp(r) ;
	rp.add(1,rhs_in);


	Vec big_r;
	VecCreate(PETSC_COMM_WORLD,&big_r);
  PetscObjectSetName((PetscObject) big_r, "big_rhs");
  VecSetSizes(big_r,PETSC_DECIDE,size_fem+size_tree);
  VecSetFromOptions(big_r);
	
	// add system.rhs into big_rp 
	for (int i=0; i<size_fem; i++) {
			VecSetValue(big_r, i , rp(i) ,INSERT_VALUES); 
	}
	
	PetscVector<Number> big_rp(big_r) ;
	
	Vec big_x;
	VecCreate(PETSC_COMM_WORLD,&big_x);
  PetscObjectSetName((PetscObject) big_x, "big_Solution");
  VecSetSizes(big_x,PETSC_DECIDE,size_fem+size_tree);
  VecSetFromOptions(big_x);
	PetscVector<Number> big_xp(big_x) ;
	///End of assembling rhs

