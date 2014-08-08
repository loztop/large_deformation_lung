// std::cout<< "Updating big coupled stiffness matrix " <<std::endl;

#if !mats
PetscViewer viewer;
#endif

#if !mats
Mat big_A;            
MatCreate(PETSC_COMM_WORLD,&big_A);
MatSetSizes(big_A,PETSC_DECIDE,PETSC_DECIDE,size_fem+size_tree,size_fem+size_tree);
MatSetFromOptions(big_A);
MatSetUp(big_A);		
MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);
#endif

clock_t begin_assemble_fem=clock();
assemble_solid(equation_systems,"Newton-update");
clock_t end_assemble_fem=clock();

clock_t begin_assemble_coupling_fast=clock();	

#if !mats
  TransientLinearImplicitSystem& system = equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
  system.update();
#endif
  
  
#if !mats
  SparseMatrix< Number > &matrix_in=*(system.matrix);
  NumericVector< Number > &solution_in= *(system.solution);
  NumericVector< Number > &rhs_in = *(system.rhs);
  Number size_mat = system.rhs->size();
#endif
  
  //Create the matrix for fem
/*
  Mat A;
  PetscMPIInt size;
  PetscBool nonzeroguess = PETSC_FALSE;

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MatCreate(PETSC_COMM_WORLD,&A);
  MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,size_mat,size_mat);
  MatSetFromOptions(A);
  MatSetUp(A);
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  */
  
 
#if !mats

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
 
  MatSetFromOptions(T);
  MatSetUp(T);	
  MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY);
  
  PetscMatrix<Number> ATP(T) ;
  Real size_t=ATP.m();
#endif

#if !mats
  //Sort out rhs - put into big vector
  Vec x;
  VecCreate(PETSC_COMM_WORLD,&x);
  PetscObjectSetName((PetscObject) x, "Solution");
  VecSetSizes(x,PETSC_DECIDE,size_fem);
  VecSetFromOptions(x);
  PetscVector<Number> xp(x) ;
#endif

  
xp.zero();
xp.add(1,solution_in);

#if !mats
Vec r;
VecCreate(PETSC_COMM_WORLD,&r);
PetscObjectSetName((PetscObject) r, "rhs");
VecSetSizes(r,PETSC_DECIDE,size_fem);
VecSetFromOptions(r);
PetscVector<Number> rp(r) ;
#endif

rp.zero();
rp.add(1,rhs_in);

#if !mats
Vec big_r;
VecCreate(PETSC_COMM_WORLD,&big_r);
PetscObjectSetName((PetscObject) big_r, "big_rhs");
VecSetSizes(big_r,PETSC_DECIDE,size_fem+size_tree);
VecSetFromOptions(big_r);
#endif

// add system.rhs into big_rp
for (int i=0; i<size_fem; i++) {
VecSetValue(big_r, i , rp(i) ,INSERT_VALUES);
}

#include "create_include_part1.cpp"

#if !fvec2
#include "add_coupling_fast.cpp"
#include "create_include_add_coupling.cpp"
#endif

#if fvec2
#include "add_coupling_fast_fvec.cpp"
#include "create_include_add_coupling_fvec.cpp"
#endif

#include "create_include_part2.cpp"

clock_t end_assemble_coupling_fast=clock();

PetscMatrix<Number> big_AP(big_A) ;
big_AP.close();

PetscVector<Number> big_rp(big_r) ;
big_rp.close();

#if !mats
Vec big_x;
VecCreate(PETSC_COMM_WORLD,&big_x);
  PetscObjectSetName((PetscObject) big_x, "big_Solution");
  VecSetSizes(big_x,PETSC_DECIDE,size_fem+size_tree);
  VecSetFromOptions(big_x);
PetscVector<Number> big_xp(big_x) ;
#endif




#if fvec2 
#include "create_destroyfvec.cpp"
#endif












































