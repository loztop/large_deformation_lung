//Create the big matrix to hold both systems
Mat big_A;            
MatCreate(PETSC_COMM_WORLD,&big_A);
MatSetSizes(big_A,PETSC_DECIDE,PETSC_DECIDE,size_fem+size_tree,size_fem+size_tree);
MatSetFromOptions(big_A);
MatSetUp(big_A);		
MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);


std::cout<<"Assemble at start to set up sparsity pattern "<<std::endl;
equation_systems.parameters.set<Real>("progress") = 0;
equation_systems.parameters.set<Real> ("time") = 0;
newton_update.current_local_solution->zero();
newton_update.solution->zero();
newton_update.update();
clock_t begin_big_assemble_b4=clock();
assemble_solid(equation_systems,"Newton-update");
clock_t end_big_assemble_b4=clock();
newton_update.current_local_solution->zero();
newton_update.solution->zero();
newton_update.solution->close();
newton_update.current_local_solution->close();
newton_update.update();   

std::cout<<"AssemblyB4, " << double(diffclock(end_big_assemble_b4,begin_big_assemble_b4)) <<std::endl;


PetscViewer    viewer;


TransientLinearImplicitSystem&  system = equation_systems.get_system<TransientLinearImplicitSystem>("Newton-update");
system.update();


SparseMatrix< Number > &matrix_in=*(system.matrix);
NumericVector< Number > &solution_in= *(system.solution);
NumericVector< Number > &rhs_in = *(system.rhs);
Number size_mat = system.rhs->size();



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
 
clock_t begin_assemble_tree=clock();
T=tree.make_tree_matrix( );
clock_t end_assemble_tree=clock();
  
MatSetFromOptions(T);
MatSetUp(T);	
MatAssemblyBegin(T,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(T,MAT_FINAL_ASSEMBLY);
  
PetscMatrix<Number> ATP(T) ;
Real size_t=ATP.m();

 

//Sort out rhs - put into big vector
Vec x;
VecCreate(PETSC_COMM_WORLD,&x);
PetscObjectSetName((PetscObject) x, "Solution");
VecSetSizes(x,PETSC_DECIDE,size_fem);
VecSetFromOptions(x);
PetscVector<Number> xp(x) ;
 
Vec r;
VecCreate(PETSC_COMM_WORLD,&r);
PetscObjectSetName((PetscObject) r, "rhs");
VecSetSizes(r,PETSC_DECIDE,size_fem);
VecSetFromOptions(r);
PetscVector<Number> rp(r) ;

Vec big_r;
VecCreate(PETSC_COMM_WORLD,&big_r);
PetscObjectSetName((PetscObject) big_r, "big_rhs");
VecSetSizes(big_r,PETSC_DECIDE,size_fem+size_tree);
VecSetFromOptions(big_r);
 

////#include "create_include_part1.cpp"

 
MPI_Comm comm;
comm = MPI_COMM_SELF;
PetscInt a_nrow,*a_rows,*a_cols;
PetscTruth done;
PetscFunctionBegin;

PetscScalar *a_array;

PetscInt t_nrow,*t_rows,*t_cols;
t_rows=(PetscInt *)malloc((4*size_tree)*sizeof(PetscInt));
t_cols=(PetscInt *)malloc((4*size_tree)*sizeof(PetscInt));
PetscScalar *t_array;
t_array=(PetscScalar *)malloc((4*size_tree)*sizeof(PetscScalar));

PetscInt *big_rows_a;
big_rows_a=(PetscInt *)malloc((a_nrow+1)*sizeof(PetscInt));

PetscInt *big_rows_t;
big_rows_t=(PetscInt *)malloc((t_nrow+a_nrow+1)*sizeof(PetscInt));

PetscInt ext=a_nrow;
PetscInt big_nrows=t_nrow+a_nrow;
 

 ////#include "add_coupling_fast.cpp"

PetscInt *end_j_f;
end_j_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

PetscReal *omega_end_j_f;
   omega_end_j_f=(PetscReal *)malloc((tree.number_nodes)*sizeof(PetscReal));
   
PetscInt *end_j_zero_f;
end_j_zero_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

PetscInt *found_end_j_f;
found_end_j_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

PetscInt *already_set_f;
already_set_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

//Create some arrays to hold inserts
PetscInt *rows_coupling_f;
rows_coupling_f=(PetscInt *)malloc((a_nrow+t_nrow)*sizeof(PetscInt));

#if !fvec
PetscInt *cols_coupling_f;
PetscMalloc((3*a_nrow+t_nrow)*sizeof(PetscInt),&cols_coupling_f); 
#endif

PetscReal *vals_coupling_f;
vals_coupling_f=(PetscReal *)malloc((3*a_nrow+t_nrow)*sizeof(PetscReal));

PetscReal *rhs_new;
rhs_new=(PetscReal *)malloc((a_nrow+t_nrow)*sizeof(PetscReal));



//////#include "create_include_add_coupling.cpp"

PetscInt *sort_idx;
sort_idx=(PetscInt *)malloc((el_counter*4)*sizeof(PetscInt));

PetscInt *cols_coupling_test_cpy;
cols_coupling_test_cpy=(PetscInt *)malloc((el_counter*4)*sizeof(PetscInt));

PetscReal *vals_coupling_test_cpy;
vals_coupling_test_cpy=(PetscReal *)malloc((el_counter*4)*sizeof(PetscReal));

PetscInt *idx_insert;
idx_insert=(PetscInt *)malloc((el_counter*4)*sizeof(PetscInt));

PetscInt *rows_new;
rows_new=(PetscInt *)malloc((a_nrow+t_nrow+1)*sizeof(PetscInt));

//#include "update_bg_matrix_again.cpp"

Vec big_x;
VecCreate(PETSC_COMM_WORLD,&big_x);
PetscObjectSetName((PetscObject) big_x, "big_Solution");
VecSetSizes(big_x,PETSC_DECIDE,size_fem+size_tree);
VecSetFromOptions(big_x);
PetscVector<Number> big_xp(big_x) ;


