#include "poro.h"

using namespace std;

void test(int a) {
  // std::cout << "ex "<< a << std::endl;
}


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
} 




Mat create_big_matrix(Mat A, Mat T){
	
  MPI_Comm           comm;
  comm = MPI_COMM_SELF;
  PetscInt       a_nrow,*a_rows,*a_cols;
  PetscBool      done;

  PetscFunctionBegin;
  MatGetRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&a_nrow,&a_rows,&a_cols,&done);

  PetscScalar *a_array;
  MatGetArray(A, &a_array);
   
  PetscInt       t_nrow,*t_rows,*t_cols;
 
  PetscFunctionBegin;
  MatGetRowIJ(T,0,PETSC_FALSE,PETSC_FALSE,&t_nrow,&t_rows,&t_cols,&done);
  
  PetscScalar *t_array;
  MatGetArray(T, &t_array);
	

	//Build bigger matrix	
	//PetscInt      big_rows_a[a_nrow+1];
	PetscInt      big_rows_t[t_nrow+a_nrow+1];
	PetscInt      ext=a_nrow;

	PetscInt      big_nrows=t_nrow+ext;


		//MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows,big_cols,0,&Big);

	
	for ( int i = 0; i < t_nrow+1 + a_nrow; i++ )
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
	

	for (int i = 0; i < num_vals_a; i++ )
   {
      big_cols_t[ i ] = a_cols[ i ];
   }
   
	for ( int i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
   {
      big_cols_t[ i ] = t_cols[ i ] + a_nrow;
   }
   
   
   //Create b_array
   	PetscScalar b_array[num_vals_t+num_vals_a];
   for ( int i = 0; i < num_vals_a+1; i++ )
   {
      b_array[ i ] = a_array[ i ];
   }
   
	for ( int i = num_vals_a+1; i < num_vals_a+num_vals_t+1; i++ )
   {
      b_array[ i ] = t_array[ i-num_vals_a-1 ];
   }
   
  
  //Insert the tree into the big matrix
	Mat   Big;
	MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows_t,big_cols_t,b_array,&Big);

	
	MatSetFromOptions(Big);
	MatSetUp(Big);	
	
	MatAssemblyBegin(Big,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(Big,MAT_FINAL_ASSEMBLY);
	
	
	
	
  
	MatView(Big,PETSC_VIEWER_STDOUT_WORLD);
	
	
	
	
	
  
	return Big;
} 
































void calculate_numeric_jacobian(EquationSystems& es, SparseMatrix< Number >& num_jac_matrix) {

	
	
			/*
			 * 
			 * Put this in poro_main.C
	//Calculate the numerical jack
	//	SparseMatrix< Number > &num_jac_matrix=*(last_non_linear_soln.matrix);
	//  num_jac_matrix.zero();
	//	calculate_numeric_jacobian(equation_systems,num_jac_matrix);
	//	num_jac_matrix.close();
	//
		//use the numerical jacobian
		newton_update.assemble_before_solve=false;
		newton_update.update();
		assemble_solid(equation_systems,"Newton-update");
		newton_update.matrix->zero();
		newton_update.matrix->add(1, num_jac_matrix);
		newton_update.matrix->close();
 			//newton_update.matrix->print(std::cout);
		*/
	
	
	TransientLinearImplicitSystem & newton_update =
  es.get_system<TransientLinearImplicitSystem> ("Newton-update");
	
	
	TransientLinearImplicitSystem&  last_non_linear_soln = es.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
	
	//probably should clone this
	
	
	int N_dofs=  newton_update.rhs->size();
	Real eps= 0.001;

	//Get the unperturbed res.
	assemble_solid(es,"Newton-update");
	AutoPtr<NumericVector<Number> > fx (newton_update.rhs->clone());
	newton_update.rhs->zero();
	newton_update.matrix->zero();
	newton_update.matrix->close();

	std::cout<<"N_dofs  "<< N_dofs << std::endl; 

	for (unsigned int i=0; i<N_dofs; ++i)
	{
		//std::cout<<"ith dof "<< i << std::endl; 

	
		//Change one of the res dofs by eps
		last_non_linear_soln.solution->add(i, eps);
		last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(i,eps);
		last_non_linear_soln.current_local_solution->close();
	
		//Calculate the new rhs (at the moment are also calculatating lhs !)
		assemble_solid(es,"Newton-update");
		newton_update.matrix->zero();
		newton_update.rhs->zero();
		newton_update.matrix->close();

		//Get the difference between the peturbed res and the original res
		AutoPtr<NumericVector<Number> > f_diff (newton_update.rhs->clone());
		f_diff->add (-1., *fx);

		AutoPtr<NumericVector<Number> > jac_col (newton_update.rhs->clone());
		jac_col->zero();
		jac_col->add (1./eps, *f_diff);
		jac_col->close();

		//Assemble into the numerical jacobian matrix
		for (unsigned int t=1; t<=N_dofs; ++t)
		{
			num_jac_matrix.set(i,t,jac_col->el(t));
		}

		//Zero out for the next test
		f_diff->close();
		f_diff->zero();

		newton_update.rhs->zero();
		//subtract again eps
		
		last_non_linear_soln.solution->add(i,-eps);
		last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(i,-eps);
		last_non_linear_soln.current_local_solution->close();
		
		
	}

 	num_jac_matrix.close();
 	num_jac_matrix.print(std::cout);
	
}
