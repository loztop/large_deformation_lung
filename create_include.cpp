
  A=AP.mat();
  MPI_Comm        comm;
  comm = MPI_COMM_SELF;
  PetscInt        a_nrow,*a_rows,*a_cols;
  PetscTruth      done;

  PetscFunctionBegin;
  MatGetRowIJ(A,0,PETSC_FALSE,PETSC_FALSE,&a_nrow,&a_rows,&a_cols,&done);

  PetscScalar *a_array;
  MatGetArray(A, &a_array);
   
  PetscInt       t_nrow,*t_rows,*t_cols;
 
  PetscFunctionBegin;
  MatGetRowIJ(T,0,PETSC_FALSE,PETSC_FALSE,&t_nrow,&t_rows,&t_cols,&done);
  													
  PetscScalar *t_array;
  MatGetArray(T, &t_array);
	
  PetscInt   *big_rows_a;
  big_rows_a=(PetscInt *)malloc((a_nrow+1)*sizeof(PetscInt));
	
  PetscInt   *big_rows_t;
  big_rows_t=(PetscInt *)malloc((t_nrow+a_nrow+1)*sizeof(PetscInt));
		
  PetscInt      ext=a_nrow;
  PetscInt      big_nrows=t_nrow+a_nrow;


  for ( int i = 0; i < t_nrow+1 + a_nrow; i++ )
  {
		 if(i<a_nrow+1){
			big_rows_t[ i ] = a_rows[i];
		 }else{
			big_rows_t[ i ] = t_rows[i-a_nrow]+big_rows_t[ a_nrow ];			
		 }
  }

  PetscInt num_vals_a=a_rows[a_nrow];
  PetscInt num_vals_t=t_rows[t_nrow];
  PetscInt   *big_cols_t;
  big_cols_t=(PetscInt *)malloc((num_vals_t+num_vals_a)*sizeof(PetscInt));
	
  for (int i = 0; i < num_vals_a; i++ )
  {
      big_cols_t[ i ] = a_cols[ i ];  
  }
    
  for ( int i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
  {
      big_cols_t[ i ] = t_cols[ i -num_vals_a] + a_nrow;
  }
   
  PetscScalar   *b_array;
  b_array=(PetscScalar *)malloc((num_vals_t+num_vals_a)*sizeof(PetscScalar));
	
   for ( int i = 0; i < num_vals_a; i++ )
   {
      b_array[ i ] = a_array[ i ];
   }
   
  for ( int i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
  {
      b_array[ i ] = t_array[ i-num_vals_a ];
  }
   
  
  MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows_t,big_cols_t,b_array,&big_A);

  MatSetFromOptions(big_A);
  MatSetUp(big_A);	

  MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);
	

  //Free arrays ???
  //free(b_array);
  //free(big_cols_t);
