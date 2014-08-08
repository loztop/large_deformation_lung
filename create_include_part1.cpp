   //A=AP.mat();
   
#if !mats
MPI_Comm comm;
comm = MPI_COMM_SELF;
PetscInt a_nrow,*a_rows,*a_cols;
PetscTruth done;
PetscFunctionBegin;
#endif
  

MatGetRowIJ(AP.mat(),0,PETSC_FALSE,PETSC_FALSE,&a_nrow,&a_rows,&a_cols,&done);

#if !mats
PetscScalar *a_array;
#endif
  
MatGetArray(AP.mat(), &a_array);
  
#if !mat
PetscInt t_nrow,*t_rows,*t_cols;
t_rows=(PetscInt *)malloc((4*size_tree)*sizeof(PetscInt));
t_cols=(PetscInt *)malloc((4*size_tree)*sizeof(PetscInt));
PetscScalar *t_array;
t_array=(PetscScalar *)malloc((4*size_tree)*sizeof(PetscScalar));
#endif 

MatGetRowIJ(T,0,PETSC_FALSE,PETSC_FALSE,&t_nrow,&t_rows,&t_cols,&done);
MatGetArray(T, &t_array);

 #if !mat
  PetscInt *big_rows_a;
  big_rows_a=(PetscInt *)malloc((a_nrow+1)*sizeof(PetscInt));

  PetscInt *big_rows_t;
  big_rows_t=(PetscInt *)malloc((t_nrow+a_nrow+1)*sizeof(PetscInt));

  PetscInt ext=a_nrow;
  PetscInt big_nrows=t_nrow+a_nrow;
#endif
  
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
  
#if !mats
  PetscInt *big_cols_t;
  big_cols_t=(PetscInt *)malloc((num_vals_t+num_vals_a)*sizeof(PetscInt));
#endif
  
  for (int i = 0; i < num_vals_a; i++ )
  {
      big_cols_t[ i ] = a_cols[ i ];
  }
    
  for ( int i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
  {
      big_cols_t[ i ] = t_cols[ i -num_vals_a] + a_nrow;
  }
 
#if !mats
  PetscScalar *b_array;
  b_array=(PetscScalar *)malloc((num_vals_t+num_vals_a)*sizeof(PetscScalar));
#endif
  
   for ( int i = 0; i < num_vals_a; i++ )
   {
      b_array[ i ] = a_array[ i ];
   }
   
  for ( int i = num_vals_a; i < num_vals_a+num_vals_t; i++ )
  {
      b_array[ i ] = t_array[ i-num_vals_a ];
  }
   