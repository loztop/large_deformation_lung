 int rowmax= a_nrow;

  //MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,big_rows_t,big_cols_t,b_array,&big_A);

  
  MatCreateSeqAIJWithArrays(comm,big_nrows,big_nrows,rows_new,cols_new,vals_new,&big_A);


  MatSetFromOptions(big_A);
  MatSetUp(big_A);	

  MatAssemblyBegin(big_A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(big_A,MAT_FINAL_ASSEMBLY);

