
//First make a copy


	Mat            T_orig;            
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MatCreate(PETSC_COMM_WORLD,&T_orig);
  MatSetSizes(T_orig,PETSC_DECIDE,PETSC_DECIDE,tree.number_nodes+tree.number_edges,tree.number_nodes+tree.number_edges);
  MatSetFromOptions(T_orig);
  MatSetUp(T_orig);
  MatAssemblyBegin(T_orig,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(T_orig,MAT_FINAL_ASSEMBLY);
	
	//copy over
	T_orig=T;
	
	//Print original
	 	 
	MatView(T_orig,viewer);
	
	//Print out CSR format
   std::cout<< "t_rows[i]" <<std::endl;
	 for ( int i = 0; i < t_nrow; i++ )
  {
			std::cout<< t_rows[i] <<std::endl;
  }
  
   std::cout<< "t_cols[i]" <<std::endl;
  for ( int i = 0; i < num_vals_t; i++ )
  {
			std::cout<< t_cols[i] <<std::endl;
  }
  
  std::cout<< "t_array[i]" <<std::endl;
  for ( int i = 0; i < num_vals_t; i++ )
  {
			std::cout<< t_array[i] <<std::endl;
  }
   

   //Define coupling additions
   
   	int size_coup=5;
   	PetscInt   *rows_coupling_test;
	rows_coupling_test=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
	
	//PetscInt   *rows_coupling_test_copy;
	//rows_coupling_test_copy=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
	
	
	PetscInt   *sort_idx;
	sort_idx=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
	
	PetscInt   *cols_coupling_test;
	cols_coupling_test=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
	PetscInt   *cols_coupling_test_cpy;
	cols_coupling_test_cpy=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
	
	PetscReal   *vals_coupling_test;
	vals_coupling_test=(PetscReal *)malloc((size_coup)*sizeof(PetscReal));
	PetscReal   *vals_coupling_test_cpy;
	vals_coupling_test_cpy=(PetscReal *)malloc((size_coup)*sizeof(PetscReal));
	
	
	PetscInt   *idx_insert;
	idx_insert=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
  
	rows_coupling_test[0]=2;
	rows_coupling_test[1]=4;
	rows_coupling_test[2]=4;
	rows_coupling_test[3]=4;
	rows_coupling_test[4]=4;

	cols_coupling_test[0]=1;
	cols_coupling_test[1]=0;
	cols_coupling_test[2]=1;
	cols_coupling_test[3]=2;
	cols_coupling_test[4]=3;

	  vals_coupling_test[0]=33;
    vals_coupling_test[1]=99;
    vals_coupling_test[2]=777;
    vals_coupling_test[3]=555;
    vals_coupling_test[4]=2222;

	  sort_idx[0]=0;
    sort_idx[1]=1;
    sort_idx[2]=2;
	  sort_idx[3]=3;
	  sort_idx[4]=4;

	
	//First sort the arrays
	PetscSortIntWithArray(size_coup,rows_coupling_test,sort_idx);
	//make copies
//	vals_coupling_test_cpy=vals_coupling_test;
//	cols_coupling_test_cpy=cols_coupling_test;

	
	  for ( int i = 0; i < size_coup; i++ )
   {
	 vals_coupling_test_cpy[i]=vals_coupling_test[i];
	 cols_coupling_test_cpy[i]=cols_coupling_test[i];
  }
	
   for ( int i = 0; i < size_coup; i++ )
   {
	 cols_coupling_test[i]=cols_coupling_test_cpy[sort_idx[i]];
	 vals_coupling_test[i]=vals_coupling_test_cpy[sort_idx[i]];
  }
	
  //Create T_new using CSR format
 
	 
Mat T_new;            
MatCreate(PETSC_COMM_WORLD,&T_new);
MatSetSizes(T_new,PETSC_DECIDE,PETSC_DECIDE,t_nrow,t_nrow);
MatSetFromOptions(T_new);
MatSetUp(T_new);		
MatAssemblyBegin(T_new,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(T_new,MAT_FINAL_ASSEMBLY);


PetscInt   *rows_new;
rows_new=(PetscInt *)malloc((t_nrow)*sizeof(PetscInt));
for ( int i = 0; i < t_nrow; i++ )
{
  rows_new[i]=0;
}
  
PetscInt   *cols_new;
cols_new=(PetscInt *)malloc((num_vals_t+size_coup)*sizeof(PetscInt));
for ( int i = 0; i < num_vals_t+size_coup; i++ )
{
  cols_new[i]=0;
}

PetscReal   *vals_new;
vals_new=(PetscReal *)malloc((num_vals_t+size_coup)*sizeof(PetscReal));
for ( int i = 0; i < num_vals_t+size_coup; i++ )
{
  vals_new[i]=0;
}

			std::cout<< " t_nrow  "<<t_nrow<<std::endl;


//Create rows_new
int idx=0;

  for ( int i = 1; i < t_nrow+1; i++ )
  {

		//Count how many rows_coupling_test[idx+j]=i
		//then idx=idx+1+j at endl
	
		if(rows_coupling_test[idx]==i-1 && idx<size_coup){
			
			int stop=0;
			int j=1;

				while ( stop<1 )
				{
					
					if(rows_coupling_test[idx+j]==rows_coupling_test[idx]){
						j=j+1;
					}else{
						stop=1;
					}
			
				}
				
				idx=idx+j;
		
		}
	  
	rows_new[i]=t_rows[i]+idx;
	
  }
//
	std::cout<< "rows_new[i]" <<std::endl;
	 for ( int i = 0; i < t_nrow+1; i++ )
  {
			std::cout<< rows_new[i] <<std::endl;
  }
  
  
//Create cols_new
for ( int i = 0; i < size_coup; i++ )
  {
	idx_insert[i]=rows_new[rows_coupling_test[i]];


		if(i>0){
			if(rows_coupling_test[i]==rows_coupling_test[i-1]&& i<size_coup){
			idx_insert[i]=idx_insert[i-1]+1;
			}
		}
		
  }
  
  idx=0;
  for ( int i = 0; i < num_vals_t+size_coup; i++ )
  {
	  if(idx_insert[idx]==i && idx<size_coup){
		cols_new[i]=cols_coupling_test[idx];
		vals_new[i]=vals_coupling_test[idx];

		idx=idx+1;
	  }else{
		cols_new[i]=t_cols[i-idx];
		vals_new[i]=t_array[i-idx];

	  }
  }
  
    	std::cout<< "sort_idx[i]" <<std::endl;
	 for ( int i = 0; i < size_coup; i++ )
  {
			std::cout<< sort_idx[i] <<std::endl;
  }
  
   	std::cout<< "idx_insert[i]" <<std::endl;
	 for ( int i = 0; i < size_coup; i++ )
  {
			std::cout<< idx_insert[i] <<std::endl;
  }
  
  	std::cout<< "cols_new[i]" <<std::endl;
	 for ( int i = 0; i < num_vals_t+size_coup; i++ )
  {
			std::cout<< cols_new[i] <<std::endl;
  }

  
    	std::cout<< "vals_new[i]" <<std::endl;
	 for ( int i = 0; i < num_vals_t+size_coup; i++ )
  {
			std::cout<< vals_new[i] <<std::endl;
  }
   
  
  
//Need to create this using the updated csr format
MatCreateSeqAIJWithArrays(comm,t_nrow,t_nrow,rows_new,cols_new,vals_new,&T_new);
  MatSetFromOptions(T_new);
  MatSetUp(T_new);	

  MatAssemblyBegin(T_new,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(T_new,MAT_FINAL_ASSEMBLY);
	
	std::cout<<"T_new "<<std::endl;
	MatView(T_new,viewer);

	 
	