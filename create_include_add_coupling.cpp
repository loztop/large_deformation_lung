int size_coup=idx_count_f;
 
#if !mats
PetscInt *sort_idx;
sort_idx=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));

PetscInt *cols_coupling_test_cpy;
cols_coupling_test_cpy=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));

PetscReal *vals_coupling_test_cpy;
vals_coupling_test_cpy=(PetscReal *)malloc((size_coup)*sizeof(PetscReal));

PetscInt *idx_insert;
idx_insert=(PetscInt *)malloc((size_coup)*sizeof(PetscInt));
#endif

  for ( int i = 0; i < size_coup; i++ )
   {
sort_idx[i]=i;
  }


//First sort the arrays
PetscSortIntWithArray(size_coup,rows_coupling_f,sort_idx);
//make copies

//THIS IS BAD
// vals_coupling_test_cpy=vals_coupling_f;
// cols_coupling_test_cpy=cols_coupling_f;

#if fvec

PetscVector<Number> cols_coupling_fvecp(cols_coupling_fvec) ;
//cols_coupling_fvecp.add(1,rhs_in);
#endif


for ( int i = 0; i < size_coup; i++ )
{
  vals_coupling_test_cpy[i]=vals_coupling_f[i];
	
#if !fvec
    cols_coupling_test_cpy[i]=cols_coupling_f[i];
#endif
	
#if fvec
		     cols_coupling_test_cpy[i]=cols_coupling_fvecp(i);
#endif

}


   for ( int i = 0; i < size_coup; i++ )
   {
vals_coupling_f[i]=vals_coupling_test_cpy[sort_idx[i]];


	#if !fvec
cols_coupling_f[i]=cols_coupling_test_cpy[sort_idx[i]];
#endif
	
#if fvec
cols_coupling_fvecp.set(i,cols_coupling_test_cpy[sort_idx[i]]);
#endif
				 
  }
  
  
#if !mats
PetscInt *rows_new;
rows_new=(PetscInt *)malloc((a_nrow+t_nrow+1)*sizeof(PetscInt));
#endif

for ( int i = 0; i < a_nrow+t_nrow; i++ )
{
  rows_new[i]=0;
}

 PetscInt *cols_new;
cols_new=(PetscInt *)malloc((num_vals_t+num_vals_a+size_coup)*sizeof(PetscInt));
 

for ( int i = 0; i < num_vals_t+num_vals_a+size_coup; i++ )
{
  cols_new[i]=0;
}

 PetscReal *vals_new;
vals_new=(PetscReal *)malloc((num_vals_t+num_vals_a+size_coup)*sizeof(PetscReal));
 

for ( int i = 0; i < num_vals_t+num_vals_a+size_coup; i++ )
{
  vals_new[i]=0;
}


//Create rows_new
int idx=0;

  for ( int i = 1; i < a_nrow+t_nrow+1; i++ )
  {

/*
if(rows_coupling_f[idx]==i-1 && idx<size_coup){
idx=idx+1;
}
//Also check the previous was not also the same (i.e. more than one entry in same row)
if(idx>0){
if(rows_coupling_f[idx]==rows_coupling_f[idx-1]&& idx<size_coup){
idx=idx+1;
}
}
*/

//Count how many rows_coupling_test[idx+j]=i
//then idx=idx+1+j at endl

if(rows_coupling_f[idx]==i-1 && idx<size_coup){

int stop=0;
int j=1;

while ( stop<1 )
{

if(rows_coupling_f[idx+j]==rows_coupling_f[idx]){
j=j+1;
}else{
stop=1;
}

}

idx=idx+j;

}


rows_new[i]=big_rows_t[i]+idx;
  }
    
//Create cols_new and vals_new
for ( int i = 0; i < size_coup; i++ )
  {
idx_insert[i]=rows_new[rows_coupling_f[i]];

if(i>0){
if(rows_coupling_f[i]==rows_coupling_f[i-1]&& i<size_coup){
idx_insert[i]=idx_insert[i-1]+1;
}
}

  }
  
  idx=0;
  for ( int i = 0; i < num_vals_t+num_vals_a+size_coup; i++ )
  {
if(idx_insert[idx]==i && idx<size_coup){



	#if !fvec
	cols_new[i]=cols_coupling_f[idx];
#endif
	
#if fvec
	cols_new[i]=cols_coupling_fvecp(idx);
#endif
				 

vals_new[i]=vals_coupling_f[idx];

idx=idx+1;
}else{
cols_new[i]=big_cols_t[i-idx];
vals_new[i]=b_array[i-idx];

}
  }
  
