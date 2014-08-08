	//Need this to repress PETSC error about inserting(adding) to zero slot
MatSetOption(big_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;
//MatSetOption(big_A,MAT_ROW_ORIENTED,PETSC_FALSE);

#if !mats
PetscInt *end_j_f;
end_j_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
#endif

   int c_f=0;
for (int j=0; j < tree.number_nodes; j++) {	
if(tree.nodes_type(j)==0){	
end_j_f[c_f]=j;
c_f=c_f+1;
}
}

 
Real res_sum_f=0;

 #if !mats
   PetscReal *omega_end_j_f;
   omega_end_j_f=(PetscReal *)malloc((tree.number_nodes)*sizeof(PetscReal));
#endif
   
for (int j=0; j < tree.number_nodes; j++) {	
omega_end_j_f[j]=0;
}

MeshBase::const_element_iterator el_omega_f = mesh.active_local_elements_begin();
const MeshBase::const_element_iterator end_el_omega_f = mesh.active_local_elements_end();
for ( ; el_omega_f != end_el_omega_f; ++el_omega_f)
{
const Elem* elem = *el_omega_f;

//Put mesh back to reference
for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = reference.current_local_solution->el(source_dof);
            (*node)(d)= value;
          }
      }
      

Point elem_pos=elem->centroid();

//Find the closest distal tree_node to q_point[qp]
int closest_j=0;
Real dist_j=0;
Real closest_dist=99999999;
int closest_end_j=0;
int closest_end_edge_j=0;
int closest_end_node=0;


for (int j=0; j < c_f; j++) {	
Real dist_j=pow(elem_pos(0)- tree.nodes(end_j_f[j])(0),2)+pow(elem_pos(1)- tree.nodes(end_j_f[j])(1),2)+pow(elem_pos(2)- tree.nodes(end_j_f[j])(2),2);

if(dist_j<closest_dist)
{
closest_dist=dist_j;

//Relies on assumption that lower node is +1 of its parent edge.
closest_end_edge_j=end_j_f[j]-1;	
}	
}


//Put mesh back to current
for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)= value;
          }
      }
      
      omega_end_j_f[closest_end_edge_j+1]=omega_end_j_f[closest_end_edge_j+1]+elem->volume();

}

 
 ///End of finding omega j

  std::vector<unsigned int> dof_indices_p_f;
const unsigned int p_var_f = system.variable_number ("p_nu");
const DofMap & dof_map_f = system.get_dof_map();

MeshBase::const_element_iterator el_coup_f = mesh.active_local_elements_begin();
const MeshBase::const_element_iterator end_el_coup_f = mesh.active_local_elements_end();

//First find all the nodes that are end nodes, later also deal with all the end nodes that didn't get chosen, display how many of them and the co-ords

#if !mats
PetscInt *end_j_zero_f;
end_j_zero_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

PetscInt *found_end_j_f;
found_end_j_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

PetscInt *already_set_f;
already_set_f=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));

//Create some arrays to hold inserts
PetscInt *rows_coupling_f;
rows_coupling_f=(PetscInt *)malloc((a_nrow+t_nrow)*sizeof(PetscInt));
#endif


#if !mats

#if !fvec
PetscInt *cols_coupling_f;
PetscMalloc((3*a_nrow+t_nrow)*sizeof(PetscInt),&cols_coupling_f); 
#endif

#endif

#if fvec
	Vec cols_coupling_fvec;
	VecCreate(PETSC_COMM_WORLD,&cols_coupling_fvec);
  PetscObjectSetName((PetscObject) cols_coupling_fvec, "Solution");
  VecSetSizes(cols_coupling_fvec,PETSC_DECIDE,3*a_nrow+t_nrow);
  VecSetFromOptions(cols_coupling_fvec);
  VecSetFromOptions(cols_coupling_fvec);
  VecAssemblyBegin(cols_coupling_fvec);
  VecAssemblyEnd(cols_coupling_fvec);
#endif
  		
   

#if !mats

PetscReal *vals_coupling_f;
vals_coupling_f=(PetscReal *)malloc((3*a_nrow+t_nrow)*sizeof(PetscReal));


PetscReal *rhs_new;
rhs_new=(PetscReal *)malloc((a_nrow+t_nrow)*sizeof(PetscReal));

#endif


 
for ( int i = 0; i < a_nrow+t_nrow; i++ )
  {
rows_coupling_f[i]=0;
rhs_new[i]=0;
}

for ( int i = 0; i < 3*a_nrow+t_nrow; i++ )
  {
vals_coupling_f[i]=0;
  }
  
 #if !fvec
   for ( int i = 0; i < 3*a_nrow+t_nrow; i++ )
  {
cols_coupling_f[i]=0;
  }
#endif

int idx_count_f=0;


for (int j=0; j < tree.number_nodes; j++) {	
found_end_j_f[j]=0;
already_set_f[j]=0;
}




for ( ; el_coup_f != end_el_coup_f; ++el_coup_f)
{
const Elem* elem = *el_coup_f;


//Put mesh back to reference
for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = reference.current_local_solution->el(source_dof);
            (*node)(d)= value;
          }
      }
      
      
      
Point elem_pos=elem->centroid();

//Find the closest distal tree_node to q_point[qp]
int closest_j=0;
Real dist_j=0;
Real closest_dist=99999999;
int closest_end_j=0;
int closest_end_edge_j=0;
int closest_end_node=0;


for (int j=0; j < c_f; j++) {	
//Real dist_j=pow(elem_pos(0)- tree.nodes_deformed(end_j[j])(0),2)+pow(elem_pos(1)- tree.nodes_deformed(end_j[j])(1),2)+pow(elem_pos(2)- tree.nodes_deformed(end_j[j])(2),2);
Real dist_j=pow(elem_pos(0)- tree.nodes(end_j_f[j])(0),2)+pow(elem_pos(1)- tree.nodes(end_j_f[j])(1),2)+pow(elem_pos(2)- tree.nodes(end_j_f[j])(2),2);

if(dist_j<closest_dist)
{
closest_dist=dist_j;

//Relies on assumption that lower node is +1 of its parent edge.
closest_end_edge_j=end_j_f[j]-1;	
}	

}

found_end_j_f[closest_end_edge_j+1]=found_end_j_f[closest_end_edge_j+1]+1;

dof_map.dof_indices (elem, dof_indices_p_f, p_var);

//Put mesh back to current
for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)= value;
          }
      }
      


Real constant=elem->volume()/(omega_end_j_f[closest_end_edge_j+1]);	


//Poro Mass conservation coupling term (source) , distal end flow is coupled to elemen volume

/************************/
//MatSetValue(big_A, dof_indices_p_f[0], AP.m()+tree.number_nodes+closest_end_edge_j, - constant,ADD_VALUES);


rows_coupling_f[idx_count_f]=dof_indices_p_f[0];

#if !fvec

    cols_coupling_f[idx_count_f]=AP.m()+tree.number_nodes+closest_end_edge_j;
	#endif
	
#if fvec
		VecSetValue(cols_coupling_fvec, idx_count_f, AP.m()+tree.number_nodes+closest_end_edge_j,ADD_VALUES); 
#endif
		
		
vals_coupling_f[idx_count_f]= - constant;
idx_count_f=idx_count_f+1;



//Need to also take care of the residual
VecSetValue(big_r, dof_indices_p_f[0], - constant*tree.edges_flowrate(closest_end_edge_j),ADD_VALUES);
rhs_new[dof_indices_p_f[0]]= rhs_new[dof_indices_p_f[0]] - constant*tree.edges_flowrate(closest_end_edge_j);



//tree bit - Set Pporo=Pdistl - this seems to work
/******************************************/
//MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , dof_indices_p_f[0] ,-elem->volume()/(omega_end_j_f[closest_end_edge_j+1]) ,ADD_VALUES);

rows_coupling_f[idx_count_f]=AP.m()+tree.number_nodes+ closest_end_edge_j;

 #if !fvec
cols_coupling_f[idx_count_f]=dof_indices_p_f[0];
#endif
#if fvec
		VecSetValue(cols_coupling_fvec, idx_count_f,dof_indices_p_f[0],ADD_VALUES); 
#endif
		
vals_coupling_f[idx_count_f]= -elem->volume()/(omega_end_j_f[closest_end_edge_j+1]);
idx_count_f=idx_count_f+1;	


//Only do this once ! so eqn is actually integrated properly
if(already_set_f[closest_end_edge_j+1]==0){

/*************************/
// MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , AP.m()+closest_end_edge_j +1 , 1 ,ADD_VALUES);

rows_coupling_f[idx_count_f]=AP.m()+tree.number_nodes+ closest_end_edge_j;
 #if !fvec
cols_coupling_f[idx_count_f]=AP.m()+closest_end_edge_j +1;
#endif
#if fvec
		VecSetValue(cols_coupling_fvec, idx_count_f,AP.m()+closest_end_edge_j +1,ADD_VALUES); 
#endif
vals_coupling_f[idx_count_f]=1;
idx_count_f=idx_count_f+1;	


already_set_f[closest_end_edge_j+1]=1;

//Do same for residual
VecSetValue(big_r, AP.m()+tree.number_nodes+ closest_end_edge_j , tree.nodes_pressure(closest_end_edge_j +1),ADD_VALUES);
rhs_new[(int)AP.m()+(int)tree.number_nodes+ (int)closest_end_edge_j]=rhs_new[(int)AP.m()+(int)tree.number_nodes+ (int)closest_end_edge_j]+tree.nodes_pressure(closest_end_edge_j +1);

res_sum_f=tree.nodes_pressure(closest_end_edge_j +1)+res_sum_f;
}



//Assemble residual for tree !!
//get current pressure
Real p_poro=last_non_linear_soln.current_local_solution->el(dof_indices_p_f[0]);
res_sum_f=-(p_poro*elem->volume())/(omega_end_j_f[closest_end_edge_j+1])+res_sum_f;

VecSetValue(big_r, AP.m()+tree.number_nodes+ closest_end_edge_j , -(p_poro*elem->volume())/(omega_end_j_f[closest_end_edge_j+1]),ADD_VALUES);
rhs_new[(int)AP.m()+(int)tree.number_nodes+ (int)closest_end_edge_j] =rhs_new[(int)AP.m()+(int)tree.number_nodes+ (int)closest_end_edge_j]-(p_poro*elem->volume())/(omega_end_j_f[closest_end_edge_j+1]);

} //end of element loop

 
//Make note of any that did get not assigned.
int z_f=0;
for (int j=0; j < c_f; j++) {	
if(found_end_j_f[end_j_f[j]]==0){	
end_j_zero_f[z_f]=j;
z_f=z_f+1;
}
}
std::cout<< "z count= " << z_f <<std::endl;

//Set outflow of nodes that didnt get assigned to zero
for (int j=0; j < z_f; j++) {	
//Set outflow to zero
//MatSetValue(big_A, AP.m()+tree.number_nodes+ end_j_f[end_j_zero_f[j]] -1, AP.m()+tree.number_nodes+ end_j_f[end_j_zero_f[j]] -1, 1 ,ADD_VALUES);

//Still need to do this
rows_coupling_f[idx_count_f]=AP.m()+tree.number_nodes+ end_j_f[end_j_zero_f[j]] -1;
 #if !fvec
    cols_coupling_f[idx_count_f]=AP.m()+tree.number_nodes+ end_j_f[end_j_zero_f[j]] -1;
#endif
		#if fvec
		VecSetValue(cols_coupling_fvec, idx_count_f,AP.m()+tree.number_nodes+ end_j_f[end_j_zero_f[j]] -1,ADD_VALUES); 
#endif	
vals_coupling_f[idx_count_f]= 1;
idx_count_f=idx_count_f+1;

}


std::cout<< "Coupling press res: "<< res_sum_f <<std::endl;