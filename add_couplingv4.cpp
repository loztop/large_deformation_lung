  std::vector<unsigned int> dof_indices_p;
	const unsigned int p_var = system.variable_number ("p_nu");
	const DofMap & dof_map = system.get_dof_map();

	MeshBase::const_element_iterator       el_coup     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el_coup = mesh.active_local_elements_end(); 
	
	//First find all the nodes that are end nodes, later also deal with all the end nodes that didn't get chosen, display how many of them and the co-ords
	
	PetscInt   *end_j;
	end_j=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
	PetscInt   *end_j_zero;
	end_j_zero=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
	PetscInt   *found_end_j;
	found_end_j=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
	int c=0;
	for (int j=0; j < tree.number_nodes; j++) {	
		  if(tree.nodes_type(j)==0){	
				end_j[c]=j;
				c=c+1;
			}
	}
		
	for (int j=0; j < tree.number_nodes; j++) {	
		found_end_j[j]=0;
	}
	
	for ( ; el_coup != end_el_coup; ++el_coup)
	{    	
	  const Elem* elem = *el_coup;
		Point elem_pos=elem->centroid();
		
		//Find the closest distal tree_node to q_point[qp]
		int closest_j=0;
		Real dist_j=0;
		Real closest_dist=99999999;
		int closest_end_j=0;
		int closest_end_edge_j=0;
		int closest_end_node=0;

		
		for (int j=0; j < c; j++) {			
				Real dist_j=pow(elem_pos(0)-  tree.nodes_deformed(end_j[j])(0),2)+pow(elem_pos(1)-  tree.nodes_deformed(end_j[j])(1),2)+pow(elem_pos(2)-  tree.nodes_deformed(end_j[j])(2),2);
											
				if(dist_j<closest_dist)
				{
					closest_dist=dist_j; 
					
					//Relies on assumption that lower node is +1 of its parent edge.
					closest_end_edge_j=end_j[j]-1;			
				}			
		}
		
		
	  	found_end_j[closest_end_edge_j+1]=found_end_j[closest_end_edge_j+1]+1;
					
		dof_map.dof_indices (elem, dof_indices_p, p_var);

		//Approx factor of distal airway resistance - 32
		Real distal_fac=1;
		
		//Need this to repress PETSC error about inserting(adding) to zero slot
		MatSetOption(big_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;
				
		//Poro Mass conservation coupling term (source)  , distal end flow is coupled to elemen volume
	 	MatSetValue(big_A, dof_indices_p[0], AP.m()+tree.number_nodes+closest_end_edge_j, - elem->volume(),ADD_VALUES); 

		//tree bit - Set Pporo=Pdistl
		MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , dof_indices_p[0]  ,-1 ,ADD_VALUES); 
	 
	  MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , AP.m()+closest_end_edge_j +1 , 1 ,ADD_VALUES); 
				
	}
	
	//Makke note of any that did get not assigned.
		int z=0;
		for (int j=0; j < c; j++) {		
			if(found_end_j[end_j[j]]==0){	
				end_j_zero[z]=j;
				z=z+1;
			}
		}
	//std::cout<< "Number of unassigned nodes, z count= " << z <<std::endl;

	//Set outflow of nodes that didnt get assigned to zero
	for (int j=0; j < z; j++) {	
	  //Set outflow to zero
	  MatSetValue(big_A, AP.m()+tree.number_nodes+ end_j[end_j_zero[j]] -1, AP.m()+tree.number_nodes+ end_j[end_j_zero[j]] -1, 1 ,ADD_VALUES);		
	}
	
  