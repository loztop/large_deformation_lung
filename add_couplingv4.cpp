 /*
 //////////////////
  	//First put solution on reference mesh then put back to deformed mesh
	 MeshBase::const_element_iterator       el_a     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el_a = mesh.local_elements_end(); 
    for ( ; el_a != end_el_a; ++el_a)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el_a;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = reference.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
    }
    equation_systems.update();
    equation_systems.allgather();
		
  ////////////////
  */
 
 	clock_t find_start=clock();	

	 
  PetscInt   *end_j;
	end_j=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
  	
  	int c=0;
	for (int j=0; j < tree.number_nodes; j++) {	
		  if(tree.nodes_type(j)==0){	
				end_j[c]=j;
				c=c+1;
			}
	}

	
	
  Real res_sum=0;
  
 	PetscReal   *omega_end_j;
	omega_end_j=(PetscReal *)malloc((tree.number_nodes)*sizeof(PetscReal));
	
   for (int j=0; j < tree.number_nodes; j++) {	
		omega_end_j[j]=0;
	}
	
	

	
	MeshBase::const_element_iterator       el_omega     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el_omega = mesh.active_local_elements_end(); 
	for ( ; el_omega != end_el_omega; ++el_omega)
	{    	
	  const Elem* elem = *el_omega;
		
		
		
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

		
		for (int j=0; j < c; j++) {			
				Real dist_j=pow(elem_pos(0)-  tree.nodes_deformed(end_j[j])(0),2)+pow(elem_pos(1)-  tree.nodes_deformed(end_j[j])(1),2)+pow(elem_pos(2)-  tree.nodes_deformed(end_j[j])(2),2);
											
				if(dist_j<closest_dist)
				{
					closest_dist=dist_j; 
					
					//Relies on assumption that lower node is +1 of its parent edge.
					closest_end_edge_j=end_j[j]-1;			
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
      
      omega_end_j[closest_end_edge_j+1]=omega_end_j[closest_end_edge_j+1]+elem->volume();
						
}
	
 
 ///End of finding omega j
  	clock_t find_end=clock();	

  

		
  std::vector<unsigned int> dof_indices_p;
	const unsigned int p_var = system.variable_number ("p_nu");
	const DofMap & dof_map = system.get_dof_map();

	MeshBase::const_element_iterator       el_coup     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el_coup = mesh.active_local_elements_end(); 
	
	//First find all the nodes that are end nodes, later also deal with all the end nodes that didn't get chosen, display how many of them and the co-ords
	
	
	PetscInt   *end_j_zero;
	end_j_zero=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
	PetscInt   *found_end_j;
	found_end_j=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
	
	PetscInt   *already_set;
	already_set=(PetscInt *)malloc((tree.number_nodes)*sizeof(PetscInt));
		
		
	for (int j=0; j < tree.number_nodes; j++) {	
		found_end_j[j]=0;
		already_set[j]=0;

	}
	

	for ( ; el_coup != end_el_coup; ++el_coup)
	{    	
	  const Elem* elem = *el_coup;
		
		
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

		
		for (int j=0; j < c; j++) {			
				//Real dist_j=pow(elem_pos(0)-  tree.nodes_deformed(end_j[j])(0),2)+pow(elem_pos(1)-  tree.nodes_deformed(end_j[j])(1),2)+pow(elem_pos(2)-  tree.nodes_deformed(end_j[j])(2),2);
				Real dist_j=pow(elem_pos(0)-  tree.nodes(end_j[j])(0),2)+pow(elem_pos(1)-  tree.nodes(end_j[j])(1),2)+pow(elem_pos(2)-  tree.nodes(end_j[j])(2),2);
					
				if(dist_j<closest_dist)
				{
					closest_dist=dist_j; 
					
					//Relies on assumption that lower node is +1 of its parent edge.
					closest_end_edge_j=end_j[j]-1;			
				}			
				
			//			std::cout<<tree.nodes(end_j[j]) <<std::endl;
			//			std::cout<< elem_pos <<std::endl;

		}
		
	  	found_end_j[closest_end_edge_j+1]=found_end_j[closest_end_edge_j+1]+1;
					
		dof_map.dof_indices (elem, dof_indices_p, p_var);

		
		
		
		
			//Put mesh back to current
		 for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)= value;
          }
      }
      
	
		//Need this to repress PETSC error about inserting(adding) to zero slot
		MatSetOption(big_A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) ;
	
		Real constant=elem->volume()/(omega_end_j[closest_end_edge_j+1]);	
		
		//Poro Mass conservation coupling term (source)  , distal end flow is coupled to elemen volume
	 	MatSetValue(big_A, dof_indices_p[0], AP.m()+tree.number_nodes+closest_end_edge_j, - constant,ADD_VALUES); 
	   //Need to also take care of the residual
	 
		VecSetValue(big_r, dof_indices_p[0], - constant*tree.edges_flowrate(closest_end_edge_j),ADD_VALUES); 

		
		
		
		//tree bit - Set Pporo=Pdistl - this seems to work
		MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , dof_indices_p[0]  ,-elem->volume()/(omega_end_j[closest_end_edge_j+1]) ,ADD_VALUES); 
				
	
		//Only do this once !  so eqn is actually integrated properly
		if(already_set[closest_end_edge_j+1]==0){
		  MatSetValue(big_A, AP.m()+tree.number_nodes+ closest_end_edge_j , AP.m()+closest_end_edge_j +1 , 1 ,ADD_VALUES); 		  
		  
		  already_set[closest_end_edge_j+1]=1;

			//Do same for residual
			VecSetValue(big_r, AP.m()+tree.number_nodes+ closest_end_edge_j , tree.nodes_pressure(closest_end_edge_j +1),ADD_VALUES); 

			res_sum=tree.nodes_pressure(closest_end_edge_j +1)+res_sum; 
		}
		
		
			
		//Assemble residual for tree !!
		//get current pressure	
	 	Real p_poro=last_non_linear_soln.current_local_solution->el(dof_indices_p[0]);
	 	res_sum=-(p_poro*elem->volume())/(omega_end_j[closest_end_edge_j+1])+res_sum;
	 
		VecSetValue(big_r, AP.m()+tree.number_nodes+ closest_end_edge_j , -(p_poro*elem->volume())/(omega_end_j[closest_end_edge_j+1]),ADD_VALUES); 

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
	
	
	std::cout<<  "tree coupling pressure residual "<< res_sum <<std::endl;
	
	
	
	
	/*
	/////////////////////////////////
	    MeshBase::const_element_iterator       el_b     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el_b = mesh.local_elements_end(); 
    for ( ; el_b != end_el_b; ++el_b)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el_b;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
    }
    equation_systems.update();
    equation_systems.allgather();
	
	//////////////////////////////////
	*/
  
	
