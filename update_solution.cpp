	std::cout<<"Update newton " <<std::endl;
	for (int i=0; i<size_mat; i++) {
		newton_update.solution->set(i,big_xp(i));
		newton_update.current_local_solution->set(i,big_xp(i));
	}
	newton_update.solution->close();
	newton_update.current_local_solution->close();
	newton_update.update();   
	
	
	
	//Copy solution back to tree
	//Update the distal pressures

	std::cout<<"Update distal pressures " <<std::endl;

	for (int i=0; i < tree.number_nodes; i++) {
		tree.nodes_pressure(i)=tree.nodes_pressure(i) -big_xp(size_mat+i);	
		
		//tree.nodes_pressure(i)= big_xp(size_mat+i);	
	}
	
	//Update the flowrates
	std::cout<<"Update flow rates " <<std::endl;
	for (int i=0; i < tree.number_edges; i++) {
		tree.edges_flowrate(i)=tree.edges_flowrate(i)-big_xp(size_mat+tree.number_nodes+i);
	
	//		tree.edges_flowrate(i)=big_xp(size_mat+tree.number_nodes+i);

	}
	
	 
	
		clock_t end_reinit=clock();
		std::cout<<"Reinit " <<std::endl;

		equation_systems.reinit();

		clock_t begin_reinit=clock();
		
		std::cout<<"reinit: " << double(diffclock(end_reinit,begin_reinit)) <<  " ms"<<std::endl;
		
	
		//update the final solution
		// xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
		//Apply a full Newton-step
		Real K=1; //Newton step size

		last_non_linear_soln.solution->add(-1*K,*newton_update.solution);
	 		//	last_non_linear_soln.reinit();
	 //newton_update.reinit();   

			last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(-1,*newton_update.current_local_solution);
		last_non_linear_soln.current_local_solution->close();
		last_non_linear_soln.update();

		
		//const DofMap & dof_map = postvars .get_dof_map();
		std::vector<unsigned int> dof_indices_j;
		const unsigned int J_var = postvars.variable_number ("J");

	
		MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
		const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 
		Real total_volume=0;
    for ( ; el != end_el; ++el)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
      
      //Write out new element volume
			
			Real elem_vol=elem->volume();
			total_volume=total_volume+elem_vol;
			dof_map.dof_indices (elem, dof_indices_j, J_var);
			Real elem_vol_ref = reference.current_local_solution->el(dof_indices_j[0]);
	//		postvars.current_local_solution->set(dof_indices_j[0], elem_vol/elem_vol_ref);
	//		postvars.solution->set(dof_indices_j[0], elem_vol/elem_vol_ref);

      
      
    }
    
    std::cout<<"total_volume_change "<< total_volume -total_volume_ref<< std::endl;
    Real global_J=total_volume / total_volume_ref;
	  std::cout<<"total_outflow  "<< total_outflow+tree.edges_flowrate(0)*dt << std::endl;
 
		
    equation_systems.update();
    equation_systems.allgather();

	
	
	clock_t begin_pos_update=clock();
    tree.update_positions(equation_systems);    
    clock_t end_pos_update=clock();

	std::cout<<"pos_update: " << double(diffclock(end_pos_update,begin_pos_update)) <<  " ms"<<std::endl;