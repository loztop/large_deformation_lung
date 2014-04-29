	for (int i=0; i<size_mat; i++) {
		newton_update.solution->set(i,big_xp(i));
		newton_update.current_local_solution->set(i,big_xp(i));
	}
	newton_update.solution->close();
	newton_update.current_local_solution->close();
	newton_update.update();   
	
	
	//Copy solution back to tree
	//Update the distal pressures
	for (int i=0; i < tree.number_nodes; i++) {
		//tree.nodes_pressure(i)=tree.nodes_pressure(i) -big_xp(size_mat+i);	
		
		tree.nodes_pressure(i)= big_xp(size_mat+i);	
	}
	
	//Update the flowrates
	for (int i=0; i < tree.number_edges; i++) {
//		tree.edges_flowrate(i)=tree.edges_flowrate(i)-big_xp(size_mat+tree.number_nodes+i);
	
		tree.edges_flowrate(i)=big_xp(size_mat+tree.number_nodes+i);

	}
	
	
	
	equation_systems.reinit();

		//update the final solution
		// xn+1 = xn + delta xn  (note that -* is from K(delatxn) = - R, ie K(-delatxn)=R )
		//Apply a full Newton-step
		Real K=1; //Newton step size

		last_non_linear_soln.solution->add(-1*K,*newton_update.solution);
		last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(-1,*newton_update.current_local_solution);
    last_non_linear_soln.current_local_solution->close();
    last_non_linear_soln.update();

		MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.local_elements_end(); 
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
    }
    equation_systems.update();
    equation_systems.allgather();
    
    tree.update_positions(equation_systems);    
    
   