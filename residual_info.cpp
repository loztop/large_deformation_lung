	  //Get the norm of the residual to check for convergence.
		TransientLinearImplicitSystem & newton_update =equation_systems.get_system<TransientLinearImplicitSystem> ("Newton-update");
		
		
		//Real solid_residual=newton_update.rhs->l2_norm ();
       Real solid_residual=big_rp.l2_norm ();
		
		
		Real n_elements=mesh.n_elem( 	);	
		Real n_nodes=mesh.n_nodes( 	);	

		Real s_start=0;
		Real s_end=3*n_nodes -1;
		Real p_start=3*n_nodes;
		Real p_end=3*n_nodes + n_elements -1;
		Real z_start=3*n_nodes + n_elements;
		Real z_end=6*n_nodes + n_elements-1;

		Real momentum_res=0;
		for (unsigned int i=s_start; i<s_end+1; ++i)
		{
			momentum_res=momentum_res+big_rp(i)*big_rp(i);
		}
		std::cout<< "momentum_res " << pow(momentum_res,0.5) << std::endl;
		
		Real mass_res=0;
		for (unsigned int i=p_start; i<p_end+1; ++i)
		{
			mass_res=mass_res+big_rp(i)*big_rp(i);
		}
		std::cout<< "mass_res " << pow(mass_res,0.5) << std::endl;
		
		Real fluid_res=0;
		for (unsigned int i=z_start; i<z_end+1; ++i)
		{
			fluid_res=fluid_res+big_rp(i)*big_rp(i);
		}
		std::cout<< "fluid_res " << pow(fluid_res,0.5) << std::endl;
		
		
		Real t_start=size_fem;
		Real t_end=size_fem+size_tree-1;

		Real tree_res=0;
		for (unsigned int i=t_start; i<t_end+1; ++i)
		{
			tree_res=tree_res+big_rp(i)*big_rp(i);
		}
		std::cout<< "tree_res " << pow(tree_res,0.5) << std::endl;
				 
		 
    /*****///Convergence and Norm Computing Stuff///***********/
    change_in_newton_update->add (-1., *newton_update.solution);
		change_in_newton_update->close();
		Real norm_delta = change_in_newton_update->l2_norm();
    change_in_newton_update->add (-1., *newton_update.solution);
    change_in_newton_update->close();
    norm_delta = change_in_newton_update->l2_norm();
    Real l2_soln=last_non_linear_soln.current_local_solution->l2_norm ();

    const unsigned int n_linear_iterations = newton_update.n_linear_iterations();
    const Real final_linear_residual = newton_update.final_linear_residual();   
    std::cout << "Poro: Linear conv at step: "
                    << n_linear_iterations
                    << ", resid: "
                    << final_linear_residual/l2_soln 
                    << ", time: " << double(diffclock(end_solid_solve,begin_solid_solve)) << " ms"
                    <<std::endl;             
    std::cout   << "Poro Nonlinear convergence: ||u - u_old|| = "
                    << norm_delta/l2_soln << ", poro res:  " << solid_residual/l2_soln
                    << std::endl;
  
   
   
   
