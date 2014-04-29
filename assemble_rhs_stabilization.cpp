 	std::vector<unsigned int> stab_dofs_cols2;
	std::vector<Real> stab_dofs_vals2;


	std::vector<Real> stab_dofs_vals_rhs;

    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
				//Only do something on the interior edges.
      } //if (elem->neighbor(s) == NULL)
			else{

      const Elem* neighbor = elem->neighbor(s);

			std::vector<unsigned int> neighbor_dof_indices_p;
      dof_map.dof_indices(neighbor, neighbor_dof_indices_p,p_var);
      const unsigned int n_neighbor_dofs_p = neighbor_dof_indices_p.size();

			Real hmax=(*elem).hmax();
			Real hmin=(*elem).hmin();

			const Real DELTA    = es.parameters.get<Real>("DELTA");
			Real delta=DELTA/dt;
			Real factor=-delta*(hmax*hmax*hmax);
			factor=factor*1;
			
			//Accumulate RHS
			std::vector< Real > p_old;
			last_non_linear_soln.old_local_solution->get(dof_indices_p,p_old);
			std::vector< Real > p_old_n;
			last_non_linear_soln.old_local_solution->get(neighbor_dof_indices_p,p_old_n);
      
			std::vector< Real > p_curr;
			last_non_linear_soln.current_local_solution->get(dof_indices_p,p_curr);
			std::vector< Real > p_curr_n;
			last_non_linear_soln.current_local_solution->get(neighbor_dof_indices_p,p_curr_n);
			
			stab_dofs_vals_rhs.push_back(-(p_curr[0]-p_old[0])*factor);
			stab_dofs_vals_rhs.push_back((p_curr_n[0]-p_old_n[0])*factor);
			
      //perf_log.pop("push back");

			}
			
		}


//add rhs
DenseVector<Number> rhs_stab;
rhs_stab.resize(dof_indices_p.size());
 	
for (int i=0; i < stab_dofs_vals_rhs.size(); i++) {
	 	rhs_stab(0)+=stab_dofs_vals_rhs[i];
	}

std::vector<unsigned int> stab_dofs_rows_rhs;
stab_dofs_rows_rhs.push_back(dof_indices_p[0]);
newton_update.rhs->add_vector(rhs_stab, stab_dofs_rows_rhs);