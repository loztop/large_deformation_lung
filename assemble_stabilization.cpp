//Pressure jump stabilisation.

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

      //perf_log.push("push back");
      stab_dofs_cols2.push_back(dof_indices_p[0]);
      stab_dofs_vals2.push_back(-factor);
      stab_dofs_cols2.push_back(neighbor_dof_indices_p[0]);
      stab_dofs_vals2.push_back(factor);
			}
		}

	//Add stabilisation contr
	//Lots of mallocs happen during the first assemble due to the unexpected sparsity pattern;
	//The implicit_neighbor_dof does not seem to work ?
	//perf_log.push("kstab");
  DenseMatrix<Number> Kstab2;
  Kstab2.resize(1, stab_dofs_vals2.size());
 	for (int i=0; i < stab_dofs_vals2.size(); i++) {
	 	Kstab2(0,i)=stab_dofs_vals2[i];
	}
 	std::vector<unsigned int> stab_dofs_rows2;
	stab_dofs_rows2.push_back(dof_indices_p[0]);
	newton_update.matrix->add_matrix(Kstab2,stab_dofs_rows2,stab_dofs_cols2);
	//perf_log.pop("kstab");
	
	