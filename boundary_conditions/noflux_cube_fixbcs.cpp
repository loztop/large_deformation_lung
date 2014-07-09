    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
        AutoPtr<Elem> side (elem->build_side(s));

        for (unsigned int ns=0; ns<side->n_nodes(); ns++)
        {

          for (unsigned int n=0; n<elem->n_nodes(); n++)
          {
            Node *node = elem->get_node(n);
			Point p;
			for (unsigned int d = 0; d < 3; ++d) {
			  unsigned int source_dof = node->dof_number(1, d, 0);
			  Real value = ref_sys.current_local_solution->el(source_dof);
			  p(d)=value;
			}

			
            int source_dof = node->dof_number(last_non_linear_soln.number(), x_var, 0);
            if((source_dof<12345678) && (source_dof>-1) && (p(0)<0.0001 || p(0)>0.999)){
									Real value = last_non_linear_soln.current_local_solution->el(source_dof) ;

								rows_values.push_back(value);	   
								rows.push_back(source_dof);
			
            }

            
            source_dof = node->dof_number(last_non_linear_soln.number(), y_var, 0);
            if((source_dof<12345678) && (source_dof>-1) && (p(1)<0.0001 || p(1)>0.999)){
							 Real value = last_non_linear_soln.current_local_solution->el(source_dof) ;
							
								rows_values.push_back(value);
								rows.push_back(source_dof);
						
            }
            
            
            
           
                    
      
      
            source_dof = node->dof_number(last_non_linear_soln.number(), z_var, 0);
            if((source_dof<12345678) && (source_dof>-1) && (p(2)<0.001 ||p(2)>0.999 ) ){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof);
							
								rows_values.push_back(value );
								rows.push_back(source_dof);
				
            }
           
            
            
            
     





          }





        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


