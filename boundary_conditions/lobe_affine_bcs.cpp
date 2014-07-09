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
            const Real xf = (*node)(0);
            Real yf = (*node)(1);
            Real zf = (*node)(2);
            
						Real pi=3.14159;
						Real fac=0.1*(1+sin(1*pi*time+(3.0/2.0)*pi));

						// fac=0;

						Point diagAt(fac*0.3,fac*0.3,fac*0.5);						
						Point bt(fac*10,fac*11,fac*45);
						//Point bt(0,0,0);


					//		Point diagAt(fac*0.3,fac*0.3,fac*0.5);						
						
		  Point p;
			for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
      	Real value = ref_sys.current_local_solution->el(source_dof);
      	p(d)=value;
			}
						
				//Dilate boundary
       //  if(p(2)>100)
					{



            int source_dof = node->dof_number(last_non_linear_soln.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
									Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);

								rows_values.push_back(value-xf*diagAt(0)+bt(0));
											   
								/*
								std::cout<<  fac << std::endl;

			      std::cout<<  xf << std::endl;
			    std::cout<<  diagAt(0) << std::endl;
			    std::cout<<  bt(0) << std::endl;

																
								std::cout<<  xf*diagAt(0)+bt(0) << std::endl;
								*/
								
								rows.push_back(source_dof);
			
            }

            source_dof = node->dof_number(last_non_linear_soln.number(), v_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							 Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
							
								rows_values.push_back(value-yf*diagAt(1)+bt(1));
								rows.push_back(source_dof);
						
            }

            source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
							
								rows_values.push_back(value-zf*diagAt(2)+bt(2));
								rows.push_back(source_dof);
				
            }
            
            
            
          }  






          }





        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


