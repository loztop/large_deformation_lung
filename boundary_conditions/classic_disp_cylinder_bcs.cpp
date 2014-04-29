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
					
					
					Real xf = p(0);
          Real yf = p(1);
          Real zf = p(2);
            
//constrain bottom in z direction - let rest slide
          if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) ))
          {
            int source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            

            source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);
             rows.push_back(source_dof);
            }
            
          }  //end if


//Fix one node
//          if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) && (xf<0.1) && (xf>-0.1) && (yf<0.1) && (yf>-0.1) ))
       if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) &&  (xf<-0.99) &&  (yf<0.1)  &&  (yf>-0.1) ))
        {

            int source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
						
            if((source_dof<12345678) && (source_dof>-1)){
						  Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);
             rows.push_back(source_dof);
            }

            source_dof = node->dof_number(last_non_linear_soln.number(), v_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
							rows_values.push_back(value);            
							rows.push_back(source_dof);
            }

             source_dof = node->dof_number(last_non_linear_soln.number(), u_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);
						 rows.push_back(source_dof);
            }
          }  //end if
 

 
//Control another node
//          if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) && (xf<0.1) && (xf>-0.1) && (yf<0.1) && (yf>-0.1) ))
       if ((elem->node(n) == side->node(ns)) && ( (zf<0.001) &&  (xf>0.99) &&  (yf<0.1)  &&  (yf>-0.1) ))
        {

           //std::cout<< (*node) <<std::endl;

            int source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);
						 rows.push_back(source_dof);
            }
            
            source_dof = node->dof_number(last_non_linear_soln.number(), v_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);             
						 rows.push_back(source_dof);
            }
         
          }  //end if
          
   

          //Compress the top
          if ((elem->node(n) == side->node(ns)) && (  (zf>0.999) ))
          {
                  
            int source_dof = node->dof_number(last_non_linear_soln.number(), w_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof)+0.01;
             rows_values.push_back(value);
						 rows.push_back(source_dof);
            }
           }  //end if


           
          //Restrict inflow at top and bottom
          if ((elem->node(n) == side->node(ns)) && (  (zf<0.001) || (zf>0.999) ))
          {
                  
            int source_dof = node->dof_number(last_non_linear_soln.number(), z_var, 0);
            if((source_dof<12345678) && (source_dof>-1)){
							
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
             rows_values.push_back(value);
						 rows.push_back(source_dof);
            }
           }  //end if
           
           
          }
        }
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


