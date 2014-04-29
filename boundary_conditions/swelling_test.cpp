for (unsigned int s=0; s<elem->n_sides(); s++){
   if (elem->neighbor(s) == NULL)
   {		
   AutoPtr<Elem> side (elem->build_side(s));

    for (unsigned int ns=0; ns<side->n_nodes(); ns++)
    {

     for (unsigned int n=0; n<elem->n_nodes(); n++){
       Node *node = elem->get_node(n);
	Point p;
	for (unsigned int d = 0; d < 3; ++d) {
      	unsigned int source_dof = node->dof_number(1, d, 0);
      	Real value = ref_sys.current_local_solution->el(source_dof);
      	p(d)=value;
   	}
	
	//Constrain normal displacements at X=0, Y=0 and Z= 0

	if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 )  )
    {
		unsigned int source_dof = node->dof_number(1, 0, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
     } 
     
    if ((elem->node(n) == side->node(ns)) && (  (p(1)<0.001) || (p(1)>0.999) ) )
    {
		unsigned int source_dof = node->dof_number(1, 1, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
     } 
     
    if ((elem->node(n) == side->node(ns)) && (p(2)<0.001 )  )
    {
		unsigned int source_dof = node->dof_number(1, 2, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
     } 

     
     /*
     //Null flux condition at Y=0,Y=1,Z=0,Z=1
     if ((elem->node(n) == side->node(ns)) && (p(1)<0.001 || p(1)>0.999 )  )
     {
		for (unsigned int d = 0; d < 3; ++d) {

		unsigned int source_dof = node->dof_number(1, 4+d, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
		
		}
     } 
     
     */
     
     if ((elem->node(n) == side->node(ns)) && (p(2)<0.001 || p(2)>0.999 )  )
     {
	   	for (unsigned int d = 0; d < 3; ++d) {

		unsigned int source_dof = node->dof_number(1, 4+d, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
		
		}
     } 
     
     
     /*
     if ((elem->node(n) == side->node(ns)) && (p(0)<0.001 || p(0)>0.999 )  )
     {
		unsigned int source_dof = node->dof_number(1, 4, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
		rows.push_back(source_dof);
		rows_values.push_back(value);
     } 
     */
     
	 /*
     //p = 0 at X=1  
     if ((elem->node(n) == side->node(ns)) && ( p(0)>0.999) )
     {
		int source_dof = elem->dof_number(1, p_var, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - 0;
        if((source_dof<12345678) && (source_dof>-1) ){
		  pressure_rows_values.push_back(value);
          pressure_rows.push_back(source_dof);
        }
      }  //end if
     */
     
     /*
	  //The pressure p bc on the inlet face X = 0 increases from zero to a limit value of 1kPa (p bc = 10 3 (1 − exp(−t 2 /0.25)) Pa)
     if ((elem->node(n) == side->node(ns)) && ( p(0)<0.001) )
     {
	   
		//Real p_value = 2000*(1-exp(-pow(((time+0.00000001)/1.0),2.0)/0.25)) ;
	   	Real p_value = 500*progress;
	   		//Real p_value = 1000 ;

		int source_dof = elem->dof_number(1, p_var, 0);
		Real value = last_non_linear_soln.current_local_solution->el(source_dof) - p_value;
        if((source_dof<12345678) && (source_dof>-1) ){
		  pressure_rows_values.push_back(value);
          pressure_rows.push_back(source_dof);
        }
      }  //end if
      */
      
    } //end nodes in element lopp
} // end nodes on side loop
  
 }//end elem->neighbor(s) == NULL 
 }// end s=0; s<elem->n_sides(); s++

