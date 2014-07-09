    for (unsigned int s=0; s<elem->n_sides(); s++)
    {
      if (elem->neighbor(s) == NULL)
      {   
        AutoPtr<Elem> side (elem->build_side(s));

        Real hmax=(*elem).hmax();
		
			
  		fe_face_f->reinit(elem,s);  

      //fluid face 
  		const std::vector<std::vector<Real> >&  phi_face_f =  fe_face_f->get_phi();
  		const std::vector<std::vector<RealGradient> >& dphi_face_f = fe_face_f->get_dphi();
  		const std::vector<Real>& JxW_face_f = fe_face_f->get_JxW();
  		const std::vector<Point>& qface_point_f = fe_face_f->get_xyz();
  		const std::vector<Point>& face_normals_f = fe_face_f->get_normals();

			
			//pressure face 
		  const std::vector<std::vector<Real> >&  phi_face_p =  fe_face_p->get_phi();
  		const std::vector<std::vector<RealGradient> >& dphi_face_p = fe_face_p->get_dphi();
  		const std::vector<Real>& JxW_face_p = fe_face_p->get_JxW();
  		const std::vector<Point>& qface_point_p = fe_face_p->get_xyz();
  		const std::vector<Point>& face_normals_p = fe_face_p->get_normals();
  		fe_face_p->reinit(elem,s); 
			
			
			for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
			{

				

	//Calculate old fluid velocity
				Point fluid_vel;
				for (unsigned int l=0; l<phi_face_f.size(); l++)			
				{
					fluid_vel(0) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_x[l]);
					fluid_vel(1) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_y[l]);
					fluid_vel(2) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_z[l]);
				}

									
						
						
					//		int source_dof_x = s->dof_number(last_non_linear_soln.number(), 4, 0);
		//	int source_dof_y = s->dof_number(last_non_linear_soln.number(), 5, 0);
 		//	int source_dof_z = s->dof_number(last_non_linear_soln.number(), 6, 0);

		
		//	fluid_vel(0)=last_non_linear_soln.current_local_solution->el(source_dof_x) - ref_sys.current_local_solution->el(source_dof_x); 
		//	fluid_vel(1)=last_non_linear_soln.current_local_solution->el(source_dof_y) - ref_sys.current_local_solution->el(source_dof_y);
		//	fluid_vel(2)=last_non_linear_soln.current_local_solution->el(source_dof_z) - ref_sys.current_local_solution->el(source_dof_z);
					
								Point fluid_vel;
	
	
		
				Point normal;
				normal(0)=face_normals_f[qp](0);
				normal(1)=face_normals_f[qp](1);
				normal(2)=face_normals_f[qp](2);
				
				normal=normal.unit();
	
	
				Real pen_bc=DELTA_BC/hmax;
					
				for (unsigned int i=0; i<phi_face_f.size(); i++){			
		 		

					Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(0);
				  Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(0);
					Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(0);
 						
				  Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(1);
				  Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(1);
					Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(1);
					
					
					Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(2);
				  Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(2);
					Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(2);

					
				}
				
				/*
		for (unsigned int i=0; i<phi_face_p.size(); i++){			
		 		

					Fp(i) +=  -JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_p[i][qp] ;
				 	Fp(i) +=  -JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_p[i][qp] ;
					Fp(i) +=  -JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_p[i][qp] ;

					
				}
				
				*/
				
		} //end qp
      } //if (elem->neighbor(s) == NULL)
    } // end boundary condition section  


