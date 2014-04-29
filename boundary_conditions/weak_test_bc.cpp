//The test probem is:
//Fix cube at the bottom in all directions.
//Then apply a force w.r.t the normal in the current configuration from the left.
//Fluid can flow out everywa

for (unsigned int s=0; s<elem->n_sides(); s++)
{
   if (elem->neighbor(s) == NULL)
   {   
			AutoPtr<Elem> side (elem->build_side(s));
			
			

			//pressure face 
		  const std::vector<std::vector<Real> >&  phi_face_p =  fe_face_p->get_phi();
  		const std::vector<std::vector<RealGradient> >& dphi_face_p = fe_face_p->get_dphi();
  		const std::vector<Real>& JxW_face_p = fe_face_p->get_JxW();
  		const std::vector<Point>& qface_point_p = fe_face_p->get_xyz();
  		const std::vector<Point>& face_normals_p = fe_face_p->get_normals();
  		fe_face_p->reinit(elem,s); 
		
			
  		fe_face_f->reinit(elem,s);  

      //fluid face 
  		const std::vector<std::vector<Real> >&  phi_face_f =  fe_face_f->get_phi();
  		const std::vector<std::vector<RealGradient> >& dphi_face_f = fe_face_f->get_dphi();
  		const std::vector<Real>& JxW_face_f = fe_face_f->get_JxW();
  		const std::vector<Point>& qface_point_f = fe_face_f->get_xyz();
  		const std::vector<Point>& face_normals_f = fe_face_f->get_normals();

			//Reference element basis functions ?
			const std::vector<RealGradient>& dxyzdxi = fe_face_f->get_dxyzdxi();
			const std::vector<RealGradient>& dxyzdeta = fe_face_f->get_dxyzdeta();
			const std::vector<RealGradient>& dxyzdzeta = fe_face_f->get_dxyzdzeta();

 
			const std::vector<std::vector<Real> >&  dphidxi =  fe_face_f->get_dphidxi();
			const std::vector<std::vector<Real> >&  dphideta =  fe_face_f->get_dphideta();
			const std::vector<std::vector<Real> >&  dphidzeta =  fe_face_f->get_dphidzeta();
 
						//fe_face_map->compute_edge_map(3,qp, &(*side) );				

			fe_face_map->init_face_shape_functions<3>(qface_point_f, &(*side) );	
			fe_face_map->compute_edge_map(3,JxW_face_f, &(*side) );				
			fe_face_map->compute_face_map(3,JxW_face_f, &(*side) );						 	
			const std::vector<std::vector<Real> >&  face_psi =  fe_face_map->get_psi();
			const std::vector<std::vector<Real> >&  dpsidxi =  fe_face_map->get_dpsidxi();
			const std::vector<std::vector<Real> >&  dpsideta =  fe_face_map->get_dpsideta();
			const std::vector<Point >& qface_point_map = fe_face_map->get_xyz();
  		const std::vector<Real>& JxW_face_map = fe_face_map->get_JxW();
			const std::vector<Point>& face_map_normals = fe_face_map->get_normals();
			
			for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
			{

				Point rX;
				for (unsigned int l=0; l<phi_face_f.size(); l++)			
				{
					rX(0) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
					rX(1) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_v[l]);
					rX(2) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_w[l]);
				}
				
				//Calculate old fluid velocity
				Point fluid_vel;
				for (unsigned int l=0; l<phi_face_f.size(); l++)			
				{
					fluid_vel(0) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_x[l]);
					fluid_vel(1) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_y[l]);
					fluid_vel(2) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_z[l]);
				}

					
					
					
				if(  rX(2)>0.9999  ){
					
					Point dxyzdxi_qp;
					Point dxyzdeta_qp;
					
					Point phidxyzdxi_qp;
					Point phidxyzdeta_qp;
					
					for (unsigned int i=0; i<dpsidxi.size(); i++) // sum over the nodes
					{
						const Point& side_point = side->point(i);
						
						dxyzdxi_qp+=side_point*dpsidxi[i][qp];
						dxyzdeta_qp+=side_point*dpsideta[i][qp];
						
						phidxyzdxi_qp+=side_point*dphidxi[i][qp];
						phidxyzdeta_qp+=side_point*dphideta[i][qp];

					}
					

					Point nhat_qp  = dxyzdxi_qp.cross(dxyzdeta_qp);		
					Point phinhat_qp  = phidxyzdxi_qp.cross(phidxyzdeta_qp);

				
						
					Point normal;
					normal(0)=face_normals_f[qp](0);
					normal(1)=face_normals_f[qp](1);
					normal(2)=face_normals_f[qp](2);
					normal=normal.unit();
					
					//std::cout<< normal <<std::endl;

					Real pen_bc=100;
					
				//std::cout<< face_normals_f[qp] <<std::endl;
				
					for (unsigned int i=0; i<phi_face_f.size(); i++){			

						/*
						//Works
						Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(0);
				//		Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(1);
				//		Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(2);
						
				//		Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(0);
						Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(1);
				//		Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(2);
						
				//		Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(0);
				//		Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(1);
						Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(2);
						*/
		
						
						Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(0)*phi_face_f[i][qp]*normal(0);
						Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(0)*phi_face_f[i][qp]*normal(1);
						Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(0)*phi_face_f[i][qp]*normal(2);
						
						Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(1)*phi_face_f[i][qp]*normal(0);
						Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(1)*phi_face_f[i][qp]*normal(1);
						Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(1)*phi_face_f[i][qp]*normal(2);
						
						Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*normal(2)*phi_face_f[i][qp]*normal(0);
						Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*normal(2)*phi_face_f[i][qp]*normal(1);
						Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*normal(2)*phi_face_f[i][qp]*normal(2);
						
						for (unsigned int j=0; j<phi_face_f.size(); j++){				
									
												
							
							Kxx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(0); 
							Kxy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(1); 
							Kxz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(0)*phi_face_f[j][qp]*normal(2); 
						

							Kyx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(0); 
							Kyy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(1); 
							Kyz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(1)*phi_face_f[j][qp]*normal(2); 
						
							
							Kzx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(0); 
							Kzy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(1); 
							Kzz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*normal(2)*phi_face_f[j][qp]*normal(2);
							
						}
					
						
					}
					
			}    		  
			
		} //end qp
		
		
	} //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


