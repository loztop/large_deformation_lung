//The test probem is:
//Fix cube at the bottom in all directions.
//Then apply a force w.r.t the normal in the current configuration from the left.
//Fluid can flow out everywa

for (unsigned int s=0; s<elem->n_sides(); s++)
{
   if (elem->neighbor(s) == NULL)
   {   
			AutoPtr<Elem> side (elem->build_side(s));
			
			//Sort out fixed(Dirichlet) Bcs
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
	
					//Constrain all displacements at  Z= 0 
					if ((elem->node(n) == side->node(ns)) && ( (p(2)<0.001 ) ) )
					{
						for (unsigned int d = 0; d < 3; ++d) {
							unsigned int source_dof = node->dof_number(1, d, 0);
							Real value = last_non_linear_soln.current_local_solution->el(source_dof) - ref_sys.current_local_solution->el(source_dof);
							
							rows.push_back(source_dof);
							rows_values.push_back(value);
							
						}
					} 
					
				} //end nodes in element lopp
			} // end nodes on side loop

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
 
					
			fe_face_map->init_face_shape_functions<3>(qface_point_f, &(*side) );	
			//fe_face_map->compute_edge_map(3,JxW_face_f, &(*side) );				
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
				
					
			if(  rX(2)>0.99  ){
					
				Point dxyzdxi_qp;
				Point dxyzdeta_qp;
				for (unsigned int i=0; i<dpsidxi.size(); i++) // sum over the nodes
				{
					const Point& side_point = side->point(i);
					dxyzdxi_qp+=side_point*dpsidxi[i][qp];
					dxyzdeta_qp+=side_point*dpsideta[i][qp];
				}
				
				Point nhat_qp  = dxyzdxi_qp.cross(dxyzdeta_qp);
				
				
				/*
				Point dphidxi_p;
				Point dphideta_p;
				dphidxi_p(0)=dxyzdxi[qp](0);
				dphidxi_p(1)=dxyzdxi[qp](1);
				dphidxi_p(2)=dxyzdxi[qp](2);

			 	dphideta_p(0)=dxyzdeta[qp](0);
				dphideta_p(1)=dxyzdeta[qp](1);
				dphideta_p(2)=dxyzdeta[qp](2);
				
				Point phinhat_qp  = dphidxi_p.cross(dphideta_p);
				std::cout<< "phinhat_qp.unit() "<< phinhat_qp.unit() <<std::endl;
				*/
		     
				side->print_info();

				std::cout<< "nhat_qp.unit() "<< nhat_qp.unit() <<std::endl;
				std::cout<< "face_map_normals[qp] "<< face_map_normals[qp] <<std::endl;
			 	std::cout<< "face_normals_f[qp] "<< face_normals_f[qp] <<std::endl;
				
				
				//Forcing traction (on solid)
				Real tn=0;
				
				Real pn=0*KPERM*2000;

				//Add the linearization of the traction
				int m=0;
				Real Fuphi;
				Real Fupsi;
							
				for (unsigned int i=0; i<phi_face_f.size(); i++){				
					

					if ( phi_face_f[i][qp]*phi_face_f[i][qp] < 0.0000001){
							//do nothing
						}
						else{
							
					
						 
							Fu(i) +=  tn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0);
							Fv(i) +=  tn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1);				
							Fw(i) +=  tn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2);
						 

							Fx(i) +=  pn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0);
							Fy(i) +=  pn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1);		
							Fz(i) +=  pn*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2);
							 

							int n=0;
							for (unsigned int j=0; j<phi_face_f.size(); j++){		
								if ( phi_face_f[j][qp]*phi_face_f[j][qp] < 0.0000001){
									//do nothing
								}
								else{
 
							//Linearization "blocks" from inside the brackets in (4.140 in Wriggers)						
						  Real lin_fac= 0*tn;

							Kuu(i,j) += lin_fac*phi_face_f[i][qp]*0; 
							Kuv(i,j) += lin_fac*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](2)-dphideta[n][qp]*dxyzdxi[qp](2));
							Kuw(i,j) += lin_fac*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](1)+dphideta[n][qp]*dxyzdxi[qp](1)); 
							
							Kvu(i,j) += lin_fac*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](2)+dphideta[n][qp]*dxyzdxi[qp](2)); 
							Kvv(i,j) += lin_fac*phi_face_f[i][qp]*0; 
							Kvw(i,j) += lin_fac*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](0)-dphideta[n][qp]*dxyzdxi[qp](0)); 
							
							Kwu(i,j) += lin_fac*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](1)-dphideta[n][qp]*dxyzdxi[qp](1)); 
							Kwv(i,j) += lin_fac*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](0)+dphideta[n][qp]*dxyzdxi[qp](0)); 
							Kww(i,j) += lin_fac*phi_face_f[i][qp]*0;
							
							
							Real lin_fac_p=1*pn;
							
							Kxu(i,j) += lin_fac_p*phi_face_f[i][qp]*0; 
							Kxv(i,j) += lin_fac_p*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](2)-dphideta[n][qp]*dxyzdxi[qp](2));
							Kxw(i,j) += lin_fac_p*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](1)+dphideta[n][qp]*dxyzdxi[qp](1)); 
							
							Kyu(i,j) += lin_fac_p*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](2)+dphideta[n][qp]*dxyzdxi[qp](2)); 
							Kyv(i,j) += lin_fac_p*phi_face_f[i][qp]*0; 
							Kyw(i,j) += lin_fac_p*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](0)-dphideta[n][qp]*dxyzdxi[qp](0)); 
							
							Kzu(i,j) += lin_fac_p*phi_face_f[i][qp]*(dphidxi[n][qp]*dxyzdeta[qp](1)-dphideta[n][qp]*dxyzdxi[qp](1)); 
							Kzv(i,j) += lin_fac_p*phi_face_f[i][qp]*(-dphidxi[n][qp]*dxyzdeta[qp](0)+dphideta[n][qp]*dxyzdxi[qp](0)); 
							Kzw(i,j) += lin_fac_p*phi_face_f[i][qp]*0;
							
							//Darcy traction linearization
							
							
								
		         	n=n+1;
						}
		       }
		        m=m+1;
		      } //end of i if
				}
			}    		  
			
		} //end qp
		
		
	} //if (elem->neighbor(s) == NULL)
}// end boundary condition section  


