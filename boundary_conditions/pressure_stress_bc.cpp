for (unsigned int s=0; s<elem->n_sides(); s++)
{
   if (elem->neighbor(s) == NULL)
    {   
        AutoPtr<Elem> side (elem->build_side(s));

        //fluid face 
  			const std::vector<std::vector<Real> >&  phi_face_f =  fe_face_f->get_phi();
  			const std::vector<std::vector<RealGradient> >& dphi_face_f = fe_face_f->get_dphi();
  			const std::vector<Real>& JxW_face_f = fe_face_f->get_JxW();
  			const std::vector<Point>& qface_point_f = fe_face_f->get_xyz();
  			const std::vector<Point>& face_normals_f = fe_face_f->get_normals();
  			fe_face_f->reinit(elem,s);  

			//pressure face 
		    const std::vector<std::vector<Real> >&  phi_face_p =  fe_face_p->get_phi();
  			const std::vector<std::vector<RealGradient> >& dphi_face_p = fe_face_p->get_dphi();
  			const std::vector<Real>& JxW_face_p = fe_face_p->get_JxW();
  			const std::vector<Point>& qface_point_p = fe_face_p->get_xyz();
  			const std::vector<Point>& face_normals_p = fe_face_p->get_normals();
  			fe_face_p->reinit(elem,s); 
				
				
				//Reference element basis functions
				
				const std::vector<RealGradient>& dxyzdxi = fe_face_f->get_dxyzdxi();
				const std::vector<RealGradient>& dxyzdeta = fe_face_f->get_dxyzdeta();
				const std::vector<RealGradient>& dxyzdzeta = fe_face_f->get_dxyzdzeta();

    for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
	  {

	  
			//Calculate normal w.r.t the reference element (Wriggers eqn. 4.136)
			//(We use notation (x,y,z)=(x_1,x_2,x_3)
			//n = varphi_{xi} x varphi_{eta} = [y_{xi}.z_{eta} - z_{xi}.y_{eta}]
			//                                 [z_{xi}.x_{eta} - x_{xi}.z_{eta}]
			//                                 [x_{xi}.y_{eta} - y_{xi}.x_{eta}]	
			
			
			
			Point normal_ref;		
			normal_ref(0)=dxyzdxi[qp](1)*dxyzdeta[qp](2) - dxyzdxi[qp](2)*dxyzdeta[qp](1);
			normal_ref(1)=dxyzdxi[qp](2)*dxyzdeta[qp](0) - dxyzdxi[qp](0)*dxyzdeta[qp](2);
			normal_ref(2)=dxyzdxi[qp](0)*dxyzdeta[qp](1) - dxyzdxi[qp](1)*dxyzdeta[qp](0);

			/*
			Real norm_normal;
			norm_normal= pow(normal_ref(0)*normal_ref(0) + normal_ref(1)*normal_ref(1) +normal_ref(2)*normal_ref(2) ,0.5);

			std::cout<< "JxW_face_f[qp] " << JxW_face_f[qp] <<std::endl;
			std::cout<< "JxW_face_f[qp]*face_normals_f[qp] " << JxW_face_f[qp]*face_normals_f[qp] <<std::endl;
			std::cout<< "normal_ref " << normal_ref <<std::endl;
					
			//Calculate linearization
			//Building the submatrix kab 4.140 in Wriggers
			//N,alpha (where a={xi,eta})  =  [0      z_{a}  -y_{a} ]
			//															 [-z_{a} 0      x_{a}  ]
			//															 [y_{a}  -x_{a} 0      ]
					
			RealTensor Nxi;
			Nxi(0,0)= 0;
			Nxi(0,1)= dxyzdxi[qp](2);
			Nxi(0,2)= -dxyzdxi[qp](1);
			
			Nxi(1,0)= -dxyzdxi[qp](2);
			Nxi(1,1)= 0;
			Nxi(1,2)= dxyzdxi[qp](0);

			Nxi(2,0)= dxyzdxi[qp](1);
			Nxi(2,1)= -dxyzdxi[qp](0);
			Nxi(2,2)= 0;
						
			RealTensor Neta;
			Neta(0,0)= 0;
			Neta(0,1)= dxyzdeta[qp](2);
			Neta(0,2)= -dxyzdeta[qp](1);
			
			Neta(1,0)= -dxyzdeta[qp](2);
			Neta(1,1)= 0;
			Neta(1,2)= dxyzdeta[qp](0);

			Neta(2,0)= dxyzdeta[qp](1);
			Neta(2,1)= -dxyzdeta[qp](0);
			Neta(2,2)= 0;
			*/
			
			
			Point lin_A;		
			lin_A(0)=-dxyzdxi[qp](1)*dxyzdeta[qp](2) + dxyzdxi[qp](2)*dxyzdeta[qp](1);
			lin_A(1)= dxyzdxi[qp](0)*dxyzdeta[qp](2) - dxyzdxi[qp](2)*dxyzdeta[qp](0) ;
			lin_A(2)= -dxyzdxi[qp](0)*dxyzdeta[qp](1) + dxyzdxi[qp](1)*dxyzdeta[qp](0);
			
			Point lin_B;		
			lin_B(0)=-dxyzdeta[qp](1)*dxyzdxi[qp](2) + dxyzdeta[qp](2)*dxyzdxi[qp](1);
			lin_B(1)= dxyzdeta[qp](0)*dxyzdxi[qp](2) - dxyzdeta[qp](2)*dxyzdxi[qp](0) ;
			lin_B(2)= -dxyzdeta[qp](0)*dxyzdxi[qp](1) + dxyzdeta[qp](1)*dxyzdxi[qp](0);
			
	  //Get previous solution to build the residual
	  Point fluid_vel;
		for (unsigned int l=0; l<phi_face_f.size(); l++)			
		{
		  fluid_vel(0) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_x[l]);
		  fluid_vel(1) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_y[l]);
		  fluid_vel(2) += phi_face_f[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_z[l]);
		}
	  
		Real p_fluid;
		for (unsigned int l=0; l<phi_face_p.size(); l++)			
		{
		  p_fluid += phi_face_p[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
		}
	  
	  Point rX;
	  for (unsigned int l=0; l<phi_face_f.size(); l++)			
		{
		rX(0) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
		rX(1) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_v[l]);
		rX(2) += phi_face_f[l][qp]*ref_sys.current_local_solution->el(dof_indices_w[l]);
		}

	  
	  Real pen_bc=1000;
	  
		Real xf = qface_point_f[qp](0);
    Real yf = qface_point_f[qp](1);
    Real zf = qface_point_f[qp](2);

            // RHS contributions
			if(  rX(0)<0.01  ){

			Real p_d=0.5*progress;

				for (unsigned int i=0; i<phi_face_f.size(); i++){				
					
					Fx(i) += - p_d*JxW_face_f[qp]*phi_face_f[i][qp]; 	
					

				}		
     
			}
		      
		      //Null flux, this is Gamma_D
		      if(  (rX(1)<0.001) ||  (rX(1)>0.999)  ){

				  for (unsigned int i=0; i<phi_face_f.size(); i++){				
					for (unsigned int j=0; j<phi_face_f.size(); j++){	
						
						
						
						
					  
/*
					  Kxx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0)*phi_face_f[j][qp]*face_normals_f[qp](0); 
					  Kxy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0)*phi_face_f[j][qp]*face_normals_f[qp](1); 
					  Kxz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](0)*phi_face_f[j][qp]*face_normals_f[qp](2); 

					  Kyx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](0); 
					  Kyy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](1); 
					  Kyz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](2); 
					  
					  Kzx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2)*phi_face_f[j][qp]*face_normals_f[qp](0); 
					  Kzy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2)*phi_face_f[j][qp]*face_normals_f[qp](1); 
					  Kzz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](2)*phi_face_f[j][qp]*face_normals_f[qp](2); 
					
					  
					// Kyy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 

					//  Kyx(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](0); 
					//  Kyy(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](1); 
					//  Kyz(i,j) += pen_bc*JxW_face_f[qp]*phi_face_f[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]*face_normals_f[qp](2); 
					}
					
					//Build residual
					
				
					Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*face_normals_f[qp](0)*phi_face_f[i][qp]*face_normals_f[qp](0); 	
					Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*face_normals_f[qp](0)*phi_face_f[i][qp]*face_normals_f[qp](1); 	
					Fx(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(0)*face_normals_f[qp](0)*phi_face_f[i][qp]*face_normals_f[qp](2); 	

					Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](0); 	
					Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](1); 	
					Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](2); 	
					
					Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*face_normals_f[qp](2)*phi_face_f[i][qp]*face_normals_f[qp](0); 	
					Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*face_normals_f[qp](2)*phi_face_f[i][qp]*face_normals_f[qp](1); 	
					Fz(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(2)*face_normals_f[qp](2)*phi_face_f[i][qp]*face_normals_f[qp](2); 	
						
					

				//	Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*phi_face_f[i][qp]; 	
				//	Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](0); 	
				//	Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](1); 	
				//	Fy(i) +=  pen_bc*JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_f[i][qp]*face_normals_f[qp](2);
				*/
				
				  }
			
			
			
			  /*
				  // (p,v.n)=(v.n,p)
				  for (unsigned int i=0; i<phi_face_f.size(); i++){				
					for (unsigned int j=0; j<phi_face_p.size(); j++){	
					  
					  Kxp(i,j) += JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_p[j][qp]*face_normals_f[qp](0); 
					  Kyp(i,j) += JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_p[j][qp]*face_normals_f[qp](1); 
					  Kzp(i,j) += JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_p[j][qp]*face_normals_f[qp](2); 
 
					}
				  }
				  
				for (unsigned int i=0; i<phi_face_f.size(); i++){	

				  	//Residual
					 Fx(i) +=  JxW_face_f[qp]*p_fluid*phi_face_f[i][qp]*face_normals_f[qp](0); 	
					 Fy(i) +=  JxW_face_f[qp]*p_fluid*phi_face_f[i][qp]*face_normals_f[qp](1); 	
					 Fz(i) +=  JxW_face_f[qp]*p_fluid*phi_face_f[i][qp]*face_normals_f[qp](2); 	
				}
				  */
				  
				  
			  /*
				  //(u.n,q)=(q,u.n)
					for (unsigned int i=0; i<phi_face_p.size(); i++){	
					  	for (unsigned int j=0; j<phi_face_f.size(); j++){				
						  Kpx(i,j) += JxW_face_p[qp]*phi_face_p[i][qp]*face_normals_f[qp](0)*phi_face_f[j][qp]; 
						  Kpy(i,j) += JxW_face_f[qp]*phi_face_p[i][qp]*face_normals_f[qp](1)*phi_face_f[j][qp]; 
						  Kpz(i,j) += JxW_face_f[qp]*phi_face_p[i][qp]*face_normals_f[qp](2)*phi_face_f[j][qp]; 			
						}							
					 }
				  
				  for (unsigned int i=0; i<phi_face_p.size(); i++){				
				  //Residual
					Fp(i) +=  JxW_face_f[qp]*fluid_vel(0)*face_normals_f[qp](0)*phi_face_p[i][qp]; 	
					Fp(i) +=  JxW_face_f[qp]*fluid_vel(1)*face_normals_f[qp](1)*phi_face_p[i][qp]; 	
					Fp(i) +=  JxW_face_f[qp]*fluid_vel(2)*face_normals_f[qp](2)*phi_face_p[i][qp];
				  }
			  */
			  
			  
			  					//  std::cout<< JxW_face_p[qp]  <<std::endl;
			  					//  std::cout<< JxW_face_f[qp]  <<std::endl;

			   }
			    		  
		
			
	} //end qp
   } //if (elem->neighbor(s) == NULL)
}// end boundary condition section  

}
