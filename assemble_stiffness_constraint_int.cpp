for (unsigned int s=0; s<elem->n_sides(); s++)
{
	//Should be integrating instead !!! to avoid problem with normals !!
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
			
			
		//	(*side).print_info(std::cout);
		//  (*side)->print_info();
	    // 	std::cout<<	side->n_nodes() <<std::endl;

			
			for (unsigned int qp=0; qp<qface_f->n_points(); qp++)
			{

					Real xf=qface_point_f[qp](0);
					Real yf=qface_point_f[qp](1);
					
					Point normal;
					normal(0)=face_normals_f[qp](0);
					normal(1)=face_normals_f[qp](1);

					Real zf=qface_point_f[qp](2);
					normal(2)=face_normals_f[qp](2);
										
					normal=normal.unit();
					
					Real pen_bc=DELTA_BC/hmax;
	
		 // if(yf > 0.999 || xf > 0.999   ){
			//std::cout<<	normal <<std::endl;
			for (unsigned int i=0; i<phi_face_f.size(); i++){			
			  for (unsigned int j=0; j<phi_face_f.size(); j++){																	
				Kcx(i,j) += normal(0)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 
				Kcy(i,j) += normal(1)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 		
				Kcz(i,j) += normal(2)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 	
			  
				Kxc(i,j) += normal(0)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 
				Kyc(i,j) += normal(1)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 			
			    Kzc(i,j) += normal(2)*JxW_face_f[qp]*phi_face_f[i][qp]*phi_face_f[j][qp]; 				
			  }					
			}
			//Enforce z.n =0

						
		//	}

		
		//assemble residal
		for (unsigned int i=0; i<phi_face_f.size(); i++){			
		  Fc(i)
		  
		}								
		
		
	} //end qp
		
		
	} 
}// end boundary condition section  

	 // system.matrix->close();
   // system.rhs->close();
