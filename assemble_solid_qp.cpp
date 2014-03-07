for (unsigned int qp=0; qp<qrule.n_points(); qp++)
  {

  Point rX;
  for (unsigned int l=0; l<n_u_dofs; l++)
  {
    rX(0) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_u[l]);
    rX(1) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_v[l]);
    rX(2) += phi[l][qp]*ref_sys.current_local_solution->el(dof_indices_w[l]);
  }

  Number   p_solid = 0.;

  grad_u_mat(0) = grad_u_mat(1) = grad_u_mat(2) = 0;
  for (unsigned int d = 0; d < dim; ++d) {
    std::vector<Number> u_undefo;
    //Fills the vector di with the global degree of freedom indices for the element. :dof_indicies
    aux_system.get_dof_map().dof_indices(elem, undefo_index,d);
    aux_system.current_local_solution->get(undefo_index, u_undefo);
    for (unsigned int l = 0; l != n_u_dofs; l++){
      grad_u_mat(d).add_scaled(dphi[l][qp], u_undefo[l]); 
    }
  }
          
  for (unsigned int l=0; l<n_p_dofs; l++)
  {
    p_solid += psi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
  }
    
  material.init_for_qp(rX,grad_u_mat, p_solid, qp,0, p_solid);

  Real Kperm=material.Kperm;
  Real J=material.J;
  Real phi_zero=material.Phi_zero;
  Real rho_s=material.Rho_s;
  Real PHI=1-((1-phi_zero)/J);
// 
  Point fluid_vel;
  for (unsigned int l=0; l<n_u_dofs; l++)
  {
    fluid_vel(0) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_x[l]);
    fluid_vel(1) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_y[l]);
    fluid_vel(2) += phi[l][qp]*last_non_linear_soln.current_local_solution->el(dof_indices_z[l]);
  }

  Point grad_p;
  for (unsigned int l=0; l<n_p_dofs; l++)
  {
    grad_p(0) += dpsi[l][qp](0)*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
    grad_p(1) += dpsi[l][qp](1)*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
    grad_p(2) += dpsi[l][qp](2)*last_non_linear_soln.current_local_solution->el(dof_indices_p[l]);
  }

  Point div_f_vel;
  for (unsigned int l=0; l<n_u_dofs; l++)
  {
    div_f_vel(0) += dphi[l][qp](0)*last_non_linear_soln.current_local_solution->el(dof_indices_x[l]);
    div_f_vel(1) += dphi[l][qp](1)*last_non_linear_soln.current_local_solution->el(dof_indices_y[l]);
    div_f_vel(2) += dphi[l][qp](2)*last_non_linear_soln.current_local_solution->el(dof_indices_z[l]);
  }

  Point div_s_vel;
  
  for (unsigned int l=0; l<n_u_dofs; l++)
  {
    div_s_vel(0) += dphi[l][qp](0)*(last_non_linear_soln.current_local_solution->el(dof_indices_u[l])-last_non_linear_soln.old_local_solution->el(dof_indices_u[l]))/dt;
    div_s_vel(1) += dphi[l][qp](1)*(last_non_linear_soln.current_local_solution->el(dof_indices_v[l])-last_non_linear_soln.old_local_solution->el(dof_indices_v[l]))/dt;
    div_s_vel(2) += dphi[l][qp](2)*(last_non_linear_soln.current_local_solution->el(dof_indices_w[l])-last_non_linear_soln.old_local_solution->el(dof_indices_w[l]))/dt;
  }

  for (unsigned int i=0; i<n_u_dofs; i++)
  {
    res.resize(dim);
    material.get_residual(res, i);
    res.scale(JxW[qp]);
    
    Fu(i) += res(0);              
    Fv(i) += res(1) ; 
    Fw(i) += res(2);  
 
    Real grav=9.8*rho_s;
 
      
    // Matrix contributions for the uu and vv couplings.
    for (unsigned int j=0; j<n_u_dofs; j++)
    {
      material.get_linearized_stiffness(stiff, i, j);
      stiff.scale(JxW[qp]);

      Kuu(i,j)+=  stiff(u_var, u_var);
      Kuv(i,j)+=  stiff(u_var, v_var);
      Kuw(i,j)+=  stiff(u_var, w_var);        
      Kvu(i,j)+=  stiff(v_var, u_var);
      Kvv(i,j)+=  stiff(v_var, v_var);
      Kvw(i,j)+=  stiff(v_var, w_var);
      Kwu(i,j)+=  stiff(w_var, u_var);
      Kwv(i,j)+=  stiff(w_var, v_var);
      Kww(i,j)+=  stiff(w_var, w_var); 

      #if GRAVITY
      Kuu(i,j)+= 1*JxW[qp]*phi[i][qp]*phi[j][qp];
      #endif
    }
  }

 
   for (unsigned int i = 0; i < n_u_dofs; i++) {
    for (unsigned int j = 0; j < n_p_dofs; j++) {
      Kup(i, j) += -JxW[qp]*dphi[i][qp](0)*psi[j][qp];
      Kvp(i, j) += -JxW[qp]*dphi[i][qp](1)*psi[j][qp];
      Kwp(i, j) += -JxW[qp]*dphi[i][qp](2)*psi[j][qp];
    }
  }
  
  Real factor=Kperm;

  //Mass Matrix needed for darcy flow
  for (unsigned int i=0; i<n_u_dofs; i++){
    for (unsigned int j=0; j<n_u_dofs; j++){
      //w.v term (u)
      Kxx(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
      Kyy(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
      Kzz(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);

	  /*
       // Laplacian for stokes flow /Brinkman
      Kxx(i,j) += brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
      Kyy(i,j) += brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
      Kzz(i,j) += brink_mu*JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
      */
    }

    for (unsigned int j=0; j<n_p_dofs; j++){
      // Weak form of grad p term
      Kxp(i,j) += -factor*JxW[qp]*dphi[i][qp](0)*psi[j][qp];
      Kyp(i,j) += -factor*JxW[qp]*dphi[i][qp](1)*psi[j][qp];
      Kzp(i,j) += -factor*JxW[qp]*dphi[i][qp](2)*psi[j][qp];
    }
  }

  // the fluid mass residual:
  for (unsigned int i = 0; i < n_p_dofs; i++) {
    Fp(i) += JxW[qp]*psi[i][qp]*(div_f_vel(0)+div_f_vel(1)+div_f_vel(2)+  div_s_vel(0)+div_s_vel(1)+div_s_vel(2) );
  }
 
  //This has the boundary integral missing + brinkmann term
  for (unsigned int i = 0; i < n_u_dofs; i++) {

    //The Fluid residual
    //This has the Fp_soli boundary term missing ?

    Fx(i) += JxW[qp]*phi[i][qp]*fluid_vel(0);
    Fx(i) += -JxW[qp]*Kperm*p_solid*dphi[i][qp](0);

    Fy(i) += JxW[qp]*phi[i][qp]*fluid_vel(1);
    Fy(i) += -JxW[qp]*Kperm*p_solid*dphi[i][qp](1);

    Fz(i) += JxW[qp]*phi[i][qp]*fluid_vel(2);
    Fz(i) += -JxW[qp]*Kperm*p_solid*dphi[i][qp](2);
			
		//std::cout<< "i " << i << std::endl;
		//		std::cout<< "phi[i][qp] " <<phi[i][qp] << std::endl;
    }

  //Mass conservation
  for (unsigned int i=0; i<n_u_dofs; i++){
    for (unsigned int j=0; j<n_p_dofs; j++){
      //the fluid velocity divergence term 
      Kpx(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
      Kpy(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
      Kpz(j,i) += JxW[qp]*psi[j][qp]*dphi[i][qp](2);   
    }
  }

  for (unsigned int i=0; i<n_p_dofs; i++){
    for (unsigned int j=0; j<n_u_dofs; j++){
      //the solid velocity divergence term 
      Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0)/dt;
      Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1)/dt;
      Kpw(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](2)/dt;   
    }
  }

}//end of qp loop