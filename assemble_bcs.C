#include "poro.h"

using namespace std;
#define PI 3.14159265

void assemble_bcs (EquationSystems& es)
{

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  TransientLinearImplicitSystem & newton_update =
  es.get_system<TransientLinearImplicitSystem> ("Newton-update");
  TransientLinearImplicitSystem & last_non_linear_soln =
  es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

  const System & ref_sys = es.get_system("Reference-Configuration"); 
  
	
	const Real DELTA_BC    = es.parameters.get<Real>("DELTA_BC");

  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = last_non_linear_soln.variable_number ("s_u");
  const unsigned int v_var = last_non_linear_soln.variable_number ("s_v");
  const unsigned int w_var = last_non_linear_soln.variable_number ("s_w");
  const unsigned int p_var = last_non_linear_soln.variable_number ("s_p");
 
  FEType fe_vel_type = last_non_linear_soln.variable_type(u_var);
  FEType fe_vel_type_ref = ref_sys.variable_type(u_var);
	

  const unsigned int dim = mesh.mesh_dimension();
  FEType fe_pres_type = last_non_linear_soln.variable_type(p_var);
  AutoPtr<FEBase> fe_vel  (FEBase::build(3, fe_vel_type));
  AutoPtr<FEBase> fe_pres (FEBase::build(3, fe_pres_type));
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
  
const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();

	
	AutoPtr<FEMap> fe_face_map (FEMap::build(fe_vel_type));

		
  std::vector< unsigned int > rows;
  std::vector< unsigned int > pressure_rows;
  std::vector< Real > rows_values;
  std::vector< Real > pressure_rows_values;


  const Real dt    = es.parameters.get<Real>("dt");
  const Real progress    = es.parameters.get<Real>("progress");
  const Real time    = es.parameters.get<Real>("time");

  //Build face

  AutoPtr<FEBase> fe_face_f (FEBase::build(3, fe_vel_type));      
	
	//AutoPtr<FEMap> fe_face_f_map (FEBase::build(3, fe_vel_type));          
		
  AutoPtr<QBase> qface_f(fe_vel_type.default_quadrature_rule(3-1));
  fe_face_f->attach_quadrature_rule (qface_f.get());

	
  AutoPtr<FEBase> fe_face_p (FEBase::build(3, fe_pres_type));          
  AutoPtr<QBase> qface_p(fe_pres_type.default_quadrature_rule(3-1));
  fe_face_p->attach_quadrature_rule (qface_p.get());

  AutoPtr<FEBase> fe_face_ref (FEBase::build(3, fe_vel_type_ref));          
  AutoPtr<QBase> qface_ref(fe_vel_type_ref.default_quadrature_rule(3-1));
  fe_face_ref->attach_quadrature_rule (qface_ref.get());
 
  const DofMap & dof_map = last_non_linear_soln.get_dof_map();

  // K will be the jacobian
  // F will be the Residual
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
  Kuu(Ke), Kuv(Ke), Kuw(Ke), 
  Kvu(Ke), Kvv(Ke), Kvw(Ke), 
  Kwu(Ke), Kwv(Ke), Kww(Ke); 
    
  DenseSubMatrix<Number>  Kup(Ke),Kvp(Ke),Kwp(Ke), Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke);
  
  DenseSubVector<Number>
  Fu(Fe),
  Fv(Fe),
  Fw(Fe);
  DenseSubVector<Number>    Fp(Fe);
  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_p;
 
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

  for ( ; el != end_el; ++el)
  {    
		const Elem* elem = *el;
    
		dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);
    dof_map.dof_indices (elem, dof_indices_p, p_var);

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size(); 
    const unsigned int n_v_dofs = dof_indices_v.size(); 
    const unsigned int n_w_dofs = dof_indices_w.size(); 
    const unsigned int n_p_dofs = dof_indices_p.size();

    fe_vel->reinit  (elem);

    Ke.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
    Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
    Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
    Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
    Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
    Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
          
    Kwu.reposition (w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
    Kwv.reposition (w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
    Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);
    Kwp.reposition (w_var*n_w_dofs, p_var*n_w_dofs, n_w_dofs, n_p_dofs);
    Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
    Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
    Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
    Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);

	Fu.reposition (u_var*n_u_dofs, n_u_dofs);
    Fv.reposition (v_var*n_u_dofs, n_v_dofs);
    Fw.reposition (w_var*n_u_dofs, n_w_dofs);
	Fp.reposition (p_var*n_u_dofs, n_p_dofs);

	const unsigned int x_var = last_non_linear_soln.variable_number ("mono_f_vel_u");
	const unsigned int y_var = last_non_linear_soln.variable_number ("mono_f_vel_v");
	const unsigned int z_var = last_non_linear_soln.variable_number ("mono_f_vel_w");

	//PF block
	DenseSubMatrix<Number>
	Kpx(Ke), Kpy(Ke), Kpz(Ke); 

	//L block
	DenseSubMatrix<Number>
	Kxu(Ke), Kxv(Ke), Kxw(Ke),
	Kyu(Ke), Kyv(Ke), Kyw(Ke),
	Kzu(Ke), Kzv(Ke), Kzw(Ke);

	//M block
	DenseSubMatrix<Number>
	Kxp(Ke), Kyp(Ke), Kzp(Ke); 

	//O block
	DenseSubMatrix<Number>
	Kxx(Ke), Kxy(Ke), Kxz(Ke),
	Kyx(Ke), Kyy(Ke), Kyz(Ke),
	Kzx(Ke), Kzy(Ke), Kzz(Ke);
	
	
	//Just creating this block for symmetry
	DenseSubMatrix<Number>
	Kux(Ke), Kvx(Ke), Kwx(Ke),
	Kuy(Ke), Kvy(Ke), Kwy(Ke),
	Kuz(Ke), Kvz(Ke), Kwz(Ke);

	DenseSubVector<Number>
	Fx(Fe), Fy(Fe), Fz(Fe);

	std::vector<unsigned int> dof_indices_x;
	std::vector<unsigned int> dof_indices_y;
	std::vector<unsigned int> dof_indices_z;

	//PF block
	Kpx.reposition (p_var*n_u_dofs, 3*n_u_dofs + n_p_dofs, n_p_dofs, n_u_dofs);
	Kpy.reposition (p_var*n_u_dofs, 4*n_u_dofs + n_p_dofs, n_p_dofs, n_v_dofs);
	Kpz.reposition (p_var*n_u_dofs, 5*n_u_dofs + n_p_dofs, n_p_dofs, n_w_dofs);

	//L block
	Kxu.reposition (3*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
	Kxv.reposition (3*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
	Kxw.reposition (3*n_u_dofs + n_p_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
	Kyu.reposition (4*n_u_dofs + n_p_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
	Kyv.reposition (4*n_u_dofs + n_p_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
	Kyw.reposition (4*n_u_dofs + n_p_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
	Kzu.reposition (5*n_u_dofs + n_p_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
	Kzv.reposition (5*n_u_dofs + n_p_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
	Kzw.reposition (5*n_u_dofs + n_p_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);

	
	
	//Symmetric to L block
	Kux.reposition ( u_var*n_u_dofs, 3*n_u_dofs + n_p_dofs, n_u_dofs, n_u_dofs);
	Kvx.reposition ( v_var*n_u_dofs,3*n_u_dofs + n_p_dofs, n_u_dofs, n_v_dofs);
	Kwx.reposition ( w_var*n_u_dofs, 3*n_u_dofs + n_p_dofs, n_u_dofs, n_w_dofs);
	Kuy.reposition ( u_var*n_v_dofs,4*n_u_dofs + n_p_dofs, n_v_dofs, n_u_dofs);
	Kvy.reposition (v_var*n_v_dofs, 4*n_u_dofs + n_p_dofs,  n_v_dofs, n_v_dofs);
	Kwy.reposition ( w_var*n_v_dofs, 4*n_u_dofs + n_p_dofs, n_v_dofs, n_w_dofs);
	Kuz.reposition ( u_var*n_w_dofs,5*n_u_dofs + n_p_dofs, n_w_dofs, n_u_dofs);
	Kvz.reposition ( v_var*n_w_dofs, 5*n_u_dofs + n_p_dofs, n_w_dofs, n_v_dofs);
	Kwz.reposition ( w_var*n_w_dofs,5*n_u_dofs + n_p_dofs,  n_w_dofs, n_w_dofs);

	
	//M block
	Kxp.reposition (3*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
	Kyp.reposition (4*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_v_dofs, n_p_dofs);
	Kzp.reposition (5*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_w_dofs, n_p_dofs);

	//O block
	Kxx.reposition (3*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_u_dofs, n_u_dofs);
	Kxy.reposition (3*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_u_dofs, n_v_dofs);
	Kxz.reposition (3*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_u_dofs, n_w_dofs);
	Kyx.reposition (4*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_v_dofs, n_u_dofs);
	Kyy.reposition (4*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_v_dofs, n_v_dofs);
	Kyz.reposition (4*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_v_dofs, n_w_dofs);
	Kzx.reposition (5*n_u_dofs + n_p_dofs, 3*n_u_dofs + n_p_dofs, n_w_dofs, n_u_dofs);
	Kzy.reposition (5*n_u_dofs + n_p_dofs, 4*n_u_dofs + n_p_dofs, n_w_dofs, n_v_dofs);
	Kzz.reposition (5*n_u_dofs + n_p_dofs, 5*n_u_dofs + n_p_dofs, n_w_dofs, n_w_dofs);

	Fx.reposition (3*n_u_dofs + n_p_dofs, n_u_dofs);
	Fy.reposition (4*n_u_dofs + n_p_dofs, n_v_dofs);
	Fz.reposition (5*n_u_dofs + n_p_dofs, n_w_dofs);

	dof_map.dof_indices (elem, dof_indices_x, x_var);
	dof_map.dof_indices (elem, dof_indices_y, y_var);
	dof_map.dof_indices (elem, dof_indices_z, z_var);

	
	///////////////////////////////////
	//THE BOUNDARY CONDITIONS
	///////////////////////////////////
	
	/*
	if(!es.parameters.get<std::string>("problem").compare("cylinder")){
		//Unconfined Compression
		#include "boundary_conditions/classic_disp_cylinder_bcs.cpp" 
	}
	
	//Swelling test problem (still no weak bcs implemented)
	if(!es.parameters.get<std::string>("problem").compare("cube")){
		  #include "boundary_conditions/swelling_test.cpp" 
		  #include "boundary_conditions/pressure_stress_bc.cpp" 
	}
	
	*/
	
	if(!es.parameters.get<std::string>("problem").compare("lung")){
	
	   #include "boundary_conditions/lobe_affine_bcs.cpp"
	
		//	  #include "boundary_conditions/weak_test_bc.cpp"	
		
	   #include "boundary_conditions/lobe_fixflux_bcs.cpp"
		
	}
	
	
		if(!es.parameters.get<std::string>("problem").compare("cube")){
	
			  #include "boundary_conditions/lobe_affine_bcs.cpp"

				//	 #include "boundary_conditions/sliding_cube_bcs.cpp"

//		  #include "boundary_conditions/sliding_cube_bcs_free.cpp"
			
		 		 #include "boundary_conditions/lobe_fixflux_bcs.cpp"

		// #include "boundary_conditions/sliding_cube_bcs.cpp"
	//	 #include "boundary_conditions/weak_test_bc.cpp"	
		
		  //#include "boundary_conditions/noflux_cube_fixbcs.cpp"	

	}
	
	
	if(!es.parameters.get<std::string>("problem").compare("cylinder")){
	
		 #include "boundary_conditions/pull_cylinder_fixall_bcs.cpp"
  // #include "boundary_conditions/weak_test_bc.cpp"	
 #include "boundary_conditions/lobe_fixflux_bcs.cpp"

	}
		
	/////////////////////////////////////////
	
	//Add the rhs contribution of the stabilization
 	#include "assemble_rhs_stabilization.cpp" 
 
  dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  newton_update.matrix->add_matrix (Ke, dof_indices);
  newton_update.rhs->add_vector    (Fe, dof_indices);
} // end of element loop

  newton_update.matrix->close();
  newton_update.matrix->zero_rows(rows,1.0);
  newton_update.rhs->close();

  for (int i=0; i < rows.size(); i++) {
	newton_update.rhs->set(rows[i],rows_values[i]);
  }

  newton_update.matrix->zero_rows(pressure_rows, 1.0);
  for (int i=0; i < pressure_rows.size(); i++) {
	newton_update.rhs->set(pressure_rows[i],pressure_rows_values[i]);
  }
  newton_update.matrix->close();
  newton_update.rhs->close();
	
  
  newton_update.matrix->close();
  newton_update.matrix->zero_rows(pressure_rows, 1.0);
  newton_update.rhs->close();
  newton_update.matrix->close();
  newton_update.update(); 

  std::cout<<"newton_update.rhs->l1_norm () "<<newton_update.rhs->l1_norm ()<<std::endl;
	
	return;
}
 

