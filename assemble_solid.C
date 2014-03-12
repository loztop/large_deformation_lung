#include "poro.h"
#include "poro_elastic_cc.h"

void assemble_solid (EquationSystems& es,
                      const std::string& system_name)
{

libmesh_assert (system_name == "Newton-update");
  
// Get a constant reference to the mesh object.
const MeshBase& mesh = es.get_mesh();
  
// The dimension that we are running
const unsigned int dim = mesh.mesh_dimension();
  
// Get a reference to the Stokes system object.
TransientLinearImplicitSystem & newton_update =
   es.get_system<TransientLinearImplicitSystem> ("Newton-update");

TransientLinearImplicitSystem & last_non_linear_soln =
    es.get_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");

const System & ref_sys = es.get_system("Reference-Configuration"); 

// Numeric ids corresponding to each variable in the system
const unsigned int u_var = last_non_linear_soln .variable_number ("s_u");
const unsigned int v_var = last_non_linear_soln .variable_number ("s_v");
const unsigned int w_var = last_non_linear_soln .variable_number ("s_w");
const unsigned int p_var = last_non_linear_soln .variable_number ("s_p");

FEType fe_vel_type = last_non_linear_soln.variable_type(u_var);
FEType fe_pres_type = last_non_linear_soln .variable_type(p_var);
AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));

QGauss qrule (dim, fe_vel_type.default_quadrature_order());
fe_vel->attach_quadrature_rule (&qrule);
test(65);
fe_pres->attach_quadrature_rule (&qrule);
const std::vector<Real>& JxW = fe_vel->get_JxW();
const std::vector<std::vector<Real> >& phi = fe_vel->get_phi();
const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();

const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
const std::vector<Point>& coords = fe_vel->get_xyz();
const DofMap & dof_map = last_non_linear_soln .get_dof_map();

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
Fu(Fe), Fv(Fe), Fw(Fe);
DenseSubVector<Number>    Fp(Fe);

std::vector<unsigned int> dof_indices;
std::vector<unsigned int> dof_indices_u;
std::vector<unsigned int> dof_indices_v;
std::vector<unsigned int> dof_indices_w;
test(661);

std::vector<unsigned int> dof_indices_p;

DenseMatrix<Real> stiff;
DenseVector<Real> res;
VectorValue<Gradient> grad_u_mat;
VectorValue<Gradient> grad_w_mat;
const Real dt    = es.parameters.get<Real>("dt");
const Real progress    = es.parameters.get<Real>("progress");

DenseVector<Real> p_stiff;
DenseVector<Real> p_res;

PoroelasticConfig material(dphi,psi);

// Just calculate jacobian contribution when we need to
material.calculate_linearized_stiffness = true;
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
  fe_pres->reinit (elem);
  
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
  test(77);

  System& aux_system = es.get_system<System>("Reference-Configuration");
  std::vector<unsigned int> undefo_index;
  std::vector<unsigned int> vel_index;      
  

 #include "assemble_solid_qp.cpp"

 #include "assemble_stabilization.cpp" 

	
  newton_update.matrix->add_matrix (Ke, dof_indices);
  newton_update.rhs->add_vector    (Fe, dof_indices);

} // end of element loop

//   dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
newton_update.rhs->close();
newton_update.matrix->close();


//Assemble the boundary conditions.
assemble_bcs(es);
 
    // newton_update.matrix->print_matlab("matrices/Kp_pnsym2000.dat");
    // newton_update.rhs->print_matlab("matrices/rp.dat");

return;
}





