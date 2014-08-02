#include "poro_elastic_cc.h"
#include "poro.h"





void PoroelasticConfig::calculate_stress_poro() {

  double mu = E / (2.0 * (1.0 + nu));
  double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  Real detF = F.det();
  RealTensor Ft = F.transpose();
  RealTensor C = Ft * F;
  RealTensor b = F * Ft;
  RealTensor identity;
  identity(0, 0) = 1.0; identity(1, 1) = 1.0; identity(2, 2) = 1.0;
  RealTensor invC = inv(C);

  //Whitely et al Neo-Hookean law (I think)
  S = 0.5 * lambda * (detF * detF - 1) * invC + mu * (identity - invC);
  S += - p_solid*J*invC;

  //convert to current configuration using a push forward operation
  tau = (F * S) * Ft;
  sigma = 1.0/detF * tau;
}

void PoroelasticConfig::init_for_qp(Point & rX,VectorValue<Gradient> & grad_u, Number & p_current, unsigned int qp, Real m, Real p_fluid) {
  
  this->current_qp = qp;
  this->p_solid = p_current;
  this->E = E_mod;
  this->nu = NU_mod;
  this->Phi_zero = PHI_ZERO;
  this->Rho_s = RHO_S;

  F.zero();
  S.zero();

  {
    RealTensor invF;
    invF.zero();
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j) {
        invF(i, j) += grad_u(i)(j);
      }
      F.add(inv(invF));
    }

    if (F.det() < -TOLERANCE) {
      std::cout << "detF < 0" << std::endl;
      libmesh_error();
    }

    this->Ft = F.transpose();
    this->C = Ft*F;
    this->c_update(C);
	this->calculate_permeability(rX);
	
    if (this->calculate_linearized_stiffness) {
      this->calculate_tangent();
    }
    this->calculate_stress_poro();
}


void PoroelasticConfig::get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) {
  Real detF = F.det();
  p_residuum.resize(1);
  p_residuum(0)=(1.0/detF) *psi[i][current_qp]*(detF-1-m); 
}

void PoroelasticConfig::c_update(RealTensor C) {     
  this-> C = C;
  this->invC = inv(C);
  this->b = F*Ft;
  this->Identity(0,0)=1.0;
  this->Identity(1,1)=1.0;
  this->Identity(2,2)=1.0;
  this->I_3 = C.det();
  this->I_1 = C(0,0)+C(1,1)+C(2,2);

  RealTensor Csqd = C*C;
  this->I_2 = 0.5*(pow(I_1,2) - (Csqd(0,0)+Csqd(1,1)+Csqd(2,2))) ;
  this->J=pow(I_3,(1.0/2.0));
  
  this->porosity=1-((1-Phi_zero)/J);
}

void PoroelasticConfig::calculate_tangent() {
  Real mu = E / (2 * (1 + nu));
  Real lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
  Real detF = F.det();
  C_mat.resize(6, 6);

  for (unsigned int i = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j) {
      if (i == j) {
        C_mat(i, j) = 2 * mu + lambda;
        C_mat(i + 3, j + 3) = mu - 0.5 * lambda * (detF * detF - 1);
      } else {
        C_mat(i, j) = lambda * detF * detF;
      }
    }
  }
     
  Real factor1=J;
  Real delta_c=-1*p_solid*factor1;
  DenseMatrix<Real> Z_mat;
  Z_mat.resize(6, 6);
  z_ref_to_voigt(Identity,Identity,Z_mat);
  Z_mat.scale((-1.0)*delta_c);
  C_mat+=Z_mat;

  Real fac=1*p_solid*factor1;
  Real delta_a=fac*(1.0/2.0);
  DenseMatrix<Real> invCinvC_mat;
  invCinvC_mat.resize(6, 6);
  tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
  invCinvC_mat.scale(delta_a);
  C_mat+=invCinvC_mat;

  /*
  //This is the linearization for S = 2*0.1 * (identity) - p_solid*invC; It converges beautifully.
  C_mat.resize(6, 6);
  Real delta_c=-1*p_solid;
  DenseMatrix<Real> Z_mat;
  Z_mat.resize(6, 6);
  z_ref_to_voigt(Identity,Identity,Z_mat);
  Z_mat.scale((-1.0)*delta_c);
  C_mat+=Z_mat;
  */

  /*
  Real factor2=0*J*(-1-m);
  Real delta_d=-1*p_solid*factor2;
  //DenseMatrix<Real> Z_mat;
  Z_mat.resize(6, 6);
  z_ref_to_voigt(Identity,Identity,Z_mat);
  Z_mat.scale((-1.0)*delta_d);
  C_mat+=Z_mat;

  fac=1*p_solid*factor2;
  Real delta_e=fac*(1.0/2.0);
  //DenseMatrix<Real> invCinvC_mat;
  invCinvC_mat.resize(6, 6);
  tensorOtensor_to_voigt(Identity,Identity,invCinvC_mat);
  invCinvC_mat.scale(delta_e);
  C_mat+=invCinvC_mat;
  */
}

void PoroelasticConfig::calculate_permeability(Point & rX) {
  
  
  //Calculate for specific co-ordiante so we can paramterise in space
  
  //Isotropic "Boring law"
  	//  Kperm=KPERM;

	/*
  if(rX(2)<0.5){
	  Kperm=KPERM;
  }else{
	  Kperm=KPERM/1000;
  }
  */
  
  //Isotropic "Lung law"
  Kperm=KPERM*pow((1-((1-J)/porosity)),2.0/3.0);
  
}

void PoroelasticConfig::get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
  // Find and write down the mathematics for this section.
  // eulerian tangent matrix at 6.15 - Bonet

  RealTensor Ft = F.transpose();
  Real detF = F.det();
  RealTensor C = Ft * F;
  RealTensor invC = inv(C);
  RealTensor A = invC*F;
  p_stiffness.resize(3);
  RealTensor invF = inv(F);

  //Build K_B (differentiate u eqn in p direction).
  RealTensor DPHIi;

  for (unsigned int z = 0; z < 3; ++z) {
    DPHIi(z,0)=dphi[i][current_qp](z); 
  }
  RealTensor invFt = invF.transpose();
  RealTensor uvw_p_stiff=-psi[j][current_qp]*invFt*DPHIi;
  p_stiffness(0)=uvw_p_stiff(0,0); 
  p_stiffness(1)=uvw_p_stiff(1,0); 
  p_stiffness(2)=uvw_p_stiff(2,0); 

  RealTensor a;
  DenseVector<Real> tmp;
  tmp.resize(3);
  RealTensor uvw_p_ref;
  uvw_p_ref =-J*invC ;
  a = 1/detF * (F * ( uvw_p_ref ) ) * Ft;

  RealTensor Identity;
  Identity(0,0)=1.0;
  Identity(1,1)=1.0;
  Identity(2,2)=1.0;

  a=-Identity;

  B_L.resize(3, 6);
  this->build_b_0_mat(i, B_L);
  DenseVector<Real> sigma_voigt(6);
  tensor_to_voigt(a, sigma_voigt);
  B_L.vector_mult(tmp, sigma_voigt);

  p_stiffness(0)=psi[j][current_qp]*tmp(0);
  p_stiffness(1)=psi[j][current_qp]*tmp(1);
  p_stiffness(2)=psi[j][current_qp]*tmp(2);
}

void PoroelasticConfig::get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
  //Build K_C (differentiate p eqn in u direction).
  RealTensor Ft = F.transpose();
  Real detF = F.det();
  RealTensor C = Ft * F;
  RealTensor invC = inv(C);
  RealTensor A = invC*F;
  p_stiffness.resize(3);
  RealTensor invF = inv(F);
  RealTensor invFt = invF.transpose();

  RealTensor DPHIj;
  for (unsigned int z = 0; z < 3; ++z) {
    DPHIj(z,0)=dphi[j][current_qp](z); 
  }
  RealTensor p_uvw_stiff= psi[i][current_qp]*invFt*DPHIj;
  p_stiffness(0)=p_uvw_stiff(0,0);
  p_stiffness(1)=p_uvw_stiff(1,0); 
  p_stiffness(2)=p_uvw_stiff(2,0); 

  p_stiffness.resize(3);
  p_stiffness(0)=psi[i][current_qp]*dphi[j][current_qp](0);
  p_stiffness(1)=psi[i][current_qp]*dphi[j][current_qp](1);
  p_stiffness(2)=psi[i][current_qp]*dphi[j][current_qp](2);

  RealTensor a;
  DenseVector<Real> tmp;
  tmp.resize(3);

  a = 1/detF * (F * ( J*invF ) ) * Ft;

  B_L.resize(3, 6);
  DenseVector<Real> a_voigt(6);
  this->build_b_0_mat(j, B_L);
  tensor_to_voigt(a, a_voigt);
  B_L.vector_mult(tmp, a_voigt);

  p_stiffness(0)=psi[i][current_qp]*tmp(0);
  p_stiffness(1)=psi[i][current_qp]*tmp(1); 
  p_stiffness(2)=psi[i][current_qp]*tmp(2);
}

