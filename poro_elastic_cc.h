#ifndef poro_elastic_CC_H_
#define poro_elastic_CC_H_

#include "poro.h"
#include "general_material_cc.h"

class PoroelasticConfig : public GeneralMaterialConfig {

public:

  PoroelasticConfig(const std::vector<std::vector<RealGradient> >& dphi, const std::vector<std::vector<Real> >& psi) :  GeneralMaterialConfig(dphi, psi){

  }

  Real p_fluid;
  Real f_density;
  Real Kperm;
  Real Phi_zero;
  Real Rho_s;
  Real Rho_f;

  Real porosity;

	
  RealTensor GRAD_w, GRAD_w_t;

  //Can get rid of m and p residual
  Real m;


  //Neo hookean stuff
  Real E;
  Real nu;

  void calculate_stress_poro();

  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) ;
  
  //standard init
  void init_for_qp(Point & rX,VectorValue<Gradient> & grad_u, Number & p_current,      unsigned int qp, Real m, Real p_fluid);
  void c_update(RealTensor C) ;
  void calculate_tangent();
  void calculate_permeability(Point & rX);
  
  void get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness,      unsigned int & i, unsigned int & j);
  void get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);

};

#endif /* NONLINEAR_NEOHOOKE_CC_H_ */

//#endif
