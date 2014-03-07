#ifndef GENERAL_MATERIAL_H_
#define GENERAL_MATERIAL_H_

#include "poro.h"

 

class GeneralMaterialConfig {
public:

  GeneralMaterialConfig(
      const std::vector<std::vector<RealGradient> >& dphi,   const std::vector<std::vector<Real> >& psi) :
      dphi(dphi),psi(psi) {
  }

  void get_residual(DenseVector<Real> & residuum, unsigned int & i);

  /**
   * Return the stiffness matrix for the current state.
   */
  void get_linearized_stiffness(DenseMatrix<Real> & stiffness,
      unsigned int & i, unsigned int & j);

  void get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j);
  void get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i);

  /**
   * Flag to indicate if it is necessary to calculate values for stiffness
   * matrix during initialization.
   */
  bool calculate_linearized_stiffness;
//private:
  void build_b_0_mat(int i, DenseMatrix<Real>& b_l_mat);

  static void tensor_to_voigt(const RealTensor &tensor, DenseVector<Real> &vec);
  static void tensorOtensor_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);
  static void z_ref_to_voigt(const RealTensor &tensorB, const RealTensor &tensorA, DenseMatrix<Real> &mat);
 
  void c_update(RealTensor C) ;

  unsigned int current_qp;

  Number p_current;
  Real p_solid;

  const std::vector<std::vector<RealGradient> >& dphi;
  const std::vector<std::vector<Real> >& psi;

	DenseMatrix<Real> C_mat;

  Real I_1, I_2, I_3, J;
  RealTensor F, S, tau, sigma, C, invC, b, Ft, Identity;
  DenseMatrix<Real> B_L;
  DenseMatrix<Real> B_K;
};

#endif /* GENERAL_MATERIAL_H_ */

