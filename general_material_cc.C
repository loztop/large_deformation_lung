#include "general_material_cc.h"
#include <math.h>
#include "poro_elastic_cc.h"



void GeneralMaterialConfig::get_residual(DenseVector<Real> & residuum, unsigned int & i) {
       B_L.resize(3, 6);
       DenseVector<Real> sigma_voigt(6);
       this->build_b_0_mat(i, B_L);
       tensor_to_voigt(sigma, sigma_voigt);
       B_L.vector_mult(residuum, sigma_voigt);  //should this be B_L^{T} * sigma_voigt ?

     //  std::cout<< " sigma " << sigma <<std::endl;

       }

#if INCOMPRESSIBLE
void GeneralMaterialConfig::get_p_residual(DenseVector<Real> & p_residuum, unsigned int & i) {
       Real detF = F.det();
       p_residuum.resize(1);
       //     std::cout<<"neo  "<<std::cout;
       p_residuum(0)=(1/detF) *psi[i][current_qp]*(detF-1);   ///Rp according to chaste fem *psi[i][current_qp] in main code       
}
#endif   

void GeneralMaterialConfig::tensor_to_voigt(const RealTensor &tensor, DenseVector<Real> &vec) {
  vec(0) = tensor(0, 0);
  vec(1) = tensor(1, 1);
  vec(2) = tensor(2, 2);
  vec(3) = tensor(0, 1);
  vec(4) = tensor(1, 2);
  vec(5) = tensor(0, 2);

}

void GeneralMaterialConfig::tensorOtensor_to_voigt(const RealTensor &tensorA, const RealTensor &tensorB, DenseMatrix<Real> &mat){


//Top three rows
for (unsigned int i = 0; i < 3; ++i) {
mat(i,0) = tensorA(i,i)*tensorB(0,0)+tensorA(i,i)*tensorB(0,0);          
mat(i,1) = tensorA(i,i)*tensorB(1,1)+tensorA(i,i)*tensorB(1,1);               
mat(i,2) = tensorA(i,i)*tensorB(2,2)+tensorA(i,i)*tensorB(2,2);
mat(i,3) = tensorA(i,i)*tensorB(0,1)+tensorA(i,i)*tensorB(1,0);          
mat(i,4) = tensorA(i,i)*tensorB(0,2)+tensorA(i,i)*tensorB(2,0);               
mat(i,5) = tensorA(i,i)*tensorB(1,2)+tensorA(i,i)*tensorB(2,1);                 
}



//fourth row
int i=3; int a=0; int b=1;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);                 

//fith row
 i=4;  a=0;  b=2;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);   

//sixth row
 i=5;  a=1;  b=2;
mat(i,0) = tensorA(a,b)*tensorB(0,0)+tensorA(a,b)*tensorB(0,0);          
mat(i,1) = tensorA(a,b)*tensorB(1,1)+tensorA(a,b)*tensorB(1,1);               
mat(i,2) = tensorA(a,b)*tensorB(2,2)+tensorA(a,b)*tensorB(2,2);
mat(i,3) = tensorA(a,b)*tensorB(0,1)+tensorA(a,b)*tensorB(1,0);          
mat(i,4) = tensorA(a,b)*tensorB(0,2)+tensorA(a,b)*tensorB(2,0);               
mat(i,5) = tensorA(a,b)*tensorB(1,2)+tensorA(a,b)*tensorB(2,1);   

}


void GeneralMaterialConfig::z_ref_to_voigt(const RealTensor &tensorA, const RealTensor &tensorB, DenseMatrix<Real> &mat){


//Top three rows
for (unsigned int i = 0; i < 3; ++i) {
mat(i,0) = tensorA(i,0)*tensorB(i,0)*2;          
mat(i,1) = tensorA(i,1)*tensorB(i,1)*2;             
mat(i,2) = tensorA(i,2)*tensorB(i,2)*2;
mat(i,3) = tensorA(i,0)*tensorB(i,1)*2;       
mat(i,4) = tensorA(i,0)*tensorB(i,2)*2;         
mat(i,5) = tensorA(i,1)*tensorB(i,2)*2;                 
}



//fourth row
int i=3; int a=0; int b=1;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);                 

//fith row
 i=4;  a=0;  b=2;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);    

//sixth row
 i=5;  a=1;  b=2;
mat(i,0) = tensorA(a,0)*tensorB(b,0)+tensorA(a,0)*tensorB(0,b);          
mat(i,1) = tensorA(a,1)*tensorB(b,1)+tensorA(a,1)*tensorB(1,b);               
mat(i,2) = tensorA(a,2)*tensorB(b,2)+tensorA(a,2)*tensorB(2,b);
mat(i,3) = tensorA(a,0)*tensorB(b,1)+tensorA(a,1)*tensorB(0,b);          
mat(i,4) = tensorA(a,0)*tensorB(b,2)+tensorA(a,2)*tensorB(0,b);               
mat(i,5) = tensorA(a,1)*tensorB(b,2)+tensorA(a,2)*tensorB(1,b);     

}



void GeneralMaterialConfig::get_linearized_stiffness(DenseMatrix<Real> & stiffness, unsigned int & i, unsigned int & j) {
       stiffness.resize(3, 3);

       double G_IK = (sigma * dphi[i][current_qp]) * dphi[j][current_qp];
       stiffness(0, 0) += G_IK;
       stiffness(1, 1) += G_IK;
       stiffness(2, 2) += G_IK;

       B_L.resize(3, 6);
       this->build_b_0_mat(i, B_L);
       B_K.resize(3, 6);
       this->build_b_0_mat(j, B_K);

       B_L.right_multiply(C_mat);
       B_L.right_multiply_transpose(B_K);
       B_L *= 1/F.det();

       stiffness += B_L;
}

#if INCOMPRESSIBLE
void GeneralMaterialConfig::get_linearized_uvw_p_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
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
//std::cout<< "DPHIi " << DPHIi <<std::endl;
RealTensor invFt = invF.transpose();

RealTensor uvw_p_stiff=psi[j][current_qp]*invFt*DPHIi;
//std::cout<< "DPHIi " << DPHIi <<std::endl;
//std::cout<< "invFt" << invFt <<std::endl;
//std::cout<< "uvw_p_stiff" << uvw_p_stiff <<std::endl;
p_stiffness(0)=uvw_p_stiff(0,0); 
p_stiffness(1)=uvw_p_stiff(1,0); 
p_stiffness(2)=uvw_p_stiff(2,0); 



////////
p_stiffness.resize(3);
p_stiffness(0)=psi[j][current_qp]*dphi[i][current_qp](0);
p_stiffness(1)=psi[j][current_qp]*dphi[i][current_qp](1);
p_stiffness(2)=psi[j][current_qp]*dphi[i][current_qp](2);
/////////


}

void GeneralMaterialConfig::get_linearized_p_uvw_stiffness(DenseVector<Real> & p_stiffness, unsigned int & i, unsigned int & j) {
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

}


#endif

void GeneralMaterialConfig::build_b_0_mat(int i, DenseMatrix<Real>& b_0_mat) {
       for (unsigned int ii = 0; ii < 3; ++ii) {
               b_0_mat(ii, ii) = dphi[i][current_qp](ii);
       }
       //See wriggers 4.94
       b_0_mat(0, 3) = dphi[i][current_qp](1);
       b_0_mat(1, 3) = dphi[i][current_qp](0);
       b_0_mat(1, 4) = dphi[i][current_qp](2);
       b_0_mat(2, 4) = dphi[i][current_qp](1);
       b_0_mat(0, 5) = dphi[i][current_qp](2);
       b_0_mat(2, 5) = dphi[i][current_qp](0);
}


void GeneralMaterialConfig::c_update(RealTensor C) {
     
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
}


