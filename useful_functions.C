#include "poro.h"

using namespace std;

void test(int a) {
 // std::cout << "ex "<< a << std::endl;
}


double diffclock(clock_t clock1,clock_t clock2)
{
	double diffticks=clock1-clock2;
	double diffms=(diffticks*1000)/CLOCKS_PER_SEC;
	return diffms;
} 



double average_stress(RealTensor A)
{
	 
  
   RealTensor I;
		 I(0, 0) = 1.0; I(1, 1) = 1.0; I(2, 2) = 1.0;

		 Real eig1,eig2,eig3,r,phie,q,p2,p;
		 
		 
		 //Algo to calculate evals
 		 
		 //To test
		// A(0, 0) = 8.0; A(0,1) = 1.0; A(0, 2) = 6.0;
		 //A(1, 0) = 3.0; A(1,1) = 5.0; A(1, 2) = 7.0;
		// A(2, 0) = 4.0; A(2,1) = 9.0; A(2, 2) = 2.0;

		 
		 //p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2
		 Real p1=A(0,1)*A(0,1) +A(0,2)*A(0,2) + A(1,2)*A(1,2);
			//if (p1 == 0) 
		   // % A is diagonal.
		 	if(p1==0){
			   eig1=A(1,1);
			   eig2=A(2,2);
			   eig3=A(3,3);
			}else{
			  // q = trace(A)/3
			   q = A.tr()/3;	   
		
			  // p2 = (A(1,1) - q)^2 + (A(2,2) - q)^2 + (A(3,3) - q)^2 + 2 * p1
			   p2=pow((A(0,0)-q),2) +pow((A(1,1)-q),2)  + pow((A(2,2)-q),2) + 2*p1;
			  
 
			  // p = sqrt(p2 / 6)
			  p=pow(p2/6,0.5);
			  // B = (1 / p) * (A - q * I)       % I is the identity matrix
			  RealTensor B;
			  B=(1 / p) * (A - q * I); 
			
				// r = det(B) / 2
			  r= B.det()/2;
			  
 			  
			}
 	
  // % In exact arithmetic for a symmetric matrix  -1 <= r <= 1
  // % but computation error can leave it slightly outside this range.
  // if (r <= -1) 
  if(r <= -1){
   //   phi = pi / 3
	phie=PI/3;
	
	  // elseif (r >= 1)
  }else if(r>=1){
   //   phi = 0
	phie=0;
	
	  // else
	  }else{
   //   phi = acos(r) / 3
		phie=acos(r)/3;
		  // end
		
 
	  }
  // % the eigenvalues satisfy eig3 <= eig2 <= eig1
   eig1 = q + 2 * p * cos(phie);
   eig3 = q + 2 * p * cos(phie + (2*PI/3));
   eig2 = 3 * q - eig1 - eig3 ;  //  % since trace(A) = eig1 + eig2 + eig3

   
		
 /// std::cout<<" eig1 " << eig1 << std::endl;
  //  std::cout<<" eig2 " << eig2 << std::endl;
 // std::cout<<" eig3 " << eig3 << std::endl;

Real av_stress= pow(eig1*eig1 + eig2*eig2 + eig3*eig3,0.5);

	return av_stress;
} 


void calculate_numeric_jacobian(EquationSystems& es, SparseMatrix< Number >& num_jac_matrix) {

	
	
			/*
			 * 
			 * Put this in poro_main.C
	//Calculate the numerical jack
	//	SparseMatrix< Number > &num_jac_matrix=*(last_non_linear_soln.matrix);
	//  num_jac_matrix.zero();
	//	calculate_numeric_jacobian(equation_systems,num_jac_matrix);
	//	num_jac_matrix.close();
	//
		//use the numerical jacobian
		newton_update.assemble_before_solve=false;
		newton_update.update();
		assemble_solid(equation_systems,"Newton-update");
		newton_update.matrix->zero();
		newton_update.matrix->add(1, num_jac_matrix);
		newton_update.matrix->close();
 			//newton_update.matrix->print(std::cout);
		*/
	
	
	TransientLinearImplicitSystem & newton_update =
  es.get_system<TransientLinearImplicitSystem> ("Newton-update");
	
	
	TransientLinearImplicitSystem&  last_non_linear_soln = es.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
	
	//probably should clone this
	
	
	int N_dofs=  newton_update.rhs->size();
	Real eps= 0.001;

	//Get the unperturbed res.
	assemble_solid(es,"Newton-update");
	AutoPtr<NumericVector<Number> > fx (newton_update.rhs->clone());
	newton_update.rhs->zero();
	newton_update.matrix->zero();
	newton_update.matrix->close();

	std::cout<<"N_dofs  "<< N_dofs << std::endl; 

	for (unsigned int i=0; i<N_dofs; ++i)
	{
		//std::cout<<"ith dof "<< i << std::endl; 

	
		//Change one of the res dofs by eps
		last_non_linear_soln.solution->add(i, eps);
		last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(i,eps);
		last_non_linear_soln.current_local_solution->close();
	
		//Calculate the new rhs (at the moment are also calculatating lhs !)
		assemble_solid(es,"Newton-update");
		newton_update.matrix->zero();
		newton_update.rhs->zero();
		newton_update.matrix->close();

		//Get the difference between the peturbed res and the original res
		AutoPtr<NumericVector<Number> > f_diff (newton_update.rhs->clone());
		f_diff->add (-1., *fx);

		AutoPtr<NumericVector<Number> > jac_col (newton_update.rhs->clone());
		jac_col->zero();
		jac_col->add (1./eps, *f_diff);
		jac_col->close();

		//Assemble into the numerical jacobian matrix
		for (unsigned int t=1; t<=N_dofs; ++t)
		{
			num_jac_matrix.set(i,t,jac_col->el(t));
		}

		//Zero out for the next test
		f_diff->close();
		f_diff->zero();

		newton_update.rhs->zero();
		//subtract again eps
		
		last_non_linear_soln.solution->add(i,-eps);
		last_non_linear_soln.solution->close();
		last_non_linear_soln.current_local_solution->add(i,-eps);
		last_non_linear_soln.current_local_solution->close();
		
		
	}

 	num_jac_matrix.close();
 	num_jac_matrix.print(std::cout);

	
}
