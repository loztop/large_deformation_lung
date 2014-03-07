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
