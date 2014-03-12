#include "poro.h"

PetscVector<Number> assemble_coupled_rhs (EquationSystems& es,Tree& tree, Mesh& mesh)
{
 
	std::cout<< "Assembling big rhs " <<std::endl;

	TransientLinearImplicitSystem & system = 
    es.add_system<TransientLinearImplicitSystem> ("Newton-update");
		
	NumericVector< Number > &solution_in= *(system.solution);
	NumericVector< Number > &rhs_in = *(system.rhs);
	Number size_mat = system.rhs->size();
	
	Real size_t=tree.number_nodes+tree.number_edges;

  //Sort out rhs	- put into big vector
	Vec x;
	VecCreate(PETSC_COMM_WORLD,&x);
    PetscObjectSetName((PetscObject) x, "Solution");
    VecSetSizes(x,PETSC_DECIDE,size_mat);
    VecSetFromOptions(x);
	PetscVector<Number> xp(x) ;
	xp.add(1,solution_in);
		
	Vec r;
	VecCreate(PETSC_COMM_WORLD,&r);
    PetscObjectSetName((PetscObject) r, "rhs");
    VecSetSizes(r,PETSC_DECIDE,size_mat);
    VecSetFromOptions(r);
	PetscVector<Number> rp(r) ;
	rp.add(1,rhs_in);

	Vec rt;
	rt=tree.make_tree_rhs( );
	PetscVector<Number> rtp(rt) ;

	Vec big_r;
	VecCreate(PETSC_COMM_WORLD,&big_r);
  PetscObjectSetName((PetscObject) big_r, "big_rhs");
  
	VecSetSizes(big_r,PETSC_DECIDE,size_mat+size_t);
  
	//VecSetSizes(big_r,PETSC_DECIDE,size_mat);
	VecSetFromOptions(big_r);
	
		
	// add system.rhs into big_rp 
	for (int i=0; i<size_mat; i++) {
		//if(abs(rhs_in(i))>0){
			VecSetValue(big_r,i,rp(i),INSERT_VALUES); 
		//}
	}
	
	// add extra tree rhs entries  should all be 0 if p0=0
	for (int i=size_mat; i<size_mat+size_t; i++) {
	  VecSetValue(big_r, i , rtp(i-size_mat) ,INSERT_VALUES); 
	}
	
	
	PetscVector<Number> big_rp(big_r) ;
  
  return big_rp;
  
}