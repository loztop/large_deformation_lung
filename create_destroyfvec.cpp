//VecDestroy(cols_coupling_fvec);
	
//VecDestroy(vals_coupling_fvec);





	Vec bla;
	VecCreate(PETSC_COMM_WORLD,&bla);
  PetscObjectSetName((PetscObject) bla, "Solution");
  VecSetSizes(bla,PETSC_DECIDE,12312312);
  VecSetFromOptions(bla);
  VecSetFromOptions(bla);
  VecAssemblyBegin(bla);
  VecAssemblyEnd(bla);
	
	VecDestroy(bla);
	
	
	
	PetscReal *lol;
	lol=(PetscReal *)malloc((234)*sizeof(PetscReal));

PetscFree(lol); 
