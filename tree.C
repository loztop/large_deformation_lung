#include "tree.h"

#define PIE 3.1415926535897932384626433832
#define MU_F 1.82e-05
void Tree::read_tree() {

	std::cout<< "Reading tree data" <<std::endl;
	
	
	//std::ifstream infile_node("lozout7.node");
	//std::ifstream infile_edge("lozout7.edge");
	
	std::ifstream infile_node("lozout.node");
	std::ifstream infile_edge("lozout.edge");
		
	//std::ifstream infile_node("cube.node");
	//std::ifstream infile_edge("cube.edge");
	
	
	//Read in node file

	int ed_count=0;
	int node_count=0;

	int lc=1;
	std::string line_node;
	std::cout<<  "Read in node file ... " <<   std::endl;
	while (std::getline(infile_node, line_node))
	{
		if(lc==1){
		  
			std::istringstream iss(line_node);
						
			if (!(iss >> number_nodes)) { break; } // error

			number_edges=number_nodes-1;
			nodes.resize(number_nodes);
			nodes_type.resize(number_nodes);
			nodes_parent_node.resize(number_nodes);
			nodes_parent_edge.resize(number_nodes);

			edges_radius.resize(number_edges);
			edges_upper_node.resize(number_edges);
			edges_lower_node.resize(number_edges);
			edges_type.resize(number_edges);
			edges_resistance.resize(number_edges);
			edges_child1.resize(number_edges);
      edges_child2.resize(number_edges);
			
		}
		
		if(lc>1){
			std::istringstream iss(line_node);							
			double num_node, x, y , z , rad, type;
			if (!(iss >>num_node>> x >> y >> z >> rad >> type)) { break; } // error
			//	std::cout<< x << " " << y << " " << z << " " << rad << " " << type <<  std::endl;
				
				Point point(x, y, z);
				nodes(node_count) = point ;
				nodes_type(node_count) = type ;

				node_count=node_count+1;
				
				//Read in resistances
				if(lc>2) {
					edges_radius(ed_count)=rad/1000; //Put into meters (from mm)
					edges_type(ed_count)=type;
					
					Real length_pipe=4.0; //Assume l=4*r.
					//edges_resistance(ed_count)=(8.0*MU_F*length_pipe)/(PIE*pow(edges_radius(ed_count),3));

					edges_resistance(ed_count)=edges_radius(ed_count);

					
					//edges_resistance(ed_count)=1;
					
					ed_count=ed_count+1;
					
					
				}
		}
		
		lc=lc+1;
	}
	
	

	

	
  lc=1;
	std::string line_edge;
	std::cout<<  "Read in edge file ... " <<   std::endl;
	while (std::getline(infile_edge, line_edge))
	{
		if(lc==1){
		  //Already have taken care of this when reading the node file.	
		}
		
		if(lc>1){
			std::istringstream iss(line_edge);							
			double edge_num, p , c;
			if (!(iss >> edge_num >> p >> c)) { break; } // error
			//	std::cout<< edge_num << " " << p << " " << c << std::endl;
				edges_upper_node(edge_num)=p;	
				edges_lower_node(edge_num)=c;	
				
					//	std::cout<< "edges_upper_node(edge_num) "<< edges_upper_node(edge_num) << std::endl;
				
		}
		
		lc=lc+1;
	}
	
	//Figure out which branches are connected.
	
	for (int i=0; i < number_edges; i++) {
	  
	   //edges that have a no branches (end edge)
	  if(edges_type(i) ==0){
			edges_child1(i)=0;
			edges_child2(i)=0;		
	  }
	  
	  //edges that have a single branch
	  if(edges_type(i) ==1){
			//This assumes that branches after each other appear in order in the file.
			edges_child1(i)=edges_lower_node(i);
			edges_child2(i)=0;		
	  }
	  
	  //edges that have a branch (2 edges leaving)
	  if(edges_type(i) ==2){
		Real lower_node=edges_lower_node(i);

		//find all edges that have this as a upper node
		int found_node=0;
		for (int k=0; k < number_edges; k++) {

		  if(edges_upper_node(k)==lower_node && found_node==0){
				edges_child1(i)=k;
				//This assumes that the two branches follow each other in the file.
				edges_child2(i)=k+1;
				found_node=1;			
		  }
		  
		}
		
		
	  }
	  
	}
	
	//Initialise node pressures
	nodes_pressure.resize(number_nodes);
	for (int i=0; i < number_nodes; i++) {
		nodes_pressure(i) = 0;
	}
	
	//Initialise distal edges_flowrate
	edges_flowrate.resize(number_edges);
	for (int i=0; i < number_edges; i++) {
		edges_flowrate(i) = 0;
	}

	
    //Initialise deformed node position
	nodes_deformed.resize(number_nodes);
	for (int i=0; i < number_nodes; i++) {
		nodes_deformed(i) = nodes(i);
	}


	//Find the upper node for each node
	for (int i=0; i < number_nodes; i++) {
	  //Starting node (it has no parent) 
	  if(i==0){
			nodes_parent_node(i)=0;		
	  }else{
			//Find the edge that has this node as its lower node
			int found_node=0;
			for (int k=0; k < number_edges; k++) {
				if(edges_lower_node(k)==i && found_node==0){
					 nodes_parent_edge(i)=k;
					 nodes_parent_node(i)=edges_upper_node(k);
					 found_node=1;			
				}  
			}
	  }		
	}
	
	/*
	  std::cout<< "nodes_parent_edge "<< nodes_parent_edge << std::endl;
	  std::cout<< "nodes_parent_node "<< nodes_parent_edge << std::endl;

	  //std::cout<< "nodes  "<< nodes_parent_edge << std::endl;
		std::cout<< "nodes_type "<< nodes_type << std::endl;
  	std::cout<< "edges_radius "<< edges_radius << std::endl;
		std::cout<< "edges_type "<< edges_type << std::endl;
		std::cout<< "edges_upper_node "<< edges_upper_node << std::endl;
		std::cout<< "edges_lower_node "<< edges_lower_node << std::endl;
		std::cout<< "edges_child1 "<< edges_child1 << std::endl;
		std::cout<< "edges_child2 "<< edges_child2 << std::endl;
 */
	
}



Vec Tree::make_tree_rhs ()
{

 // std::cout<<"Make tree rhs "<<std::endl;

  //Petsc linear system example
  int NumberOfEntries = number_nodes+number_edges;

  //  Matrices and vectors in the linear system Au=b
  Mat A;
  Vec u, b;

  //  Create a PETSc Krylov space linear solver
  KSP ksp;

	//  Create the matrix A with 3 non-zeros per row
  MatCreateSeqAIJ(PETSC_COMM_SELF, NumberOfEntries, NumberOfEntries, NumberOfEntries,   PETSC_NULL, &A);
  MatSetFromOptions(A);

	//  Create the vector u
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, NumberOfEntries);
  VecSetFromOptions(u);
  //  Create the vector b to have identical size to u
  VecDuplicate(u, &b);

  VecSet(b, 0.0);
  // Set up array to access entries of b directly
  PetscScalar *bvec;
  VecGetArray(b, &bvec);

	//Not doing anything !!
  
	//Set outflow BCs (prescribe flow at outlets)
	for (double j=0; j < number_edges ; j++) {
		
		//Apply outflow BCS for end branches
		if(edges_type(j)==0){
			VecSetValue(b,j+number_nodes, edges_flowrate(j) ,ADD_VALUES);
		}
		
	}	
	
  VecAssemblyBegin(b);VecAssemblyEnd(b);


  return b;
}

Mat Tree::make_tree_matrix ()
{
	
	std::cout<<"Assembling Tree matrix "<<std::endl;

  //Petsc linear system example
  int NumberOfEntries = number_nodes+number_edges;

  //  Matrices and vectors in the linear system Au=b
  Mat A;

  //  Create the matrix A with 3 non-zeros per row
  MatCreateSeqAIJ(PETSC_COMM_SELF, NumberOfEntries, NumberOfEntries, NumberOfEntries,   PETSC_NULL, &A);
  MatSetFromOptions(A);

  //Initialise the entries of A, b to zero
  MatZeroEntries(A);	

	
  //Deal with the nodes (pressures)
  //P_0 this will be a boundary condition (atm P_0=0)

  for (double j=0; j < number_nodes ; j++) {	
	//P_j_parent-P_j-Q_{j-1}*resistance=0
			
		if(j==0){
			  MatSetValue(A, 0, 0, 1 ,ADD_VALUES); 			
		}else{
				MatSetValue(A, j, nodes_parent_node(j), 1 ,ADD_VALUES); 
				MatSetValue(A, j, j, -1 ,ADD_VALUES); 
				MatSetValue(A, j, number_nodes+ nodes_parent_edge(j), -edges_resistance(j-1),ADD_VALUES);
		}
		
  }	
  
  //Deal with flow rates 
  for (double j=0; j < number_edges ; j++) {
	
		//If edge has branch comming off it
		if((edges_type(j)==2) ||(edges_type(j)==1) ){
			//Q_{j}-Q_{j*2}-Q_{j*2 +1}=0	
			MatSetValue(A, j+number_nodes , j +number_nodes  , 1 ,ADD_VALUES); 
			
			
			if(edges_child1(j)>0){
				MatSetValue(A, j+number_nodes , number_nodes + edges_child1(j) , -1 ,ADD_VALUES); 
			}
			if(edges_child2(j)>0){
			MatSetValue(A, j+number_nodes , number_nodes + edges_child2(j) , -1 ,ADD_VALUES); 	
			}
			
		}
	
		//If edge is the final branch
		if(edges_type(j)==0){
			//Q_j - (1/rd) P_{j} + ( (1/r_d)*(1/omega_j)*p_poro this in intro.C )
			MatSetValue(A, j+number_nodes, j+number_nodes ,1,ADD_VALUES);
			MatSetValue(A, j+number_nodes, edges_lower_node(j) , -1/edges_resistance(j-1)   ,ADD_VALUES); 			
		}
	
  }
  
  
	
  //  Update the matrix A to reflect the changes made above
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	//Print matrix
	//MatView(A,PETSC_VIEWER_STDOUT_WORLD);
  
	
	return A;
}

void Tree::update_resistances (EquationSystems& es)
{

  const Real dt    = es.parameters.get<Real>("dt");
  const Real time    = es.parameters.get<Real>("time");
  const Real progress    = es.parameters.get<Real>("progress");

  return;
}


void Tree::calculate_omega_j (EquationSystems& es)
{

	
  const Real dt    = es.parameters.get<Real>("dt");
  const Real progress    = es.parameters.get<Real>("progress");
  const Real time    = es.parameters.get<Real>("time");
 
  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();
  
  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();
  
  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Last_non_linear_soln");

 //TransientLinearImplicitSystem & system =    es.get_system<TransientLinearImplicitSystem> ("Stokes");
  // Numeric ids corresponding to each variable in the system
  const unsigned int u_var = system.variable_number ("s_u");
  const unsigned int v_var = system.variable_number ("s_v");
  #if THREED
  const unsigned int w_var = system.variable_number ("s_w");
  #endif
  const unsigned int p_var = system.variable_number ("s_p");
  const unsigned int x_var = system.variable_number ("x");
  const unsigned int y_var = system.variable_number ("y");
   #if THREED
  const unsigned int z_var = system.variable_number ("z");
  #endif
  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_disp_type = system.variable_type(u_var);
  FEType fe_vel_type = system.variable_type(x_var);

  // Get the Finite Element type for "p".
  FEType fe_pres_type = system.variable_type(p_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  AutoPtr<FEBase> fe_disp  (FEBase::build(dim, fe_disp_type));
  AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
    
  // Build a Finite Element object of the specified type for
  // the pressure variables.
  AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
  
  // A Gauss quadrature rule for numerical integration.
  // Let the \p FEType object decide what order rule is appropriate.
  QGauss qrule (dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_disp->attach_quadrature_rule (&qrule);
  fe_vel->attach_quadrature_rule (&qrule);
  fe_pres->attach_quadrature_rule (&qrule);
  
  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.   
  const std::vector<Real>& JxW = fe_vel->get_JxW();
  
  const std::vector<Point>& q_point = fe_vel->get_xyz();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> >& dphi = fe_disp->get_dphi();
  const std::vector<std::vector<Real> >& phi = fe_disp->get_phi();
  const std::vector<std::vector<RealGradient> >& f_dphi = fe_vel->get_dphi();
  const std::vector<std::vector<Real> >& f_phi = fe_vel->get_phi();

  // The element shape functions for the pressure variable
  // evaluated at the quadrature points.
  const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
  const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseMatrix<Number> Kstab;

#if !THREED
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke), Kux(Ke), Kuy(Ke),
    Kvu(Ke), Kvv(Ke), Kvp(Ke), Kvx(Ke), Kvy(Ke),
    Kpu(Ke), Kpv(Ke), Kpp(Ke), Kpx(Ke), Kpy(Ke),
    Kxu(Ke), Kxv(Ke), Kxp(Ke), Kxx(Ke), Kxy(Ke),
    Kyu(Ke), Kyv(Ke), Kyp(Ke), Kyx(Ke), Kyy(Ke);
#endif

#if THREED
  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke), Kup(Ke), Kux(Ke), Kuy(Ke), Kuz(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke), Kvp(Ke), Kvx(Ke), Kvy(Ke), Kvz(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke), Kwp(Ke), Kwx(Ke), Kwy(Ke), Kwz(Ke),
    Kpu(Ke), Kpv(Ke), Kpw(Ke), Kpp(Ke), Kpx(Ke), Kpy(Ke), Kpz(Ke),
    Kxu(Ke), Kxv(Ke), Kxw(Ke), Kxp(Ke), Kxx(Ke), Kxy(Ke), Kxz(Ke),
    Kyu(Ke), Kyv(Ke), Kyw(Ke), Kyp(Ke), Kyx(Ke), Kyy(Ke), Kyz(Ke),
    Kzu(Ke), Kzv(Ke), Kzw(Ke), Kzp(Ke), Kzx(Ke), Kzy(Ke), Kzz(Ke);
#endif

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe),
    Fx(Fe),
    Fy(Fe);

  #if THREED
  DenseSubVector<Number>
    Fw(Fe),
    Fz(Fe);
  #endif

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_p;
  std::vector<unsigned int> dof_indices_x;
  std::vector<unsigned int> dof_indices_y;
  #if THREED
  std::vector<unsigned int> dof_indices_w;
  std::vector<unsigned int> dof_indices_z;
  #endif
  
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  

  //vectors needed to store bcs. (wha happened to condense ?)
  std::vector<  int > rows;
  std::vector<  int > pressure_rows;
  std::vector< Real > rows_values;
  std::vector< Real > pressure_rows_values;

  std::vector<unsigned int> stab_dofs_rows;
  std::vector<unsigned int> stab_dofs_cols;
  std::vector<Real> stab_dofs_vals;

  //DenseVector<Real> omega_j;
  omega_j.resize(pow(2,N_level-1));
	  
  for ( ; el != end_el; ++el)
    {    

      const Elem* elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);
      dof_map.dof_indices (elem, dof_indices_p, p_var);
      dof_map.dof_indices (elem, dof_indices_x, x_var);
      dof_map.dof_indices (elem, dof_indices_y, y_var);
      #if THREED
      dof_map.dof_indices (elem, dof_indices_w, w_var);
      dof_map.dof_indices (elem, dof_indices_z, z_var);
      #endif

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size(); 
      const unsigned int n_v_dofs = dof_indices_v.size();
      const unsigned int n_p_dofs = dof_indices_p.size();
      const unsigned int n_x_dofs = dof_indices_x.size(); 
      const unsigned int n_y_dofs = dof_indices_y.size();
      #if THREED
      const unsigned int n_w_dofs = dof_indices_w.size();
      const unsigned int n_z_dofs = dof_indices_z.size();
      #endif
      
      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe_disp->reinit  (elem);
      fe_vel->reinit  (elem);
      fe_pres->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
      Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
      Kux.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs , n_u_dofs, n_x_dofs);
      Kuy.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_u_dofs, n_y_dofs);
      #if THREED
      Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
      Kuz.reposition (u_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kvx.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_v_dofs, n_x_dofs);
      Kvy.reposition (v_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_v_dofs, n_y_dofs);
      #if THREED
      Kvw.reposition (v_var*n_u_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
      Kuz.reposition (v_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      #if THREED
      Kwu.reposition (w_var*n_w_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kwv.reposition (w_var*n_w_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
      Kwp.reposition (w_var*n_w_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
      Kwx.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs , n_v_dofs, n_x_dofs);
      Kwy.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_v_dofs, n_y_dofs);
      Kww.reposition (w_var*n_w_dofs, w_var*n_u_dofs, n_v_dofs, n_w_dofs);
      Kwz.reposition (w_var*n_w_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_u_dofs, n_z_dofs);
      #endif

      Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
      Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
      Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
      Kpx.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs , n_p_dofs, n_x_dofs);
      Kpy.reposition (p_var*n_v_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_p_dofs, n_y_dofs);
      #if THREED
      Kpw.reposition (p_var*n_u_dofs, w_var*n_u_dofs, n_p_dofs, n_w_dofs);
      Kpz.reposition (p_var*n_u_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_p_dofs, n_z_dofs);
      #endif

      Kxu.reposition (p_var*n_u_dofs + n_p_dofs, u_var*n_u_dofs, n_x_dofs, n_u_dofs);
      Kxv.reposition (p_var*n_u_dofs + n_p_dofs, v_var*n_u_dofs, n_x_dofs, n_v_dofs);
      Kxp.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs, n_x_dofs, n_p_dofs);
      Kxx.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs , n_x_dofs, n_x_dofs);
      Kxy.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_x_dofs, n_y_dofs);
      #if THREED
      Kxw.reposition (p_var*n_u_dofs + n_p_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kxz.reposition (p_var*n_u_dofs + n_p_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif


      Kyu.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, u_var*n_u_dofs, n_y_dofs, n_u_dofs);
      Kyv.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, v_var*n_u_dofs, n_y_dofs, n_v_dofs);
      Kyp.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs, n_y_dofs, n_p_dofs);
      Kyx.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs , n_y_dofs, n_x_dofs);
      Kyy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_y_dofs, n_y_dofs);
      #if THREED
      Kyw.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kyz.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif

      #if THREED
      Kzu.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, u_var*n_u_dofs, n_y_dofs, n_u_dofs);
      Kzv.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, v_var*n_u_dofs, n_y_dofs, n_v_dofs);
      Kzp.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs, n_y_dofs, n_p_dofs);
      Kzx.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs , n_y_dofs, n_x_dofs);
      Kzy.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs+n_x_dofs , n_y_dofs, n_y_dofs);
      Kzw.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, w_var*n_u_dofs, n_x_dofs, n_w_dofs);
      Kzz.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, p_var*n_u_dofs + n_p_dofs+2*n_x_dofs , n_x_dofs, n_z_dofs);
      #endif



      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);
      Fp.reposition (p_var*n_u_dofs, n_p_dofs);
      Fx.reposition (p_var*n_u_dofs + n_p_dofs, n_x_dofs);
      Fy.reposition (p_var*n_u_dofs + n_p_dofs+n_x_dofs, n_y_dofs);
      #if THREED
      Fw.reposition (w_var*n_u_dofs, n_w_dofs);
      Fz.reposition (p_var*n_u_dofs + n_p_dofs+2*n_x_dofs, n_y_dofs);
      #endif
    
	  
		
      // Now we will build the element matrix.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {


		  
		  	//Find the closest distal tree_node to q_point[qp]
				int closest_j=0;
				Real dist_j=0;
				Real closest_dist=99999999999999999;
				int closest_j_idx=0;

				for (int j=pow(2,N_level-1)+1; j <= pow(2,N_level); j++) {	
					Real dist_j=pow(q_point[qp](0)-  nodes(j-1)(0),2)+pow(q_point[qp](1)-  nodes(j-1)(1),2)+pow(q_point[qp](2)-  nodes(j-1)(2),2);
											
					if(dist_j<closest_dist)
					{
						closest_dist=dist_j; 
						closest_j_idx=j;
					}
				}
				
				   
				for (unsigned int l=0; l<n_p_dofs; l++)
				{
				  				   
                  omega_j(closest_j_idx-pow(2,N_level-1) -1 ) += psi[l][qp]*JxW[qp];
				}
					

} // end qp



} // end of element loop
  

  return;
}


