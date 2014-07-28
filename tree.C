#include "tree.h"
#include <math.h>
 #include <sstream>
#include <string>
#include <petscksp.h>
 

#define PIE 3.1415926535897932384626433832
#define MU_F 1.92e-05



void Tree::read_tree(EquationSystems& es) {

	std::cout<< "Reading tree data" <<std::endl;
		
	std::string node_file_name;
  node_file_name = es.parameters.get<std::string>("tree_input") + ".node";
  std::string edge_file_name;
  edge_file_name = es.parameters.get<std::string>("tree_input") + ".edge";
	
	std::ifstream infile_node(node_file_name.c_str());
	std::ifstream infile_edge(edge_file_name.c_str());

	//Read in node file

	//int max_line=20000;

	
	int ed_count=0;
	int node_count=0;
	int distal_count=0;

	int lc=1;
	std::string line_node;
	std::cout<<  "Read in node file ... " <<   std::endl;
	while (std::getline(infile_node, line_node) )
	{
		if(lc==1){
		  
			std::istringstream iss(line_node);
						
			if (!(iss >> number_nodes)) { break; } // error

			//number_nodes=max_line-1;
			
			number_edges=number_nodes-1;
			nodes.resize(number_nodes);
			nodes_type.resize(number_nodes);
			nodes_parent_node.resize(number_nodes);
			nodes_parent_edge.resize(number_nodes);

			edges_diseased.resize(number_edges);
			edges_radius.resize(number_edges);
			edges_length.resize(number_edges);		
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

				//Change sign because tree has been previosuly flipped
				//Only do this for particular meshes as outlined below...
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/whole_lung_246.msh")){
			//	nodes(node_count) = -point ;
				}
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/whole_lung_751.msh")){
			//	nodes(node_count) = -point ;
				}
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048r_fine1447.msh")){
				nodes(node_count) = -point ;
				}
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048r_fine2881.msh")){
				nodes(node_count) = -point ;
				}			
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/N048_node6598.msh")){
				nodes(node_count) = -point ;
				}
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/N051_node870.msh")){
				nodes(node_count) = -point ;
				}
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/APLE_36266rmerge_insp_765.msh")){
				nodes(node_count) = -point ;
				}
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/APLE_36266rmerge_insp_2010.msh")){
				nodes(node_count) = -point ;
				}
				
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/APLE_36266rmerge_insp_2967.msh")){
				nodes(node_count) = -point ;
				}
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/A65rmerge_insp2233.msh")){
				nodes(node_count) = -point ;
				}
				if(!es.parameters.get<std::string>("mesh_input").compare("meshes/lung/A65rmerge_inspfine3174.msh")){
				nodes(node_count) = -point ;
				}
				
				nodes_type(node_count) = type ;

				//coount distal edges
				if( type<1){ distal_count=distal_count+1;}
				
				node_count=node_count+1;
				
				//Read in resistances
				if(lc>2) {
					edges_radius(ed_count)=rad; //Put into meters (from mm)
					edges_type(ed_count)=type;
					
					Real length_pipe=0.1; //Assume l=4*r.
					
					ed_count=ed_count+1;
					
					
				}
		}
		
		lc=lc+1;
	}
		
		
	std::cout<<  "distal_count " << distal_count<<  std::endl;

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
			}
		
		lc=lc+1;
	}

	Real disease_count=0;
	//Figure out which branches are connected.
	for (int i=0; i < number_edges; i++) {
	   //edges that have a no branches (end edge)
	  if(edges_type(i) ==0){
			edges_child1(i)=0;
			edges_child2(i)=0;		
			
			
			
			//Deal with disease here 

		 /*
 			if(	  abs(nodes(i+1)(1))> 160    ){

			  edges_diseased(i)=1;
		    disease_count=disease_count+1;
			  
				//Instead reduce the radius
				edges_radius(i)=edges_radius(i)/2;

			}
		 */
			
			 
			
			
	  }
	  
	  //edges that have a single branch !!!!!!!!!!!!! fix below
	  if(edges_type(i) ==1){
		
	     edges_child2(i)=0;		

					
		Real lower_node=edges_lower_node(i);
		//find all edges that have this as a upper node
		int found_node=0;
		 int found_node_second=0;

		for (int k=0; k < number_edges; k++) {

		  if(edges_upper_node(k)==lower_node && found_node==0){
				edges_child1(i)=k;

				found_node=1;		
				k=k+1;
		  }
		  
		}
	  }
	  
	  //edges that have a branch (2 edges leaving)
	  if(edges_type(i) ==2){
		Real lower_node=edges_lower_node(i);

		//find all edges that have this as a upper node
		int found_node=0;
		 int found_node_second=0;

		for (int k=0; k < number_edges; k++) {

		  if(edges_upper_node(k)==lower_node && found_node==0){
				edges_child1(i)=k;

				found_node=1;		
				k=k+1;
		  }
		  
		   if(edges_upper_node(k)==lower_node && found_node_second==0&& found_node==1){
				edges_child2(i)=k;

				found_node_second=1;			
		  }
		  
		}
		
		
	  }
	  
	}
	
	std::cout<<"disease_count "<<disease_count<<std::endl;
	
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
	
	//calculate the length and then resistance of an edge
	edges_flowrate.resize(number_edges);
	for (int i=0; i < number_edges; i++) {
		
	  //Wrong !!
		Point u_node=nodes(edges_upper_node(i));
		Point l_node=nodes(edges_lower_node(i)); //Bad assumpion - wrong !!!!!!!!!!!!
		Real dist =pow( pow(u_node(0)-l_node(0),2)+pow(u_node(1)-l_node(1),2)+pow(u_node(2)-l_node(2),2) , 0.5);
		
		//std::cout<< "dist  "<< dist << std::endl;
		//std::cout<< "u_node  "<< u_node << std::endl;
		//std::cout<< "l_node  "<< l_node << std::endl;

		edges_length(i) = dist;	
		
		//Make tree a bit more healthy
		
		if( dist> 5){
			
		 // std::cout<< "dist  "<< dist << std::endl;
		//	std::cout<< "u_node  "<< u_node << std::endl;
		//	std::cout<< "l_node  "<< l_node << std::endl;			
		}
		
		
		edges_resistance(i)=1*(8.0*MU_F*(edges_length(i)))/(PIE*pow(1*edges_radius(i),4));
		
		
		//edges_resistance(ed_count)=1;
	}
	
	
/*
	  std::cout<< "nodes_parent_edge "<< nodes_parent_edge << std::endl;
	  std::cout<< "nodes_parent_node "<< nodes_parent_edge << std::endl;
		std::cout<< "nodes_type "<< nodes_type << std::endl;
  	std::cout<< "edges_radius "<< edges_radius << std::endl;
		std::cout<< "edges_length "<< edges_child2 << std::endl;
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

  Vec u, b;
  
  //  Create a PETSc Krylov space linear solver
  KSP ksp;

	//  Create the vector u
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, NumberOfEntries);
  VecSetFromOptions(u);
  VecCreate(PETSC_COMM_WORLD, &b);
  VecSetSizes(b, PETSC_DECIDE, NumberOfEntries);
  VecSetFromOptions(b);
  
  //  Create the vector b to have identical size to u
  VecDuplicate(u, &b);

  VecSet(b, 0.0);
  // Set up array to access entries of b directly
  PetscScalar *bvec;
  VecGetArray(b, &bvec);

	//Not doing anything !!!!!!!!!!!!!!
  
/*
	//Set outflow BCs (prescribe flow at outlets)
	for (double j=0; j < number_edges ; j++) {
		
		//Apply outflow BCS for end branches
		if(edges_type(j)==0){
			VecSetValue(b,j+number_nodes, edges_flowrate(j) ,ADD_VALUES);
		}
		
	}	
	*/


  VecAssemblyBegin(b);VecAssemblyEnd(b);


  return b;
}

Mat Tree::make_tree_matrix ()
{
	
	//std::cout<<"Assembling Tree matrix "<<std::endl;

  //Petsc linear system example
  int NumberOfEntries = number_nodes+number_edges;

  //  Matrices and vectors in the linear system Au=b
  Mat A;

  //  Create the matrix A with 3 non-zeros per row
  MatCreateSeqAIJ(PETSC_COMM_SELF, NumberOfEntries, NumberOfEntries, 15,   PETSC_NULL, &A);
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
		
				if(edges_resistance(j-1)<0){ std::cout<< "OMG"<<std::endl;}
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
		
		
	 /*
		//If edge is the final branch
		if(edges_type(j)==0){
			//Q_j - (1/rd) P_{j} + ( (1/r_d)*(1/omega_j)*p_poro this in intro.C )
			
			MatSetValue(A, j+number_nodes, j+number_nodes ,1,ADD_VALUES);
		
			MatSetValue(A, j+number_nodes, edges_lower_node(j) , -(1.0)/edges_resistance(j)   ,ADD_VALUES); 		

			
			//changed from edges_resistance(j-1) to edges_resistance(j) -- might need to debug this !!
			
		}
		*/
		 
	 
	
  }
  
  
	
  //  Update the matrix A to reflect the changes made above
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	//Print matrix
//	MatView(A,PETSC_VIEWER_STDOUT_WORLD);
  
	
	return A;
}

void Tree::update_resistances (EquationSystems& es)
{

  const Real dt    = es.parameters.get<Real>("dt");
  const Real time    = es.parameters.get<Real>("time");
  const Real progress    = es.parameters.get<Real>("progress");

  return;
}



void Tree::update_positions (EquationSystems& es)
{

	TransientLinearImplicitSystem&  last_non_linear_soln = es.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
	const MeshBase& mesh = es.get_mesh();
	const unsigned int u_var = last_non_linear_soln.variable_number ("s_u");

	const Real time    = es.parameters.get<Real>("time");

	for (double j=0; j <number_nodes ; j++) {
		for (unsigned int d = 0; d < 3; ++d) {
			Real def_value = 0;
			try
			{
				def_value = last_non_linear_soln.point_value(u_var+d,nodes_deformed(j));
				

				
				throw 20;
			}catch (int e){ }
				
			if(def_value>0)
			{
			nodes_deformed(j)(d)=def_value;
			}	 
		}
    }
      
 return;
 
	/*
		TransientLinearImplicitSystem&  last_non_linear_soln = es.get_system<TransientLinearImplicitSystem>("Last-non-linear-soln");
			
			TransientLinearImplicitSystem&  reference =   es.get_system<TransientLinearImplicitSystem>("Reference-Configuration");

			
			//Update the position of the airway tree
			//std::cout<<"Updating the position of the tree " <<std::endl;
    	
	  const MeshBase& mesh = es.get_mesh();

	
	//First put solution on reference mesh then put back to deformed mesh
	 MeshBase::const_element_iterator       el_a     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el_a = mesh.local_elements_end(); 
    for ( ; el_a != end_el_a; ++el_a)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el_a;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = reference.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
    }
    

    es.update();
    es.allgather();
		
		
	
	
			

			
			const Real time    = es.parameters.get<Real>("time");

			for (double j=0; j <number_nodes ; j++) {
		      for (unsigned int d = 0; d < 3; ++d) {
				Real def_value = 0;
				try
				{
				  const unsigned int u_var = last_non_linear_soln.variable_number ("s_u");
				  def_value = last_non_linear_soln.point_value(u_var+d,nodes(j)) - reference.point_value(u_var+d,nodes(j));
				  throw 20;
				}catch (int e){ }
				
				if(def_value==0)
				{	
					
				 
					//  do this if we are using a lung
				  if(!es.parameters.get<std::string>("problem").compare("lung")){
				  //If point is outside then apply the known registration !
				  Real pi=3.14159;
				  Real fac=0.25*(1+sin(1*pi*time+(3.0/2.0)*pi));
				  Point diagAt(fac*0.3,fac*0.3,fac*0.5);						
				  Point bt(fac*10,fac*11,fac*45);
			    nodes_deformed(j)(d)=nodes(j)(d)+ nodes(j)(d)*diagAt(d)+bt(d) ;
					
					//Ignore the above
					nodes_deformed(j)(d)=nodes(j)(d) ;

					}
					
					 
				}else{	
				  nodes_deformed(j)(d)=def_value+nodes(j)(d) ;
				}
			  }
      }
      
     
     
     //Put the mesh back (doesn't quite seem to work!!)
     
    MeshBase::const_element_iterator       el_b     = mesh.local_elements_begin();
    const MeshBase::const_element_iterator end_el_b = mesh.local_elements_end(); 
    for ( ; el_b != end_el_b; ++el_b)
    {    
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem* elem = *el_b;
      for (unsigned int n=0; n<elem->n_nodes(); n++){
        Node *node = elem->get_node(n);
          for (unsigned int d = 0; d < 3; ++d) {
            unsigned int source_dof = node->dof_number(1, d, 0);
            Real value = last_non_linear_soln.current_local_solution->el(source_dof);
            (*node)(d)=value;
          }
      }
    }
    es.update();
    es.allgather();
 
 return;
 */
	
}
