#include "poro.h"

using namespace std;

void read_parameters(EquationSystems& es,  int& argc, char**& argv){

  
 es.parameters.set<std::string> ("problem") = "lung";

		
	
 if(!es.parameters.get<std::string>("problem").compare("cube")){
     es.parameters.set<Real> ("end_time") =1;
		es.parameters.set<Real> ("n_timesteps")=5;
		es.parameters.set<Real> ("N_eles") = 4;	
		es.parameters.set<Real> ("DELTA") = 0.001;
		es.parameters.set<Real> ("DELTA_BC") = 1000000000;

		es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	
		es.parameters.set<std::string> ("result_file_name") = "data/cubeX";
		es.parameters.set<std::string> ("tree_file_name") = "data/cubeX";
	
	 //es.parameters.set<std::string> ("tree_input") = "meshes/tree/cube_tree_single";

	es.parameters.set<std::string> ("tree_input") = "meshes/tree/cube_tree";

		
 }

	
  
   if(!es.parameters.get<std::string>("problem").compare("lung")){
    es.parameters.set<Real> ("end_time") =2;
		es.parameters.set<Real> ("n_timesteps") =100;
		es.parameters.set<Real> ("N_eles") = 4;	
		es.parameters.set<Real> ("DELTA") = 0.001;
		es.parameters.set<Real> ("DELTA_BC") = 10000000000;

		es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	
		es.parameters.set<std::string> ("result_file_name") = "data/751LUNGFULLlong";
		es.parameters.set<std::string> ("tree_file_name") = "data/751LUNGFULLlong";
	
	
    //Old patient //Whole tree pruned down to imaging data	

	//		es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_246.msh";
		
				es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_751.msh";

				
		// es.parameters.set<std::string> ("tree_input") = "meshes/tree/single_tree";
	
		
	//	es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048r_fine1447.msh";

		
	 //es.parameters.set<std::string> ("tree_input") = "meshes/tree/simple_tree";
		 
	//	es.parameters.set<std::string> ("tree_input") = "meshes/tree/IM4branch";
		
		 //es.parameters.set<std::string> ("mesh_input") ="meshes/lung/whole_lung_751.msh";
	
	 es.parameters.set<std::string> ("tree_input") = "meshes/tree/OTHERminIM";
	
		
		
		//Coarse setup "Broken tree"
 //    es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N051_node272.msh";
	//	es.parameters.set<std::string> ("tree_input") = "meshes/tree/lozout_n051_g6";

		
		
 }
 
 
 
  
   if(!es.parameters.get<std::string>("problem").compare("cylinder")){
    es.parameters.set<Real> ("end_time") =1;
		es.parameters.set<Real> ("n_timesteps") =7;
		es.parameters.set<Real> ("N_eles") = 4;	
		es.parameters.set<Real> ("DELTA") = 0.001;
		es.parameters.set<Real> ("DELTA_BC") = 100000000;

		es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	
		es.parameters.set<std::string> ("result_file_name") = "data/CYLfixF";
		es.parameters.set<std::string> ("tree_file_name") = "data/CYLfixF";
	
	
    //Old patient //Whole tree pruned down to imaging data	
	 	es.parameters.set<std::string> ("mesh_input") = "meshes/cylinder_sym133.msh";
	 	
		es.parameters.set<std::string> ("tree_input") = "meshes/tree/single_tree_cyl";
	 
		//es.parameters.set<std::string> ("tree_input") = "meshes/tree/branch_tree_cyl";

		
 }
 
  
  

}


