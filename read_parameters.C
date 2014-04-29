#include "poro.h"

using namespace std;

void read_parameters(EquationSystems& es,  int& argc, char**& argv){

  
 es.parameters.set<std::string> ("problem") = "lung";

		
		
 if(!es.parameters.get<std::string>("problem").compare("cylinder")){
   	es.parameters.set<Real> ("end_time") =1000;
	es.parameters.set<Real> ("n_timesteps") =100;
	es.parameters.set<Real> ("N_eles") = 4;
	es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	es.parameters.set<std::string> ("result_file_name") = "cylinder_sym5567";
	es.parameters.set<Real> ("DELTA") = 0.001;
 }
 

	
	
 if(!es.parameters.get<std::string>("problem").compare("cube")){
    es.parameters.set<Real> ("end_time") =100;
	es.parameters.set<Real> ("n_timesteps") =5;
	es.parameters.set<Real> ("N_eles") = 4;
	es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	es.parameters.set<std::string> ("result_file_name") = "data/cube";
	es.parameters.set<Real> ("DELTA") = 0.001;
 }

	
  
   if(!es.parameters.get<std::string>("problem").compare("lung")){
    es.parameters.set<Real> ("end_time") =10;
		es.parameters.set<Real> ("n_timesteps") =5;
		es.parameters.set<Real> ("N_eles") = 3;	
		es.parameters.set<Real> ("DELTA") = 0.01;
		es.parameters.set<Real> ("DELTA_BC") = 1000000;

		es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	
		es.parameters.set<std::string> ("result_file_name") = "data/BTig";
		es.parameters.set<std::string> ("tree_file_name") = "data/BTig";
	
	
    //Old patient //Whole tree pruned down to imaging data	
		// es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_751.msh";
		//	es.parameters.set<std::string> ("tree_input") = "meshes/tree/OTHERminIM";
	
		//Coarse setup
    es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N051_node272.msh";
		es.parameters.set<std::string> ("tree_input") = "meshes/tree/lozout_n051_g6";

		
			//Cube test problem
		es.parameters.set<std::string> ("tree_input") = "meshes/tree/cube_tree";

		
		
		
 }

  
  

}


