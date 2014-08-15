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
		es.parameters.set<Real> ("E") =8000;
    es.parameters.set<Real> ("NU") = 0.3;
    es.parameters.set<Real> ("KPERM") =1.e-5;
		es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	
		es.parameters.set<std::string> ("result_file_name") = "data/cubeX";
		es.parameters.set<std::string> ("tree_file_name") = "data/cubeX";
	
		es.parameters.set<std::string> ("mesh_input") = "cube";
			
	 // es.parameters.set<std::string> ("tree_input") = "meshes/tree/cube_tree_single";

	es.parameters.set<std::string> ("tree_input") = "meshes/tree/cube_tree";

		
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
 
 
	
  
   if(!es.parameters.get<std::string>("problem").compare("lung")){
	 
	 
	 
	 if(argc>8){
	   std::cout<<"Loading aguments "<<std::endl;
	   es.parameters.set<Real> ("airway_disease") =  atof( argv[1] );
		es.parameters.set<Real> ("tissue_disease") = atof( argv[2] );
	es.parameters.set<std::string> ("result_file_name") =argv[3] ;
	es.parameters.set<std::string> ("tree_file_name") = argv[3] ;
	
	   
	 }else{
	   
	   es.parameters.set<Real> ("airway_disease") = 0;
		es.parameters.set<Real> ("tissue_disease") = 0;
		es.parameters.set<std::string> ("output_file_name") = "data/effective_test.mat";
	es.parameters.set<std::string> ("result_file_name") = "data/effective_test";
	es.parameters.set<std::string> ("tree_file_name") = "data/effective_test";
	
	   }
	
	
    es.parameters.set<Real> ("end_time") =2;
	es.parameters.set<Real> ("n_timesteps") =10;
	es.parameters.set<Real> ("N_eles") = 4;	
	es.parameters.set<Real> ("DELTA") = 0.00001;

	//Add the registration information N048 exp -> insp
	Point A( 0.1559, 0.0393, 0.2603);
	es.parameters.set<Point> ("A") =A;
	Point b(-4.3526,-3.2864,-12.8055);
	es.parameters.set<Point> ("b") =b;

	
		

	
	Point disease_cent(-60,-170,-90);
	es.parameters.set<Point>("disease_cent")=disease_cent;
	es.parameters.set<Real>("rad_dis")=25;
		
	es.parameters.set<Real> ("E") =8000;
    es.parameters.set<Real> ("NU") = 0.3;
    es.parameters.set<Real> ("KPERM") =1.e-5;
  
	
	
	
 	//Trees for for APLE36(266)
  //es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_tree_APLE36P3";
  //es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_tree_APLE36P10";
	 
///////////////////////////////////////////////////////////////		
//Lobe meshes for N048
 //  es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_246.msh";
//    es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_751.msh";
 es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N48rmerge_290.msh";

 // es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048r_fine1447.msh";
  //  es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048r_fine2881.msh";
  // es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048_node5036.msh";

	 
	  //Trees for N048
  	///////	es.parameters.set<std::string> ("tree_input") = "meshes/tree/single_tree";
//	 es.parameters.set<std::string> ("tree_input") = "meshes/tree/IM4branch";
  es.parameters.set<std::string> ("tree_input") = "meshes/tree/OTHERminIM";
// es.parameters.set<std::string> ("tree_input") = "meshes/tree/simple_tree";
  // es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_treeP5";
 // es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_treeP3";
	 //  es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_tree";

	 
 }
 
 
 
  

 
  
  

}