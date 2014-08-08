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
  
		
		if (argc >8){
		  
		std::cout<<"Readig in arguments "<<std::endl;
		std::cout<<"argc   "<<argc << std::endl;

		es.parameters.set<Real> ("airway_disease") = atof( argv[1] );
		es.parameters.set<Real> ("tissue_disease")= atof( argv[2] );
		es.parameters.set<std::string> ("output_file_name") = argv[3] ;
		es.parameters.set<std::string> ("result_file_name") = argv[3] ;
		es.parameters.set<std::string> ("tree_file_name") = argv[3] ;
	
	
		}else{

	es.parameters.set<Real> ("tissue_disease")=0;
	es.parameters.set<Real> ("airway_disease")=0;

	es.parameters.set<std::string> ("output_file_name") = "data/Const3174_0.00001fix_12NT2T_";
	es.parameters.set<std::string> ("result_file_name") = "data/Const3174_0.00001fix_12NT2T_";
	es.parameters.set<std::string> ("tree_file_name") = "data/Const3174_0.00001fix_12NT2T_";
	
	
}


std::cout<<"AD = "<< es.parameters.get<Real> ("airway_disease") <<std::endl;
std::cout<<"TD = "<< es.parameters.get<Real> ("tissue_disease") <<std::endl;

	es.parameters.set<Real> ("end_time") =2;
	es.parameters.set<Real> ("n_timesteps") =12;
	es.parameters.set<Real> ("N_eles") = 4;	
	es.parameters.set<Real> ("DELTA") = 0.00001;
	
	
	
	
	
	

	
	Point disease_cent(-60,-100,120);
	es.parameters.set<Point> ("disease_cent") =disease_cent;
	es.parameters.set<Real> ("rad_dis") =30;

	es.parameters.set<Real> ("E") =8000;
	es.parameters.set<Real> ("NU") = 0.3;
	es.parameters.set<Real> ("KPERM") =1.e-5;
  
	

	
	

	
	//Lobe meshes for N048 /////////////
 //es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_246.msh";
  //  es.parameters.set<std::string> ("mesh_input") = "meshes/lung/whole_lung_751.msh";
 // es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048r_fine1447.msh";
 //   es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048r_fine2881.msh";
 //  es.parameters.set<std::string> ("mesh_input") = "meshes/lung/N048_node6598.msh";

	 
	  //Trees for N048 /////////////////
	// es.parameters.set<std::string> ("tree_input") = "meshes/tree/IM4branch";
  // es.parameters.set<std::string> ("tree_input") = "meshes/tree/OTHERminIM";
	//  es.parameters.set<std::string> ("tree_input") = "meshes/tree/simple_tree";
  // es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_treeP5";
  //es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_treeP3";
	//   es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_tree";
	  
	//Point A( 0.0743, 0.0561, 0.2782);
	//es.parameters.set<Point> ("A") =A;
	//Point b(-0.776,4.5526,-56.2963);
	//es.parameters.set<Point> ("b") =b;
	
	
	
	
	  
	/////////////////////////////////////////////////// A65
	////Lobe meshes for A65
	//	es.parameters.set<std::string> ("mesh_input") = "meshes/lung/A65rmerge_insp_198.msh";
		es.parameters.set<std::string> ("mesh_input") = "meshes/lung/A65rmerge_insp_356.msh";
	//es.parameters.set<std::string> ("mesh_input") = "meshes/lung/A65rmerge_insp_1339.msh";

	//es.parameters.set<std::string> ("mesh_input") = "meshes/lung/A65rmerge_insp2233.msh";
	//es.parameters.set<std::string> ("mesh_input") = "meshes/lung/A65rmerge_inspfine3174.msh";
		 
	//Trees for A65
es.parameters.set<std::string> ("tree_input") = "meshes/tree/half_tree_A65P3";
	


//Add the registration information A65 exp -> insp
	Point A( 0.0673, 0.0511, 0.2782);
	es.parameters.set<Point> ("A") =A;
	Point b(-0.7076,4.1526,-56.2963);
	es.parameters.set<Point> ("b") =b;
	
 }
 
 
 
  

 
  
  

}


