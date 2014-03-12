#include "poro.h"

using namespace std;

void read_parameters(EquationSystems& es,  int& argc, char**& argv){

 
	es.parameters.set<Real> ("n_timesteps") = 4;
	es.parameters.set<Real> ("N_eles") = 3;
	es.parameters.set<std::string> ("output_file_name") = "data/test_2D.mat";
	es.parameters.set<std::string> ("result_file_name") = "traction_test";
	es.parameters.set<Real> ("DELTA") = 0.001;
 

  std::cout<<"n_timesteps "<< es.parameters.get<Real>("n_timesteps") <<" \n";
  std::cout<<"N_eles "<< es.parameters.get<Real>("N_eles") <<" \n";
  std::cout<<"output_file_name "<< es.parameters.get<std::string>("output_file_name") <<" \n";
  std::cout<<"result_file_name "<< es.parameters.get<std::string>("result_file_name") <<" \n";
  std::cout<<"DELTA "<< es.parameters.get<Real>("DELTA") <<" \n";

}


