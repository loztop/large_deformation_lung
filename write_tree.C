//#include "main_header.cpp"
#include "tree.h"
#include <math.h>
//#include "assemble.h"
#include <sstream>
#include <string>
#include <petscksp.h>
//#include <petsc_linear_solver.h>
#include <iostream>
#include <fstream>


void Tree::write_tree(EquationSystems& es) {

  std::string output_file_name;
  std::ostringstream convert;   // stream used for the conversion
  
  //timestep
  convert << es.parameters.get<unsigned int>("step");  
    
  output_file_name = es.parameters.get<std::string>("tree_file_name") + convert.str() + '.' + "vtk";
  
	  
  std::ofstream outFile;
  outFile.open (output_file_name.c_str());
  std::cout<<"Writing "<<output_file_name<<std::endl;
  
  outFile << "# vtk DataFile Version 3.0" << "\n";
  outFile << "vtk output" << "\n";
  outFile << "ASCII" << "\n";
  outFile << "DATASET POLYDATA" << "\n";
  
 
  outFile << "POINTS " << number_nodes<<" float"<<"\n"; 
  for (int i=0; i < number_nodes; i++) {
	  outFile <<nodes_deformed(i)(0) << " "<<nodes_deformed(i)(1)<< " "<<nodes_deformed(i)(2)<<"\n";
  }

	
  outFile << "LINES " << number_edges<<" "<<  number_edges*3  <<"\n"; 
  for (int i=0; i < number_edges; i++) {
	  outFile <<2 << " "<<edges_upper_node(i)<< " "<<edges_lower_node(i)<<"\n";
  }
	
	
  outFile << "CELL_DATA " <<  number_edges  <<"\n"; 
  outFile << "scalars Radius float"  <<"\n"; 
  outFile << "LOOKUP_TABLE default"  <<"\n"; 
  for (int i=0; i < number_edges; i++) {
	  outFile <<edges_radius(i)<< " ";
  }
 	
  outFile <<"\n"; 

  outFile << "scalars Flowrate float"  <<"\n"; 
  outFile << "LOOKUP_TABLE default"  <<"\n"; 
  for (int i=0; i < number_edges; i++) {
	  outFile <<edges_flowrate(i)<< " ";
  }
  outFile <<"\n"; 
  
  outFile << "scalars Airway_resistance float"  <<"\n"; 
  outFile << "LOOKUP_TABLE default"  <<"\n"; 
  for (int i=0; i < number_edges; i++) {
	  outFile <<edges_total_resistance(i)<< " ";
  }
  outFile <<"\n"; 
	
  outFile << "POINT_DATA " << number_nodes  <<"\n"; 
  outFile << "scalars Pressure float"  <<"\n"; 
  outFile << "LOOKUP_TABLE default"  <<"\n"; 
  for (int i=0; i < number_nodes; i++) {
	  outFile <<nodes_pressure(i)<< " ";
  }
  outFile <<"\n"; 
 
  
  outFile.close();
  	
}


