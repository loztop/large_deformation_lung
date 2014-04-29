/*
tree.write_tree(equation_systems);

#ifdef LIBMESH_HAVE_EXODUS_API
 if ((t_step+1)%write_interval == 0)
 { 
  #if WRITE_MESH
  GmshIO(mesh).write(mesh_out_file_name);
  #endif

  std::stringstream file_name;
  file_name << result_file_name;
  file_name << std::setw(2) << std::setfill('0') << t_step;
  file_name << ".e-s.";
  file_name << std::setw(3) << std::setfill('0') << t_step+1;
  
  exo.write_timestep(file_name.str(), equation_systems,t_step+1,time);
	exo.write_element_data(equation_systems);
  std::cout<<"Wrote "<< file_name.str() <<std::endl;


  #if WRITE_TEC
  std::stringstream file_name_tec;
  file_name_tec <<result_file_name << "_"<< t_step<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
  std::cout<<"Wrote "<< file_name_tec.str() <<std::endl;
  #endif



 }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
*/


    std::stringstream file_name;
    file_name << equation_systems.parameters.get<std::string>("result_file_name");
    file_name << std::setw(2) << std::setfill('0') << t_step;
    file_name << ".e-s.";
    file_name << std::setw(3) << std::setfill('0') << t_step;
		//exo.write_timestep(file_name.str(), equation_systems,t_step+1,time); 
		exo.write_timestep(file_name.str(), equation_systems,t_step+1,t_step); 
    std::cout<<"exodus "<< file_name.str() <<std::endl;
		exo.write_element_data(equation_systems);
   
  tree.write_tree(equation_systems);
