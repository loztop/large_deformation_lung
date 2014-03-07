#ifdef LIBMESH_HAVE_EXODUS_API
 if ((t_step+1)%write_interval == 0)
 { 
  #if WRITE_MESH
  GmshIO(mesh).write(mesh_out_file_name);
  #endif

  std::stringstream file_name;
  file_name << "data/S_";  
  file_name << result_file_name;
  file_name << std::setw(2) << std::setfill('0') << t_step;
  file_name << ".e-s.";
  file_name << std::setw(3) << std::setfill('0') << t_step+1;
  
  exo.write_timestep(file_name.str(), equation_systems,t_step+1,time);
	exo.write_element_data(equation_systems);
  std::cout<<"Wrote "<< file_name.str() <<std::endl;


  #if WRITE_TEC
  std::stringstream file_name_tec;
  file_name_tec << "data/"<<result_file_name << "_"<< t_step<< ".tec" ;
  tec.write_equation_systems (file_name_tec.str(),equation_systems);
  std::cout<<"Wrote "<< file_name_tec.str() <<std::endl;
  #endif

  #if WRITE_VTK
  std::stringstream file_name_vtk;
  file_name_vtk << "VTK_SIM"<< t_step<< ".nem" ;;
  //vtk.write_equation_systems (file_name_vtk.str(),equation_systems);
  #endif


 }
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
