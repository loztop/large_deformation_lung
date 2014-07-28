	
    std::cout<< "renaming_rows" <<std::endl;
  	Real sum_entries=0;
	 rows_coupling[0] =0; 
	 for ( int i = 1; i < a_nrow; i++ )
  {
		 sum_entries=sum_entries+rows_coupling[i];
		 rows_coupling[i] = sum_entries ;
	}
  

  /*
  	//Print out
  	std::cout<<  "after  "  <<std::endl;
  	 for ( int i = 0; i < a_nrow; i++ )
  {
  
		std::cout<<  rows_coupling[i] <<std::endl;

	}
  */