#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

#example_name=systems_of_equations_ex4
example_name=nonlinear_poroelasticity

options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu --implicit_neighbor_dofs"

run_example "$example_name" "$options"


#options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu"

#options="-ksp_view"



data_dir="plot_data/"


#####################

#fname="5036_Ce"
#AD=1
#TD=0
#loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
#run_example "$example_name" "$loz_opt" "$options"

AD=0.5
TD=0
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0.25
TD=0
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0.1
TD=0
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0.01
TD=0
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"


fname="5036_We"
AD=0
TD=1
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0
TD=0.5
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0
TD=0.25
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0
TD=0.1
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"

AD=0
TD=0.01
loz_opt="$AD $TD $data_dir$fname"_"$AD"_"$TD"_""
run_example "$example_name" "$loz_opt" "$options"