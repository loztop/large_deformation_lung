#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

#example_name=systems_of_equations_ex4
example_name=whole_lung_large_deformation

options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu"

#options="-ksp_view"

run_example "$example_name" "$options"
