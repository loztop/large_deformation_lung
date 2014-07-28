#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

#example_name=systems_of_equations_ex4
example_name=nonlinear_poroelasticity

options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu --implicit_neighbor_dofs"

#options="-ksp_type preonly -pc_factor_mat_solver_package mumps -pc_type lu"

#options="-ksp_view"

run_example "$example_name" "$options"
