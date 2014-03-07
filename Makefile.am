example_name  = systems_of_equations_ex4
check_SCRIPTS = run.sh
install_dir   = $(examples_install_path)/systems_of_equations/ex4
data          = systems_of_equations_ex4.C poro_main.C read_parameters.C  general_material_cc.C poro_elastic_cc.C assemble_solid.C assemble_bcs.C setup_equationsystem.C useful_functions.C poro.h general_material_cc.h poro_elastic_cc.h run.sh
sources       = $(data) run.sh

CLEANFILES = displacement.e

##############################################
# include common example environment
include $(top_srcdir)/examples/Make.common
