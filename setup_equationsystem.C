#include "poro.h"

using namespace std;

#define ELEMENT_TYPE LAGRANGE 
#define ELEMENT_TYPE_PRESS MONOMIAL

#define ORDER_HIGH FIRST
#define ORDER_LOW CONSTANT 


void setup_equationsystem(EquationSystems& equation_systems) {

TransientLinearImplicitSystem & last_non_linear_soln_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Last-non-linear-soln");
last_non_linear_soln_system.add_variable ("s_u", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("s_v", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("s_w", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("s_p", ORDER_LOW,ELEMENT_TYPE_PRESS);

 TransientLinearImplicitSystem & ref_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Reference-Configuration");
ref_system.add_variable ("u_ref", ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("v_ref", ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("w_ref", ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("p_ref",ORDER_LOW,ELEMENT_TYPE_PRESS);

TransientLinearImplicitSystem & newton_update_system = equation_systems.add_system<TransientLinearImplicitSystem> ("Newton-update");
newton_update_system.add_variable ("u_nu", ORDER_HIGH,ELEMENT_TYPE);
newton_update_system.add_variable ("v_nu", ORDER_HIGH,ELEMENT_TYPE);
newton_update_system.add_variable ("w_nu", ORDER_HIGH,ELEMENT_TYPE);
newton_update_system.add_variable ("p_nu", ORDER_LOW,ELEMENT_TYPE_PRESS);
newton_update_system.attach_assemble_function (assemble_solid);

newton_update_system.add_variable ("x_nu", ORDER_HIGH,ELEMENT_TYPE);
newton_update_system.add_variable ("y_nu", ORDER_HIGH,ELEMENT_TYPE);
newton_update_system.add_variable ("z_nu", ORDER_HIGH,ELEMENT_TYPE);

last_non_linear_soln_system.add_variable ("mono_f_vel_u", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("mono_f_vel_v", ORDER_HIGH,ELEMENT_TYPE);
last_non_linear_soln_system.add_variable ("mono_f_vel_w", ORDER_HIGH,ELEMENT_TYPE);

ref_system.add_variable ("x_ref",ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("y_ref",ORDER_HIGH,ELEMENT_TYPE);
ref_system.add_variable ("z_ref",ORDER_HIGH,ELEMENT_TYPE);

}

