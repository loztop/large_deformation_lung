#ifndef TREE_H_
#define TREE_H_

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/transient_system.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/fe_macro.h"
#include "libmesh/petsc_matrix.h"

// The definition of a geometric element
#include "libmesh/elem.h"
#include <petscksp.h>
#include <petsctime.h>


//#include "poro.h"

#include "libmesh/equation_systems.h"

class Tree{
public:

	void  read_tree( EquationSystems& es);
	
	
	void  write_tree(EquationSystems& es);

	void calculate_omega_j (EquationSystems& es);
	Mat make_tree_matrix( );
	Vec make_tree_rhs( );

 	void  update_resistances(EquationSystems& es);

	
 	void  update_positions(EquationSystems& es);

	
	Number number_nodes;
	Number number_edges;
	Number number_distal;
	Number N_level;

	DenseVector<Point> nodes;
	DenseVector<Point> nodes_deformed;
	DenseVector<Real> nodes_pressure;
	DenseVector<Real> nodes_type;
	DenseVector<Real> nodes_parent_node;
	DenseVector<Real> nodes_parent_edge;

	DenseVector<Point> edges;
	DenseVector<Real> edges_type;
	DenseVector<Real> edges_resistance;
	DenseVector<Real> edges_upper_node;
	DenseVector<Real> edges_lower_node;
	DenseVector<Real> edges_child1;
	DenseVector<Real> edges_child2;

	DenseVector<Real> edges_flowrate;
	DenseVector<Real> distal_resistance;
	DenseVector<Real> edges_radius;
	DenseVector<Real> omega_j;
	
	
		DenseVector<Real> edges_diseased;
	DenseVector<Real> edges_length;


	Real p_out;

};



#endif 
