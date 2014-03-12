#ifndef ASSEMBLE_H_
#define ASSEMBLE_H_


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

#include <petscksp.h>
#include <petsctime.h>

#include <cmath>
#include <time.h>


#include "tree.h"


//Material parameters
#define KPERM 0.001

#define E_mod 1000

#define NU_mod 0.3

#define PHI_ZERO 0.9

#define RHO_S 1


// Bring in everything from the libMesh namespace
using namespace libMesh;



void assemble_solid (EquationSystems& es,
                      const std::string& system_name);

void assemble_bcs (EquationSystems& es);

PetscMatrix<Number> assemble_coupled_stiffness (EquationSystems& es,Tree& tree, Mesh& mesh);

PetscVector<Number> assemble_coupled_rhs (EquationSystems& es,Tree& tree, Mesh& mesh);


void read_parameters(EquationSystems& es, int& argc, char**& argv) ;

void test(int a);

void setup_equationsystem(EquationSystems& equation_systems);

void calculate_numeric_jacobian(EquationSystems& es, SparseMatrix< Number >& num_jac_matrix );

double diffclock(clock_t clock1,clock_t clock2);



Mat create_big_matrix(Mat A, Mat T);



//void setup_projection_es(EquationSystems& projection_es,EquationSystems& equation_systems);

//void copy_es_to_projection_es(EquationSystems& projection_es,EquationSystems& equation_systems);

//void read_options(unsigned int &  n_timesteps, unsigned int &  N_eles, std::string& result_file_name,int& argc, char**& argv) ;


//Point get_expanding_sphere_bcs(EquationSystems& es, const Elem* elem, int n,double scale);
//Point constrain_tet_nodes(EquationSystems& es, const Elem* elem, int n);

//void verify_jack(EquationSystems& es);


  //void get_traction(DenseVector<Real> & traction, Point rX, const Real progress);

  //void get_traction_current(DenseVector<Real> & traction, Point rX, const Real progress);

 // void get_traction_test(DenseVector<Real> & traction, Point rX, const Real progress);


 // void get_bodyforce(DenseVector<Real> & body_force, Point rX, const Real progress);

  //void tensor_mult_vector(DenseVector<Real> & ans, RealTensor tens, Point normal);

template <typename T> TypeTensor<T> inv(const TypeTensor<T> &A ) {
  double Sub11, Sub12, Sub13;
  Sub11 = A._coords[4]*A._coords[8] - A._coords[5]*A._coords[7];
  Sub12 = A._coords[3]*A._coords[8] - A._coords[6]*A._coords[5];
  Sub13 = A._coords[3]*A._coords[7] - A._coords[6]*A._coords[4];
  double detA = A._coords[0]*Sub11 - A._coords[1]*Sub12 + A._coords[2]*Sub13;
  libmesh_assert( std::fabs(detA)>1.e-15 );
  TypeTensor<T> Ainv(A);
  Ainv._coords[0] =  Sub11/detA;
  Ainv._coords[1] = (-A._coords[1]*A._coords[8]+A._coords[2]*A._coords[7])/detA;
  Ainv._coords[2] = ( A._coords[1]*A._coords[5]-A._coords[2]*A._coords[4])/detA;
  Ainv._coords[3] = -Sub12/detA;
  Ainv._coords[4] = ( A._coords[0]*A._coords[8]-A._coords[2]*A._coords[6])/detA;
  Ainv._coords[5] = (-A._coords[0]*A._coords[5]+A._coords[2]*A._coords[3])/detA;
  Ainv._coords[6] =  Sub13/detA;
  Ainv._coords[7] = (-A._coords[0]*A._coords[7]+A._coords[1]*A._coords[6])/detA;
  Ainv._coords[8] = ( A._coords[0]*A._coords[4]-A._coords[1]*A._coords[3])/detA;

  return Ainv;
}

#endif 
