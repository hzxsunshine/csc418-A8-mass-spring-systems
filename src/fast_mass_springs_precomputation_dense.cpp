#include "fast_mass_springs_precomputation_dense.h"
#include "signed_incidence_matrix_dense.h"
#include <Eigen/Dense>
// Inputs:
//   V  #V by 3 list of vertex positions
//   E  #E by 2 list of edge indices into rows of V
//   k  spring stiffness
//   m  #V list of masses
//   b  #b list of "pinned"/fixed vertices as indices into rows of V
//   delta_t  time step in seconds
// Outputs:
//   r  #E list of edge lengths
//   M  #V by #V mass matrix
//   A  #E by #V signed incidence matrix
//   C  #b by #V selection matrix
//   prefactorization  LLT prefactorization of energy's quadratic matrix
bool fast_mass_springs_precomputation_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::MatrixXd & M,
  Eigen::MatrixXd & A,
  Eigen::MatrixXd & C,
  Eigen::LLT<Eigen::MatrixXd> & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  int n = V.rows(); // #V
  double w = 1e10;

  // output r:
  r.resize(E.rows());
  Eigen::RowVector3d V1, V2;
  for (int j = 0; j < E.rows(); ++j) {
      V1 = V.row(E(j, 0));
      V2 = V.row(E(j, 1));
      r[j] = (V1 - V2).norm();
  }

  // output M:
  M = Eigen::MatrixXd::Zero(n, n);
  for (int l = 0; l < V.rows(); ++l) {
      M(l, l) = m[l];
  }

  // output A:
  signed_incidence_matrix_dense(n, E, A);

  // output C:
  C = Eigen::MatrixXd::Zero(b.rows(), n);
  for (int i = 0; i < b.rows(); ++i) {
      C(i, b(i)) = 1;
  }

  // prefactorization
  Eigen::MatrixXd Q;
  // penalty energy for fixed vertices
  Eigen::MatrixXd W = w * C.transpose() * C;
  /////////////////////////////////////////////////////////////////////////////
  Q = (k * A.transpose() * A) + (M / pow(delta_t, 2)) + W;
  prefactorization.compute(Q);

  return prefactorization.info() != Eigen::NumericalIssue;
}
