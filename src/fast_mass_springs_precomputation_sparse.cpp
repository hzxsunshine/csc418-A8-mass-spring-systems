#include "fast_mass_springs_precomputation_sparse.h"
#include "signed_incidence_matrix_sparse.h"
#include <vector>

bool fast_mass_springs_precomputation_sparse(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXd & m,
  const Eigen::VectorXi & b,
  const double delta_t,
  Eigen::VectorXd & r,
  Eigen::SparseMatrix<double>  & M,
  Eigen::SparseMatrix<double>  & A,
  Eigen::SparseMatrix<double>  & C,
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > & prefactorization)
{
  /////////////////////////////////////////////////////////////////////////////
  // Replace with your code

  int n = V.rows();
  double w = 1e10;

  // output r
  r.resize(E.rows());
  Eigen::RowVector3d V1, V2;
  for (int j = 0; j < E.rows(); ++j) {
      V1 = V.row(E(j, 0));
      V2 = V.row(E(j, 1));
      r[j] = (V1 - V2).norm();
  }

  // output M
  M.resize(n, n);
  std::vector<Eigen::Triplet<double> > tri_M;
  tri_M.reserve(V.rows());
  for (int l = 0; l < V.rows(); ++l) {
      tri_M.emplace_back(l, l, m[l]);
  }
  M.setFromTriplets(tri_M.begin(), tri_M.end());


  // output A:
  signed_incidence_matrix_sparse(n, E, A);

  // output C:
  C.resize(b.rows(), n);
  std::vector<Eigen::Triplet<double> > tri_C;
  tri_C.reserve(b.rows());
  for (int i1 = 0; i1 < b.rows(); ++i1) {
      tri_C.emplace_back(i1, b(i1), 1);
  }
  C.setFromTriplets(tri_C.begin(), tri_C.end());

  Eigen::SparseMatrix<double> W(n, n);
  W = w * C.transpose() * C;

  Eigen::SparseMatrix<double> Q(n,n);
  Q = (k * A.transpose() * A) + (M / pow(delta_t, 2)) + W;
  /////////////////////////////////////////////////////////////////////////////
  prefactorization.compute(Q);

  return prefactorization.info() != Eigen::NumericalIssue;
}
