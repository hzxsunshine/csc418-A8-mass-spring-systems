#include "signed_incidence_matrix_dense.h"

void signed_incidence_matrix_dense(
  const int n,
  const Eigen::MatrixXi & E,
  Eigen::MatrixXd & A)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  // each row is a spring
  A = Eigen::MatrixXd::Zero(E.rows(),n);

  //   n  number of vertices (#V)
  //   E  #E by 2 list of edge indices into rows of V
  //   A  #E by n signed incidence matrix

  // Eigen::Matrix3d Id = Eigen::Matrix3d::Identity(); // might be needed in the future.
  for (int i = 0; i < E.rows(); ++i){
      A(i, E(i, 0)) += 1;
      A(i, E(i, 1)) += -1;
  }


  //////////////////////////////////////////////////////////////////////////////
}
