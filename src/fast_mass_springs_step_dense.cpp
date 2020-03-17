#include "fast_mass_springs_step_dense.h"
#include <igl/matlab_format.h>

// Conduct a single step of the "Fast Simulation of Mass-Spring Systems" method.
//
// Inputs:
//   V  #V by 3 list of **rest** vertex positions
//   E  #E by 2 list of edge indices into rows of V
//   k  spring stiffness
//   b  #b list of indices of fixed vertices as indices into rows of V
//   delta_t  time step in seconds
//   fext  #V by 3 list of external forces
//   r  #E list of edge lengths
//   M  #V by #V mass matrix
//   A  #E by #V signed incidence matrix
//   C  #b by #V selection matrix
//   prefactorization  LLT prefactorization of energy's quadratic matrix
//   Uprev  #V by 3 list of previous vertex positions (at time t-∆t)
//   Ucur  #V by 3 list of current vertex positions (at time t)
// Outputs:
//   Unext #V by 3 list of next vertex positions (at time t+∆t)

void fast_mass_springs_step_dense(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & E,
  const double k,
  const Eigen::VectorXi & b,
  const double delta_t,
  const Eigen::MatrixXd & fext,
  const Eigen::VectorXd & r,
  const Eigen::MatrixXd & M,
  const Eigen::MatrixXd & A,
  const Eigen::MatrixXd & C,
  const Eigen::LLT<Eigen::MatrixXd> & prefactorization,
  const Eigen::MatrixXd & Uprev,
  const Eigen::MatrixXd & Ucur,
  Eigen::MatrixXd & Unext)
{
  //////////////////////////////////////////////////////////////////////////////
  // Replace with your code
  const Eigen::MatrixXd& cur = Ucur;
  const Eigen::MatrixXd& pre = Uprev;

  double w = 1e10;

  Eigen::MatrixXd y;
  y.resize(V.rows(), 3);
  y = M * (2 * cur - pre) / pow(delta_t, 2) + fext;

  Eigen::MatrixXd d;
  Eigen::RowVector3d V1;
  Eigen::RowVector3d V2;
  d.resize(E.rows(), 3);

  // w * C^T * C * p_rest for pinned vertices
  Eigen::MatrixXd f_b;
  f_b.resize(V.rows(), 3);

  Unext = cur;

  for(int iter = 0;iter < 50;iter++)
  {
      // Step 1 (local): Given current values of p determine dij for each spring
      // this is the dij at t+1
      // for each spring;
      for (int i = 0; i < E.rows(); ++i)
      {
          double length  = r[i];
          V1 = Unext.row(E(i, 0));
          V2 = Unext.row(E(i, 1));
          d.row(i) = (V1 - V2).normalized() * length;
      }
      f_b = k * A.transpose() * d + y + w * C.transpose() * C * V;

      // Step 2 (global): Given all dij vectors, find positions p that minimize quadratic energy $\tilde{E}$.
      Unext = prefactorization.solve(f_b);
  }
  //////////////////////////////////////////////////////////////////////////////
}
