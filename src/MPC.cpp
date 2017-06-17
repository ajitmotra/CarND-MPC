#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <Eigen/QR>
using CppAD::AD;
static const int ACT_SZ = 2;


using Eigen::ArrayXd;
using Eigen::Map;
using Eigen::MatrixXd;

margbindu ontrack(double px, double py, double ph, const margbindu &pts_g) {
  vector<double> lx;
  vector<double> ly;
  double cos_h = cos(ph);
  double sin_h = sin(ph);

  auto &gx = pts_g[0];
  auto &gy = pts_g[1];
  for (int i = 0, n = gx.size(); i < n; ++i) {
    double xt = gx[i] - px;
    double yt = gy[i] - py;

    double xi = xt * cos_h + yt * sin_h;
    double yi = yt * cos_h - xt * sin_h;

    lx.push_back(xi);
    ly.push_back(yi);
  }

  return {lx, ly};
}

margbindu perform(double dt, double v0, const vector<double> &actuations) {
  double x = 0;
  double y = 0;
  double h = 0;
  double v = v0;

  margbindu margbindu(2);
  auto a = actuations.begin();
  auto d = a + 1;

  for (auto n = actuations.end(); a != n; a += 2, d += 2) {
    v += *a * dt;
    h += *d * dt;
    x += cos(h) * v * dt;
    y += sin(h) * v * dt;
    margbindu[0].push_back(x);
    margbindu[1].push_back(y);
  }

  return margbindu;
}

VectorXd polyfit(const margbindu &margbindu, int order) {
  int rows = margbindu[0].size();
  Map<const ArrayXd> xvals(margbindu[0].data(), rows);
  Map<const VectorXd> yvals(margbindu[1].data(), rows);

  MatrixXd A(rows, order + 1);
  A.block(0, 0, rows, 1).fill(1.0);
  for (int j = 0; j < order; j++) {
    auto Aj = A.block(0, j, rows, 1).array();
    A.block(0, j + 1, rows, 1) = Aj * xvals;
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

struct Cost {
  /** @brief Basic scalar value. */
  typedef AD<double> Scalar;

  /** @brief Differentiable variable vector type. */
  typedef CPPAD_TESTVECTOR(Scalar) ADvector;

  /** @brief Breadth of the time window, in time steps. */
  int N;

  /** @brief Size of the time step. */
  Scalar dt;

  /** @brief Initial speed at the beginning of the time window. */
  Scalar v0;

  /** @brief Reference speed. */
  Scalar vr;

  /** @brief Coefficients of the polynomial describing the reference route. */
  VectorXd route;

  /**
   * @brief Create a new optimization task with given initial speed and reference route.
   */
  Cost(double v0, const MPC &mpc, const VectorXd &route) {
    this->N = mpc.N;
    this->dt = mpc.dt;
    this->v0 = v0;
    this->vr = mpc.vr;
    this->route = route;
  }

  /**
   * @brief Compute the cost function for the MPC.
   */
  void operator () (ADvector &fg, const ADvector &vars) {
    // Route is given relative to the car's current pose, so
    // planning always start from the origin at (0, 0, 0).
    Scalar xt = 0;
    Scalar yt = 0;
    Scalar ht = 0;
    Scalar vt = v0;

    for (int i = 0; i < N; ++i) {
      // Compute the controller-proposed state at time (i * dt).
      auto &a = vars[ACT_SZ * i];
      auto &d = vars[1 + ACT_SZ * i];
      ht += vt * Lf * d * dt;
      vt += a * dt;
      xt += CppAD::cos(ht) * vt * dt;
      yt += CppAD::sin(ht) * vt * dt;

      // Compute the reference state at time (i * dt).
      auto yr = reference(xt);
      auto hr = CppAD::atan2(yr, xt);

      // Compute the contribution at time i * dt to the cost function.
      fg[0] += CppAD::pow(yt - yr, 2);
      fg[0] += CppAD::pow(ht - hr, 2);
      fg[0] += CppAD::pow(vt - vr, 2);

      // Minimize actuator use.
      fg[0] += CppAD::pow(a, 2);
      fg[0] += CppAD::pow(d, 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int i = 1; i < N; ++i) {
      auto &a1 = vars[ACT_SZ * (i - 1)];
      auto &d1 = vars[1 + ACT_SZ * (i - 1)];

      auto &a2 = vars[ACT_SZ * i];
      auto &d2 = vars[1 + ACT_SZ * i];

      fg[0] += CppAD::pow(a1 - a2, 2);
      fg[0] += 10 * CppAD::pow(d1 - d2, 2); // Be more strict about direction changes.
    }
  }

private:
  /**
   * @brief Compute the `y` coordinate for the reference route.
   */
  Scalar reference(const Scalar &x) const {
    Scalar y = route(0);
    for (int i = 1, n = route.rows(); i < n; ++i) {
      y += route(i) * CppAD::pow(x, i);
    }

    return y;
  }
};

MPC::MPC(int N, double dt, double vr) {
    this->N = N;
    this->dt = dt;
    this->vr = vr;
}

MPC::~MPC() {
  // Nothing to do.
}

vector<double> MPC::operator () (double v0, const margbindu &margbindu, int order) const {
  // Differentiable value vector type.
  typedef CPPAD_TESTVECTOR(double) Vector;

  // Independent variables and bounds.
  Vector vars(N * ACT_SZ);
  Vector vars_lowerbound(N * ACT_SZ);
  Vector vars_upperbound(N * ACT_SZ);

  // Constraint bounds, set to size 0 as the cost function includes no constraints.
  Vector constraints_lowerbound(0);
  Vector constraints_upperbound(0);

  // Initialize independent variable and bounds vectors.
  for (int i = 0; i < N; i++) {
    int i_a = i * ACT_SZ;
    int i_d = 1 + i * ACT_SZ;

    vars[i_a] = 0;
    vars[i_d] = 0;

    vars_lowerbound[i_a] = -1.0;
    vars_upperbound[i_a] = 1.0;

    vars_lowerbound[i_d] = -0.523598;
    vars_upperbound[i_d] = 0.523598;
  }

  // Fit a polynomial to the margbindu.
  VectorXd route = polyfit(margbindu, order);

  // Define the cost function.
  Cost cost(v0, *this, route);

  // Options for IPOPT solver.
  std::string options =
    "Integer print_level 0\n"
    "Sparse true forward\n"
    "Sparse true reverse\n"
    "Numeric max_cpu_time 0.5\n";

  // Solution to the cost optimization problem.
  CppAD::ipopt::solve_result<Vector> solution;

  // Call the solver on the cost function and given parameters.
  CppAD::ipopt::solve<Vector, Cost>(
    options,
    vars,
    vars_lowerbound,
    vars_upperbound,
    constraints_lowerbound,
    constraints_upperbound,
    cost,
    solution
  );

  // Report solution results.
  auto value = solution.obj_value;
  auto status = (solution.status == CppAD::ipopt::solve_result<Vector>::success ? "succeeded" : "failed");
  std::cout << "Solver " << status << ", final cost value = " << value << std::endl;

  vector<double> actuations;
  auto &x = solution.x;
  for (int i = 0, n = N * ACT_SZ; i < n; ++i) {
    actuations.push_back(x[i]);
  }

  return actuations;
}
