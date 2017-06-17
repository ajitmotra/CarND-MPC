#include <Eigen/Dense>
#include <vector>
#ifndef MPC_H
#define MPC_H

using namespace std;
using Eigen::VectorXd;
typedef vector<vector<double>> margbindu;
const double Lf = 0.1;
// Our favourite function for 2nd term as well
VectorXd polyfit(const margbindu &margbindu, int order);
// convert global --> local
margbindu ontrack(double px, double py, double ph, const margbindu &pts_g);

// compute way point from the actuations
margbindu perform(double dt, double v0, const vector<double> &actuations);


struct MPC {
  // Time window width
  int N;

  // Delta time
  double dt;

  // Speed
  double vr;

// New MPC controller over time window, step size and rference
  MPC(int N, double dt, double vr);

  virtual ~MPC();

//Calc sequence of activation to reach  
vector<double> operator () (double v0, const margbindu &margbindu, int order = 2) const;
};

#endif /* MPC_H */
