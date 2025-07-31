#include <Rcpp.h>
using namespace Rcpp;

// Helper function equivalent to R's smooth_max0
inline double smooth_max0(double x, double eps = 1e-9) {
  return 0.5 * (x + sqrt(x*x + eps*eps));
}

// Helper function equivalent to R's smooth_switch
inline double smooth_switch(double x, double k = 1e-3) {
  return 1.0 / (1.0 + exp(-x / k));
}

// [[Rcpp::export]]
List yang_rhs_cpp(double t, NumericVector y, NumericVector p) {
  // State variables
  double NL = y[0];
  double ND = y[1];
  double G  = y[2];
  
  // Parameters
  double kp    = p[0];
  double Kp    = p[1];
  double kdmax = p[2];
  double Kd50  = p[3];
  double nd    = p[4];
  double kbys  = p[5];
  double vglc  = p[6];
  
  // Model equations
  double Gp = smooth_max0(G);
  double mu = kp * Gp / (Kp + Gp);
  double kS = kdmax * pow(Kd50, nd) / (pow(Gp, nd) + pow(Kd50, nd));
  double kB = kbys * ND;
  
  double dNL = (mu - kS - kB) * NL;
  double dND = (kS + kB) * NL;
  double gate = smooth_switch(G, 1e-3);
  double dG = -vglc * NL * gate;
  
  return List::create(NumericVector::create(dNL, dND, dG));
}
