#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>

Eigen::ArrayXd tridiag_solve(Eigen::ArrayXd &a, Eigen::ArrayXd &b, Eigen::ArrayXd &c, Eigen::ArrayXd &d) {
  int n = d.size();
  Eigen::ArrayXd x(n);
  for(int i = 1; i < n; i++) {
    double w = a[i]/b[i - 1];
    b[i] = b[i] - w*c[i - 1];
    d[i] = d[i] - w*d[i - 1];
  }
  x[n - 1] = d[n - 1]/b[n - 1];
  for(int i = n - 2; i >= 0; i--) {
    x[i] = (d[i] - c[i]*x[i + 1])/b[i];
  } 
  return x;
}

int main() {
  //grid constants
  int Tm = 10*24;
  double tau = 0.5;
  int N = Tm/tau, M = 100;
  double L = 500, h = L/M;

  //physical parameters
  double pi = 150e5, pp = 50e5;
  double k = 1e-14, mu = 1e-5/36, phi = 0.2;
  double cf = 1e-9;
  double r0 = 1e3, p0 = 120e5;
  double pinit = 100e5;

  std::vector<Eigen::ArrayXd> p(N + 1, Eigen::ArrayXd::Zero(M + 1));

  //initial state
  p[0] += pinit;

  auto r = [cf, p0, r0](double p) -> double {
    return r0*(1 + cf*(p - p0));
  };

  std::ofstream ofi("p.dat");

  //main cycle
  for(int n = 0; n < N; n++) {
    double t = tau*n, dt = 2.4, tn = dt;
    Eigen::ArrayXd a(M + 1), b(M + 1), c(M + 1), d(M + 1);
    for(int i = 1; i < M; i++) {
      //weighting r against the current
      
      double rp = (p[n](i) >= p[n](i + 1)) ? r(p[n](i)) : r(p[n](i + 1));
      double rm = (p[n](i - 1) >= p[n](i)) ? r(p[n](i - 1)) : r(p[n](i));
      std::cout << rp << " " << rm << std::endl;

      //tridiagonal system coefficient
      c(i) = k*rm/mu/h/h;
      b(i) = k*rp/mu/h/h;
      a(i) = -c(i) - b(i) - phi*cf*r0/tau;
      d(i) = -phi*cf*r0*p[n](i)/tau;
    }
    a(0) = 1; b(0) = 0; a(M) = 1; c(M) = 0;
    c(0) = 0; b(M) = 0;
    d(0) = pi; d(M) = pp;

    p[n + 1] = tridiag_solve(c, a, b, d);

    //write output 
    if(t >= tn) {
      ofi << "\"time = " << std::setprecision(2) << t/24 << " days\"" << std::endl << std::setprecision(5);
      for(int i = 0; i <= M; i++) {
        ofi << i*h << " " << p[n](i) << std::endl;
      }
      ofi << std::endl << std::endl;
      tn += dt;
    }
  }

  (void)execl("/usr/bin/gnuplot", "/usr/bin/gnuplot", "plot.gp", (char *)0);

  return 0;
}