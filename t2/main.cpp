#include <iostream>
#include <fstream>
#include <array>
#include <eigen3/Eigen/Dense>
#include <float.h>
#include <cmath>
#include <unistd.h>

int main() 
{
  //gas parameters
  double g = 5./3.;
  double vl = 0., vr = 0.;
  double rl = 13., rr = 1.3;
  double pl = 1e6, pr = 1e5;
  double L = 10., T = 0.02;

  //grid parameters
  double CFL_max = 0.01;
  const int N = 100;
  double h = 2*L/N;
  double t = 0, tau_init = 1e-6, tau = tau_init;

  std::array<Eigen::Vector3d, N> wn, wnn;
  std::array<Eigen::Vector3d, N> vn, vnn;
  std::array<double, N> pn, pnn;
  std::array<Eigen::Matrix3d, N> omtn, An, omtmn, Ln;

  //initial conditions
  for(int i = 0; i < N/2; i++) {
    vnn[i] << rl, vl, pl/rl/(g - 1.0);
    wnn[i] << rl, 0, pl/(g - 1.);
    pnn[i] = pl;
  }
  for(int i = N/2; i < N; i++) {
    vnn[i] << rr, vr, pr/rr/(g - 1.);
    wnn[i] << rr, 0, pr/(g - 1.);
    pnn[i] = pr;
  }

  std::ofstream anim("anim.dat");
  double nt = 0.001;

  while(t < T - tau) 
  {
    //compute matrices
    for(int i = 0; i < N; i++) {
      double c = sqrt(g*(g - 1)*vnn[i](2));
      An[i] << 0,                      1,           0,
               -vnn[i](1)*vnn[i](1),   2*vnn[i](1), g - 1,
               -g*vnn[i](1)*vnn[i](2), g*vnn[i](2), vnn[i](1);
      Ln[i] << vnn[i](1) + c, 0,         0,
               0,             vnn[i](1), 0,
               0,             0,         vnn[i](1) - c;
      omtn[i] << -vnn[i](1)*c, c,  g - 1,
                 -c*c,         0,  g - 1,
                  vnn[i](1)*c, -c, g - 1;
      omtmn[i] << 1./2./c/c,              -1./c/c,        1./2./c/c,
                  (vnn[i](1) + c)/2./c/c, -vnn[i](1)/c/c, (vnn[i](1) - c)/2./c/c,
                  1./2./(g - 1),          0,              1./2./(g - 1);
    }

    //check CFL value
    //if too big, make time step smaller
    tau = tau_init;
    bool chflag = false, fineflag = false;
    double CFL = 0;
    while(!fineflag) {
      for(int i = 0; i < N; i++) {
        CFL = tau / h * Ln[i].cwiseAbs().maxCoeff();
        if(CFL > CFL_max) {
          chflag = true;
          break;
        }
        fineflag = true;
      }
      if(chflag) tau /= 2;
    }

    //solve next layer
    t += tau;
    for(int i = 1; i < N - 1; i++) {
      wn[i] = wnn[i] - tau/(2*h)*(An[i]*(wnn[i + 1] - wnn[i - 1])) + tau/(2*h)*(omtmn[i]*(Ln[i].cwiseAbs())*omtn[i])*(wnn[i + 1] - 2*wnn[i] + wnn[i - 1]);

    //compute pressure and primitive variables vector
      pn[i] = (g - 1)*wn[i](2);
      vn[i](0) = wn[i](0);
      vn[i](1) = wn[i](1) / wn[i](0);
      vn[i](2) = wn[i](2) / wn[i](0);

    //rewrite previous layer
      wnn[i] = wn[i];
      vnn[i] = vn[i];
      pnn[i] = pn[i];
    }

    //boundary conditions
    wn[0] = wn[1];   wn[N - 1] = wn[N - 2];
    vn[0] = vn[1];   vn[N - 1] = vn[N - 2];
    pn[0] = pn[1];   pn[N - 1] = pn[N - 2];
    wnn[0] = wnn[1]; wnn[N - 1] = wnn[N - 2];
    vnn[0] = vnn[1]; vnn[N - 1] = vnn[N - 2];
    pnn[0] = pnn[1]; pnn[N - 1] = pnn[N - 2];

    //output if time is at .015 
    if(t > 0.015 - tau/2 && t < 0.015 + tau/2) {
      std::ofstream ofi("t15.dat");
      for(int i = 0; i < N; i++) {
        ofi << -L + i*h << " " << pn[i] << " " << vn[i](0) << " " << vn[i](1) << " " << vn[i](2) << std::endl;
      }
      std::cout << "written t=" << t << std::endl;
    }    

    //write for animation
    if(t > nt - tau/2 && t < nt + tau/2) {
      anim << "\"t = " << t << "\"" << std::endl;
      for(int i = 0; i < N; i++) {
        anim << -L + i*h << " " << pn[i] << " " << vn[i](0) << " " << vn[i](1) << " " << vn[i](2) << std::endl;
      } 
      anim << std::endl << std::endl;
      nt += 0.001;
      std::cout << "written t=" << t << std::endl;
    }
  }


  //write final result
  std::ofstream ofi("t2.dat");
  for(int i = 0; i < N; i++) {
    ofi << -L + i*h << " " << pn[i] << " " << vn[i](0) << " " << vn[i](1) << " " << vn[i](2) << std::endl;
  }
  std::cout << "written t=" << t << std::endl;

  //call plotting scripts
  (void)execl("/usr/bin/gnuplot", "/usr/bin/gnuplot", "plot.gp", (char *)0);
  (void)execl("/usr/bin/gnuplot", "/usr/bin/gnuplot", "anim.gp", (char *)0);

  return 0;
}

/*
i drew an owl for fun in the process, don't mind her
    ^______^
  /          \
 /  <.>  <.>  \
/ \    __    \ \
\ \    \/    \ /
 \\          \/
  \          /
     ^    ^

*/