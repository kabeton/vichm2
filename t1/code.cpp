#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>

int main(int argc, char *argv[], char *envp[]) {
  assert(argc == 6);
  //initial conditions
  char *scheme = argv[1];
  double x0 = atof(argv[2]);
  double L = atof(argv[3]);
  double T = atof(argv[4]);
  double CFL = atof(argv[5]);

  double h = 0.5;
  double tau = h * CFL;
  
  //grid generation
  int M = (int)(L/h);
  int N = (int)(T/tau);
  std::cout << N << std::endl;
  double *grid = new double[N*M];
  
  //u(x, 0) = sin(4pi*x/L)
  for(int i = 0; i < M; i++) {
    grid[i] = sin(4*3.1415926*(x0 + i*h)/L);
  }

  //right angle scheme
  if(!strcmp(scheme, "ra")) {
    for(int n = 0; n < N - 1; n++) {
      for(int m = 1; m < M; m++) {
        grid[(n + 1)*M + m] = CFL*(grid[n*M + (m - 1)] - grid[n*M + m]) + grid[n*M + m];
      }
      //periodic boundary condition
      grid[(n + 1)*M] = grid[(n + 2)*M - 1];
    }
  }
  
  //lax-wendroff scheme
  if(!strcmp(scheme, "lv")) {
  }

  //writing answer to file
  char *oname;
  asprintf(&oname, "%s_%.1f_%.1f_%.1f.dat", scheme, L, T, CFL);
  std::ofstream ofi(oname);
  
  for(int i = 0; i < M; i++) {
    ofi << i*h << " ";
    for(int j = 0; j < N; j += (int)T/36/tau + 1) {
      ofi << grid[j*M + i] << " ";
    }
    ofi << std::endl;
  }
  
  delete [] grid;
  free(oname);
  return 0;
}
