#include <iostream>
#include <vector>
#include <array>
#include <fstream>

int main() {
  //constants
  double a = 5.0, b = 1.0;
  double tau = 0.001;
  double h = a * tau / 0.8;

  //generating grid with initial conditions
  const int size_x = (int)(2 / h);
  const int size_v = (int)(10 / tau);
  std::vector<std::vector<double>> grid(size_v);
  for(auto it = grid.begin(); it != grid.end(); it++) {
    (*it).reserve(size_x);
  }

  for(int i = 0; i < size_x; i++) {
    grid[0][i] = 0;
  }
  for(int i = 0; i < size_v; i++) {
    grid[i][0] = 1;
  }
  grid[0][0] = 0;

  //computation
  for(int m = 0; m < size_v - 1; m++) {
    for(int n = 1; n < size_x; n++) {
      grid[m+1][n] = b*tau + (a*tau/h)*grid[m][n-1] - (a*tau/h - 1)*grid[m][n];
    }
  }

  /*
  for(int i = 0; i < size_v; i++) {
    for(int j = 0; j < size_x; j++) {
      std::cout << grid[i][j] << " ";
    }
    std::cout << std::endl;
  }
  */

  //output
  std::vector<double> tf01, tf1, tf10;
  tf01.reserve(size_x);
  tf1.reserve(size_x);
  tf10.reserve(size_x);

  for(int i = 0; i < size_x; i++) {
    tf01[i] = grid[100 - 1][i];
    tf1[i] = grid[200 - 1][i];
    tf10[i] = grid[300 - 1][i];
  }

  std::ofstream ofi("output.dat");
  for(int i = 0; i < size_x; i++) {
    ofi << i*h << " " << tf01[i] << " " << tf1[i] << " " << tf10[i] << std::endl;
  }
  
  return 0;
}