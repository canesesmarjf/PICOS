#ifndef VRANIC_H
#define VRANIC_H

#include <vector>
#include <armadillo>
#include <cassert>
using namespace std;

class merge_cell_TYP
{
public:
arma::vec wi;
arma::vec xi;
arma::vec yi;
arma::vec zi;
int n_elem;

// Constructor:
merge_cell_TYP(int I)
{
  n_elem = I;
  wi.zeros(I);
  xi.zeros(I);
  yi.zeros(I);
  zi.zeros(I);
}

// Methods:

};

class vranic_TYP
{
public:

  // Constructor:
  vranic_TYP(){};

  // Methods:
  void down_sample_node_3D(merge_cell_TYP * set_M, merge_cell_TYP * set_N);
  void down_sample_node_2D(merge_cell_TYP * set_M, merge_cell_TYP * set_N);
  void print_stats(merge_cell_TYP * set);
  double get_sigma(merge_cell_TYP * set);
};

#endif
