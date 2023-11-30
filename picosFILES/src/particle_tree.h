#ifndef PARTICLETREE_H
#define PARTICLETREE_H

#include <vector>
#include <armadillo>

#include "binary_tree.h"
#include "quad_tree.h"
#include "vranic_downsampling.h"

using namespace std;
using namespace arma;

class particle_tree_TYP
{
// private:


public:
  // Variables:
  bt_TYP bt; // Binary tree object that hold x-space ROOT nodes:
  vector<qt_TYP> qt; // vector of quadtrees to hold v-space ROOT nodes
  bt_params_TYP * bt_params;  // Pointer to tree parameters
  qt_params_TYP * qt_params; // Pointer to tree parameters
  vec xq; // x query grid

  // Variables:
  vector<node_TYP *> leaf_x; // Vector of x-space LEAF nodes
  vector<vector<q_node_TYP *>> leaf_v; // Vector of v-space LEAF nodes
  ivec ip_count; // Hold the particle count for each leaf_x node

  // Constructor:
  particle_tree_TYP(){};
  particle_tree_TYP(bt_params_TYP * bt_params, qt_params_TYP * qt_params, vec * x_p, mat * v_p, vec * a_p);

  // Methods:
  void populate_tree(string calculation_type);
  void resample_distribution();
  void clear_all_contents();
  void assess_conservation(string output_dir, string suffix);
  void save_leaf_v_structure(string output_dir);

private:
  void downsample_surplus_nodes(vector<uint> * ip_free);
  void upsample_deficit_nodes(vector<uint> * ip_free);
  void populate_binary_tree();
  void populate_quad_tree();
  void get_mean_ip_count();
  void create_x_query_grid();
  void calculate_leaf_v();
  void calculate_leaf_x();
  void assemble_qt_vector();

  int mean_ip_count;
  vec * x_p; // Pointer to particle x data
  mat * v_p; // Pointer to particle v data
  mat * a_p; // Pointer to particle weight data
};

#endif
