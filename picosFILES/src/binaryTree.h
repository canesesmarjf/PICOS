#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <armadillo>
#include <random>
#include "types.h"
using namespace std;

// This header file contains three main classes:
// - tree_params_TYP
// - binaryTree_TYP
// - node_TYP

// Is important to node that a tree is composed of nodes; thus a tree as a whole will have parameters such as:
// - Maximum depth
// - Left and right boundaries
// - Number of nodes in the final layer
// Moreover, each node will also have its own parameters:
// - current depth
// - Left and right boundaries
// - Data stored in node

// =====================================================================================
class tree_params_TYP
{
public:
  // Interface variables:
  double x_left;  // Left boundary of tree
  double x_right; // Right boundary of tree
  int depth_max; // Maximum depth of tree
  int num_elem;   // Total number of elements to hold in tree

  // Derived variables:
  double length; // Length of entire domain
  double dx;
  int num_nodes; // Number of nodes based on 2^depth_max
  int mean_elems_per_node; // number of elements per node based on a uniform distribution
  arma::vec node_centers; // Vector containing the center locations of the nodes

  // Constructor:
  tree_params_TYP();
  tree_params_TYP(double x_left, double x_right, int depth_max, int num_elem);
};

// =====================================================================================
class node_TYP
{
public:
  // Data variables:
  int x_count;            // Number of points indexed in the current node
  std::vector<uint> ix;   // Indices of data appended to this node

  // Node parameters:
  double x_center;        // Center position of the node
  double x_left;          // Left boundary of node
  double x_right;         // Right boundary of node
  int depth;

  // Constructor:
  node_TYP();
  node_TYP(double x_left, double x_right, int depth, tree_params_TYP * tree_params);

  // Methods:
  void insert(uint i, arma::vec * r, bool write_data); // Insert the ith element of vector *r
  void insert_all(arma::vec * r); // Insert all elements of vector *r
  node_TYP * find(double xq); // Find and return pointer of node corresponding to position xq
  node_TYP * get_subnode(int index); // ?

private:
  // Variables:
  tree_params_TYP * tree_params; // Pointer to tree parameters

  // Subnodes within this node:
  // subnode[0] : right_node
  // subnode[1] : left_node
  //   +------------------+------------------+
  //   |  left_node = 1   |  right_node = 0  |
  std::vector<node_TYP *> subnode;

  // Methods:
  bool IsPointInsideBoundary(double p);
  bool HasNodeReachMaxDepth();
  int WhichSubNodeDoesItBelongTo(double p);
  bool DoesSubNodeExist(int subNode);
  void CreateSubNode(int subNode);
};

// =====================================================================================
class binaryTree_TYP
{
public:
  // Constructor:
  binaryTree_TYP(){};
  binaryTree_TYP(double x_left, double x_right, int depth_max, int num_elem);

  // Tree parameters:
  tree_params_TYP * tree_params; // Pointer to tree attributes

  // Tree node list:
  std::vector<node_TYP *> node_list; // List of pointers to nodes at maxmimum depth

  // Data variables:
  arma::ivec count_profile; // Profile of x_count
  arma::ivec delta_profile; // Record nodes that have +ve (surplus) or -ve (deficit) number of particles relative to mean_elems_per_node
  uint total_deficit_elems; // Total number of elements in deficit. Represents the number of elemnts from the surplues region to repurpose.
  std::vector<uint> ix_repurpose; // List of indices to repurpose from surplus nodes

  // Methods:
  void insert_all(arma::vec * r);
  int get_num_nodes();
  arma::vec get_node_centers();
  int get_max_depth();
  void clear_all();
  void print_info(int ii);
  void save_data_all(string prefix);
  void calculate_delta_profile();
  void gather_all_surplus_indices();
  void renormalize_surplus_nodes(ions_TYP * IONS);
  void renormalize_deficit_nodes(ions_TYP * IONS);

private:
  // Variables:
  node_TYP * root; // Root node of tree

  // Methods:
  void assemble_empty_tree();
  void assemble_node_list();
  void assemble_count_profile();
  void save_data(int ii, string prefix);
};

#endif
