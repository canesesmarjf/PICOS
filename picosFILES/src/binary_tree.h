#ifndef BINARYTREE_H
#define BINARYTREE_H

#include <vector>
#include <armadillo>
#include <cassert>
using namespace std;
using namespace arma;

// This header file contains three main classes:
// - bt_params_TYP
// - bt_TYP
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
class bt_params_TYP
{
public:
  uint dimensionality;
  vec min;
  vec max;
  uvec max_depth;
  uvec dim_levels; // Determines when the tree needs to change to the next dimension
  uint num_nodes; // Number of nodes based on 2^depth_max

  // Derived variables:
  // double length; // Length of entire domain
  // double dx;
  // int num_elem; // Total number of elements to hold in tree

  // int elem_per_node; // number of elements per node based on a uniform distribution
  // vec node_centers; // Vector containing the center locations of the nodes
};

// =====================================================================================
class node_TYP
{
public:
  // Data variables:
  int ip_count;            // Number of points indexed in the current node
  vector<uint> ip;   // Indices of data appended to this node

  // Node parameters:
  uint depth;
  vec min;
  vec max;
  vec center;

  // Constructor:
  node_TYP(){};
  node_TYP(vec min, vec max, uint depth, bt_params_TYP * bt_params);

  // Methods:
  void insert(uint i, vector<vec *> data, bool write_data); // Insert the ith element of data
  void insert_all(vector<vec *> data); // Insert all elements of the data
  node_TYP * find(uint i, vector<vec *> data, int search_dimensionality);
  void clear_node();
  node_TYP * find(double xq); // Find and return pointer of node corresponding to position xq
  node_TYP * find(double xq, int dim); // Find and return pointer of node corresponding to position xq searched along dimension "dim"
  int count_leaf_nodes(int k);
  int count_leaf_points(int k);
  void delete_nodes();

private:
  // Variables:
  bt_params_TYP * bt_params; // Pointer to tree parameters
  vector<node_TYP *> subnode;

  // Subnodes within this node:
  // subnode[0] : right_node
  // subnode[1] : left_node
  //   +------------------+------------------+
  //   |  left_node = 1   |  right_node = 0  |
  //   +------------------+------------------+

  // Methods:
  bool IsPointInsideBoundary(double p, int dim);
  bool HasNodeReachMaxDepth(int dim);
  int WhichSubNodeDoesItBelongTo(double p, int dim);
  bool DoesSubNodeExist(int subNode);
  void CreateSubNode(int subNode, int dim);
};

// =====================================================================================
class bt_TYP
{
public:
  // Constructor:
  bt_TYP();
  bt_TYP(bt_params_TYP * bt_params);

  // Variables:
  node_TYP * root; // Root node of tree
  bt_params_TYP * bt_params;  // Pointer to tree attributes

  // Methods:
  void insert_all(vector<vec *> data);
  node_TYP * find(int i,vector<vec *> data,int search_dimensionality);
  node_TYP * find(double xq);
  void clear_all();
  int count_leaf_nodes();
  int count_leaf_points();
  void delete_nodes();

private:

};

#endif
