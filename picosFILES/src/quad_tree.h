#ifndef QUADTREE_H
#define QUADTREE_H

#include <vector>
#include <armadillo>
using namespace std;
using namespace arma;

// =====================================================================================
class qt_params_TYP
{
public:
  vec min; // Minimum lengths that bound domain covered by tree
  vec max; // Maximum lengths that bound domain covered by tree
  uint min_depth;
  uint max_depth; // Maximum depth to be reached of min_count is not reached first
  uint min_count; // Minimum number of particles in cell allowed.

  // Tree structure information:
  int subnodes_created;

  // Constructor:
  qt_params_TYP()
  {
    subnodes_created = 0;
  }
};

// =====================================================================================
class q_node_TYP
{
public:
  // Node "data" attributes:
  int ip_count;            // Number of points indexed in the current node
  vector<uint> ip;   // Indices of data appended to this node
  mat *v;           // Pointer to data to be indexed

  // Node "natural" attributes:
  uint depth;
  vec min;
  vec max;
  vec center;

  // Constructor:
  q_node_TYP(){};
  q_node_TYP(vec min, vec max, uint depth, qt_params_TYP * qt_params, vector<uint> ip,mat * v);

  // Methods:
  void populate_subnodes();
  void organze_ip_into_proposed_subnodes();
  void clear_node();
  void delete_nodes();
  int count_leaf_points(int k);
  void get_subnode_bounds(int node_index, vec * min_l, vec * max_l);
  void get_leaf_nodes(vector<q_node_TYP *> * leafs);
  void create_subnode(int n, vector<uint> ip);
  void update_subnode(int n,vector<uint> ip);
  void clear_ip_subnode();

private:
  // Variables:
  qt_params_TYP * qt_params; // Pointer to tree parameters
  vector<q_node_TYP *> subnode;
  vector<vector<uint>> ip_subnode;
  bool is_leaf;

  // subnodes within this node:
  //   +------------------+------------------+
  //   |    subnode[1]    |    subnode[0]    |
  //   +------------------+------------------+
  //   |    subnode[2]    |    subnode[3]    |
  //   +------------------+------------------+

  // Methods:
  bool IsPointInsideBoundary(vec r);
  int WhichSubNodeDoesItBelongTo(vec r);
  bool DoesSubNodeExist(int subNode);
  void CreateSubNode(int subNode);
  bool is_node_leaf(int method);
  int apply_conditionals_ip_subnode();

};

// =====================================================================================
class qt_TYP
{
public:
  // Constructor:
  qt_TYP();
  qt_TYP(qt_params_TYP * qt_params, vector<uint> ip, mat * v);

  // Variables:
  q_node_TYP * root; // Root node of tree
  qt_params_TYP * qt_params;  // Pointer to tree attributes

  // Methods:
  void populate_tree();
  void clear_tree();
  void delete_tree();
  vector<q_node_TYP *> get_leaf_nodes();
  int count_leaf_points();

  // QUESTION:
  // Everytime we add data to the quadtree, do we need to delete all leaf nodes? This means releasing memory and deleting dangling pointers. We need to consider this carefully as we will be using the quadtrees potetially millions of times.

private:
  // Variables:
  // node_TYP * root; // Root node of tree

  // Methods:
  // void assemble_node_list();
  // void save_data(int ii, string prefix);
};

#endif
