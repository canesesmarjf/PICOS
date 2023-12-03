#include "quad_tree.h"
#include <iostream>

// =======================================================================================
qt_TYP::qt_TYP()
{
  // All pointers must be made NULL:
  this->qt_params = NULL;
  root = NULL;
}

// =======================================================================================
qt_TYP::qt_TYP(qt_params_TYP * qt_params, vector<uint> ip, mat * v)
{
  // Store pointer to tree parameters
  this->qt_params = qt_params;

  // Create root node:
  uint depth_root = 0;
  vec min = qt_params->min;
  vec max = qt_params->max;
  root = new q_node_TYP(min,max,depth_root,qt_params,ip,v);
}

// =======================================================================================
void qt_TYP::populate_tree()
{
  this->root->populate_subnodes();
}

// =======================================================================================
void qt_TYP::clear_tree()
{
  if (NULL == root)
  {
    //cout << "Tree has NULL root node" << endl;
    return;
  }

  root->clear_node();
}

// =======================================================================================
void qt_TYP::delete_tree()
{
  this->root->delete_nodes();
}

// =======================================================================================
int qt_TYP::count_leaf_points()
{
  int k = 0;
  if (this->root != NULL)
    return this->root->count_leaf_points(k);
  else
    return k;
}

// =======================================================================================
vector<q_node_TYP *> qt_TYP::get_leaf_nodes()
{
  vector<q_node_TYP *> leafs;
  root->get_leaf_nodes(&leafs);
  return leafs;
}

// =======================================================================================
void q_node_TYP::get_leaf_nodes(vector<q_node_TYP *> * leafs)
{
  // Check if present node is a leaf:
  int method = 2;
  if (is_node_leaf(method) == true)
  {
    // Add data to output variable:
    leafs->push_back(this);

    // If a leaf has been found, then no need to drill in deeper:
    return;
  }

  // Traverse the tree:
  for(int nn = 0; nn < 4; nn++)
  {
    // If subnode[nn] exists, drill into it
    if (this->subnode[nn] != NULL)
    {
      this->subnode[nn]->get_leaf_nodes(leafs);
    }
  }
}

// =======================================================================================
void q_node_TYP::delete_nodes()
{
  // This method loops over ALL subnodes of the present node and identifies which subnodes have pointers to other nodes.
  // If a subnode has a valid pointer, it then drills into it (recursion).
  // The recursion finishes whenever a given node has NULL for ALL its subnodes.
  // This condition is the definition of a leaf node.
  // When that condition is net, delete_nodes does nothing and just returns to the calling stack and deletes the subnode[nn] from which we just returned.
  // In other words, it deletes the leaf subnode and them moves to the next nn

  for (int nn = 0; nn < 4; nn++)
  {
    if (this->subnode[nn] != NULL) // if subnode[nn] exists
    {
      this->subnode[nn]->delete_nodes();
      // This section can only be reached if all elements of this->subnode == NULL

      // Release memory pointed by subnode[nn]:
      delete this->subnode[nn];

      // Prevent dangling pointer:
      this->subnode[nn] = NULL;

    } // if
  } // for
}

// =======================================================================================
// Overloaded constructor:
q_node_TYP::q_node_TYP(vec min, vec max, uint depth, qt_params_TYP * qt_params,vector<uint> ip, mat * v)
{
  // Node "natural" attributes:
  this->min     = min;
  this->max     = max;
  this->depth   = depth;
  this->qt_params = qt_params;

  // Node "data" attributes:
  this->ip = ip;
  this->v  = v;

  // Derived node "natural" attributes:
  this->center  = (min + max)/2;

  // Derived "data" attributes:
  this->ip_count = ip.size();
  this->is_leaf = false;

  // Allocate memory for subnodes:
  this->subnode.resize(4, NULL);
  this->ip_subnode.resize(4);
}

// =======================================================================================
int q_node_TYP::apply_conditionals_ip_subnode()
{
  // Tree parameters:
  int min_depth = qt_params->min_depth;
  int max_depth = qt_params->max_depth;
  int min_count = qt_params->min_count;

  // Node depths
  int parent_node_depth = this->depth;
  int subnode_depth = parent_node_depth + 1;

  // Reject new subnodes by default:
  int accept_new_subnodes = 0;

  // Determine in which of 3 state the subnode will be in:
  int state;
  if (subnode_depth <= min_depth)
  {
    state = 1;
  }
  else if (subnode_depth > max_depth)
  {
    state = 3;
  }
  else
  {
    state = 2;
  }

  // Compute result based on state:
  switch(state)
  {
    case 1: // subnode_depth <= min_depth
      accept_new_subnodes = 1;
      break;
    case 3: // subnode_depth > max_depth
      accept_new_subnodes = 0;
      break;
    case 2: // subnode_depth is between [min_depth, max_depth].
      // Default:
      accept_new_subnodes = 0;

      // Accept new subnodes if at least ONE of them has enough particles (> min_count)
      for (int n = 0; n < 4; n++)
      {
        if (ip_subnode[n].size() > min_count)
        {
          accept_new_subnodes = 1;
          break; // Break out of loop
        }
      }
      break;
  }

  return accept_new_subnodes;
}

// =======================================================================================
void q_node_TYP::create_subnode(int n, vector<uint> ip)
{
  // Bounds of new subnode:
  vec min_local(2);
  vec max_local(2);
  get_subnode_bounds(n,&min_local,&max_local);

  // Create subnode:
  uint depth = this->depth + 1;
  subnode[n] = new q_node_TYP(min_local,max_local,depth,qt_params,ip,v);
}

// =======================================================================================
void q_node_TYP::update_subnode(int n,vector<uint> ip)
{
  // Update node "data" attributes of the subnode which alread exists:
  subnode[n]->ip = ip;
  subnode[n]->ip_count = ip.size();
  subnode[n]->is_leaf = false;
}

// =======================================================================================
void q_node_TYP::clear_ip_subnode()
{
  for (int n = 0; n < 4; n++)
  {
    ip_subnode[n].clear();
  }
}

// =======================================================================================
void q_node_TYP::populate_subnodes()
{
  // When each q_node is created (constructor) it is automatically populated with an ip
  // vector. This corresponds to the data associated with this node.
  // If we want to proceed (drill deeper), we need to populate the subnodes.

  // Organize ip data into ip_subnode to determine if new subnodes need to be created:
  // -------------------------------------------------------------------------------------
  organze_ip_into_proposed_subnodes();

  // Run tests to determine if new subnodes are to be accepted:
  // -------------------------------------------------------------------------------------
  // Do we accept new proposed subnodes?
  int accept_new_subnodes = apply_conditionals_ip_subnode();

  // If not accepted, declare parent node a leaf node
  // -------------------------------------------------------------------------------------
  if (accept_new_subnodes == 0)
  {
    // Declare parent node as a leaf node:
    this->is_leaf = true;
    // cout << "leaf_node" << endl;
  }

  // If accepted, populate subnodes
  // If they dont exist yet, create them:
  // -------------------------------------------------------------------------------------
  if (accept_new_subnodes == 1)
  {
    // Declare parent node as NOT a leaf node:
    this->is_leaf = false;

    // Loop over proposed subnodes.
    // if ip_subnode[n] is not empty, do the following
    // If they already exist, just update subnode[n].ip variable
    // If they dont exist, create them using the new operator

    // Loop over proposed subnodes:
    for (int n = 0; n < 4; n++)
    {
      // Check that proposed subnode has data:
      if (ip_subnode[n].size() > 0)
      {
        // Create new subnode if it doesnt exist:
        if (subnode[n] == NULL)
        {
          create_subnode(n,ip_subnode[n]);
        }
        else // subnode[n] already exists:
        {
          update_subnode(n,ip_subnode[n]);
        }

        // Move in deeper:
        subnode[n]->populate_subnodes();
      }
    } // end, for loop

    // Clear ip on parent node since they have now being distributed amongsnt new subnodes:
    this->ip.clear();
    this->ip_count = 0;
  }

  // Clear ip_subnode. Data is either not needed OR has been inserted into new subnodes:
  // -------------------------------------------------------------------------------------
  clear_ip_subnode();

  return;
}

// =======================================================================================
void q_node_TYP::get_subnode_bounds(int node_index, vec * min_local, vec * max_local)
{
  switch (node_index)
  {
  case 0: // subnode 0:
      {
        // Attributes for new subnode:
        *min_local = this->center;
        *max_local = this->max;

        // Exit:
        break;
      }
  case 1: // subnode 1:
      {
        // Attributes for new subnode:
        (*min_local)[0] = this->min[0];
        (*min_local)[1] = this->center[1];
        (*max_local)[0] = this->center[0];
        (*max_local)[1] = this->max[1];

        // Exit:
        break;
      }
  case 2: // subnode 2:
      {
        // Attributes for new subnode:
        *min_local = this->min;
        *max_local = this->center;

        // Exit:
        break;
      }
  case 3: // subnode 3:
      {
        // Attributes for new subnode:
        (*min_local)[0] = this->center[0];
        (*min_local)[1] = this->min[1];
        (*max_local)[0] = this->max[0];
        (*max_local)[1] = this->center[1];

        // Exit:
        break;
      }
  }
}

// =======================================================================================
void q_node_TYP::organze_ip_into_proposed_subnodes()
{
  for (int ii = 0; ii < this->ip_count; ii++)
  {
    // Get global index:
    uint jj = this->ip[ii];

    // Current data point:
    double y = (*v)(jj,0);
    double z = (*v)(jj,1);
    vec r = {y,z};

    // Check if data is within node's boundaries:
    // ========================================
    if (!IsPointInsideBoundary(r))
    {
      cout << "point (" << r[0] << "," << r[1] << ") is out of bounds" << endl;
      return;
    }

    // Determine which subnode to insert point:
    // ========================================
    int node_index = WhichSubNodeDoesItBelongTo(r);

    // Store index jj in corresponding ip_subnode:
    // ========================================
    ip_subnode[node_index].push_back(jj);
  }
}

// =======================================================================================
bool q_node_TYP::is_node_leaf(int method)
{
  // This method determines if current node is a leaf node.
  // There are at least two ways to identify leaf_nodes.
  // 1- if ALL subnodes == NULL.
  // 2- if is_leaf == true where the flag is set during the formation of the quad tree in "populate_subnodes()".
  // The second method is favoured because the quad tree heap memory infrastructure can be reutilized many times; whereas the NULL method relies on realeasing the quad-tree memory every time it is to be refreshed because the data *v has changed.

  bool leaf_flag = true;

  switch (method)
  {
    case 1: // NULL method
    {
      // Check all subnodes:
      // If at LEAST ONE subode exist, current node is NOT a leaf:
      for (int nn = 0; nn < 4; nn++)
      {
        if (this->subnode[nn] != NULL)
        {
          leaf_flag = false;
          break;
        }
      }
      break;
    }
    case 2: // Read is_leaf variable
    {
      leaf_flag = this->is_leaf;
      break;
    }
    default:
    {
      // Error: Invalid method
      cerr << "Invalid method. Choose either method 1 or method 2." << endl;
      break;
    }
  }

  return leaf_flag;
}

// =======================================================================================
int q_node_TYP::count_leaf_points(int k)
{
  // This method counts all the points present in a quad tree
  // First, it deterines if present node is a leaf.
  // If so, it accumulates ip_count
  // If not so, it proceeds to traverse tree using recursion
  // Every time a new node is accessed, leaf status is checked.

  // If current is node is a leaf node, accumulate ip_count and return to calling stack:
  int method = 2;
  if (is_node_leaf(method) == true)
  {
    return k + this->ip_count;
  }

  // Traverse the tree:
  for(int nn = 0; nn < 4; nn++)
  {
    // If subnode[nn] exists, drill into it
    if (this->subnode[nn] != NULL)
    {
      k = this->subnode[nn]->count_leaf_points(k);
    }
  }

  // Return to calling stack:
  return k;
}

// =======================================================================================
bool q_node_TYP::IsPointInsideBoundary(vec r)
{
    // Objective:
    // if r is inside the boundaries of the node, return true, otherwise false

    // coordinates of point:
    double y = r[0];
    double z = r[1];

    // Define boundaries of node:
    double y_min = this->min[0];
    double y_max = this->max[0];
    double z_min = this->min[1];
    double z_max = this->max[1];

    // Create boolean result:
    bool flag_y = (y >= y_min) && (y <= y_max);
    bool flag_z = (z >= z_min) && (z <= z_max);

    return flag_y && flag_z;
}

// =======================================================================================
int q_node_TYP::WhichSubNodeDoesItBelongTo(vec r)
{
    //   |-----------------+------------------+
    //   |  node_left = 1  |  node_right = 0  |
    //   +-----------------o------------------+
    //   |  node_left = 2  |  node_right = 3  |
    //   +-----------------+------------------+

    // Origin of parent node:
    vec r0 = center;

    //  Vector pointing in the direction of point r relative to center of node:
    vec d = r - r0;

    // Number associated with subnode:
    int node_index;

    // Select node index based on vector d
    if (d[0] >= 0) // Quadrants 0 or 3
    {
      if (d[1] >=0) // Quadrant 0
      {
        node_index = 0;
      }
      else // Quadrant 3
      {
        node_index = 3;
      }
    }
    else // Quadrants 1 or 2
    {
      if (d[1] >=0)
      {
        node_index = 1;
      }
      else
      {
        node_index = 2;
      }
    }

    return node_index;
}

// =======================================================================================
bool q_node_TYP::DoesSubNodeExist(int node_index)
{
    if (NULL == subnode[node_index])
    {
        // Does not exist:
        return 0;
    }
    else
    {
        // It already exists:
        return 1;
    }
}

// =======================================================================================
void q_node_TYP::clear_node()
{
  // This method just removes the particle count and indices stored in this node. It does not remove the depth, min, max, center information as these are the node parameters and will not change if we use another dataset in this->v;

  this->ip_count = 0;
  this->ip.clear();
  this->is_leaf = false;

  for (int nn = 0; nn < 4; nn ++)
  {
    if (NULL != this->subnode[nn])
    {
      this->subnode[nn]->clear_node();
    }
  }
}
