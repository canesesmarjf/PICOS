#include "binary_tree.h"
#include <iostream>

// =======================================================================================
bt_TYP::bt_TYP()
{
  // All pointers must be made NULL:
  this->bt_params = NULL;
  root = NULL;
}

// =======================================================================================
bt_TYP::bt_TYP(bt_params_TYP * bt_params, vec * x_p, vec * a_p)
{
  // Update bt_params:
  bt_params->num_nodes = pow(2,sum(bt_params->max_depth));
  bt_params->dim_levels = cumsum(bt_params->max_depth);

  // Store pointer to tree parameters
  this->bt_params = bt_params;

  // Create root node:
  uint depth_root = 0;
  vec min = bt_params->min;
  vec max = bt_params->max;
  root = new node_TYP(min,max,depth_root,bt_params,x_p,a_p);
}

// =======================================================================================
void bt_TYP::insert_all(vector<vec *> data)
{
  // Insert points into nodes:
  this->root->insert_all(data);
}

// =======================================================================================
void bt_TYP::clear_all()
{
  if (NULL == root)
  {
    // cout << "Tree has NULL root node" << endl;
    return;
  }

  root->clear_node();
}

// =======================================================================================
int bt_TYP::count_leaf_nodes()
{
  int k = 0;
  return this->root->count_leaf_nodes(k);
}

// =======================================================================================
int bt_TYP::count_leaf_points()
{
  int k = 0;
  return this->root->count_leaf_points(k);
}

// =======================================================================================
void bt_TYP::delete_nodes()
{
  this->root->delete_nodes();
}

// =======================================================================================
int node_TYP::count_leaf_nodes(int k)
{
  if (this->subnode[0] != NULL)
  {
    k = this->subnode[0]->count_leaf_nodes(k);
  }

  if (this->subnode[1] != NULL)
  {
    k = this->subnode[1]->count_leaf_nodes(k);
  }

  // If both subnodes are NULL, the this must be a leaf_node:
  if (this->subnode[0] == NULL && this->subnode[1] == NULL)
  {
    return k + 1;
  }
  else
  {
    return k;
  }

}

// =======================================================================================
int node_TYP::count_leaf_points(int k)
{
  // If subnode[0] exists, drill into it
  if (this->subnode[0] != NULL)
  {
    k = this->subnode[0]->count_leaf_points(k);
  }

  // If subnode[1] exists, drill into it
  if (this->subnode[1] != NULL)
  {
    k = this->subnode[1]->count_leaf_points(k);
  }

  // If both subnodes are NULL, the this must be a leaf_node:
  if (this->subnode[0] == NULL && this->subnode[1] == NULL)
  {
    return k + this->ip_count;
  }
  else
  {
    return k;
  }
}

// =======================================================================================
void node_TYP::delete_nodes()
{
  if (this->subnode[0] != NULL)
  {
    this->subnode[0]->delete_nodes();
    delete this->subnode[0]; // Release memory pointed by subnode[0]
    this->subnode[0] = NULL; // Prevent dangling pointer
  }

  if (this->subnode[1] != NULL)
  {
    this->subnode[1]->delete_nodes();
    delete this->subnode[1]; // Release memory pointed by subnode[1]
    this->subnode[1] = NULL; // Prevent dangling pointer
  }
}

// =======================================================================================
node_TYP::node_TYP(vec min, vec max, uint depth, bt_params_TYP * bt_params, vec * x_p, vec * a_p)
{
  // Node attributes:
  this->center  = (min + max)/2;
  this->min     = min;
  this->max     = max;
  this->depth       = depth;
  this->bt_params = bt_params;

  // Allocate memory for subnodes:
  this->subnode.reserve(2);
  this->subnode[0] = NULL;
  this->subnode[1] = NULL;
  this->ip_count    = 0;

  // Pointers to data:
  this->x_p = x_p;
  this->a_p = a_p;
}

// insert method:
// =======================================================================================
void node_TYP::insert(uint i, vector<vec *> data, bool write_data)
{
    // Objective:
    // insert point into a subnode of current node. When maximum depth is reached, append point to node.

    // Determine dimension to operate on:
    int dim = sum(depth > bt_params->dim_levels);

    // Current data point:
    double p = as_scalar(data[dim]->at(i)); // make it a function that takes in i, data,dim

    // Get the particle weight:
    double w = (*a_p)[i];
    if (w < 1e-2)
      return;

    // Check if data is within node's boundaries:
    // ===============================================
    if (!IsPointInsideBoundary(p,dim))
    {
        // Warning message:
        cout << "point " << p <<" is outside domain" << endl;
        return;
    }

    // Check if maximum tree depth has been reached:
    // ===============================================
    // If yes, append point's index to leaf and RETURN up the stack
    if (HasNodeReachMaxDepth(dim))
    {
        // Append point:
        if (write_data)
        {
          this->ip.push_back(i);
          this->ip_count++;
        }

        // Return control to calling stack IF max_depth for all dimensions has been reached:
        uint total_dims = bt_params->dimensionality;
        if (depth >= bt_params->dim_levels(total_dims-1))
        {
          return;
        }
        else
        {
          // Evaluate data at next dimension:
          dim++;
          p = as_scalar(data[dim]->at(i));
        }
    }

    // Determine which subnode to insert point:
    // ========================================
    int node_index = WhichSubNodeDoesItBelongTo(p,dim);

    // Check if subnode needs to be created:
    // ====================================
    if ( !DoesSubNodeExist(node_index) )
    {
        CreateSubNode(node_index,dim);
    }

    // Insert point into subnode:
    // ==========================
    this->subnode[node_index]->insert(i,data,write_data);

} // node_TYP::insert

// insert_all method:
// ================================================================================================================
void node_TYP::insert_all(vector<vec *> data)
{
  for (int i = 0; i < data[0]->size(); i++)
  {
      // if ( i == 9500 - 1 )
      // {
      //   cout << i << endl;
      // }
      this->insert(i,data,true);
  }
}


// =================================================================================================================
bool node_TYP::IsPointInsideBoundary(double p, int dim)
{
    // Objective:
    // if p is inside the boundaries of the node, return true, otherwise false

    // Define boundaries of node:
    double x_left  = this->min[dim];
    double x_right = this->max[dim];

    // Create boolean result:
    bool flag = ((p >= x_left) && (p <= x_right));

    return flag;
}

// ================================================================================================================
bool node_TYP::HasNodeReachMaxDepth(int dim)
{
    int depth = this->depth;

    if (depth >= bt_params->dim_levels[dim]) // Greater only? equal only?
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// ================================================================================================================
int node_TYP::WhichSubNodeDoesItBelongTo(double p, int dim)
{
    //   +------------------+------------------+
    //   |  node_left = 1  |  node_right = 0  |

    // Center location of present node:
    double x_center = this->center[dim];

    // Number associated with subnode:
    int node_index;

    // Determine which location point belongs to:
    if ( p < x_center)
    {
        node_index = 1; // left
    }
    else
    {
        node_index = 0; // right
    }

    return node_index;
}

// find method:
// =================================================================================================================
node_TYP * bt_TYP::find(int i,vector<vec *> data,int search_dimensionality)
{
  return this->root->find(i,data,search_dimensionality);
}

// find method:
// =================================================================================================================
node_TYP * bt_TYP::find(double xq)
{
  return this->root->find(xq);
}

// find method:
// =================================================================================================================
node_TYP * node_TYP::find(double xq, int dim)
{
    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq,dim))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth(dim))
    {
        return this;
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq,dim);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(xq,dim);

}

// find method:
// =================================================================================================================
node_TYP * node_TYP::find(double xq)
{
    // 1D data input implies the following:
    int search_dimensionality = 1;

    // Calculate the relevant dimension based on the depth of the present node:
    int dim = sum(depth > bt_params->dim_levels);

    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq,dim))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth(dim))
    {
        return this;
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq,dim);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(xq);

}

// find method:
// =================================================================================================================
node_TYP * node_TYP::find(uint i, vector<vec *> data, int search_dimensionality)
{
    // Check the value of search _dimensionality:
    assert(search_dimensionality <= bt_params->dimensionality && "search_dimensionality needs to be less or equal to dimensionality of tree");

    assert(search_dimensionality > 0 && "search_dimensionality needs to be greater than zero");

    // Calculate the relevant dimension based on the depth of the present node:
    int dim = sum(depth > bt_params->dim_levels);

    // Extract the data based on the dimension being considered:
    double xq = as_scalar(data[dim]->at(i));

    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq,dim))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth(dim))
    {
      // Return control to calling stack IF maximum search dimension has been reached:
      if (depth >= bt_params->dim_levels(search_dimensionality-1))
      {
        return this;
      }
      else
      {
        // Evaluate data at next dimension:
        dim++;
        xq = as_scalar(data[dim]->at(i));
      }
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq,dim);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(i,data,search_dimensionality);
}

// ================================================================================================================
bool node_TYP::DoesSubNodeExist(int node_index)
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

// ================================================================================================================
void node_TYP::CreateSubNode(int node_index, int dim)
{
    // Attributes for new subnode:
    uint depth = this->depth + 1;
    vec min = this->min; //x_left;
    vec max = this->max; //x_right;
    ip.reserve(500);

    switch (node_index)
    {
    case 0: // Right subnode:
        {
            // Attributes for new subnode:
            min[dim] = this->center[dim];
            max[dim] = this->max[dim];

            // Exit:
            break;
        }
    case 1: // Left subnode:
        {
            // Attributes for new subnode:
            min[dim] = this->min[dim];
            max[dim] = this->center[dim];

            // Exit:
            break;
        }
    }

    // Create new subnode:
    this->subnode[node_index] = new node_TYP(min, max, depth, bt_params,x_p,a_p);
}

// ================================================================================================================
void node_TYP::clear_node()
{
  this->ip_count = 0;
  this->ip.clear();

  if (NULL != this->subnode[0])
  {
    this->subnode[0]->clear_node();
  }

  if (NULL != this->subnode[1])
  {
    this->subnode[1]->clear_node();
  }
}
