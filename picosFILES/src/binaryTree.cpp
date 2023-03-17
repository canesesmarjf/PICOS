#include "binaryTree.h"
#include <iostream>

using namespace std;

// Constructors:
// ================================================================================================================
tree_params_TYP::tree_params_TYP(double x_left, double x_right, int depth_max, int num_elem)
{
  this->x_left = x_left;
  this->x_right = x_right;
  this->depth_max = depth_max;
  this->num_elem = num_elem;

  this->num_nodes = pow(2,depth_max);
  this->length = x_right - x_left;
  this->dx = (x_right - x_left)/this->num_nodes;
  this->elem_per_node = round(num_elem/this->num_nodes);
  this->node_centers = arma::linspace(x_left,(x_right-dx),this->num_nodes) + dx/2;
}

// ================================================================================================================
binaryTree_TYP::binaryTree_TYP(double x_left, double x_right, int depth_max, int num_elem)
{
  // Interface variables:
  tree_params = new tree_params_TYP(x_left, x_right, depth_max, num_elem);

  // Allocate memory to node list, which is a vector of pointers:
  node_list.resize(tree_params->num_nodes);

  // Create root node:
  int depth_root = 0;
  root = new node_TYP(tree_params->x_left,tree_params->x_right,depth_root,tree_params);

  // Traverse entire tree:
  // This is done by inserting node_centers so as to create all the node
  // at the final depth while not writting any data to the nodes.
  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    bool write_data = false;
    root->insert(nn,&tree_params->node_centers,write_data);
  }

  // Assemble list of nodes into a vector of pointers:
  assemble_node_list();
}

// ================================================================================================================
void binaryTree_TYP::insert_all(arma::vec * r)
{
  // Insert points into nodes:
  this->root->insert_all(r);
}

// ================================================================================================================
int binaryTree_TYP::get_num_nodes()
{
  return tree_params->num_nodes;
}

// ================================================================================================================
arma::vec binaryTree_TYP::get_node_centers()
{
  return tree_params->node_centers;
}

// ================================================================================================================
int binaryTree_TYP::get_max_depth()
{
  return tree_params->depth_max;
}

// ================================================================================================================
void binaryTree_TYP::assemble_node_list()
{
  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    // Get the location of the nth node center:
    double xq = tree_params->node_centers.at(nn);

    // Copy the pointer of the node containing xq into the node list:
    node_list.at(nn) = root->find(xq);
  }
}

// ================================================================================================================
void binaryTree_TYP::clear_all()
{
  for (int nn = 0; nn < node_list.size(); nn++)
  {
    node_list[nn]->x_count = 0;
    node_list[nn]->ix.clear();
  }
}

// ================================================================================================================
void binaryTree_TYP::print_info(int ii)
{
  cout << " " << endl;
  cout << "tree.node_list[ii]->x_left        = " << node_list[ii]->x_left << endl;
  cout << "tree.node_list[ii]->x_center      = " << node_list[ii]->x_center << endl;
  cout << "tree.node_list[ii]->x_right       = " << node_list[ii]->x_right << endl;
  cout << "tree.node_list[ii]->depth         = " << node_list[ii]->depth << endl;
  cout << "tree.node_list[ii]->ix.capacity() = " << node_list[ii]->ix.capacity() << endl;
  cout << "tree.node_list[ii]->x_count       = " << node_list[ii]->x_count << endl;
  cout << "tree.node_list[ii]->ix.size()     = " << node_list[ii]->ix.size() << endl;
  cout << " " << endl;
}

// ================================================================================================================
void binaryTree_TYP::save_data(int ii, string prefix)
{
  // Put data into armadillo containers:
  arma::uvec ix = arma::conv_to<arma::uvec>::from(node_list[ii]->ix);

  // Save data to file:
  string fileName = "output_files/" + prefix + to_string(ii) + ".txt";
  cout << "Saving data with fileName = " << fileName << endl;
  ix.save(fileName,arma::raw_ascii);
}

// ================================================================================================================
void binaryTree_TYP::save_data_all(string prefix)
{
  for (int ii = 0; ii < get_num_nodes(); ii++)
  {
      if (node_list[ii]->ix.empty() == false)
      {
          // Save data for node ii:
          save_data(ii, prefix);
      }
      else
      {
          cout << "no points in node i =  " << ii << endl;
      }
  } // for
}

// ================================================================================================================
node_TYP::node_TYP()
{
  cout << "default constructor" << endl;
}

// ================================================================================================================
node_TYP::node_TYP(double x_left, double x_right, int depth, tree_params_TYP * tree_params)
{
  // Node attributes:
  this->x_center    = (x_left + x_right)/2;
  this->x_left      = x_left;
  this->x_right     = x_right;
  this->depth       = depth;
  this->tree_params = tree_params;

  // Allocate memory for subnodes:
  this->subnode.reserve(2);
  this->subnode[0] = NULL;
  this->subnode[1] = NULL;
  this->x_count    = 0;

  if (depth == tree_params->depth_max)
  {
    // When depth_max == depth, we have reached the final layer of nodes. At this point. we can reserve
    // memory to this node as it will be a data holding node.

    // Reserve memory:
    this->ix.reserve(2*tree_params->elem_per_node);
  }
}

// insert method:
// ================================================================================================================
void node_TYP::insert(uint i, arma::vec * r, bool write_data)
{
    // Objective:
    // insert point into a subnode of current node. When maximum depth is reached, append point to node.

    // Current data point:
    double p = arma::as_scalar(r->at(i));

    // Check if data is within node's boundaries:
    // ===============================================
    if (!IsPointInsideBoundary(p))
    {
        // Warning message:
        cout << "point " << p <<" is outside domain" << endl;
        return;
    }

    // Check if maximum tree depth has been reached:
    // ===============================================
    // If yes, append point's index to leaf and RETURN up the stack
    if (HasNodeReachMaxDepth())
    {
        // Append point:
        if (write_data)
        {
          this->ix.push_back(i);
          this->x_count++;
        }

        // Diagnostics:
        if (0)
        {
            cout << "i = " << i << endl;
            cout << "size of ix = " << this->ix.size() << endl;
            cout << "this->x_count = " << this->x_count << endl;
            cout << "x_left = " << this->x_left << endl;
            cout << "x_right = " << this->x_right << endl;
            for (int j = 0; j < this->ix.size(); j++)
            {
                cout << "j = " << j << ", ix[j] = " << this->ix[j] <<", value = " << r->at(this->ix[j]) <<endl;
            }
        }

        // Return control to calling stack:
        return;
    }

    // Determine which subnode to insert point:
    // ========================================
    int node_index = WhichSubNodeDoesItBelongTo(p);

    // Check if subnode needs to be created:
    // ====================================
    if ( !DoesSubNodeExist(node_index) )
    {
        CreateSubNode(node_index);
    }

    // Insert point into subnode:
    // ==========================
    this->subnode[node_index]->insert(i,r,write_data);

} // node_TYP::insert

// insert_all method:
// ================================================================================================================
void node_TYP::insert_all(arma::vec * r)
{
  for (int i = 0; i < r->size(); i++)
  {
      this->insert(i,r,true);
  }
}

// =================================================================================================================
bool node_TYP::IsPointInsideBoundary(double p)
{
    // Objective:
    // if p is inside the boundaries of the node, return true, otherwise false

    // Define boundaries of node:
    double x_left  = this->x_left;
    double x_right = this->x_right;

    // Create boolean result:
    bool flag = ((p >= x_left) && (p <= x_right));

    return flag;
}

// ================================================================================================================
bool node_TYP::HasNodeReachMaxDepth()
{
    int depth     = this->depth;
    // int depth_max = this->depth_max;

    if (depth >= tree_params->depth_max)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

// ================================================================================================================
int node_TYP::WhichSubNodeDoesItBelongTo(double p)
{
    //   +------------------+------------------+
    //   |  node_left = 1  |  node_right = 0  |

    // Center location of present node:
    double x_center = this->x_center;

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
node_TYP * node_TYP::find(double xq)
{

    // Check if data is within node's boundaries:
    if (!IsPointInsideBoundary(xq))
    {
        // Warning message:
        cout << "point " << xq <<" is outside domain" << endl;
        return NULL;
    }

    // Check if we have reached maximum depth:
    if (HasNodeReachMaxDepth())
    {
        return this;
    }

    // Determine which subnode to move into:
    int node_index = WhichSubNodeDoesItBelongTo(xq);

    // Check if subnode exists:
    if (!DoesSubNodeExist(node_index))
    {
        return NULL;
    }

    // Drill further into subnode:
    return this->subnode[node_index]->find(xq);

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
void node_TYP::CreateSubNode(int node_index)
{
    // Attributes for new subnode:
    int depth     = this->depth + 1;
    double x_left;
    double x_right;

    switch (node_index)
    {
    case 0: // Right subnode:
        {
            // Attributes for new subnode:
            x_left  = this->x_center;
            x_right = this->x_right;

            // Exit:
            break;
        }
    case 1: // Left subnode:
        {
            // Attributes for new subnode:
            x_left  = this->x_left;
            x_right = this->x_center;

            // Exit:
            break;
        }
    }

    // Create new subnode:
    this->subnode[node_index] = new node_TYP(x_left, x_right, depth, tree_params);
}

// ================================================================================================================
node_TYP * node_TYP::get_subnode(int index)
{
    return this->subnode[index];
}
