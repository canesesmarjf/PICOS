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
  this->mean_elems_per_node = round(num_elem/this->num_nodes);
  this->node_centers = arma::linspace(x_left,(x_right-dx),this->num_nodes) + dx/2;
}

// ================================================================================================================
binaryTree_TYP::binaryTree_TYP(double x_left, double x_right, int depth_max, int num_elem)
{
  // Interface variables:
  tree_params = new tree_params_TYP(x_left, x_right, depth_max, num_elem);

  // Allocate memory to node list, which is a vector of pointers:
  node_list.resize(tree_params->num_nodes);

  // Allocate memory to count profile, which is an arma::vec created ti store x_count for all nodes
  count_profile.set_size(tree_params->num_nodes);
  delta_profile.set_size(tree_params->num_nodes);
  total_deficit_elems = 0;

  // Create root node:
  int depth_root = 0;
  root = new node_TYP(tree_params->x_left,tree_params->x_right,depth_root,tree_params);

  // Assemble empty binary tree:
  assemble_empty_tree();

  // Assemble list of nodes into a vector of pointers:
  assemble_node_list();
}

// ================================================================================================================
void binaryTree_TYP::insert_all(arma::vec * r)
{
  // Insert points into nodes:
  this->root->insert_all(r);

  // Populate the count_profile to allow classification of deficit and surplus nodes:
  assemble_count_profile();
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
void binaryTree_TYP::assemble_count_profile()
{
  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    // Copy the value of x_value to each element of count_profile;
    count_profile.at(nn) = node_list.at(nn)->x_count;
  }
}

// ================================================================================================================
void binaryTree_TYP::assemble_empty_tree()
{
  // The objective of this method is to traverse the entire tree down to the final layer of nodes.
  // This is done by inserting positions corresponding to the node_centers so as to create all the nodes
  // at the final depth while not writting any data to the nodes.
  // This process allows us to preallocate memory to the final nodes based on the estimated
  // mean elements per node (See tree_params)

  for (int nn = 0; nn < tree_params->num_nodes; nn++)
  {
    bool write_data = false;
    root->insert(nn,&tree_params->node_centers,write_data);
  }
}

// ================================================================================================================
void binaryTree_TYP::clear_all()
{
  ix_repurpose.clear();
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
void binaryTree_TYP::calculate_delta_profile()
{
  // Based on the mean_elems_per_node, we can calculate whether a node is in deficit or surplus.
  // This calculation is done on a per MPI rank basis:
  delta_profile = count_profile - tree_params->mean_elems_per_node;

  // Calculate total number of memory location to repurpose using the deficit nodes:
  arma::uvec rng = arma::find(delta_profile < 0);
  arma::ivec deficit_elems   = -1*delta_profile.elem(rng);
  total_deficit_elems = arma::sum(deficit_elems);
}

// ================================================================================================================
void binaryTree_TYP::gather_all_surplus_indices()
{
  for (int nn = 0; nn < node_list.size(); nn++)
  {
    if (delta_profile[nn] > 0)
    {
      // Number of surplus computational particles:
      int num_surplus = delta_profile[nn];

      // Determine indices of all surplus computational particles:
      std::vector<uint>& ix = node_list[nn]->ix;
      auto it_end   = ix.end();
      auto it_start = std::prev(it_end,num_surplus);
      std::vector<uint> ix_subset(it_start,it_end);

      // Copy indices of surplus computational particles to ix_repurpose:
      for (int jj = 0; jj < num_surplus; jj++)
      {
        ix_repurpose.push_back(ix.back());
        ix.pop_back();
      }
    }
  }
}

// ================================================================================================================
void binaryTree_TYP::renormalize_surplus_nodes(ions_TYP * IONS)
{
  for (int nn = 0; nn < node_list.size(); nn++)
  {
    if (delta_profile[nn] > 0)
    {
      // Number of surplus computational particles:
      int num_surplus = delta_profile[nn];

      // New number of computational particles:
      int x_count_new = node_list[nn]->x_count - num_surplus;

      // Calculate new particle weight:
      int delta_x_count = num_surplus;
      double weight_factor = 1 + static_cast<double>(delta_x_count)/x_count_new;

      // Loop over all particles in node and rescale weight:
      for (int jj = 0; jj < node_list[nn]->ix.size(); jj++)
      {
        int ii = node_list[nn]->ix[jj];
        IONS->a_p(ii) *= weight_factor;
      }

      // Update x_count:
      node_list[nn]->x_count = x_count_new;
    }
  }

  // Update count profile:
  assemble_count_profile();
};

// ================================================================================================================
void binaryTree_TYP::renormalize_deficit_nodes(ions_TYP * IONS)
{
  // Create random number generator variables:
  arma::arma_rng::set_seed_random();

  for (int nn = 0; nn < node_list.size(); nn++)
  {
    if (delta_profile[nn] < 0)
    {
      // Number of deficit computational particles:
      int num_deficit = delta_profile[nn];

      // Number of computational particles indexed in present node:
      int x_count_start = node_list[nn]->x_count;
      int x_count_end   = x_count_start;

      // Vector containing indices of computational particles in present node:
      std::vector<uint>& ix = node_list[nn]->ix;

      /*
      // New number of computational particles:
      int x_count_new = node_list[nn]->x_count - num_deficit;

      // Calculate new particle weight:
      int delta_x_count = num_deficit;
      double weight_factor = 1 + static_cast<double>(delta_x_count)/x_count_new;

      // Loop over all particles in node, create copies and rescale weight:
      for (int jj = 0; jj < node_list[nn]->ix.size(); jj++)
      {
        int ii = node_list[nn]->ix[jj];
        IONS->a_p(ii) *= weight_factor;
      }

      // Update x_count:
      node_list[nn]->x_count = x_count_new;
      */

      for (int jj = 0; jj < -num_deficit; jj ++)
      {
        // Select a computational particle from present node at random:
        uint jj_random = arma::randi<uint>(arma::distr_param(0,x_count_start-1));
        uint ii = ix[jj_random];

        // Copy attributes of iith computational particle:
        double x_p    = IONS->x_p(ii);
        double vpar_p = IONS->v_p(ii,0);
        double vper_p = IONS->v_p(ii,1);
        double a_p    = IONS->a_p(ii);

        // Create a copy of the iith particle with the memory locations to be repurposed:
        uint ii_copy = ix_repurpose.back();
        IONS->x_p(ii_copy)   = x_p;
        IONS->v_p(ii_copy,0) = vpar_p;
        IONS->v_p(ii_copy,1) = vper_p;

        // Split the computational particle weight to preserve mass:
        IONS->a_p(ii)      = a_p/2;
        IONS->a_p(ii_copy) = a_p/2;

        // Remove repurposed index:
        ix_repurpose.pop_back();

        // Increment x_count:
        x_count_end++;
      }

      // Update x_count:
      node_list[nn]->x_count = x_count_end;
    }
  }

  // Update count profile:
  assemble_count_profile();
};

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
    this->ix.reserve(2*tree_params->mean_elems_per_node);
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
