#include "particle_tree.h"

using namespace arma;

/* The following is an example of how to use the particle tree in PICOS++:
vector<particle_tree_TYP> IONS_tree;
IONS_tree.leaf_x[xx]->ip; // Get the indices of particles in the x_node at xx
*/

// ======================================================================================
particle_tree_TYP::particle_tree_TYP(bt_params_TYP * bt_params, qt_params_TYP * qt_params, vec * x_p, mat * v_p, vec * a_p)
{
  // Store pointer to tree parameters:
  this->bt_params = bt_params;
  this->qt_params = qt_params;

  // Store pointers to data:
  this->x_p = x_p;
  this->v_p = v_p;
  this->a_p = a_p;

  // Instance of binary tree for x dimension:
  bt = bt_TYP(bt_params,x_p,a_p);

  // Instance of quad tree vector, one for each node of binary tree:
  int Nx = bt_params->num_nodes;
  qt.resize(Nx);

  // Allocate memory and resize leaf_x and leaf_v:
  leaf_v.resize(Nx);
  leaf_x.resize(Nx);

  // Create x dimension query grid for binary_tree:
  this->create_x_query_grid();

  // Allocate memory and resize ip_count:
  ip_count = zeros<ivec>(Nx);
}

// ======================================================================================
void particle_tree_TYP::create_x_query_grid()
{
  // Number of x grid nodes:
  int Nx = pow(2,bt_params->max_depth[0]);

  // Calculate total length and node length in x-space:
  double Lx = bt_params->max[0] - bt_params->min[0];
  double dx = Lx/Nx;

  // Create query grid:
  this->xq = bt_params->min[0] + dx/2 + regspace(0,Nx-1)*dx;
}

// ======================================================================================
void particle_tree_TYP::calculate_leaf_x()
{
  // Assemble leaf_x and ip_count:
  int Nx = bt_params->num_nodes;
  for (int xx = 0; xx < Nx ; xx++)
  {
   leaf_x[xx] = bt.find(xq(xx));
   ip_count[xx] = leaf_x[xx]->ip_count;
  }
}

// ======================================================================================
void particle_tree_TYP::get_mean_ip_count()
{
  // This operation requires the division between two integers; hence, care needs to be taken since this leads to loss of precision. Method used the the following:
  // 1- Cast numerator as a double so that division retains fractional part
  // 2- Round UP to the nearest LARGEST integer so as to over-estimate total number of particles. We later test for total particle usage to avoid "out of bounds" access.

  // This calculation is needed to determine the target particle count we need to have on each x-node in order to achieve a flat or near flat computational particle density along the x dimension:

  // Total number of x-nodes:
  int Nx = bt_params->num_nodes;

  // Total number of particles indexed in x-nodes:
  int p_count_sum = sum(ip_count);

  // Calculate the mean particle number along the x dimension:
  mean_ip_count = ceil((double)p_count_sum/Nx);
}

// ======================================================================================
void particle_tree_TYP::assemble_qt_vector()
{
  int Nx = bt_params->num_nodes;
  for (int xx = 0; xx < Nx ; xx++)
  {
    // Calculate surplus or deficit of ip_count profile:
    int delta_ip_count = ip_count[xx] - mean_ip_count;

    // Populate quadtree[xx] only if leaf_x[xx] is surplus
    if (delta_ip_count > 0)
    {
      // Data for quadtree:
      vector<uint> ip = leaf_x[xx]->ip;

      // Populate quadtree:
      qt[xx] = qt_TYP(qt_params,ip,this->v_p);

      // Populate quadtree:
      qt[xx].populate_tree();
    }
  }
}

// ======================================================================================
void particle_tree_TYP::calculate_leaf_v()
{
  for (int xx = 0; xx < leaf_x.size() ; xx++)
  {
    // Create quadtree[xx] only if leaf_x[xx] is surplus:
    int delta_ip_count = ip_count[xx] - mean_ip_count;
    if (delta_ip_count > 0)
    {
      // Extract all leaf nodes:
      leaf_v[xx] = qt[xx].get_leaf_nodes();
    }
    else
    {
      leaf_v[xx] = {NULL};
    }
  }
}

// ======================================================================================
void particle_tree_TYP::populate_binary_tree()
{
  // Create a vector with pointers to the data:
  vector<vec *> x_data = {x_p};

  // Insert 1D points into binary tree:
  bt.insert_all(x_data);

  // Calculate leaf_x list:
  this->calculate_leaf_x();

  // Calculating the mean number of particles per node:
  // This is to be used by quad tree
  this->get_mean_ip_count();
}

// ======================================================================================
void particle_tree_TYP::populate_quad_tree()
{
  // This method requires that we have a valid binary tree.
  // Assumes we have used this->populate_binary_tree() prior.

  // Assemble quadtrees for each leaf_x element:
  this->assemble_qt_vector();

  // Assembling leaf_v 2D vector:
  this->calculate_leaf_v();
}

// ======================================================================================
void particle_tree_TYP::populate_tree(string calculation_type)
{
  // Clear contents to prevent doubling up on data:
  this->clear_all_contents();

  // Select the calculation type:
  if (calculation_type == "binary and quad")
  {
    this->populate_binary_tree();
    this->populate_quad_tree();
  }
  else if (calculation_type == "binary only")
  {
    this->populate_binary_tree();
  }
  else
  {
    cout << "Error, calculation type not recognized" << endl;
  }
}

// =======================================================================================
void particle_tree_TYP::clear_all_contents()
{
  // Clear contents of binary tree:
  this->bt.clear_all();

  // Clear contents of quad tree:
  int Nx = this->qt.size();
  for (int xx = 0; xx < Nx; xx++)
  {
    this->qt[xx].clear_tree();
  }

}

// =======================================================================================
void particle_tree_TYP::save_leaf_v_structure(string output_dir)
{
  // This method saves data from leaf_v so that the entire tree structure (binary + quad) can be visualized in MATLAB:

  // For every leaf_x node, there is a corresponding quadtree in velocity space.
  // From each of those quadtrees, we have extracted leaf_v which contains a list of all the leaf nodes in velocity space for a given "x" location.
  // What we are doing in this section is to save the particle count, coordinates and dimensions of the leaf_v nodes on each "x" location,
  // This data can then be plotted in MATLAB as a group of nested bisected quadrants in velocity space and then we can superimpose the particle data to confirm whether of not the quadtree algorithm is working.
  // We can visually inspect of the particle count, coordinate and dimensions of the node match the observed number of particles in that region based on the particle data.

  int Nx = bt_params->num_nodes;
  for (int xx = 0; xx < Nx ; xx++)
  {
    if (leaf_v[xx][0] != NULL)
    {
      // Get the number of leaf nodes on current leaf_v:
      int num_v_nodes = leaf_v[xx].size();

      // Create variables to contain data:
      ivec node_ip_count(num_v_nodes);
      mat  node_center(num_v_nodes,2);
      mat  node_dim(num_v_nodes,2);

      // Assemble data:
      for (int vv = 0; vv < num_v_nodes; vv++)
      {
        node_ip_count(vv)   = leaf_v[xx][vv]->ip_count;
        node_center.row(vv) = leaf_v[xx][vv]->center.t();
        node_dim.row(vv)    = leaf_v[xx][vv]->max.t() - leaf_v[xx][vv]->min.t();
      }

      // Save data:
      string file_name;
      file_name = \
      output_dir + "/"  + "leaf_v_" + "ip_count" + "_xx_" + to_string(xx) + ".csv";
      node_ip_count.save(file_name, csv_ascii);

      file_name = \
      output_dir + "/"  + "leaf_v_" + "node_center" + "_xx_" + to_string(xx) + ".csv";
      node_center.save(file_name, csv_ascii);

      file_name = \
      output_dir + "/"  + "leaf_v_" + "node_dim" + "_xx_" + to_string(xx) + ".csv";
      node_dim.save(file_name, csv_ascii);
    }
  }
}

// =======================================================================================
void particle_tree_TYP::assess_conservation(string output_dir, string suffix)
{
  // Only requires that we assemble the binary tree and leaf_x:

  int Nx = bt_params->num_nodes;
  vec m_t(Nx, fill::zeros);
  vec p_x(Nx, fill::zeros);
  vec p_r(Nx, fill::zeros);
  vec KE(Nx, fill::zeros);

  for (int xx = 0; xx < Nx ; xx++)
  {
    if (leaf_x[xx] != NULL)
    {
      uvec ip = conv_to<uvec>::from(leaf_x[xx]->ip);
      m_t[xx] = sum(a_p->elem(ip));
      mat v_p_subset = v_p->rows(ip);
      p_x[xx] = dot(a_p->elem(ip),v_p_subset.col(0));
      p_r[xx] = dot(a_p->elem(ip),v_p_subset.col(1));

      KE[xx] = (dot(a_p->elem(ip),pow(v_p_subset.col(0),2)) + dot(a_p->elem(ip),pow(v_p_subset.col(1),2)));
    }
  }

  cout << "Suffix: " + suffix + ". Total KE = " + to_string(sum(KE)) << endl;

  m_t.save(output_dir + "/" + "m_profile_" + suffix + ".csv",csv_ascii);
  p_x.save(output_dir + "/" + "p_x_profile_" + suffix + ".csv",csv_ascii);
  p_r.save(output_dir + "/" + "p_r_profile_" + suffix  + ".csv",csv_ascii);
  KE.save(output_dir + "/" + "KE_profile_" + suffix + ".csv",csv_ascii);
}

// ======================================================================================
void particle_tree_TYP::downsample_surplus_nodes(vector<uint> * ip_free)
{
  // Create vranic down-sampling object:
  vranic_TYP vranic;

  // Set the min and maximum number of computational particles per node:
  int N_min = 7;
  int N_max = 300;

  // Set the particle set sizes for the input and output sets:
  int N; // Input set
  int M = 6; // Output set0
  int particle_surplus;

  // Down-sample distribution and populate ip_free vector:
  for (int xx = 0; xx < leaf_x.size() ; xx++)
  {
    // Calculate particle surplus
    particle_surplus = leaf_x[xx]->ip_count - mean_ip_count;

    // Apply Vranic method if x-node is NOT NULL:
    if (leaf_v[xx][0] != NULL)
    {
      // Total number of v-nodes for this xx position:
      int Nv = leaf_v[xx].size();

      // Assemble depth and ip_count vectors:
      vec depth(Nv);
      vec ip_count(Nv);
      vec node_metric(Nv);
      for(int vv = 0; vv < Nv; vv++){depth(vv) = leaf_v[xx][vv]->depth;}
      for(int vv = 0; vv < Nv; vv++){ip_count(vv) = leaf_v[xx][vv]->ip_count;}

      // Calculate metric for every node based on depth*ip_count (elementwise multiplication). The higher the depth AND number of particles, the higher the metric:
      node_metric = depth%ip_count;

      // Create sorted list starting from highest depth to lowest:
      // By giving priority to nodes with highest metric, we are operating on the regions that minimize the changes to the distribution function since they are more localized in velocity space:
      uvec sorted_index_list = sort_index(node_metric,"descend");

      // Loop over sorted leaf_v and apply vranic method:
      for (int vv = 0; vv < sorted_index_list.n_elem; vv++)
      {
        // Sorted index:
        int tt = sorted_index_list(vv);

        // Total number of particles in leaf_v cube:
        N = leaf_v[xx][tt]->ip_count;

        // Test if current node is suitable for resampling:
        if (N < N_min){continue;}

        // Upper limit to the number of particles to down-sample in present cell:
        if (N > N_max){N = N_max;}

        // Check that we do not go into deficit:
        if (particle_surplus - (N-M) < 0){N = particle_surplus + M;}

        // Create objects for down-sampling:
        merge_cell_TYP set_N(N);
        merge_cell_TYP set_M(M);

        // Define particle indices for set N:
        uvec ip = conv_to<uvec>::from(leaf_v[xx][tt]->ip);
        uvec ip_subset = ip.head(N);

        // Assign particle attributes for set N:
        mat v_p_subset = v_p->rows(ip_subset);
        set_N.xi = x_p->elem(ip_subset);
        set_N.yi = v_p_subset.col(0);
        set_N.zi = v_p_subset.col(1);
        set_N.wi = a_p->elem(ip_subset);

        for (int jj = 0; jj < N; jj++)
        {
          if (set_N.wi(jj) == -1)
            cout << "RS, wi == -1" << endl;

          //if (set_N.wi(jj) < 1e-6)
            //cout << "RS, wi < 1e-6" << endl;
        }

        // Calculate set M based on set N:
        // We have observed that using the 2D implementation we get conservation of energy down to the machine precision (1E-14); however, using the 3D method, we get conservation of energy not as good (1E-5):
        vranic.down_sample_node_2D(&set_N, &set_M);
        //vranic.down_sample_node_3D(&set_N, &set_M);

        // Diagnostics:
        if (false)
        {
          // Print statistics:
          cout << "Set N: " << endl;
          vranic.print_stats(&set_N);

          // Print statistics:
          cout << "Set M: " << endl;
          vranic.print_stats(&set_M);
        }

        // Apply changes to distribution function:
        // Loop over all particles that underwent down-sampling:
        for (int ii = 0; ii < N; ii++)
        {
          // Get global index:
          int jj = ip(ii);

          // Apply changes to x_p, v_p, a_p OR create new mem locations:
          if (ii < set_M.n_elem)
          {
            // Down-sample global distribution:
            // Use set_N for xi in order to remove oscillations in density:
            // Use set_M for all other quantities (v and a):

            // Set N:
            (*x_p)(jj) = set_N.xi(ii);

            if (isnan(set_M.yi(ii)))
              cout << "RS, isnan" << endl;

            if (isnan(set_M.zi(ii)))
              cout << "RS, isnan" << endl;

            if (set_M.wi(ii) == -1)
              cout << "wi == -1" << endl;

            // Set M:
            (*v_p)(jj,0) = set_M.yi(ii);
            (*v_p)(jj,1) = set_M.zi(ii);
            (*a_p)(jj)   = set_M.wi(ii);

          }
          else
          {
            // Record global indices that correspond to memory locations that are to be repurposed:
            ip_free->push_back(jj);

            // Set values to -1 to flag them as undefined values memory locations:
            (*x_p)(jj)   = -1;
            (*v_p)(jj,0) = -1;
            (*v_p)(jj,1) = -1;
            (*a_p)(jj)   = -1;

            // Decrement particle surplus only after ip_free.push_back operation:
            particle_surplus--;
          }
        } // for ip loop
      } // sorted index Loop
    } // if not NULL
  } // xx loop

}

// ======================================================================================
void particle_tree_TYP::upsample_deficit_nodes(vector<uint> * ip_free)
{
  // Initialize replication flag:
  int ip_free_flag = 0;
  int particle_deficit;

  // Vector to keep track of particle counts in x-nodes:
  int Nx = leaf_x.size();
  ivec node_counts(Nx, fill::zeros);
  for (int xx = 0; xx < Nx; xx++)
    node_counts(xx) = leaf_x[xx]->ip_count;

  // Vector to keep track how many times we need to replicate particles per node:
  ivec rep_num_vec(Nx,fill::ones);
  for (int xx = 0; xx < Nx; xx++)
    rep_num_vec(xx) = ceil((double)mean_ip_count/node_counts(xx));

  // Layers:
  // The idea is to use the ip_free memory locations to fill in the locations where particles are needed but we aim to this in a layered manner so that we topup the deficits in stages. This prevents neglecting certain x-nodes in cases where insufficient number of particles where extracted from the down sampling process
  vec layer_fraction = {1/2, 1/4, 1/8, 1/16, 1/16};
  ivec layer(5);
  layer(0) = (int)round((double)mean_ip_count/2);
  layer(1) = (int)round((double)mean_ip_count/4);
  layer(2) = (int)round((double)mean_ip_count/8);
  layer(3) = (int)round((double)mean_ip_count/16);
  layer(4) = mean_ip_count - sum(layer.subvec(0,layer.n_elem - 2));

  for (int ll = 0; ll < layer.n_elem; ll++)
  {
    // cout << "ll = " << ll << endl;
    for (int xx = 0; xx < Nx ; xx++)
    {
      // If ip_free is empty, then stop replication:
      if (ip_free_flag == 1)
        break;

      // Calculate particle deficit
      int target_counts = sum(layer.subvec(0,ll));
      particle_deficit = node_counts(xx) - target_counts;

      if (particle_deficit < 0)
      {
        // Number of particles to replicate (original particles in node):
        int num_0 = leaf_x[xx]->ip_count;

        // Replication number: represents how many times a particle needs to be replicated
        // rep_num - 1 gives you the number of new particles per parent particle
        // int rep_num = ceil((double)layer(ll)/num_0);
        int rep_num = rep_num_vec(xx);

        // Loop over all particles in node:
        for (int ii = num_0 - 1; ii >= 0; ii--)
        {
          // Get global index of particle to replicate:
          uint jj = leaf_x[xx]->ip[ii];

          // On the last layer, remove the indexes since they have been fully used:
          if (ll == layer.n_elem)
          {
            leaf_x[xx]->ip.pop_back();
            leaf_x[xx]->ip_count--;
          }

          // Diagnostics:
          if (leaf_x[xx]->ip_count < 0)
            cout << "error:" << endl;

          // Parent particle attributes:
          double xi = (*x_p)(jj);
          double yi = (*v_p)(jj,0);
          double zi = (*v_p)(jj,1);
          double wi = (*a_p)(jj);

          // Calculate number of new daughter particles to create:
          int num_free_left = ip_free->size();
          int num_deficit_left = -particle_deficit;
          int num_requested = rep_num - 1;
          ivec num_vec = {num_requested, num_free_left, num_deficit_left};
          int num_new = num_vec(num_vec.index_min());

          // Create num_new daughter particles:
          for (int rr = 0; rr < num_new; rr++)
          {
            uint jj_free = ip_free->back();
            ip_free->pop_back();
            (*x_p)(jj_free)   = xi;//- sign(xi)*0.05;
            (*v_p)(jj_free,0) = yi;
            (*v_p)(jj_free,1) = zi;
            (*a_p)(jj_free)   = wi/((double)num_new + 1.0);

            // bool res = wi/((double)num_new + 1.0) < 1e-10;
            // if (res)
            // {
            //   cout << "upsampling wi < 1e-10" << endl;
            // }

            // Modify deficit:
            particle_deficit++;
            node_counts(xx)++;
          }

          // Adjust weight of parent particle to conserve mass:
          (*a_p)(jj) = wi/((double)num_new + 1.0);

          // if size of ip_free vanishes, then stop all replication:
          int num_free = ip_free->size();
          if (num_free == 0)
          {
            ip_free_flag = 1;
            break;
          }

          // If deficit becomes +ve, then stop:
          if (particle_deficit >= 0)
          {
            break;
          }
        } // particle Loop

        if (particle_deficit != 0)
        {
          //abort();
        }

      } // deficit if
    } // xx loop
  } // ll Loop

}

// ======================================================================================
void particle_tree_TYP::resample_distribution()
{
  // Create vector to hold indices of computational particles that can be reused:
  // -------------------------------------------------------------------------------------
  // Could make this vector a private member of particle_tree_TYP. This would mean that we allocate memory to it once and over time it grows to a size that is sufficient for the problem. this means that we are not faced with multiple reallocations when the capacity is exceeded evertime we resample
  vector<uint> ip_free;

  // Down-sample surplus nodes:
  // -------------------------------------------------------------------------------------
  this->downsample_surplus_nodes(&ip_free);

  // Up-sample deficit nodes via replication:
  // -------------------------------------------------------------------------------------
  this->upsample_deficit_nodes(&ip_free);

  // Clear contents of particle tree:
  // -------------------------------------------------------------------------------------
  // Since the distribution has been changed, the particle tree is no longer current and needs to be cleared.
  this->clear_all_contents();

}
