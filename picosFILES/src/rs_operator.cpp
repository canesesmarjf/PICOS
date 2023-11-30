#include "rs_operator.h"

RS_TYP::RS_TYP()
{
  cout << "RS_TYP constructor" << endl;
}

// =======================================================================================
RS_TYP::RS_TYP(params_TYP * params, CS_TYP *CS, vector<ions_TYP> * IONS, vector<particle_tree_TYP> * particle_tree)
{
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    cout << "Overloaded constructor for RS_TYP" << endl;

    // Get geometric data:
    _L  = params->mesh_params.Lx_max - params->mesh_params.Lx_min;
    _dx = params->mesh_params.dx;
    _Nx = params->mesh_params.Nx;

    // Number of ion species:
    _numIONS = IONS->size();

    // Mean computational particle density:
    _mean_ncp_m = zeros(_numIONS);
    for (int ss = 0; ss < _numIONS; ss++)
    {
      vec ncp_m = IONS->at(ss).ncp_m.subvec(0,_Nx-2);
      _mean_ncp_m.at(ss) = sum(ncp_m)*_dx/_L;
    }
    _mean_ncp_m.print("_mean_ncp_m = ");

    // Create empty binary tree vector:
    //particle_tree->resize(_numIONS);

    // Calculate the maximum velocity to use in quad tree grid:
    double E_max = F_E*30e3/CS->energy;
    double mass = IONS->at(1).M;
    double v_p_max = sqrt(2*F_KB_DS*E_max/mass);

    // Create tree parameters:
    bt_params_TYP bt_params;
    bt_params.dimensionality = 1;
    bt_params.min       = {params->mesh_params.Lx_min};
    bt_params.max       = {params->mesh_params.Lx_max};
    bt_params.max_depth = {+6};

    qt_params_TYP qt_params;
    qt_params.min = {-v_p_max,-v_p_max};
    qt_params.max = {+v_p_max,+v_p_max};
    qt_params.min_depth = +5;
    qt_params.max_depth = +6;
    qt_params.min_count = 7;

    // Construct all elements of particle_tree vector:
    for (int ss = 0; ss < _numIONS; ss++)
    {
      vec& x_p = IONS->at(ss).x_p;
      mat& v_p = IONS->at(ss).v_p;
      vec& a_p = IONS->at(ss).a_p;
      particle_tree->push_back(particle_tree_TYP(&bt_params,&qt_params,&x_p,&v_p,&a_p));
    }
  }
}
/*
bool RS_TYP::IsResamplingNeeded(params_TYP * params, vector<ions_TYP> * IONS, int ss)
{
  // Get computational particle density:
  vec ncp_m = IONS->at(ss).ncp_m.subvec(0,_Nx-2);

  // // Mean computational particle density
  // double mean_ncp_m = sum(ncp_m)*_dx/_L;

  // Calculate metric:
  arma::vec dncp_m = ncp_m - _mean_ncp_m[ss]*0.66;
  // arma::vec dncp_m = ncp_m - _mean_ncp_m[ss]/3;

  // Check if dncp_m becomes negative:
  for (int m = 0; m < _Nx; m++)
  {
    if (dncp_m.at(m) < 0)
    {
      return true;
    }
  }

  return false;
}

void RS_TYP::ApplyResampling_AllSpecies(params_TYP * params, mesh_TYP * mesh, vector<ions_TYP> * IONS, vector<binaryTree_TYP> * tree)
{
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    for (int ss = 0; ss < _numIONS; ss++)
    {
      if (IsResamplingNeeded(params,IONS,ss))
      {
        cout << "Apply resample, species " << ss << endl;
        // Clear all previous states:
        tree->at(ss).clear_all();

        // Insert all computational particle positions into binary tree:
        tree->at(ss).insert_all(&IONS->at(ss).x_p);

        // Determine which nodes have surplus or deficit of computational particles:
        // calculate_delta_profile(&tree->at(ss));
        tree->at(ss).calculate_delta_profile();

        // Gather all computational particle indices from surplus nodes:
        tree->at(ss).gather_all_surplus_indices();

        // Renormalize distributions in all surplus nodes:
        tree->at(ss).renormalize_surplus_nodes(&IONS->at(ss));

        // Resample and renormalize distributions in all deficit nodes:
        tree->at(ss).renormalize_deficit_nodes(&IONS->at(ss));

        // Loop over surplus nodes to gather memory locations
        //    + Copy surplus indices to repurpose_list
        //    + Set -1 to repurposed indices
        //    + scale weights of current surplus cell using a factor a_new = a_old*(1 + dN/N1),
        // where dN is number of indices repurposed and N1 is new number of indices which should approach the mean
        //  end loop
        // Loop over deficit nodes:
        //    + replicate each particle until all indices_to_be_repurposed have been used
        //    + Scale weight by 1/2
        //  end loop
      }
      else
      {
        // cout << "Do not resample, species " << ss << endl;
      }
    }

  } // MPI_COLOR

}
*/
