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

    // Initialize the resamaple count:
    resample_count = zeros<uvec>(_numIONS);

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
    bt_params.dimensionality = 1;
    bt_params.min       = {params->mesh_params.Lx_min};
    bt_params.max       = {params->mesh_params.Lx_max};
    bt_params.max_depth = {+6};

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

bool RS_TYP::IsResamplingNeeded(params_TYP * params, vector<ions_TYP> * IONS, mesh_TYP * mesh, vector<particle_tree_TYP> * tree, int ss)
{
  // Get computational particle density:
  vec ncp_m = IONS->at(ss).ncp_m.subvec(0,_Nx-1);
  vec& x_m  = mesh->xm;

  // x grid for leaf_x nodes:
  vec& xq = tree->at(ss).xq;
  double dxq = xq[2] - xq[1];
  vec x_min = xq - dxq/2;
  vec x_max = xq + dxq/2;
  vec ncp_q = zeros(xq.size());

  // Average ncp_m over each x node from binary tree:
  for (int xx = 0; xx < tree->at(ss).leaf_x.size(); xx++)
  {
    uvec indices = find(x_m >= x_min[xx] && x_m < x_max[xx]);
    ncp_q[xx] = mean(ncp_m.elem(indices));
  }

  // Calculate metric:
  arma::vec dncp_q = ncp_q - _mean_ncp_m[ss]/3;

  // Check if dncp_q becomes negative:
  int count = 0;
  for (int q = 0; q < ncp_q.size(); q++)
  {
    if (dncp_q.at(q) < 0)
    {
      return true;
    }
  }

  return false;
}

void RS_TYP::check_for_nans(params_TYP * params, vector<ions_TYP> * IONS)
{
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    for (int ss = 0; ss < _numIONS; ss++)
    {
      for (int ii = 0; ii < IONS->at(ss).N_CP_MPI; ii++)
      {
        double vpar = IONS->at(ss).v_p(ii,0);
        double vper = IONS->at(ss).v_p(ii,1);
        double a    = IONS->at(ss).a_p(ii);

        if (isnan(vpar) || (vpar == -1))
          cout << "RS operator, problematic vpar" << endl;

        if (isnan(vper) || (vper == -1))
          cout << "RS operator, problematic vper" << endl;

        if (a == -1)
          cout << "RS operator, problematic a == -1" << endl;
      }
    }
  } // MPI_COLOR
}

void RS_TYP::ApplyResampling_AllSpecies(params_TYP * params, mesh_TYP * mesh, vector<ions_TYP> * IONS, vector<particle_tree_TYP> * particle_tree)
{
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    for (int ss = 0; ss < _numIONS; ss++)
    {
      if (true)//(IsResamplingNeeded(params,IONS,mesh,particle_tree,ss))
      {

        if (params->mpi.IS_PARTICLES_ROOT)
        {
          //cout << "Apply resample, species " << ss << endl;
        }
        
        particle_tree->at(ss).populate_tree("binary and quad");
        particle_tree->at(ss).resample_distribution();
        resample_count[ss] = resample_count[ss] + 1;

        // Release memory if needed:
        if (resample_count[ss] % 1 == 0)
        {
          particle_tree->at(ss).release_memory();
        }

      }
    }
  } // MPI_COLOR
}
