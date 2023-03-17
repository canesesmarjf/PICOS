#include "resamplingOperator.h"

RS_TYP::RS_TYP()
{
  cout << "Constructor for RS_TYP" << endl;
}

RS_TYP::RS_TYP(params_TYP * params, vector<ions_TYP> * IONS, vector<binaryTree_TYP> * tree)
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
  //tree->resize(_numIONS);

  // Construct all elements of binary tree vector:
  for (int ss = 0; ss < _numIONS; ss++)
  {
    double x_left  = params->mesh_params.Lx_min;
    double x_right = params->mesh_params.Lx_max;
    int depth_max  = 6;
    int num_elem   = IONS->at(ss).N_CP_MPI;
    tree->push_back(binaryTree_TYP(x_left, x_right, depth_max, num_elem));
  }

}

bool RS_TYP::IsResamplingNeeded(params_TYP * params, vector<ions_TYP> * IONS, int ss)
{
  // Get computational particle density:
  vec ncp_m = IONS->at(ss).ncp_m.subvec(0,_Nx-2);

  // // Mean computational particle density
  // double mean_ncp_m = sum(ncp_m)*_dx/_L;

  // Calculate metric:
  arma::vec dncp_m = ncp_m - _mean_ncp_m[ss]/2;

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
        tree->at(ss).insert_all(&IONS->at(ss).x_p);
        
        // Determine deficit/surplus nodes:
        // Calculate required memory locations to repurpose
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
        cout << "Do not resample, species " << ss << endl;
      }
    }

  } // MPI_COLOR

}