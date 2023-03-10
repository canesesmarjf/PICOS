#include "resamplingOperator.h"

RS_operator_TYP::RS_operator_TYP()
{
  cout << "Constructor for RS_TYP" << endl;
}

bool RS_operator_TYP::IsResamplingNeeded(params_TYP * params, ions_TYP * IONS)
{
  // Get geometric data:
  double L = params->mesh_params.Lx_max - params->mesh_params.Lx_min;
  double dx = params->mesh_params.dx;
  int Nx = params->mesh_params.Nx;

  // Get computational particle density:
  vec ncp_m = IONS->ncp_m.subvec(0,Nx-2);

  // Mean computational particle density
  double mean_ncp_m = sum(ncp_m)*dx/L;

  // Calculate metric:
  arma::vec dncp_m = ncp_m - mean_ncp_m/2;

  // Check if dncp_m becomes negative:
  for (int m = 0; m < Nx; m++)
  {
    if (dncp_m.at(m) < 0)
    {
      return true;
    }
  }

  return false;
}

void RS_operator_TYP::ApplyResampling_AllSpecies(params_TYP * params, mesh_TYP * mesh, vector<ions_TYP> * IONS)
{
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    for (int ss = 0; ss < IONS->size(); ss++)
    {
      if (IsResamplingNeeded(params,&IONS->at(ss)))
      {
        cout << "Apply resample, species " << ss << endl;
        // Assemble binary tree
        // Assemble vector of nodes:
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
