#ifndef H_RESAMPLINGOPER
#define H_RESAMPLINGOPER

// System libraries:
#include <iostream>
#include <cmath>
#include <vector>

// User-defined libraries:
#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "binaryTree.h"

// Parallelization libraries:
#include <omp.h>
#include "mpi_main.h"

using namespace std;
using namespace arma;

class RS_TYP
{
  private:
  double _L;  // Length of computational domain:
  double _dx; // spatial increment
  int _Nx;    // Number of grid cells
  int _numIONS; // Total number of ion species
  vec _mean_ncp_m; // Mean computational particle density. For each ion species

  public:
  // Default constructor:
  RS_TYP();

  // Overaloaded constructor:
  RS_TYP(params_TYP * params, vector<ions_TYP> * IONS, vector<binaryTree_TYP> * tree);

  // Methods:
  bool IsResamplingNeeded(params_TYP * params, vector<ions_TYP> * IONS, int ss);
  void ApplyResampling_AllSpecies(params_TYP * params, mesh_TYP * mesh, vector<ions_TYP> * IONS, vector<binaryTree_TYP> * tree);
  // void calculate_delta_profile(binaryTree_TYP * tree);
};


#endif
