#ifndef H_RS_OPERATOR
#define H_RS_OPERATOR

// System libraries:
#include <iostream>
#include <cmath>
#include <vector>

// User-defined libraries:
#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "particle_tree.h"

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
  RS_TYP(params_TYP * params, CS_TYP * CS, vector<ions_TYP> * IONS, vector<particle_tree_TYP> * tree);

  // Methods:
  bool IsResamplingNeeded(params_TYP * params, vector<ions_TYP> * IONS, int ss);
  void ApplyResampling_AllSpecies(params_TYP * params, mesh_TYP * mesh, vector<ions_TYP> * IONS, vector<particle_tree_TYP> * tree);

};

#endif
