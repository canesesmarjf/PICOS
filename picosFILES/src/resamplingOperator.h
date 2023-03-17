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

// What does the construct need?
// It needs to method that takes in IONS, loops over all species and first calculates the metric
// Whether or not resampling is needed.
// IF resampling is needed, then proceed to take in IONS, params and/or mesh and create a binary tree.
// we could create a binary tree inmeddiateby after the initial placement of ions, This will create all the nodes we will ever require
// Each ion species requires its own resample_operatar, so in main we can create vector <resample_TYP> resample
// We also ahve:
//     if (params.SW.resample == 1)
//    {
//
//    }

  // node binaryTree
  // std::vector<node*> node_vector;

};


#endif
