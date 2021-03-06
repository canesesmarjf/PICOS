#ifndef H_FIELDS_SOLVER
#define H_FIELDS_SOLVER

#include <iostream>
#include <cmath>
#include <vector>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"

#include "mpi_main.h"

using namespace std;
using namespace arma;

class fields_solver_TYP
{

  	int NX_S; // Number of grid cells per subdomain, including 2 ghost cells.
  	int NX_T; // Number of grid cells in entire simulation domain, including 2 ghost cells.
  	int NX_R; // Number of grid cells in entire simulation domain, not including ghost cells.

  	// Electron density:
  	arma::vec ne;    // Current time tt
  	arma::vec ne_;   // tt - 1
  	arma::vec ne__;  // tt - 2
  	arma::vec ne___; // tt - 3

    // Electron temperature:
    arma::vec Te;

    // Electron pressure:
    arma::vec Pe;

    // Electron pressure gradient:
    arma::vec dPe;

  	// Electron density gradient:
  	//arma::vec dne;

    // Electric field:
    arma::vec EX_m;

  	// Grid cell increment
  	double dx;

  	// MPI functions:
  	void MPI_Allgathervec(const params_TYP * params, arma::vec * field);

  	void MPI_SendVec(const params_TYP * params, arma::vec * v);

    // Ghost cells:
    void fillGhosts(arma::vec * C);

    void fill4Ghosts(arma::vec * v);

  	// Smoothing:
  	void smooth(arma::vec * v, double as);

	public:

  	fields_solver_TYP(){};

  	fields_solver_TYP(const params_TYP * params, CS_TYP * CS);

  	//void advanceBField(const params_TYP * params, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

  	void advanceEfield(const params_TYP * params, fields_TYP * fields, CS_TYP * CS, vector<ionSpecies_TYP> * IONS, electrons_TYP * electrons);
};

#endif
