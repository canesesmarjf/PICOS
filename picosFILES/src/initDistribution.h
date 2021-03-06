#ifndef H_INIT_DIST
#define H_INIT_DIST

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

#include "types.h"
#include "mpi_main.h"

using namespace std;
using namespace arma;

class initDist_TYP
{
    // Cartesian  unitary vectors
    //arma::vec x = {1.0, 0.0, 0.0};
    //arma::vec y = {0.0, 1.0, 0.0};
    //arma::vec z = {0.0, 0.0, 1.0};

    //arma::vec b1; // Unitary vector along B field
    //arma::vec b2; // Unitary vector perpendicular to b1
    //arma::vec b3; // Unitary vector perpendicular to b1 and b2

    arma::vec x_profile;

    double target(const params_TYP * params, ionSpecies_TYP * IONS, double X, double V3, double V2, double V1);

    public:

    initDist_TYP(const params_TYP * params);

    void uniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS);

    void nonuniform_maxwellianDistribution(const params_TYP * params, ionSpecies_TYP * IONS);
};

#endif
