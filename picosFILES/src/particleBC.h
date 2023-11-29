#ifndef H_PARTICLEBOUNDARYCONDITIONS
#define H_PARTICLEBOUNDARYCONDITIONS

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

#define ARMA_ALLOW_FAKE_GCC
#include "armadillo"
#include "types.h"
#include "mpi_main.h"
#include "omp.h"

using namespace std;
using namespace arma;

class particleBC_TYP
{

private:

    // void particleReinjection(int ii, params_TYP * params, const CS_TYP * CS, fields_TYP * fields, ions_TYP * IONS);
    void particleReinjection(int ii, params_TYP * params, const CS_TYP * CS, fields_TYP * fields, ions_TYP * IONS, int ss);

    void MPI_AllreduceDouble(params_TYP * params, double * v);

    template <typename vec_TYP> void MPI_OMP_AllreduceVec(params_TYP * params, vec_TYP * V, double * S);

    void checkBoundaryAndFlag(params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ions_TYP> * IONS);

    void getFluxesAcrossBoundaries(params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ions_TYP> * IONS);

    void calculateParticleWeight(params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ions_TYP> * IONS);

    void getParticleInjectionRates(params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ions_TYP> * IONS);

public:

    double N1_dot;
    double N2_dot;
    double N5_dot;
    double E1_dot;
    double E2_dot;
    double E5_dot;
    double DT;

    particleBC_TYP();

    void applyParticleReinjection(params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ions_TYP> * IONS);
};

#endif
