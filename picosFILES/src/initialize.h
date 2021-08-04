#ifndef H_INITIALIZE
#define H_INITIALIZE

// Intrinsic header files:
// =======================
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>

// Armadillo header:
// =================
#include <armadillo>

// User-defined headers:
// =====================
#include "types.h"
//#include "quietStart.h"
//#include "PIC.h"
#include "mpi_main.h"

using namespace std;
using namespace arma;

class init_TYP
{
	vector<string> split(const string& str, const string& delim);

	map<string, string> ReadAndloadInputFile(string *  inputFile);

	void initializeParticlesArrays(const params_TYP * params, fields_TYP * fields, ionSpecies_TYP * IONS);

	void initializeBulkVariablesArrays(const params_TYP * params, ionSpecies_TYP * IONS);

public:

	init_TYP(params_TYP * params, int argc, char* argv[]);

	void loadInputParameters(params_TYP * params, int argc, char* argv[]);

	void loadMeshGeometry(params_TYP * params, FS_TYP * FS);

  void loadIonParameters(params_TYP * params, vector<ionSpecies_TYP> * IONS);

  void loadPlasmaProfiles(params_TYP * params, vector<ionSpecies_TYP> * IONS);

	void setupIonsInitialCondition(const params_TYP * params, const CS_TYP * CS, fields_TYP * fields, vector<ionSpecies_TYP> * IONS);

	void initializeFieldsSizeAndValue(params_TYP * params, fields_TYP * fields);

	void initializeFields(params_TYP * params, fields_TYP * fields);

};

#endif
