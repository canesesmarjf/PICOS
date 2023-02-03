#ifndef H_OUTPUT_H5
#define H_OUTPUT_H5

#include<iostream>
#include<armadillo>
#include<string>
#include<cmath>
#include "H5Cpp.h"
#include "types.h"

using namespace H5;
using namespace arma;
using namespace std;

class HDF_simple_TYP
{
public:
  void saveData(string fileName, params_TYP * params, fields_TYP * fields, vector<ions_TYP> * IONS);

  void H5_writeToFile(string fileName,string groupName,string datasetName,arma::vec * v ,int access);

  void H5_writeToFile(string fileName,string groupName,string datasetName,arma::ivec * v,int access);

  void H5_writeToFile(string fileName,string groupName,string datasetName,arma::mat * m ,int access);
};
#endif
