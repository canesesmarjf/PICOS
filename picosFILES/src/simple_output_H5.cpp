#include "simple_output_H5.h"
using namespace arma;

// ==========================================================================================
void HDF_simple_TYP::H5_writeToFile(string fileName,string groupName,string datasetName,arma::vec * v,int access)
{
  // Declare HDF5 objects:
    // ===========================================================================
    H5::H5File * file;
    H5::Group * group;
    H5::DataSpace * space;
    H5::DataSet * dataset;

    // Create file:
    // ===========================================================================
    if (access == 0) // New file
      file = new H5File(fileName,H5F_ACC_TRUNC);
    else if (access == 1) // Existing file
      file = new H5File(fileName,H5F_ACC_RDWR);

    // Create/open group:
    // ===========================================================================
    if (access == 0) // New file
      group = new H5::Group(file->createGroup(groupName));
    else if (access == 1) // Existing file
    {
      if(!file->exists(groupName))
        group = new H5::Group(file->createGroup(groupName));
      else
        group = new H5::Group(file->openGroup(groupName));
    }

    // Create dataspace:
    // ===========================================================================
    int rank = 1;
    hsize_t dims[1] = {v->n_elem};
    space = new H5::DataSpace(rank,dims);

    // Create dataset:
    // ===========================================================================
    dataset = new H5::DataSet(group->createDataSet(datasetName,PredType::NATIVE_DOUBLE,*space));


    // Write data:
    // ===========================================================================
    dataset->write(v->memptr(),PredType::NATIVE_DOUBLE);

    // Release memory:
    // ===========================================================================
    delete file;
    delete group;
    delete space;
    delete dataset;
}

// ==========================================================================================
void HDF_simple_TYP::H5_writeToFile(string fileName,string groupName,string datasetName,arma::ivec * v,int access)
{
  // Declare HDF5 objects:
    // ===========================================================================
    H5::H5File * file;
    H5::Group * group;
    H5::DataSpace * space;
    H5::DataSet * dataset;

    // Create file:
    // ===========================================================================
    if (access == 0) // New file
      file = new H5File(fileName,H5F_ACC_TRUNC);
    else if (access == 1) // Existing file
      file = new H5File(fileName,H5F_ACC_RDWR);

    // Create/open group:
    // ===========================================================================
    if (access == 0) // New file
      group = new H5::Group(file->createGroup(groupName));
    else if (access == 1) // Existing file
    {
      if(!file->exists(groupName))
        group = new H5::Group(file->createGroup(groupName));
      else
        group = new H5::Group(file->openGroup(groupName));
    }

    // Create dataspace:
    // ===========================================================================
    int rank = 1;
    hsize_t dims[1] = {v->n_elem};
    space = new H5::DataSpace(rank,dims);

    // Create dataset:
    // ===========================================================================
    H5::IntType datatype(H5T_STD_I64LE);
    // dataset = new H5::DataSet(group->createDataSet(datasetName,PredType::NATIVE_INT,*space));
    dataset = new H5::DataSet(group->createDataSet(datasetName,datatype,*space));

    // Write data:
    // ===========================================================================
    // dataset->write(v->memptr(),PredType::NATIVE_INT);
    dataset->write(v->memptr(),datatype);

    // Release memory:
    // ===========================================================================
    delete file;
    delete group;
    delete space;
    delete dataset;
}

// ==========================================================================================
void HDF_simple_TYP::H5_writeToFile(string fileName,string groupName,string datasetName,arma::mat * m,int access)
{
  // Declare HDF5 objects:
    // ===========================================================================
    H5::H5File * file;
    H5::Group * group;
    H5::DataSpace * space;
    H5::DataSet * dataset;

    // Create file:
    // ===========================================================================
    if (access == 0) // New file
      file = new H5File(fileName,H5F_ACC_TRUNC);
    else if (access == 1) // Existing file
      file = new H5File(fileName,H5F_ACC_RDWR);

    // Create/open group:
    // ===========================================================================
    if (access == 0) // New file
      group = new H5::Group(file->createGroup(groupName));
    else if (access == 1) // Existing file
    {
      if(!file->exists(groupName))
        group = new H5::Group(file->createGroup(groupName));
      else
        group = new H5::Group(file->openGroup(groupName));
    }

    // Create dataspace:
    // ===========================================================================
    int rank = 2;
    hsize_t dims[2] = {m->n_cols,m->n_rows};
    space = new H5::DataSpace(rank,dims);

    // Create dataset:
    // ===========================================================================
    dataset = new H5::DataSet(group->createDataSet(datasetName,PredType::NATIVE_DOUBLE,*space));

    // Write data:
    // ===========================================================================
    dataset->write(m->memptr(),PredType::NATIVE_DOUBLE);

    // Release memory:
    // ===========================================================================
    delete file;
    delete group;
    delete space;
    delete dataset;
}

// ==========================================================================================
void HDF_simple_TYP::saveData(string fileName, params_TYP * params, fields_TYP * fields, vector<ions_TYP> * IONS)
{
  // Location where data is to be stored:
  string targetDir = "outputFiles/HDF5_simple/";

  // fileName with full path:
  string pathName = targetDir + fileName;

  for(int s = 0; s < IONS->size(); s++)
  {
    stringstream so;
    so << s;
    string groupName = "/ions_" + so.str();
    int access_flag;

    if (s == 0)
      access_flag = 0;
    else
      access_flag = 1;

    this->H5_writeToFile(pathName,groupName,"x_p",&IONS->at(s).x_p,access_flag);
    this->H5_writeToFile(pathName,groupName,"v_p",&IONS->at(s).v_p,1);
    this->H5_writeToFile(pathName,groupName,"a_p",&IONS->at(s).a_p,1);

    arma::vec u = conv_to<vec>::from(IONS->at(s).mn);
    this->H5_writeToFile(pathName,groupName,"m_n",&u ,1);
    // this->H5_writeToFile(pathName,groupName,"m_n",&IONS->at(s).mn,1);

    this->H5_writeToFile(pathName,groupName,"n_m",&IONS->at(s).n_m,1);
    arma::vec upar_m = IONS->at(s).nv_m/IONS->at(s).n_m;
    this->H5_writeToFile(pathName,groupName,"upar_m",&upar_m,1);
    this->H5_writeToFile(pathName,groupName,"Tpar_m",&IONS->at(s).Tpar_m,1);
    this->H5_writeToFile(pathName,groupName,"Tper_m",&IONS->at(s).Tper_m,1);
    this->H5_writeToFile(pathName,groupName,"ncp_m",&IONS->at(s).ncp_m,1);
  }

  // Write field data into H5 file:
  string groupName = "/fields";
  this->H5_writeToFile(pathName,groupName,"x_m"   ,&fields->x_mg  ,1);
  this->H5_writeToFile(pathName,groupName,"Ex_m"  ,&fields->Ex_m  ,1);
  this->H5_writeToFile(pathName,groupName,"Bx_m"  ,&fields->Bx_m  ,1);
  this->H5_writeToFile(pathName,groupName,"dBx_m" ,&fields->dBx_m ,1);
  this->H5_writeToFile(pathName,groupName,"ddBx_m",&fields->ddBx_m,1);
}
