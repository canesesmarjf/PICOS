#include "initialize.h"
// #include "initDistribution.h"

void debug(int i,params_TYP * params)
{
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "debug " << i << endl;
  }
}

void debug(int i)
{
  cout << "debug " << i << endl;
}

// Function to split strings:
// =============================================================================
vector<string> init_TYP::split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);

        if (pos == string::npos)
        {
            pos = str.length();
        }

        string token = str.substr(prev, pos-prev);

        if (!token.empty())
        {
            tokens.push_back(token);
        }

        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return tokens;
}

// Function to read and load data from inputfile.input:
// =============================================================================
map<string,string> init_TYP::readTextFile(string * inputFile)
{
    // Create stream object:
    // =====================
    fstream reader;

    // Create map object:
    // ==================
    std::map<string,string> readMap;

    // Open input file using reader object:
    // ====================================
    reader.open(inputFile->data(),ifstream::in);

    // Handle error:
    // =============
    if (!reader)
    {
        MPI_Barrier(MPI_COMM_WORLD);

    	cerr << "PRO++ ERROR: The input file couldn't be opened." << endl;
    	MPI_Abort(MPI_COMM_WORLD, -101);
    }
    debug(10);

    // Parse through file:
    // ===================
    string lineContent;
    vector<string> keyValuePair;
    while ( reader.good() )
    {
        // Read entire line:
        getline(reader,lineContent);

        // Search for comment symbol:
        size_t commentCharPos = lineContent.find("//",0);

        // Check for comment symbol:
        if (commentCharPos == 0 || lineContent.empty())
        {
            // Skip line
        }
        else
        {
            // Get value pair:
            keyValuePair = split(lineContent," ");

            // Update map:
            readMap[ keyValuePair[0] ] = keyValuePair[1];
        }
    }
    debug(11);

    // Close stream object:
    // ===================
    reader.close();

    // Return map:
    // ==========
    return readMap;
}

// Constructor:
// =============================================================================
init_TYP::init_TYP(params_TYP * params, int argc, char* argv[])
{
  // Get RANK and SIZE of nodes within COMM_WORLD:
  // =============================================
  MPI_Comm_size(MPI_COMM_WORLD, &params->mpi.NUMBER_MPI_DOMAINS);
  MPI_Comm_rank(MPI_COMM_WORLD, &params->mpi.MPI_DOMAIN_NUMBER);

  // params->mpi.NUMBER_MPI_DOMAINS = 4;
  // params->mpi.MPI_DOMAIN_NUMBER = 0;

  // MPI_Barrier(MPI_COMM_WORLD);
  cout << "params->mpi.NUMBER_MPI_DOMAINS = " << params->mpi.MPI_DOMAIN_NUMBER << endl;
  // MPI_Barrier(MPI_COMM_WORLD);

  // Error codes:
  // ============
  params->errorCodes[-100] = "Odd number of MPI processes";
  params->errorCodes[-102] = "MPI's Cartesian topology could not be created";
  params->errorCodes[-103] = "Grid size violates assumptions of hybrid model for the plasma -- DX smaller than the electron skind depth can not be resolved";
  params->errorCodes[-106] = "Inconsistency in iniital ion's velocity distribution function";

  // Program information:
  // ===========================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
    cout << "* PICOS++, a 1D-2V GC hybrid PIC code for Open plasma Systems           *" << endl;
    cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
    cout << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Arguments and paths to main function:
  // =====================================
  params->PATH = argv[2];
	params->argc = argc;
	params->argv = argv;

  // Check number of MPI domains:
  // ============================
	if( fmod( (double)params->mpi.NUMBER_MPI_DOMAINS, 2.0 ) > 0.0 )
  {
    MPI_Barrier(MPI_COMM_WORLD);

		if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
			cerr << "PICOS++ ERROR: The number of MPI processes must be an even number." << endl;
		}

		//MPI_Abort(MPI_COMM_WORLD,-100);
	}

  // Stream date when simulation is started:
  // =======================================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      time_t current_time = std::time(NULL);
      cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
      cout << "STARTING " << params->argv[1] << " SIMULATION ON: " << std::ctime(&current_time) << endl;
      cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
  }
}

// Populate params with data from input file:
// =============================================================================
void init_TYP::read_inputFile(params_TYP * params)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
      cout << "READING INPUT FILE ..." << endl;
  }

    // Get name of path to input file:
    // ===============================
	string name;
	if(params->argc > 3)
    {
		string argv(params->argv[3]);
		name = "inputFiles/input_file_" + argv + ".input";
		params->PATH += "/" + argv;
	}
    else
    {
		name = "inputFiles/input_file.input";
		params->PATH += "/";
	}

  // Read input file and assemble map:
  // ================================
	std::map<string,string> parametersStringMap;
	parametersStringMap = readTextFile(&name);

    // Create HDF5 folders if they don't exist:
    // ========================================
	if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
  		string mkdir_outputs_dir = "mkdir " + params->PATH;
  		const char * sys = mkdir_outputs_dir.c_str();
  		int rsys = system(sys);

      // HDF5 directory:
  		string mkdir_outputs_dir_HDF5 = mkdir_outputs_dir + "/HDF5";
  		sys = mkdir_outputs_dir_HDF5.c_str();
  		rsys = system(sys);

      // HDF5_simple directory:
      mkdir_outputs_dir_HDF5 = mkdir_outputs_dir + "/HDF5_simple";
      sys = mkdir_outputs_dir_HDF5.c_str();
      rsys = system(sys);
	   }

  // Populate "params" with data from input file:
  // ============================================
  params->mpi.MPIS_FIELDS  = stoi( parametersStringMap["mpisForFields"] );

  params->numberOfParticleSpecies = stoi( parametersStringMap["numberOfParticleSpecies"] );
  params->numberOfTracerSpecies   = stoi( parametersStringMap["numberOfTracerSpecies"] );
  params->advanceParticleMethod   = stoi( parametersStringMap["advanceParticleMethod"] );

  // Characteristic values:
  // -------------------------------------------------------------------------
  params->CV.ne   = stod( parametersStringMap["CV_ne"] );
  params->CV.Te   = stod( parametersStringMap["CV_Te"] )*F_E/F_KB; // [K]
  params->CV.B    = stod( parametersStringMap["CV_B"] );
  params->CV.Tpar = stod( parametersStringMap["CV_Tpar"] )*F_E/F_KB; // [K]
  params->CV.Tper = stod( parametersStringMap["CV_Tper"] )*F_E/F_KB; // [K]

  // Simulation time:
  // -------------------------------------------------------------------------
  params->DTc            = stod( parametersStringMap["DTc"] );
  params->simulationTime = std::stod( parametersStringMap["simulationTime"] );

  // Switches:
  // -------------------------------------------------------------------------
  params->SW.resample      = stoi( parametersStringMap["SW_resample"] );
  params->SW.EfieldSolve   = stoi( parametersStringMap["SW_EfieldSolve"] );
  params->SW.BfieldSolve   = stoi( parametersStringMap["SW_BfieldSolve"] );
  params->SW.Collisions    = stoi( parametersStringMap["SW_Collisions"] );
  params->SW.RFheating     = stoi( parametersStringMap["SW_RFheating"] );
  params->SW.advancePos    = stoi( parametersStringMap["SW_advancePos"] );
  params->SW.linearSolve   = stoi( parametersStringMap["SW_linearSolve"] );

  // Magnetic field initial conditions:
  // -------------------------------------------------------------------------
  params->fields_IC.fileName  = parametersStringMap["IC_fields_fileName"];
  params->fields_IC.Ex_offset = stod( parametersStringMap["IC_Ex_offset"] );
  params->fields_IC.Ex_scale  = stod( parametersStringMap["IC_Ex_scale"] );
  params->fields_IC.Bx_offset = stod( parametersStringMap["IC_Bx_offset"] );
  params->fields_IC.Bx_scale  = stod( parametersStringMap["IC_Bx_scale"] );

  // Electron initial conditions:
  // -------------------------------------------------------------------------
  params->electrons_IC.fileName  = parametersStringMap["IC_electrons_fileName"];
  params->electrons_IC.n_offset = stod(parametersStringMap["IC_n_offset"]);
  params->electrons_IC.n_scale  = stod(parametersStringMap["IC_n_scale"]);
  params->electrons_IC.T_offset = stod(parametersStringMap["IC_T_offset"])*F_E/F_KB; // [K]
  params->electrons_IC.T_scale  = stod(parametersStringMap["IC_T_scale"]);

  // Mesh:
  // -------------------------------------------------------------------------
  params->mesh_params.dx_norm = stod(parametersStringMap["M_dx_norm"]);
  params->mesh_params.x0      = stod(parametersStringMap["M_x0"]);
  params->mesh_params.r0_min  = stod(parametersStringMap["M_r0_min"]);
  params->mesh_params.r0_max  = stod(parametersStringMap["M_r0_max"]);
  params->mesh_params.Lx_min  = stod(parametersStringMap["M_Lx_min"]);
  params->mesh_params.Lx_max  = stod(parametersStringMap["M_Lx_max"]);
  params->mesh_params.getA0();

  // RF parameters
  // -------------------------------------------------------------------------
  params->RF.Prf        = stod( parametersStringMap["RF_Prf"] );
  params->RF.n_harmonic = stoi( parametersStringMap["RF_n_harmonic"] );
  params->RF.freq       = stod( parametersStringMap["RF_freq"]);
  params->RF.x1         = stod( parametersStringMap["RF_x1"]  );
  params->RF.x2         = stod( parametersStringMap["RF_x2"]  );
  params->RF.t_ON       = stoi( parametersStringMap["RF_t_ON"]  );
  params->RF.t_OFF      = stoi( parametersStringMap["RF_t_OFF"]  );
  params->RF.kpar       = stod( parametersStringMap["RF_kpar"]);
  params->RF.kper       = stod( parametersStringMap["RF_kper"]);
  params->RF.handedness = stoi( parametersStringMap["RF_handedness"]);
  params->RF.Prf_NS     = stoi( parametersStringMap["RF_Prf_NS"] );
  params->RF.Prf_fileName = parametersStringMap["RF_Prf_fileName"];

  // Output variables:
  // -------------------------------------------------------------------------
  params->outputCadence           = stod( parametersStringMap["outputCadence"] );
  string nonparsed_variables_list = parametersStringMap["outputs_variables"].substr(1, parametersStringMap["outputs_variables"].length() - 2);
  params->outputs_variables       = split(nonparsed_variables_list,",");

  // Data smoothing:
  // -------------------------------------------------------------------------
  params->smoothingParameter        = stod( parametersStringMap["smoothingParameter"] );
  params->filtersPerIterationFields = stoi( parametersStringMap["filtersPerIterationFields"] );
  params->filtersPerIterationIons   = stoi( parametersStringMap["filtersPerIterationIons"] );

  MPI_Barrier(MPI_COMM_WORLD);

  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "READING INPUT FILE COMPLETED" << endl;
      cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
  }
}

// Read and populate ion parameters input file:
// =============================================================================
void init_TYP::read_ionsPropertiesFile(params_TYP * params)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * LOADING ION PARAMETERS * * * * * * * * * * * * * * * * * *\n";
  	cout << "+ Number of ion species: " << params->numberOfParticleSpecies << endl;
  	cout << "+ Number of tracer species: " << params->numberOfTracerSpecies << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Assemble path to "ion_properties.ion":
  // ======================================
  string name;
  if(params->argc > 3)
  {
  	string argv(params->argv[3]);
  	name = "inputFiles/ions_properties_" + argv + ".ion";
  }
  else
  {
  	name = "inputFiles/ions_properties.ion";
  }

  // Read data from "ion_properties.ion" into "parametersMap":
  // =========================================================
  std::map<string,string> parametersMap;
  debug(1,params);
  cout << "ion propertuies reading start, rank = " << params->mpi.MPI_DOMAIN_NUMBER << endl;
  parametersMap = readTextFile(&name);
  cout << "ion propertuies reading end, rank = " << params->mpi.MPI_DOMAIN_NUMBER << endl;
  cout << "IC_ions_fileName = " << parametersMap["IC_ions_fileName"] << endl;
  debug(2,params);

  MPI_Barrier(MPI_COMM_WORLD);

  // Determine the total number of ION species:
  // =========================================
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Get ion parameters:
  // ====================
  params->ions_params.resize(totalNumSpecies);
  for(int ss = 0; ss<totalNumSpecies; ss++)
  {
    string name;
    stringstream kk;
    kk << ss + 1;

    // General:
    name = "SPECIES" + kk.str();
    params->ions_params[ss].SPECIES = stoi(parametersMap[name]);
    //IONS->at(ss).SPECIES         = stoi(parametersMap[name]);

    name = "N_CP" + kk.str();
    params->ions_params[ss].N_CP_IN_SIM = stoi(parametersMap[name]);
    //IONS->at(ss).N_CP_IN_SIM         = stoul(parametersMap[name]);

    name = "pct_N_CP_Output" + kk.str();
    params->ions_params[ss].pct_N_CP_Output = stoi(parametersMap[name]);

    name = "Z" + kk.str();
    params->ions_params[ss].Z = stoi(parametersMap[name]);
    //IONS->at(ss).Z         = stoi(parametersMap[name]);

    name = "M" + kk.str();
    params->ions_params[ss].M = F_U*stod(parametersMap[name]);
    //IONS->at(ss).M         = stoi(parametersMap[name]);
  }

  // Get ion IC and BC:
  // ======================================================
  params->ions_IC.resize(totalNumSpecies);
  params->ions_BC.resize(totalNumSpecies);
  for(int ss = 0; ss<totalNumSpecies; ss++)
  {
    params->ions_IC[ss].fileName = parametersMap["IC_ions_fileName"];
    params->ions_IC[ss].CP_fileName = parametersMap["IC_CP_pdf_fileName"];

    string name;
    stringstream kk;
    kk << ss + 1;

    // IC:
    name = "IC_Tpar_offset" + kk.str();
    params->ions_IC[ss].Tpar_offset = stod(parametersMap[name])*F_E/F_KB; // [K]

    name = "IC_Tpar_scale" + kk.str();
    params->ions_IC[ss].Tpar_scale = stod(parametersMap[name]);

    name = "IC_Tper_offset" + kk.str();
    params->ions_IC[ss].Tper_offset = stod(parametersMap[name])*F_E/F_KB; // [K]

    name = "IC_Tper_scale" + kk.str();
    params->ions_IC[ss].Tper_scale = stod(parametersMap[name]);

    name = "IC_upar_offset" + kk.str();
    params->ions_IC[ss].upar_offset = stod(parametersMap[name]);

    name = "IC_upar_scale" + kk.str();
    params->ions_IC[ss].upar_scale = stod(parametersMap[name]);

    name = "IC_densityFraction" + kk.str();
    params->ions_IC[ss].densityFraction = stod(parametersMap[name]);

    name = "IC_mean_ai" + kk.str();
    params->ions_IC[ss].mean_ai = stod(parametersMap[name]);

    // BC:
    name = "BC_type" + kk.str();
    params->ions_BC[ss].type = stod(parametersMap[name]);

    name = "BC_T" + kk.str();
    params->ions_BC[ss].T = stod(parametersMap[name])*F_E/F_KB; // [K]

    name = "BC_E" + kk.str();
    params->ions_BC[ss].E = stod(parametersMap[name])*F_E/F_KB; // [K]

    name = "BC_eta" + kk.str();
    params->ions_BC[ss].eta = stod(parametersMap[name]);

    name = "BC_mean_x" + kk.str();
    params->ions_BC[ss].mean_x = stod(parametersMap[name]);

    name = "BC_sigma_x" + kk.str();
    params->ions_BC[ss].sigma_x = stod(parametersMap[name]);

    name = "BC_G" + kk.str();
    params->ions_BC[ss].G = stod(parametersMap[name]); // [particles/sec]

    name = "BC_G_fileName" + kk.str();
    params->ions_BC[ss].G_fileName = parametersMap[name];
  }

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * * ION PARAMETERS LOADED * * * * * * * * * * * * * * * * * *\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void init_TYP::create_mesh(params_TYP * params, mesh_TYP * mesh)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * CREATING MESH * * * * * * * * * * * * * * * * * * * * *\n";
  }

  // Prior to calculating mesh quantities, we need to compute the characteristic skin depth using
  // the characteristic value for the density, mass and charge.
  // The characteristic density is given in the inputfile.
  // For the other terms, we use the charge and mass of the majority ion.

  // Calculating Nx and dx:
  // ==================================================================
  // Calculate characteristic skin depth:
  params->getCharacteristicIonSkinDepth();
  double ionSkinDepth = params->mesh_params.ionSkinDepth;

  // Calculate Nx and Nx_PER_MPI:
  // This step makes use of the ion skin depth and x_norm from the input to
  // Produce Nx and dx which satify this requirement.
  params->get_Nx_dx(ionSkinDepth);
  int Nx    = params->mesh_params.Nx;
  double dx = params->mesh_params.dx;

  // Calculate mesh-defined cell center grid xm and xmg:
  // ==================================================================
  arma::vec xm  = arma::zeros(Nx);
  arma::vec xmg = arma::zeros(Nx+2);
  double Lx_min = params->mesh_params.Lx_min;

  // No ghost cells:
  for (int i = 0; i < Nx; i++)
  {
    xm.at(i) = Lx_min + (i + 0.5)*dx;
  }

  // With ghost cells:
  for (int i = 0; i < (Nx+2); i++)
  {
    xmg.at(i) = Lx_min + (i - 0.5)*dx;
  }

  // Assign values:
  // ==================================================================
  mesh->Nx_PER_MPI = params->mesh_params.Nx_PER_MPI;
  mesh->Nx  = Nx;
  mesh->dx  = dx;
  mesh->xm  = xm;
  mesh->xmg = xmg;

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * * MESH CREATED * * * * * * * * * * * * * * * * * * * * *\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

// Read "ion properties" file and populate IC object with profiles:
void init_TYP::read_IC_profiles(params_TYP * params, mesh_TYP * mesh, IC_TYP * IC)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * *  LOADING INITIAL CONDITION PROFILES  * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Directory where profile data is located:
  string directory = "inputFiles/";

  // Read and scale electron profiles:
  // =====================================================
  // Assemble full path to HDF5 file with profiles:
  string fileName  = params->electrons_IC.fileName;
  string fullPath = directory + fileName;

  // temp containers:
  arma::vec y;
  double offset;
  double scale;

  cout << "here I am" << endl;
  cout << "fullPath = " << fullPath << endl;

  // Load and assign data:
  y.load(arma::hdf5_name(fullPath,"x"));
  IC->electrons.x = y; // [m]

  y.load(arma::hdf5_name(fullPath,"n"));
  offset = params->electrons_IC.n_offset;
  scale  = params->electrons_IC.n_scale;
  IC->electrons.n = offset + scale*y;

  y.load(arma::hdf5_name(fullPath,"T"));
  offset = params->electrons_IC.T_offset;
  scale  = params->electrons_IC.T_scale;
  IC->electrons.T = offset + scale*y*F_E/F_KB; // [K]

  // Read and scale fields profiles:
  // ==================================================
  fileName  = params->fields_IC.fileName;
  fullPath = directory + fileName;

  // Load and assign data:
  y.load(arma::hdf5_name(fullPath,"x"));
  IC->fields.x = y;

  y.load(arma::hdf5_name(fullPath,"Bx"));
  offset = params->fields_IC.Bx_offset;
  scale  = params->fields_IC.Bx_scale;
  IC->fields.Bx = offset + scale*y;

  y.load(arma::hdf5_name(fullPath,"Ex"));
  offset = params->fields_IC.Ex_offset;
  scale  = params->fields_IC.Ex_scale;
  IC->fields.Ex = offset + scale*y;

  // Compute 1st and 2nd derivatives of Bx:
  double dx = IC->fields.x(2) - IC->fields.x(1);
  int nx    = IC->fields.Bx.n_elem;
  arma::vec B(nx,1);
  arma::vec dB(nx,1);
  arma::vec ddB(nx,1);

  B = IC->fields.Bx;
  dB.subvec(1,nx-2) = (B.subvec(2,nx-1) - B.subvec(0,nx-3))/(2*dx);
  dB(0)    = dB(1);
  dB(nx-1) = dB(nx-2);
  IC->fields.dBx = dB;

  dB = IC->fields.dBx;
  ddB.subvec(1,nx-2) = (dB.subvec(2,nx-1) - dB.subvec(0,nx-3))/(2*dx);
  ddB(0)    = ddB(1);
  ddB(nx-1) = ddB(nx-2);
  IC->fields.ddBx = ddB;

  // Read and scale ions profiles:
  // ==================================================
  // Determine the total number of ION species:
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Allocate memory:
  IC->ions.resize(totalNumSpecies);

  // Assemble full path to HDF5 file with ion profiles:
  fileName  = params->ions_IC.at(0).fileName;
  fullPath = directory + fileName;

  // Load ion profile data from H5 file into containers:
  for (int ss = 0; ss < totalNumSpecies; ss++)
  {
    string dataset;
    string group;
    stringstream kk;
    kk << ss + 1;
    group = "ions_" + kk.str();

    dataset = group + "/x";
    y.load(arma::hdf5_name(fullPath,dataset));
    IC->ions.at(ss).x = y;

    scale  = params->ions_IC.at(ss).densityFraction;
    IC->ions.at(ss).n = scale*IC->electrons.n;

    dataset = group + "/Tpar";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).Tpar_offset;
    scale  = params->ions_IC.at(ss).Tpar_scale;
    IC->ions.at(ss).Tpar = offset + scale*y*F_E/F_KB; // [K]

    dataset = group + "/Tper";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).Tper_offset;
    scale  = params->ions_IC.at(ss).Tper_scale;
    IC->ions.at(ss).Tper = offset + scale*y*F_E/F_KB; // [K]

    dataset = group + "/upar";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).upar_offset;
    scale  = params->ions_IC.at(ss).upar_scale;
    IC->ions.at(ss).upar = offset + scale*y;
  }

  // Read computational particle profiles:
  // ==================================================
  // Assemble full path to HDF5 file with computational particle (CP) profile:
  fileName  = params->ions_IC.at(0).CP_fileName;
  fullPath = directory + fileName;

  // Load CP profile data from H5 file into containers:
  for (int ss = 0; ss < totalNumSpecies; ss++)
  {
    y.load(arma::hdf5_name(fullPath,"n_pdf"));
    offset = 1;
    scale  = 0;
    IC->ions.at(ss).ncp_pdf = offset + scale*y;
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void interp_scalar(arma::vec * xv, arma::vec * yv, double * x, double * y)
{
  // Convert input doubles into 1D arma vectors:
  arma::vec xq(1); xq = {*x};
  arma::vec yq(1);

  // Interpolate:
  arma::interp1(*xv,*yv,xq,yq);

  // Return value:
  *y = arma::as_scalar(yq);
}

void init_TYP::interpolate_IC_profiles(params_TYP * params, mesh_TYP * mesh, IC_TYP * IC)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * INTERPOLATING IC PROFILES * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Query x-vector with ghost cells:
  arma::vec xq = mesh->xmg;

  // Interpolate electron profiles:
  // ==========================================================================================
  // x-vector for input profiles:
  arma::vec xv = IC->electrons.x;

  // Assign x-vector to output profiles:
  IC->electrons.x_mg = xq;

  // Interpolate profiles at cell-center grid:
  arma::interp1(xv,IC->electrons.n,xq,IC->electrons.n_mg);
  arma::interp1(xv,IC->electrons.T,xq,IC->electrons.T_mg);

  // Interpolate field profiles:
  // ==========================================================================================
  // x-vector for input profiles:
  xv = IC->fields.x;

  // Assign x-vector to output profiles:
  IC->fields.x_mg = xq;

  // Interpolate fields at cell-center grid:
  arma::interp1(xv,IC->fields.Bx  ,xq,IC->fields.Bx_mg);
  arma::interp1(xv,IC->fields.dBx ,xq,IC->fields.dBx_mg);
  arma::interp1(xv,IC->fields.ddBx,xq,IC->fields.ddBx_mg);
  arma::interp1(xv,IC->fields.Ex  ,xq,IC->fields.Ex_mg);

  // Get reference magnetic field:
  double x0 = params->mesh_params.x0;
  interp_scalar(&IC->fields.x,&IC->fields.Bx,&x0,&params->mesh_params.B0);

  // Interpolate ion profiles:
  // ==========================================================================================
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);
  for (int ss = 0; ss < totalNumSpecies; ss ++)
  {
    // x-vector for input profiles:
    xv = IC->ions.at(ss).x;

    // Assign x-vector to output profiles:
    IC->ions.at(ss).x_mg = xq;

    // Interpolate profiles at cell-center grid:
    arma::interp1(xv,IC->ions.at(ss).n      ,xq,IC->ions.at(ss).n_mg);
    arma::interp1(xv,IC->ions.at(ss).Tpar   ,xq,IC->ions.at(ss).Tpar_mg);
    arma::interp1(xv,IC->ions.at(ss).Tper   ,xq,IC->ions.at(ss).Tper_mg);
    arma::interp1(xv,IC->ions.at(ss).upar   ,xq,IC->ions.at(ss).upar_mg);
    arma::interp1(xv,IC->ions.at(ss).ncp_pdf,xq,IC->ions.at(ss).ncp_shape_mg);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * IC PROFILES INTERPOLATED  * * * * * * * * * * * * * * * * * *" << endl;
  }
}

arma::vec pdf(arma::vec y, double x)
{
  arma::vec value = y/sum(y*x);
  return value;
}

void init_TYP::calculate_IC_particleWeight(params_TYP * params, IC_TYP * IC, vector<ions_TYP> * IONS)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * CALCULATING PARTICLE WEIGHT INITIAL PROFILE * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Define total number of ions species:
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Number of particle MPIS:
  int Nk = params->mpi.MPIS_PARTICLES;

  // Allocate space for IONS vector:
  IONS->resize(totalNumSpecies);

  // Mesh parameters:
  int Nx    = params->mesh_params.Nx;
  double dx = params->mesh_params.dx;
  double A0 = params->mesh_params.A0;
  double B0 = params->mesh_params.B0;

  // Allocate memory for temporary variables (no ghost cells):
  arma::vec xm(Nx,fill::zeros);
  arma::vec nm_cp_in_SIM(Nx,fill::zeros);
  arma::vec nm_cp_MPI(Nx,fill::zeros);
  arma::vec nm_cp_shape(Nx,fill::zeros);
  arma::vec nm(Nx,fill::zeros);
  arma::vec Am(Nx,fill::zeros);
  arma::vec am(Nx,fill::zeros);

  // Loop over all ion species:
  for (int ss = 0; ss < totalNumSpecies; ss ++)
  {
    // Calculate N_CP and nm_cp:
    // ===================================================
    // The value of N_CP has to be both divisable by dx in order to produce an integer
    // number of computational particles per cell and it has to be divisable by MPIS_PARTICLES
    // so that each MPI gets a integer number of particles:

    // Initial N_CP value provided by user:
    int N_CP_IN_SIM = params->ions_params.at(ss).N_CP_IN_SIM;
    int N_CP_MPI    = N_CP_IN_SIM/Nk;

    cout << "intial N_CP_IN_SIM = " << N_CP_IN_SIM << endl;
    cout << "initial N_CP_MPI = "    << N_CP_MPI    << endl;

    // Initial estimate of nm_cp_in_SIM:
    nm_cp_shape  = IC->ions.at(ss).ncp_shape_mg.subvec(1,Nx);
    nm_cp_in_SIM = N_CP_IN_SIM*pdf(nm_cp_shape,dx);

    // Make exactly divisable by dx and multiple of Nk:
    nm_cp_in_SIM = round(nm_cp_in_SIM*dx/Nk)*Nk/dx;

    // Calculate new N_CP_IN_SIM:
    N_CP_IN_SIM = sum(nm_cp_in_SIM)*dx;

    // Calculate new N_CP_MPI:
    N_CP_MPI = N_CP_IN_SIM/Nk;

    // Produce the final  computational particle density profile:
    nm_cp_in_SIM = N_CP_IN_SIM*pdf(nm_cp_in_SIM,dx);
    nm_cp_MPI    = N_CP_MPI*pdf(nm_cp_in_SIM/Nk,dx);

    cout << "N_CP_IN_SIM = " << N_CP_IN_SIM << endl;
    cout << "N_CP_MPI = "    << N_CP_MPI    << endl;

    cout << "sum(nm_cp_in_SIM)*dx = " << sum(nm_cp_in_SIM)*dx << endl;
    cout << "sum(nm_cp_MPI)*dx    = " << sum(nm_cp_MPI)*dx << endl;

    // Assign value:
    IONS->at(ss).N_CP_IN_SIM = N_CP_IN_SIM;
    IONS->at(ss).N_CP_MPI    = N_CP_MPI;
    IC->ions.at(ss).ncp_in_SIM_mg.zeros(Nx+2);
    IC->ions.at(ss).ncp_MPI_mg.zeros(Nx+2);
    IC->ions.at(ss).ncp_in_SIM_mg.subvec(1,Nx) = nm_cp_in_SIM;
    IC->ions.at(ss).ncp_MPI_mg.subvec(1,Nx)    = nm_cp_MPI;

    // Calculate N_R:
    // =====================================================
    Am = A0*B0/IC->fields.Bx_mg.subvec(1,Nx);
    nm = IC->ions.at(ss).n_mg.subvec(1,Nx);
    double N_R = sum(nm%Am)*dx;

    // Assign value:
    IONS->at(ss).N_R = N_R;

    // Calculate N_SP:
    // =====================================================
    double mean_ai = params->ions_IC.at(ss).mean_ai;
    double N_SP = mean_ai*N_CP_IN_SIM;

    // Assign value:
    IONS->at(ss).N_SP = N_SP;

    // Calculate normalization constant K for distribution function:
    // =====================================================
    // This value remains constant over entire simulations
    IONS->at(ss).K = N_R/N_SP;

    // Calculate particle weight IC profile (only for t = 0):
    // =====================================================
    am = (N_SP/N_R)*(nm%Am/nm_cp_in_SIM);

    // Assign value:
    IC->ions.at(ss).a_mg.zeros(Nx+2);
    IC->ions.at(ss).a_mg.subvec(1,Nx) = am;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * * * * * * * COMPLETE * * * * * * * * * * * * * * * * * *" << endl;
  }
}

void init_TYP::initialize_fields(params_TYP * params, IC_TYP * IC, fields_TYP * fields)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * * INITIALIZING ELECTROMAGNETIC FIELDS * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Size of vectors wth ghost cells:
  int Nxg = params->mesh_params.Nx + 2;

  // Allocate memory to field quantities, including ghost cells:
  fields->x_mg.zeros(Nxg);
  fields->Bx_m.zeros(Nxg);
  fields->dBx_m.zeros(Nxg);
  fields->ddBx_m.zeros(Nxg);
  fields->Ex_m.zeros(Nxg);

  // Assign value to profiles:
  fields->x_mg   = IC->fields.x_mg;
  fields->Bx_m   = IC->fields.Bx_mg;
  fields->dBx_m  = IC->fields.dBx_mg;
  fields->ddBx_m = IC->fields.ddBx_mg;
  fields->Ex_m   = IC->fields.Ex_mg;

  // Calculate cross sectional area vector based on latest B field:
  double B0 = params->mesh_params.B0;
  double A0 = params->mesh_params.A0;
  fields->getAm(A0,B0);

  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * * ELECTROMAGNETIC FIELDS INITIALIZED  * * * * * * * * * * * * * * * * * *" << endl;
  }

}

void init_TYP::initialize_electrons(params_TYP * params, IC_TYP * IC, electrons_TYP * electrons)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * * INITIALIZING ELECTRON FLUID * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Size of vectors wth ghost cells:
  int Nxg = params->mesh_params.Nx + 2;

  // Allocate memory to profile quantities, including ghost cells:
  electrons->Te_m.zeros(Nxg);

  // Assign value to profiles:
  electrons->Te_m = IC->electrons.T_mg;

  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * ELECTRON FLUID INITIALIZED  * * * * * * * * * * * * * * * * * *" << endl;
  }
}

void init_TYP::initialize_ions(params_TYP * params, IC_TYP * IC, mesh_TYP * mesh, vector<ions_TYP> * IONS)
{
  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * * SETTING UP IONS INITIAL CONDITION * * * * * * * * * * * * * * * * * *" << endl;
  }

  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Define ion parameters:
  // ===========================================================
  for (int s = 0; s < totalNumSpecies; s++)
  {
    IONS->at(s).SPECIES = params->ions_params[s].SPECIES;
    IONS->at(s).M       = params->ions_params[s].M;
    IONS->at(s).Z       = params->ions_params[s].Z;
    IONS->at(s).Q       = F_E*IONS->at(s).Z;

    // Calculating number of computational particles to store EACH MPI's output:
    IONS->at(s).pct_N_CP_MPI_Output = params->ions_params[s].pct_N_CP_Output;
    IONS->at(s).N_CP_MPI_Output     = floor((IONS->at(s).pct_N_CP_MPI_Output/100.0)*IONS->at(s).N_CP_MPI);

    double M = IONS->at(s).M;
    double Z = IONS->at(s).Z;
    double Q = IONS->at(s).Q;
    double f = params->ions_IC.at(s).densityFraction;

    // Characteristic frequencies:
    IONS->at(s).Wc    = Q*params->CV.B/M;
    IONS->at(s).Wp    = sqrt(f*params->CV.ne*Q*Q/(F_EPSILON*M));

    // Characteristic thermal velocities:
    // IONS->at(s).VTper = sqrt(2.0*F_KB*params->CV.Tper/M);
    // IONS->at(s).VTpar = sqrt(2.0*F_KB*params->CV.Tpar/M);

    // Characteristic lengths and time scales:
    // IONS->at(s).LarmorRadius = IONS->at(s).VTper/IONS->at(s).Wc;
    IONS->at(s).SkinDepth    = F_C/IONS->at(s).Wp;
    IONS->at(s).GyroPeriod   = 2.0*M_PI/IONS->at(s).Wc;
  }

  // Characteristic length and time scales of simulation:
  // ===================================================
  // Assume that species 0 is the majority species
  params->ionSkinDepth  = IONS->at(0).SkinDepth;
  params->ionGyroPeriod = IONS->at(0).GyroPeriod;

  // Allocate memory for ion arrays:
  // ===========================================================
  for (int s = 0; s < totalNumSpecies; s++)
  {
    allocate_meshDefinedIonArrays(params,&IONS->at(s));
    if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
    {
        allocate_particleDefinedIonArrays(params, &IONS->at(s));
    }
  }

  // Distribute computational particles according to IC profiles:
  // ===========================================================
  if (params->mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    int Nx    = params->mesh_params.Nx;
    double dx = params->mesh_params.dx;

    // Loop over all ion species:
    for (int s = 0; s < totalNumSpecies; s++)
    {
      int i_start = 0;
      double M = IONS->at(s).M;
      uint N_CP = IONS->at(s).N_CP_MPI;
      arma::vec N_CP_cell = round(IC->ions.at(s).ncp_MPI_mg*dx);
      arma::vec x_center  = IC->ions.at(s).x_mg;
      arma::vec a_m       = IC->ions.at(s).a_mg;
      arma::vec Tpar_m    = IC->ions.at(s).Tpar_mg;
      arma::vec Tper_m    = IC->ions.at(s).Tper_mg;
      arma::vec upar_m    = IC->ions.at(s).upar_mg;

      // Temporary variables to hold particle-defined attributes:
      arma::vec x_p(N_CP);
      arma::vec a_p(N_CP);
      arma::mat v_p(N_CP,2);

      // Loop over all grid cells:
      for (int m = 0; m < Nx; m++)
      {
        // Number of computational particles in mth cell:
        int N = N_CP_cell(m+1);

        // Uniform random numbers between -1/2 to +1/2:
        arma::arma_rng::set_seed(params->mpi.MPI_DOMAIN_NUMBER);
        arma::vec R1 = arma::randu(N) - 0.5;

        // Index range:
        int i_end   = i_start + N - 1;

        // Particle position:
        // ==================
        x_p.subvec(i_start,i_end) = x_center(m+1) + R1*dx;

        // Particle weight:
        // ==================
        a_p.subvec(i_start,i_end) = a_m(m+1)*arma::ones(N);

        // Parallel velocity:
        // ==================
        double vTpar = sqrt(2*F_KB*Tpar_m(m+1)/M);
        arma::arma_rng::set_seed(params->mpi.MPI_DOMAIN_NUMBER);
        arma::vec X1 = arma::randu(N);
        arma::vec X2 = arma::randu(N);
        arma::vec xpar = sqrt(-log(X1))%cos(2*M_PI*X2);

        v_p(span(i_start,i_end),0) = upar_m(m+1) + xpar*vTpar;

        // Perpendicular velocity:
        // ==================
        double vTper = sqrt(2*F_KB*Tper_m(m+1)/M);
        arma::arma_rng::set_seed(params->mpi.MPI_DOMAIN_NUMBER);
        X1 = arma::randu(N);
        X2 = arma::randu(N);
        arma::vec xy = sqrt(-log(X1))%cos(2*M_PI*X2);
        arma::vec xz = sqrt(-log(X1))%sin(2*M_PI*X2);
        arma::vec xper = sqrt( pow(xy,2) + pow(xz,2) );

        v_p(span(i_start,i_end),1) = xper*vTper;

        // Increment particle counter:
        i_start = i_end + 1;
      }

      // Now we need to "shuffle" the position of each computational particle in memory:
      // Create random permutation of unsigned integers 1:N_CP list :
      arma::uvec random_index_list = arma::randperm(N_CP);

      // Assign to IONS->at(s).x_p, a_p, v_p:
      IONS->at(s).x_p = x_p(random_index_list);
      IONS->at(s).a_p = a_p(random_index_list);
      IONS->at(s).v_p = v_p.rows(random_index_list);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * ION INITIALIZATION FINAIZED * * * * * * * * * * * * * * * * * *" << endl;
  }
}

void init_TYP::allocate_meshDefinedIonArrays(params_TYP * params, ions_TYP * IONS)
{
  MPI_Barrier(MPI_COMM_WORLD);

  int Nx = params->mesh_params.Nx;

  // Computational particle density per MPI:
  IONS->ncp_m.zeros(Nx + 2);

  // Ion density:
  IONS->n_m.zeros(Nx + 2);
  IONS->n_m_.zeros(Nx + 2);
  IONS->n_m__.zeros(Nx + 2);
  IONS->n_m___.zeros(Nx + 2);

  // Ion flux density:
  IONS->nv_m.zeros(Nx + 2);
  IONS->nv_m_.zeros(Nx + 2);
  IONS->nv_m__.zeros(Nx + 2);

  // Pressure tensors:
  IONS->P11_m.zeros(Nx + 2);
  IONS->P22_m.zeros(Nx + 2);

  // Derived quantities:
  IONS->Tpar_m.zeros(Nx + 2);
  IONS->Tper_m.zeros(Nx + 2);

  MPI_Barrier(MPI_COMM_WORLD);
}

void init_TYP::allocate_particleDefinedIonArrays(params_TYP * params, ions_TYP * IONS)
{
  // Numnber of computational particles per MPI:
  uint N_CP = IONS->N_CP_MPI;

  // Initialize particle-defined quantities:
  // ==================================
  IONS->x_p.zeros(N_CP);
  IONS->v_p.zeros(N_CP,2);
  IONS->a_p.zeros(N_CP);
  IONS->mu_p.zeros(N_CP);
  IONS->mn.zeros(N_CP);

  IONS->Ex_p.zeros(N_CP);
  IONS->Bx_p.zeros(N_CP);
  IONS->dBx_p.zeros(N_CP);
  IONS->ddBx_p.zeros(N_CP);

  IONS->wxc.zeros(N_CP);
  IONS->wxl.zeros(N_CP);
  IONS->wxr.zeros(N_CP);

  IONS->n_p.zeros(N_CP);
  IONS->nv_p.zeros(N_CP);
  IONS->Tpar_p.zeros(N_CP);
  IONS->Tper_p.zeros(N_CP);
  IONS->Te_p.zeros(N_CP);

  // Initialize particle defined flags:
  // ==================================
  IONS->f1.zeros(N_CP);
  IONS->f2.zeros(N_CP);
  IONS->f3.zeros(N_CP);
  //IONS->f4.zeros(N_CP);
  IONS->f5.zeros(N_CP);

  // Initialize particle kinetic energy at boundaries:
  // ================================================
  IONS->dE1.zeros(N_CP);
  IONS->dE2.zeros(N_CP);
  IONS->dE3.zeros(N_CP);
  //IONS->dE4.zeros(N_CP);
  IONS->dE5.zeros(N_CP);

  // Initialize resonance number:
  // ============================
  IONS->resNum.zeros(N_CP);
  IONS->resNum_.zeros(N_CP);

  // Rf terms:
  // ========
  IONS->udErf.zeros(N_CP);
  IONS->doppler.zeros(N_CP);
  IONS->udE3.zeros(N_CP);
}
