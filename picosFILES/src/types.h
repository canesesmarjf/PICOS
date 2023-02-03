#ifndef H_TYPES
#define H_TYPES

#include <typeinfo>
#include <iostream>
#include <vector>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <string>
#include <map>

#include <omp.h>
#include "mpi.h"

using namespace std;

//  Define macros:
// =============================================================================
#define FIELDS_MPI_COLOR 0
#define PARTICLES_MPI_COLOR 1
#define FIELDS_TAG 100
#define PARTICLES_TAG 200

#define float_zero 1E-7
#define double_zero 1E-15

// Physical constants
// =============================================================================
#define PRO_ZERO 1.0E-15	// Definition of zero in PROMETHEUS
#define F_E 1.602176E-19	// Electron charge in C (absolute value)
#define F_ME 9.109382E-31	// Electron mass in kg
#define F_MP 1.672621E-27	// Proton mass in kg
#define F_U 1.660538E-27	// Atomic mass unit in kg
#define F_KB 1.380650E-23	// Boltzmann constant in Joules/Kelvin
#define F_EPSILON 8.854E-12 // Vacuum permittivity in C^2/(N*m^2)
#define F_C 299792458.0 	// Light speed in m/s
#define F_MU (4*M_PI)*1E-7 	// Vacuum permeability in N/A^2
extern double F_EPSILON_DS; // Dimensionless vacuum permittivity
extern double F_E_DS; 		// Dimensionless electron charge
extern double F_ME_DS; 		// Dimensionless electron mass
extern double F_MU_DS; 		// Dimensionless vacuum permeability
extern double F_C_DS; 		// Dimensionless speed of light

// Class to hold ion parameters:
// =============================================================================
class ions_params_TYP
{
public:
  int SPECIES;
  int N_CP_IN_SIM; // Number of computational particles over entire simulation
  int pct_N_CP_Output; // percentage of N_CP recorded in output file
  int Z;
  double M;
};

// Class to hold the initial conditions PARAMETERS for each ion species:
// =============================================================================
class ions_IC_params_TYP
{
public:
    string fileName;
    double Tpar_offset;
    double Tpar_scale;
    double Tper_offset;
    double Tper_scale;
    double upar_offset;
    double upar_scale;
    double densityFraction;
    double mean_ai; // Mean particle weight

    string CP_fileName; // Filename that indicates where the computational particle density profile is to be found.
};

// Class to hold the boundary conditions PARAMETERS for each ion species:
// =============================================================================
struct ions_BC_params_TYP
{
    int type;
    double T;
    double E;
    double eta;
    double mean_x;
    double sigma_x;

    double G;
    string G_fileName; // file that stores time dependent source rate
};

// Class to hold the initial conditions PARAMETERS for electrons:
// =============================================================================
class electrons_IC_params_TYP
{
public:
    string fileName;
    double n_offset;
    double n_scale;
    double T_offset;
    double T_scale;
};

// Class to hold the initial conditions PARAMETERS for fields:
// =============================================================================
class fields_IC_params_TYP
{
public:
    string fileName;
    double Ex_offset;
    double Ex_scale;
    double Bx_offset;
    double Bx_scale;
};

// Class to store mesh PARAMETERS:
// =============================================================================
class mesh_params_TYP
{
public:
  double dx_norm;
  double x0;
  double r0_min;
  double r0_max;
  double Lx_min;
  double Lx_max;

  double ionSkinDepth;
  int Nx;
	int Nx_PER_MPI;
  double dx;
  double A0;
  double B0;

  // Methods:
  void getA0();
};

// Class to store position of mesh cell-centers:
// =============================================================================
class mesh_TYP
{
public:
  int Nx;
	int Nx_PER_MPI;
  double dx;
  arma::vec xm;
  arma::vec xmg;
  // arma::vec Am;
  // arma::vec Bxm;

  // Overloaded constructor:
  // mesh_TYP(params_TYP * params);

};

// Class to hold initial conditon PROFILES for each ion species:
// =============================================================================
class ions_IC_TYP
{
public:
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec n;
    arma::vec Tpar;
    arma::vec Tper;
    arma::vec upar;
    arma::vec ncp_pdf; // PDF that defines initial computational particle density profile

    // Store profiles interpolated at the cell-centers of the mesh with ghost cells included:
    arma::vec x_mg;
    arma::vec n_mg;
    arma::vec Tpar_mg;
    arma::vec Tper_mg;
    arma::vec upar_mg;
    arma::vec ncp_shape_mg;

    // IC profiles calculated in the code (with ghost cells included):
    arma::vec a_mg;
    arma::vec ncp_in_SIM_mg;
		arma::vec ncp_MPI_mg;
};

// Class to hold the initial conditions PROFILES for the electrons:
// =============================================================================
struct electrons_IC_TYP
{
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec n;
    arma::vec T;

    // Store profiles interpolated at the cell-centers of the mesh:
    arma::vec x_mg;
    arma::vec n_mg;
    arma::vec T_mg;
};

// Class to hold the initial conditions PROFILES for the fields:
// =============================================================================
struct fields_IC_TYP
{
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec Bx;
    arma::vec Ex;
    arma::vec dBx;
    arma::vec ddBx;

    // Store profiles interpolated at the cell-centers of the mesh:
    arma::vec x_mg;
    arma::vec Bx_mg;
    arma::vec Ex_mg;
    arma::vec dBx_mg;
    arma::vec ddBx_mg;
};

// Class to hold all initial condition PROFILES:
// =============================================================================
class IC_TYP
{
public:
    vector<ions_IC_TYP> ions;
    electrons_IC_TYP electrons;
    fields_IC_TYP fields;
};

// Class to represent each simulated ion species:
// =============================================================================
class ions_TYP
{

	// Notes:
  // Consider double get_K();
  // Consider double get_NR(); Keep in mind that this calculation will be local to the MPI
  // To get the real N_R we would need to reduce over all MPI process and then broadcast

public:

	int SPECIES;     // 0: tracer 1: GC particles
	uint N_CP_IN_SIM; // total number of computational particles in ENTIRE simulation
	uint N_CP_MPI;    // Number of computational particles per MPI process
	double K;         // Distribution function normalization constant K = N_R/N_SP

	double N_R;    // Number of real particles in ENTIRE simulation (Dynamic)
	double N_SP;   // Number of super-particles in ENTIRE simulation (Dynamic)

	double Z;      // Atomic number
	double M;      // Mass
	double Q;      // Charge

	// Variables to control computational particle output for EACH MPI:
	double pct_N_CP_MPI_Output; // Percentage of computational particles PER MPI to record in each MPI's output file
	int N_CP_MPI_Output;        // Number of computational particles PER MPI to record in output file

	// Characteristic scales based on CV ne and CV B:
	// ======================
	// double LarmorRadius;
	double GyroPeriod;
	double SkinDepth;
	// double VTper;
	// double VTpar;
	double Wc;
	double Wp;

	// Particle-defined quantities:
	// ============================
	// Attributes:
	arma::vec x_p;	// Particle position
	arma::mat v_p;  // Velocity vector, V(0): parallel, V(1): perpendicular
	arma::vec a_p;  // Computational particle weight
	arma::vec mu_p; // Magnetic moment

	// Nearest grid point:
	arma::ivec mn;  // Ions' position in terms of the index of mesh node

	// Particle-defined fields:
	arma::vec Ex_p;
	arma::vec Bx_p;
	arma::vec dBx_p;
	arma::vec ddBx_p;

	// Particle-defined ion and electron moments:
	arma::vec n_p;
	arma::vec nv_p;
	arma::vec Tpar_p;
	arma::vec Tper_p;
	arma::vec Te_p;

	// Assignment function values for TSC shape function
	arma::vec wxl;				// Left node
	arma::vec wxc;				// Central node
	arma::vec wxr;				// Right node

	// Mesh-defined quantities:
	// ============================
	// These quantities record the values reduced over all MPI processes.
	// Each MPI process will have a partial mesh-defined quantity which then needs to be reduced over all MPI processes to get the final value for the mesh-defined quantity.
	arma::vec n_m;
	arma::vec n_m_;
	arma::vec n_m__;
	arma::vec n_m___;

	arma::vec nv_m;
	arma::vec nv_m_;
	arma::vec nv_m__;
	arma::vec nv_m___;

	arma::vec P11_m;				// Ion pressure tensor, component 1,1
	arma::vec P22_m;				// Ion pressure tensor, component 2,2
	arma::vec Tpar_m;			// Ion parallel temperature
	arma::vec Tper_m;			// Ion perpendicular temperature

	arma::vec ncp_m;     // Computational particle density in each MPI process

	// Particle defined flags:
	arma::ivec f1;             	// Flag for left boundary
	arma::ivec f2;              // Flag for Right boundary
	arma::ivec f3;              // Flag for RF operator
	arma::ivec f4;              // Flag for Collisiona operator
	arma::ivec f5; 							// Flag for injected particle

	// Particle kinetic energy at boundaries:
	arma::vec dE1;              // left boundary
	arma::vec dE2;              // Right boundary
	arma::vec dE3;              // RF operator
	arma::vec dE4;							// Collisional energy loss
	arma::vec dE5;							// Kinetic energy if injected particle

	// Resonance numnber:
	arma::vec resNum;
	arma::vec resNum_;

	// Rf terms:
	arma::vec udErf;           // Perpendicular energy gained by particle during unit electric field RF interaction
	arma::vec doppler;         // Doppler shift term
	arma::vec udE3;						 // Total energy gained by particle during unit Electric field RF interaction

	// Constructor:
	ions_TYP(){};

	// Destructor:
	~ions_TYP(){};
};

//  Define ELECTRONS DERIVED TYPE:
// =============================================================================
class electrons_TYP
{
public:

	// Mesh-defined temperature:
	arma::vec Te_m;

	// Constructor:
	electrons_TYP(){};

	// Destructor:
	~electrons_TYP(){};
};


// Class to represent electromagnetic fields in the simulation:
// =============================================================================
class fields_TYP
{
public:

    arma::vec x_mg;
    arma::vec Ex_m;
    arma::vec Bx_m;
    arma::vec dBx_m;
    arma::vec ddBx_m;

    arma::vec Am;

    fields_TYP(){};

    //fields_TYP(unsigned int N) : Ex_m(N), Bx_m(N), dBx_m(N), ddBx_m(N) {};

    ~fields_TYP(){};

    //void zeros(unsigned int N);
    //void fill(double A);
    void getAm(double A0, double B0);
};

//  Define structure to store characteristic values for the normalization:
// =============================================================================
struct CV_TYP
{
	double ne;
	double Te;
	double B;
	double Tpar;
	double Tper;

	CV_TYP()
	{
		ne   = 0;
		Te   = 0;
		B    = 0;
		Tpar = 0;
		Tper = 0;
	}
};

//  Define structure to store switches that control physics modules:
// =============================================================================
struct SW_TYP
{
	int EfieldSolve;
	int BfieldSolve;
	int Collisions;
	int RFheating;
	int linearSolve;
	int advancePos;

	SW_TYP()
	{
		EfieldSolve   = 0;
		BfieldSolve   = 0;
		Collisions    = 0;
		RFheating     = 0;
		linearSolve   = 0;
		advancePos    = 0;
	}

};

//  Define structure to hold MPI parameters:
// =============================================================================
struct mpi_params_TYP
{
	int NUMBER_MPI_DOMAINS;
	int MPI_DOMAIN_NUMBER;
	int FIELDS_ROOT_WORLD_RANK;
	int PARTICLES_ROOT_WORLD_RANK;

	int MPIS_FIELDS;
	int MPIS_PARTICLES;

	int MPI_DOMAINS_ALONG_X_AXIS;

	MPI_Comm MPI_TOPO; // Cartesian topology

	// Particle pusher and field solver communicator params
	int COMM_COLOR;
	int COMM_SIZE;
	int COMM_RANK;
	MPI_Comm COMM;

	bool IS_FIELDS_ROOT;
	bool IS_PARTICLES_ROOT;

	int MPI_CART_COORDS_1D[1];
	int MPI_CART_COORDS_2D[2];
	vector<int *> MPI_CART_COORDS;

	unsigned int iIndex;
	unsigned int fIndex;

	unsigned int irow;
	unsigned int frow;
	unsigned int icol;
	unsigned int fcol;

	int MPI_DOMAIN_NUMBER_CART;
	int LEFT_MPI_DOMAIN_NUMBER_CART;
	int RIGHT_MPI_DOMAIN_NUMBER_CART;
	//int UP_MPI_DOMAIN_NUMBER_CART;
	//int DOWN_MPI_DOMAIN_NUMBER_CART;
};

// Define structure to hold characteristic scales:
// =============================================================================
struct CS_TYP
{
	double time;
	double velocity;
	double momentum;
	double energy;
	double length;
	double volume;
	double mass;
	double charge;
	double density;
	double eField;
	double bField;
	double temperature;
	double magneticMoment;
	double vacuumPermeability;
	double vacuumPermittivity;

	// Constructor:
	CS_TYP()
	{
		time = 0.0;
		velocity = 0.0;
		momentum = 0.0;
		energy = 0.0;
		length = 0.0;
		volume = 0.0;
		mass = 0.0;
		charge = 0.0;
		density = 0.0;
		eField = 0.0;
		bField = 0.0;
		magneticMoment = 0.0;
		vacuumPermeability = 0.0;
		vacuumPermittivity = 0.0;
	}
};

// Define structure to hold RF operator parameters and global values:
// =============================================================================
struct RF_TYP
{
	// RF parameters:
	// =============
	double Prf;
	int n_harmonic;
	double freq;
	double x1;
	double x2;
	double t_ON;
	double t_OFF;
	double kpar;
	double kper;
	int handedness;

	// Name and storage time-dependent RF power trace:
	// ========================================
	string Prf_fileName;
	int Prf_NS;
	arma::vec Prf_profile;

	// Total RF power accumulated over all species:
	// ============================================
	double E3;

	// Power accumulated over all species per unit electric field:
	// ==========================================================
	double uE3;

	// Global RF electric field:
	// =========================
	double Erf;

	// Constructor:
	// ============
	RF_TYP()
	{
		Prf  = 0;
		n_harmonic = 1;
		freq = 0;
		x1   = 0;
		x2   = 0;
		kpar = 0;
		kper = 0;
		Prf_NS = 0;
		handedness = 0;
		E3   = 0;
		uE3  = 0;
		Erf  = 0;
	}
};

//  Define structure to store simulation parameters:
// =============================================================================
struct params_TYP
{
	// List of variables in the outputs
	vector<string> outputs_variables;

	// Select method for particle RK4 integrator
	int advanceParticleMethod;

	//Control parameters for the simulation:
	// Path to save the outputs
	string PATH;

	int argc;
	char **argv;

	double smoothingParameter;
	double simulationTime; // In units of the shorter ion gyro-period in the simulation
	double currentTime = 0;
	int timeIterations;

	// Consider eliminating one of the following:
	// Or call the 2nd one dt_norm
	double DT;//Time step
	double DTc;//Cyclotron period fraction.

	// Consider deleting of not neded:
	// ===============================
	//int loadFields;
	int loadGrid;
	int usingHDF5;
	int outputCadenceIterations;
	arma::file_type outputFormat;//Outputs format (raw_ascii,raw_binary).
	//types_info typesInfo;
	// ===============================

	double outputCadence;//Save variables each "outputCadence" times the background ion cycloperiod.

	// Mesh parameters:
	mesh_params_TYP mesh_params;

	// Ions properties:
	int numberOfParticleSpecies; // This species are evolved self-consistently with the fields
	int numberOfTracerSpecies; // This species are not self-consistently evolved with the fields

	// Ion parameters:
	vector<ions_params_TYP> ions_params;

	// Initial conditions parameters:
	electrons_IC_params_TYP electrons_IC;
	fields_IC_params_TYP fields_IC;
	vector<ions_IC_params_TYP> ions_IC;

	// Boundary condition parameters:
	vector<ions_BC_params_TYP> ions_BC;

	// Simulation Characterstic values:
	CV_TYP CV;

	// Simulation switches:
	SW_TYP SW;

	// RF operator conditions:
	RF_TYP RF;

	int filtersPerIterationFields;
	int filtersPerIterationIons;

	// double ionLarmorRadius;
	double ionSkinDepth;
	double ionGyroPeriod;

	//double DrL;
	// double dp;

	// int checkStability;
	// int rateOfChecking;//Each 'rateOfChecking' iterations we use the CLF criteria for the particles to stabilize the simulation

	// MPI parameters
	mpi_params_TYP mpi;

	// Error codes
	map<int,string> errorCodes;

	//  Methods:
	void getCharacteristicIonSkinDepth();
	void get_Nx_dx(double ionSkinDepth);

	// Constructor
	params_TYP(){};
};

// Define structure to hold fundamental scales:
// =============================================================================
struct FS_TYP
{
	double electronSkinDepth;
	double electronGyroPeriod;
	double electronGyroRadius;
	double * ionSkinDepth;
	double * ionGyroPeriod;
	double * ionGyroRadius;

	// Constructor:
	FS_TYP(params_TYP * params)
	{
		electronSkinDepth = 0.0;
		electronGyroPeriod = 0.0;
		electronGyroRadius = 0.0;
		ionSkinDepth = new double[params->numberOfParticleSpecies];
		ionGyroPeriod = new double[params->numberOfParticleSpecies];
		ionGyroRadius = new double[params->numberOfParticleSpecies];
	}
};

#endif
