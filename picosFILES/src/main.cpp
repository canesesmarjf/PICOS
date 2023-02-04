// PICOS++
// Include intrinsic header files:
// =============================================================================
#include <iostream>
#include <vector>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <cmath>
#include <ctime>
#include <utility>

// Include user-defined header files:
// =============================================================================
#include "types.h"
#include "initialize.h"
#include "simple_output_H5.h"
#include "units.h"
#include "outputHDF5.h"
#include "PIC.h"
#include "fieldSolve.h"
#include "particleBC.h"
#include "collisionOperator.h"
// #include "rfOperator.h"

// Include headers for parallelization:
// =============================================================================
#include <omp.h>
#include "mpi_main.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[])
{
  // Initialize MPI process:
  MPI_Init(&argc, &argv);

  // MPI object to hold topology information:
  MPI_MAIN_TYP mpi_main;

  // Input parameters for simulation:
  params_TYP params;

  // Ion species vector:
  vector<ions_TYP> IONS;

  // Electron species object:
  electrons_TYP electrons;

  // Electromagnetic fields:
  fields_TYP fields;

  // Characteristic scales:
  CS_TYP CS;

  // Mesh object:
  mesh_TYP mesh;

  // Object to hold IC condition profiles:
  IC_TYP IC;

  // UNITS object:
  units_TYP units;

  // Collision operator object:
  coll_operator_TYP coll_operator;

  // Particle boundary condition operator:
  particleBC_TYP particleBC;

  // Initialize object:
  init_TYP init(&params, argc, argv);

  // =========================================================================
  // - Read "input_file.input" into "params"
  init.read_inputFile(&params);

  // - Read "ions_properties.ion" into "params":
  init.read_ionsPropertiesFile(&params);

  // - Create mesh using input parameters:
  init.create_mesh(&params,&mesh);

  // Create MPI topology:
  mpi_main.createMPITopology(&params);

  // - Read IC profiles from external files into "IC":
  init.read_IC_profiles(&params,&mesh,&IC);

  // - Interpolate IC profiles to mesh grid:
  init.interpolate_IC_profiles(&params,&mesh,&IC);

  // Calculate particle weight initial condition profile:
  init.calculate_IC_particleWeight(&params,&IC,&IONS);

  // - Initialize fields using IC field profiles:
  init.initialize_fields(&params,&IC,&fields);

  // - Initialize electrons using IC profiles:
  init.initialize_electrons(&params,&IC,&electrons);

  // - Initialize ions using IC profiles:
  init.initialize_ions(&params,&IC,&mesh,&IONS);

  // Define characteristic scales and broadcast them to all processes in COMM_WORLD:
  units.defineCharacteristicScalesAndBcast(&params, &IONS, &CS);

  // =========================================================================
  //  CONSIDER REMOVING FS
  // This will require that FS be removed from outputHDF5 source codes
  // =========================================================================
  // Create object and allocate memory according to "params":
  // FS_TYP FS(&params);

  // Define fundamental scales and broadcast them to all processes in COMM_WORLD:
  // units.calculateFundamentalScalesAndBcast(&params, &IONS, &FS);

  // Check that mesh size is consistent with hybrid approximation:
  // units.spatialScalesSanityCheck(&params, &FS);
  // =========================================================================

  // HDF object constructor and create "main.h5"
  HDF_TYP HDF(&params, &mesh, &IONS);

  // Define time step based on ion CFL condition:
  units.defineTimeStep(&params, &IONS);

  // Normalize "params", "IONS", "electrons", "fields" using "CS"
  units.normalizeVariables(&params, &mesh, &IONS, &electrons, &fields, &CS);

  // #########################################################################
  /**************** All the quantities below are dimensionless ****************/
  // #########################################################################

  // Definition of variables for advancing in time particles and fields:
  // =========================================================================
  double t1 = 0.0;
  double t2 = 0.0;
  params.currentTime = 0.0;
  int outputIterator = 0;
  int numberOfIterationsForEstimator = 1000;

  // Create EM solver:
  // =========================================================================
  fields_solver_TYP fields_solver(&params, &CS);

  // Create PIC solver:
  // =========================================================================
  PIC_TYP PIC(&params, &mesh, &fields, &IONS, &electrons);

  HDF_simple_TYP HDF_simple;
  if (params.mpi.COMM_COLOR == PARTICLES_MPI_COLOR)
  {
    string fileName;
    stringstream kk;
    kk << params.mpi.MPI_DOMAIN_NUMBER;

    // General:
    fileName = "file_1_rank_" + kk.str() + ".h5";
    HDF_simple.saveData(fileName,&params,&fields,&IONS);
  }

  // Create RF operator object:
  // =========================================================================
  // RF_Operator_TYP RF_operator(&params,&CS,&fields,&IONS);

  // Save 1st output:
  // =========================================================================
  HDF.saveOutputs(&params, &IONS, &electrons, &fields, &CS, 0, 0);

  // Start timing simulations:
  // =========================================================================
  t1 = MPI_Wtime();

  // #########################################################################
  // Start time iterations:
  // #########################################################################
  for(int tt=0; tt<params.timeIterations; tt++)
  {

    // if (params.mpi.IS_PARTICLES_ROOT)
    // {
    //   if (fmod((double)(tt + 1), 100) == 0)
    //   {
    //       cout << "time = " << tt*params.DT*CS.time*1E3 << " [ms] "<< endl;
    //
    //
    //       //cout << "N1_dot = " << particleBC.N1_dot/CS.time << " [1/s]" << endl;
    //       //cout << "N2_dot = " << particleBC.N2_dot/CS.time << " [1/s]" << endl;
    //       //cout << "E1_dot = " << (particleBC.E1_dot*CS.energy/CS.time)/1000 << " [kW]" << endl;
    //       //cout << "E2_dot = " << (particleBC.E2_dot*CS.energy/CS.time)/1000 << " [kW]" << endl;
    //       //cout << "E5_dot = " << (particleBC.E5_dot*CS.energy/CS.time)/1000 << " [kW]" << endl;
    //       //cout << "N5_dot = " << particleBC.N5_dot/CS.time << " [1/s]" << endl;
    //
    //   }
    // }

    // Advance particles and re-inject:
    // =====================================================================
    if (params.SW.advancePos == 1)
    {
        // Advance particle position and velocity to level X^(N+1):
        PIC.advanceParticles(&params, &mesh, &fields, &IONS);

        // Re-inject particles that leave computational domain:
        particleBC.applyParticleReinjection(&params,&CS,&fields,&IONS);

        // Assign cell:
        PIC.assignCell_AllSpecies(&params,&mesh,&IONS);

        // Interpolate all fields:
        PIC.interpolateFields_AllSpecies(&params,&IONS,&fields);

        // Interpolate electron temperature:
    	  PIC.interpolateElectrons_AllSpecies(&params,&IONS,&electrons);
    }

    // Calculate ion moments:
    // =====================================================================
    PIC.extrapolateMoments_AllSpecies(&params,&fields,&IONS);

    // Apply collision operator:
    // =====================================================================
    if (params.SW.Collisions == 1)
    {
      coll_operator.ApplyCollisions_AllSpecies(&params,&CS,&IONS,&electrons);
    }

    // Apply RF operator:
    // =====================================================================
    if (params.SW.RFheating == 1)
    {
        if (params.currentTime >= params.RF.t_ON*CS.time && params.currentTime <= params.RF.t_OFF*CS.time)
        {
            // RF_operator.ApplyRfHeating_AllSpecies(&params,&CS,&fields,&IONS);
        }
    }

    // if (params.mpi.IS_PARTICLES_ROOT)
    // {
    //     if (params.currentTime >= params.RF.t_ON*CS.time & params.currentTime <= params.RF.t_OFF*CS.time)
    //     {
    //         cout << "t = " << params.currentTime*1e3 << " [ms] " << endl;
    //     }
    // }

    // Field solve:
    // =====================================================================
    // Magnetic field:
    if (params.SW.BfieldSolve == 1)
    {
        // Use Faraday's law to advance the magnetic field to level B^(N+1).
        //fields_solver.advanceBfield(&params, &fields, &IONS);
    }

    if (params.SW.EfieldSolve == 1)
    {
        // Use Ohm's law to advance the electric field:
        fields_solver.advanceEfield(&params,&fields,&CS,&IONS,&electrons);
    }

    // Advance time:
    // =====================================================================
    params.currentTime += params.DT*CS.time;

    // Save data:
    // =====================================================================
    if(fmod((double)(tt + 1), params.outputCadenceIterations) == 0)
    {
        vector<ions_TYP> IONS_OUT = IONS;

        HDF.saveOutputs(&params, &IONS_OUT, &electrons, &fields, &CS, outputIterator+1, params.currentTime);

        outputIterator++;

        if(params.mpi.MPI_DOMAIN_NUMBER == 0)
        {
          cout << "data saved at tt = " << tt << endl;
        }
    }

    // Estimate simulation time:
    // =====================================================================
    if(tt == numberOfIterationsForEstimator)
    {
        t2 = MPI_Wtime();

        double estimatedSimulationTime = ( (double)params.timeIterations*(t2 - t1)/(double)numberOfIterationsForEstimator )/60.0;

        if(params.mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cout << "ESTIMATED TIME OF COMPLETION: " << estimatedSimulationTime <<" MINUTES" << endl;
        }
    }

  }

  // #########################################################################
  // End time iterations.
  // #########################################################################

  // Finalizing MPI communications:
  // =========================================================================
  mpi_main.finalizeCommunications(&params);

  cout << "end" << endl;

	return(0);
}
