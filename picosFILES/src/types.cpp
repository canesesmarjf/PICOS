#include "types.h"

double F_EPSILON_DS = F_EPSILON;    // Dimensionless vacuum permittivity
double F_E_DS = F_E;                // Dimensionless electric charge
double F_ME_DS = F_ME;              // Dimensionless electron mass
double F_MU_DS = F_MU;              // Dimensionless vacuum permeability
double F_C_DS = F_C;                // Dimensionless speed of light
double F_KB_DS = F_KB;              // Dimensionless Boltzmann constant

// =============================================================================
void params_TYP::getCharacteristicIonSkinDepth()
{
  // Calculate characteristic skin depth:
  double n = this->CV.ne;
  double M = this->ions_params[0].M;
  double Z = this->ions_params[0].Z;
  double Q = F_E*Z;
  double omega_pi = sqrt(n*Q*Q/(F_EPSILON*M));
  double ionSkinDepth = F_C/omega_pi;

  // Save ion skin depth:
  this->mesh_params.ionSkinDepth = ionSkinDepth;
}

// =============================================================================
void params_TYP::get_Nx_dx(double ionSkinDepth)
{
  // First estimate of dx:
  double dx_norm = this->mesh_params.dx_norm;
  double dx = dx_norm*ionSkinDepth;

  // First estimate of NX:
  double L = this->mesh_params.Lx_max - this->mesh_params.Lx_min;
  int Nx = round(L/dx);

  // Define final NX by making it a multiple of 2:
  if (fmod(Nx, 2.0) > 0.0)
  {
      Nx -= 1;
  }
  this->mesh_params.Nx = Nx;

  // Get Nx per MPIs (field MPIs):
  this->mesh_params.Nx_PER_MPI = Nx/this->mpi.MPIS_FIELDS;

  // Final value of dx:
  this->mesh_params.dx = L/Nx;
}

// =============================================================================
void mesh_params_TYP::getA0()
{
  // Reference radii of flux surface that bounds the plasma:
  double rmax = this->r0_max;
  double rmin = this->r0_min;

  // Cross sectional area:
  this->A0 = M_PI*( pow(rmax,2) - pow(rmin,2) );
}

// =============================================================================
void fields_TYP::getAm(double A0, double B0)
{
  this->Am = A0*B0/this->Bx_m;
}
