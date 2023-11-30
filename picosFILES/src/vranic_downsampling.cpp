#include "vranic_downsampling.h"

// ======================================================================================
void vranic_TYP::down_sample_node_3D(merge_cell_TYP * set_N, merge_cell_TYP * set_M)
{
  // Create refereces:
  // This does NOT require copy operation
  arma::vec& wi   = set_N->wi;
  arma::vec& xi   = set_N->xi;
  arma::vec& yi   = set_N->yi;
  arma::vec& zi   = set_N->zi;
  arma::vec& wi_M = set_M->wi;
  arma::vec& xi_M = set_M->xi;
  arma::vec& yi_M = set_M->yi;
  arma::vec& zi_M = set_M->zi;

  // Calculate merge-cell statistics:
  double wt_N = sum(wi);
  arma::vec ri = wi/wt_N;

  // Expectation values:
  double E_x = dot(ri,xi);
  double E_y = dot(ri,yi);
  double E_z = dot(ri,zi);

  // Calculate deltas:
  arma::vec dx = xi - E_x;
  arma::vec dy = yi - E_y;
  arma::vec dz = zi - E_z;
  arma::vec dr = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

  // Calculate skewness vector for deltas:
  int n = 2;
  double mu3_dx = dot(ri,dx%pow(dr,n));
  double mu3_dy = dot(ri,dy%pow(dr,n));
  double mu3_dz = dot(ri,dz%pow(dr,n));
  double mu3_dr = sqrt(pow(mu3_dx,2.0) + pow(mu3_dy,2.0) + pow(mu3_dz,2.0));

  // Calculate the merge-cell coordinate system:
  arma::vec x_hat = {1,0,0};
  arma::vec x_hat_prime = {mu3_dx/mu3_dr,mu3_dy/mu3_dr,mu3_dz/mu3_dr};
  arma::vec z_hat_prime = cross(x_hat,x_hat_prime);
  z_hat_prime = z_hat_prime/norm(z_hat_prime);
  arma::vec y_hat_prime = cross(z_hat_prime,x_hat_prime);

  // Merge-cell coordinate system matrix:
  arma::mat e_prime(3,3);
  e_prime.col(0) = x_hat_prime;
  e_prime.col(1) = y_hat_prime;
  e_prime.col(2) = z_hat_prime;

  // Rotation matrix:
  arma::mat R = e_prime.t();

  // Convert deltas from standard coordinate to merge-cell coordinate system:
  arma::vec xi_prime = R(0,0)*xi + R(0,1)*yi + R(0,2)*zi;
  arma::vec yi_prime = R(1,0)*xi + R(1,1)*yi + R(1,2)*zi;
  arma::vec zi_prime = R(2,0)*xi + R(2,1)*yi + R(2,2)*zi;

  // Calculate merge-cell statistics in new coordinate system:
  // Expectation values:
  double E_x_prime = dot(ri,xi_prime);
  double E_y_prime = dot(ri,yi_prime);
  double E_z_prime = dot(ri,zi_prime);

  // Deltas:
  arma::vec dx_prime = xi_prime - E_x_prime;
  arma::vec dy_prime = yi_prime - E_y_prime;
  arma::vec dz_prime = zi_prime - E_z_prime;
  arma::vec dr_prime = sqrt(pow(dx_prime,2.0) + pow(dy_prime,2.0) + pow(dz_prime,2.0));

  // Standard deviation:
  double sigma_x_prime = sqrt(dot(ri,pow(dx_prime,2.0)));
  double sigma_y_prime = sqrt(dot(ri,pow(dy_prime,2.0)));
  double sigma_z_prime = sqrt(dot(ri,pow(dz_prime,2.0)));
  double sigma_r_prime = sqrt(dot(ri,pow(dr_prime,2.0)));

  // Calculate deltas of new set M:
  int M = set_M->n_elem;
  arma::vec dx_prime_M(M);
  arma::vec dy_prime_M(M);
  arma::vec dz_prime_M(M);

  dx_prime_M(0) = + sqrt(M/2.0)*sigma_x_prime;
  dx_prime_M(1) = - sqrt(M/2.0)*sigma_x_prime;
  dx_prime_M(2) = 0;
  dx_prime_M(3) = 0;
  dx_prime_M(4) = 0;
  dx_prime_M(5) = 0;

  dy_prime_M(0) = 0;
  dy_prime_M(1) = 0;
  dy_prime_M(2) = + sqrt(M/2.0)*sigma_y_prime;
  dy_prime_M(3) = - sqrt(M/2.0)*sigma_y_prime;
  dy_prime_M(4) = 0;
  dy_prime_M(5) = 0;

  dz_prime_M(0) = 0;
  dz_prime_M(1) = 0;
  dz_prime_M(2) = 0;
  dz_prime_M(3) = 0;
  dz_prime_M(4) = + sqrt(M/2.0)*sigma_z_prime;
  dz_prime_M(5) = - sqrt(M/2.0)*sigma_z_prime;

  // Convert deltas of M into standard frame:
  R = e_prime;
  arma::vec dx_M = R(0,0)*dx_prime_M + R(0,1)*dy_prime_M + R(0,2)*dz_prime_M;
  arma::vec dy_M = R(1,0)*dx_prime_M + R(1,1)*dy_prime_M + R(1,2)*dz_prime_M;
  arma::vec dz_M = R(2,0)*dx_prime_M + R(2,1)*dy_prime_M + R(2,2)*dz_prime_M;

  // Produce vectors for M set:
  wi_M = arma::ones<arma::vec>(M)*wt_N/M;
  xi_M = E_x + dx_M;
  yi_M = E_y + dy_M;
  zi_M = E_z + dz_M;
}

// ======================================================================================
void vranic_TYP::down_sample_node_2D(merge_cell_TYP * set_N, merge_cell_TYP * set_M)
{
  // Create refereces:
  // This does NOT require copy operation
  arma::vec& wi   = set_N->wi;
  arma::vec xi    = set_N->xi;
  arma::vec& yi   = set_N->yi;
  arma::vec& zi   = set_N->zi;
  arma::vec& wi_M = set_M->wi;
  arma::vec& xi_M = set_M->xi;
  arma::vec& yi_M = set_M->yi;
  arma::vec& zi_M = set_M->zi;

  // Since this method has been repurposed from the 3D case, to make it 2D we need to neglect the non-velocity variable which in this case is the "xi" variable (physical space):
  xi = xi*0;

  // Calculate merge-cell statistics:
  double wt_N = sum(wi);
  arma::vec ri = wi/wt_N;

  // Expectation values:
  double E_x = dot(ri,xi);
  double E_y = dot(ri,yi);
  double E_z = dot(ri,zi);

  // Calculate deltas:
  arma::vec dx = xi - E_x;
  arma::vec dy = yi - E_y;
  arma::vec dz = zi - E_z;
  arma::vec dr = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

  // Calculate skewness vector for deltas:
  int n = 2;
  double mu3_dx = dot(ri,dx%pow(dr,n));
  double mu3_dy = dot(ri,dy%pow(dr,n));
  double mu3_dz = dot(ri,dz%pow(dr,n));
  double mu3_dr = sqrt(pow(mu3_dx,2.0) + pow(mu3_dy,2.0) + pow(mu3_dz,2.0));

  // Calculate the merge-cell coordinate system:
  arma::vec x_hat = {1,0,0};
  arma::vec x_hat_prime = {mu3_dx/mu3_dr,mu3_dy/mu3_dr,mu3_dz/mu3_dr};
  arma::vec z_hat_prime = cross(x_hat,x_hat_prime);
  z_hat_prime = z_hat_prime/norm(z_hat_prime);
  arma::vec y_hat_prime = cross(z_hat_prime,x_hat_prime);

  // Merge-cell coordinate system matrix:
  arma::mat e_prime(3,3);
  e_prime.col(0) = x_hat_prime;
  e_prime.col(1) = y_hat_prime;
  e_prime.col(2) = z_hat_prime;

  // Rotation matrix:
  arma::mat R = e_prime.t();

  // Convert deltas from standard coordinate to merge-cell coordinate system:
  arma::vec xi_prime = R(0,0)*xi + R(0,1)*yi + R(0,2)*zi;
  arma::vec yi_prime = R(1,0)*xi + R(1,1)*yi + R(1,2)*zi;
  arma::vec zi_prime = R(2,0)*xi + R(2,1)*yi + R(2,2)*zi;

  // Calculate merge-cell statistics in new coordinate system:
  // Expectation values:
  double E_x_prime = dot(ri,xi_prime);
  double E_y_prime = dot(ri,yi_prime);
  double E_z_prime = dot(ri,zi_prime);

  // Deltas:
  arma::vec dx_prime = xi_prime - E_x_prime;
  arma::vec dy_prime = yi_prime - E_y_prime;
  arma::vec dz_prime = zi_prime - E_z_prime;
  arma::vec dr_prime = sqrt(pow(dx_prime,2.0) + pow(dy_prime,2.0) + pow(dz_prime,2.0));

  // Standard deviation:
  double sigma_x_prime = sqrt(dot(ri,pow(dx_prime,2.0)));
  double sigma_y_prime = sqrt(dot(ri,pow(dy_prime,2.0)));
  double sigma_z_prime = sqrt(dot(ri,pow(dz_prime,2.0)));
  double sigma_r_prime = sqrt(dot(ri,pow(dr_prime,2.0)));

  // Calculate deltas of new set M:
  int M = set_M->n_elem;
  arma::vec dx_prime_M(M);
  arma::vec dy_prime_M(M);
  arma::vec dz_prime_M(M);

  dx_prime_M(0) = + sqrt(M/2.0)*sigma_x_prime;
  dx_prime_M(1) = - sqrt(M/2.0)*sigma_x_prime;
  dx_prime_M(2) = 0;
  dx_prime_M(3) = 0;
  dx_prime_M(4) = 0;
  dx_prime_M(5) = 0;

  dy_prime_M(0) = 0;
  dy_prime_M(1) = 0;
  dy_prime_M(2) = + sqrt(M/2.0)*sigma_y_prime;
  dy_prime_M(3) = - sqrt(M/2.0)*sigma_y_prime;
  dy_prime_M(4) = 0;
  dy_prime_M(5) = 0;

  dz_prime_M(0) = 0;
  dz_prime_M(1) = 0;
  dz_prime_M(2) = 0;
  dz_prime_M(3) = 0;
  dz_prime_M(4) = + sqrt(M/2.0)*sigma_z_prime;
  dz_prime_M(5) = - sqrt(M/2.0)*sigma_z_prime;

  // Convert deltas of M into standard frame:
  R = e_prime;
  arma::vec dx_M = R(0,0)*dx_prime_M + R(0,1)*dy_prime_M + R(0,2)*dz_prime_M;
  arma::vec dy_M = R(1,0)*dx_prime_M + R(1,1)*dy_prime_M + R(1,2)*dz_prime_M;
  arma::vec dz_M = R(2,0)*dx_prime_M + R(2,1)*dy_prime_M + R(2,2)*dz_prime_M;

  // Produce vectors for M set:
  wi_M = arma::ones<arma::vec>(M)*wt_N/M;
  xi_M = E_x + dx_M;
  yi_M = E_y + dy_M;
  zi_M = E_z + dz_M;
}

// ======================================================================================
void vranic_TYP::print_stats(merge_cell_TYP * set)
{
  // Create references:
  arma::vec& wi = set->wi;
  arma::vec& xi = set->xi;
  arma::vec& yi = set->yi;
  arma::vec& zi = set->zi;

  // Expectation values:
  double wt = sum(wi);
  arma::vec ri = wi/wt;
  double E_x = dot(ri,xi);
  double E_y = dot(ri,yi);
  double E_z = dot(ri,zi);

  // Calculate deltas:
  arma::vec dx = E_x - xi;
  arma::vec dy = E_y - yi;
  arma::vec dz = E_z - zi;
  arma::vec dr = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));
  double sigma_r = sqrt(dot(ri,pow(dr,2.0)));

  cout << "wt = " << wt << endl;

  cout << "E_x = " << E_x << endl;

  cout << "E_y   = " << E_y   << endl;

  cout << "E_z   = " << E_z   << endl;

  cout << "sigma_r     = " << sigma_r     << endl;
}

double vranic_TYP::get_sigma(merge_cell_TYP * set)
{
  // Create references:
  arma::vec& wi = set->wi;
  arma::vec& yi = set->yi;
  arma::vec& zi = set->zi;

  // Expectation values:
  double wt = sum(wi);
  arma::vec ri = wi/wt;
  double E_y = dot(ri,yi);
  double E_z = dot(ri,zi);

  // Calculate deltas:
  arma::vec dy = E_y - yi;
  arma::vec dz = E_z - zi;
  arma::vec dr = sqrt(pow(dy,2.0) + pow(dz,2.0));
  double sigma_r = sqrt(dot(ri,pow(dr,2.0)));

  return sigma_r;
}
