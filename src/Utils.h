#ifndef _UTILS_HEADER
#define _UTILS_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <map>
#include <gsl/gsl_sf_bessel.h>
#ifdef _COMPLEX_BESSEL
#include <complex_bessel.h>
#endif
#include "Spline.h"
#include "ODESolver.h"

// The constants used in this code. Everything is here in SI units
extern struct ConstantsAndUnits {
  // Basic units (here we use SI)
  const double m           = 1.0;			                    // Length (in meters)
  const double s           = 1.0;                         // Time (in seconds)
  const double kg          = 1.0;                         // Kilo (in kilos)
  const double K           = 1.0;                         // Temperature (in Kelvins)

  // Derived units
  const double km          = 1e3 * m;                     // Kilometers
  const double N           = kg*m/(s*s);                  // Newton
  const double J           = N*m;                         // Joule
  const double W           = J/s;                         // Watt
  const double Mpc         = 3.08567758e22 * m;		        // Megaparsec
  const double eV          = 1.60217653e-19 * J;		      // Electronvolt
  
  // Physical constants    
  const double k_b         = 1.38064852e-23 * J/K;	      // Bolzmanns constant
  const double m_e         = 9.10938356e-31 * kg;	        // Mass of electron
  const double m_H         = 1.6735575e-27 * kg;	        // Mass of hydrogen atom
  const double c           = 2.99792458e8 * m/s;	        // Speed of light
  const double G           = 6.67430e-11 * N*m*m/(kg*kg);	// Gravitational constant
  const double hbar        = 1.054571817e-34 * J*s;		    // Reduced Plancks constant
  const double sigma_T     = 6.6524587158e-29 * m*m;		  // Thomas scattering cross-section
  const double lambda_2s1s = 8.227 / s;                   // Transition time is
  const double H0_over_h   = 100 * km/s/Mpc;              // H0 / h
  const double epsilon_0   = 13.605693122994 * eV;	      // Ionization energy for the ground state of hydrogen
  const double xhi0        = 24.587387 * eV;			        // Ionization energy for neutral Helium
  const double xhi1        = 4.0 * epsilon_0;		          // Ionization energy for singly ionized Helium
  
  // Min and max k-value
  const double k_min = 0.00005 / Mpc;
  const double k_max = 0.3     / Mpc;
  
  // Min and max x-value
  const double x_start = log(1e-8);
  const double x_end   = 0.0;

  // Include polarization and/or neutrinos?
  const bool polarization  = true;
  const bool neutrinos     = true;

  // For integration of perturbations (number of equations and positions in arrays)
  const int n_scalars           = 5;
  const int n_ell_theta         = 8;
  const int n_ell_thetap        = 8 * polarization;
  const int n_ell_neutrinos     = 8 * neutrinos;
  const int n_ell_tot_full      = n_scalars + n_ell_theta + n_ell_thetap + n_ell_neutrinos;
  const int ind_deltacdm        = 0; 
  const int ind_deltab          = 1;
  const int ind_vcdm            = 2;
  const int ind_vb              = 3;
  const int ind_Phi             = 4;
  const int ind_start_theta     = n_scalars;
  const int ind_start_thetap    = ind_start_theta  + n_ell_theta;
  const int ind_start_nu        = ind_start_thetap + n_ell_thetap;
 
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_scalars_tc        = 5;
  const int n_ell_theta_tc      = 2;
  const int n_ell_thetap_tc     = 0;
  const int n_ell_neutrinos_tc  = n_ell_neutrinos;
  const int n_ell_tot_tc        = n_scalars_tc + n_ell_theta_tc + n_ell_thetap_tc + n_ell_neutrinos;
  const int ind_deltacdm_tc     = 0; 
  const int ind_deltab_tc       = 1;
  const int ind_vcdm_tc         = 2;
  const int ind_vb_tc           = 3;
  const int ind_Phi_tc          = 4;
  const int ind_start_theta_tc  = n_scalars_tc;
  const int ind_start_thetap_tc = ind_start_theta_tc  + n_ell_theta_tc;
  const int ind_start_nu_tc     = ind_start_thetap_tc + n_ell_thetap_tc;

} Constants;

namespace Utils {
  
  // Find the x-value such that y(x) = y_value
  double binary_search_for_value(
      Spline &y, 
      double y_value, 
      std::pair<double,double> xrange = {0.0,0.0}, 
      double epsilon = 1e-7);

  // Ordinary bessel functions
  double J_n(const int n, const double x);
  
  // Spherical bessel functions
  double j_ell(const int n, const double x);
  std::vector<double> j_ell_array(const int lmax, const double x);
  
  // Generate an array with n equispaced points from xmin to xmax
  std::vector<double> linspace(double xmin, double xmax, int num);
  
  // Take the derivative of a function (simple 2pt stencil)
  std::vector<double> derivative(std::vector<double> &x, std::vector<double> &f);
  
  // For timing
  void StartTiming(std::string &name);
  void EndTiming(std::string &name);
  void StartTiming(const char *name);
  void EndTiming(const char *name);
  std::chrono::time_point<std::chrono::steady_clock> getTime();
  double timeInSeconds(
      std::chrono::time_point<std::chrono::steady_clock> & time_start, 
      std::chrono::time_point<std::chrono::steady_clock> & time_end);

}

// Macro to be able to take usual math functions of a whole vector
#define FUNS(FUN) \
  std::vector<double> FUN(const std::vector<double>& x);
FUNS(exp); FUNS(log); FUNS(cos); FUNS(sin); FUNS(tan); FUNS(fabs); FUNS(atan);
#undef FUNS

#endif
