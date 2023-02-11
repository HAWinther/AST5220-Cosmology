//====================================================================================================
//
// Method to be used with code template to allow you to do MCMC and get constraints from supernova data
// 
// Things needed to be done to use this:
// * Check that line 109-116 is how you set up your background class and do the solving
// * Edit line 124 with the call to your luminosity function method (mine is called get_luminosity_distance_of_x). 
//   As written here this assumes that the luminosity distance is returned in meters from this method
// * Include this file in Main.cpp: #include "SupernovaFitting.h"
// * Call the function below in main: mcmc_fit_to_supernova_data("supernovadata.txt", "results.txt");
//
// The input [supernovadata] is the path to the file supernovadata.txt. The output file with all the samples is [outfile]
// This runs for [maxsteps] steps before ending. You can reduce maxsteps or just kill the run if it takes too long
//
// How to analyze the resulting chains:
// * Load the chains and skip the first few hundred samples (the burnin of the chains). E.g. loadtxt(file,skiprows=200) in python
// * Find the minimum chi2 and the corresponding best-fit parameters (you can use np.argmin to get index of the minvalue in python)
// * Select all samples of OmegaM and OmegaLambda (computed from OmegaM and OmegaK) that satisfy chi2 < chi2_min + 3.53 
//   (e.g. OmegaM[chi2 < chi2min + 3.53] in python)
// * Scatterplotting these gives you the 1sigma (68.4%) confidence region
// * Find the standard deviation of the samples to get the 1sigma confidence region of the parameters (assuming the posterior is a gaussian)
// * Make and plot a histogram of the samples for the different parameters (OmegaM, OmegaK, OmegaLambda, H0)
// * You can also compute the mean and standard deviation of the chain values and use this to overplot a gaussian with the same mean and variance for comparison.
//
//====================================================================================================

#include "Utils.h"
#include <climits>
#include <random>
#include <array>
#include <vector>

void mcmc_fit_to_supernova_data(std::string supernovadata_filename, std::string outfile){

  // The number of parameters to fit
  const int nparam = 3;
  // Maximal number of samples to generate
  const int maxsteps = 10000;
  // Seed for random number generator
  const int seed = 1234;
  // How many meters in a Gpc
  const double Gpc = Constants.Mpc * 1000;

  // Read data from file
  std::vector<double> z_arr;
  std::vector<double> L_arr;
  std::vector<double> dL_arr;
  auto read_data = [&](std::string filename){
    std::string header;
    std::ifstream fp(filename.c_str());
    if(!fp)
      throw std::runtime_error("Error: cannot open file " + filename);
    std::getline(fp, header);
    std::cout << "Reading luminosity data from file:\n";
    while(1){
      // Read line by line
      double z, L, dL;
      fp >> z;
      if(fp.eof()) break;
      fp >> L;
      fp >> dL;
      z_arr.push_back(z);
      L_arr.push_back(L);
      dL_arr.push_back(dL);
      std::cout << "z: " << z << " " << L << " " << dL << "\n";
    }
    std::cout << "We found n = " << z_arr.size() << " points\n";
  };
  read_data(supernovadata_filename);

  // Set up random number generator
  std::default_random_engine gen;
  gen.seed(seed);
  std::normal_distribution<double> ndist(0.0,1.0);
  std::uniform_real_distribution<double> udist(0.0,1.0);

  // Our priors (chi^2 -> Inf if outside these ranges)
  const std::array<double, nparam> prior_high {1.5, 1.0, 1.0};
  const std::array<double, nparam> prior_low {0.5, 0.0, -1.0};
  
  // Starting point for chain and step size
  // 0: h
  // 1: OmegaM (OmegaCDM + OmegaB)
  // 2: OmegaK
  std::array<double, nparam> parameters{0.7, 0.25, 0.0};
  std::array<double, nparam> stepsize{0.007, 0.05, 0.05};
  for(int i = 0; i < nparam; i++){
    parameters[i] = prior_low[i] + (prior_high[i]-prior_low[i])*udist(gen);
  }

  // Best-fit as we go along in the chain
  std::array<double, nparam> best_parameters = parameters;
  double chi2_min = std::numeric_limits<double>::max();
  
  // The chi^2 function
  auto comp_chi2 = [&](std::array<double, nparam> & parameters){
    // Priors: if outside range return huuuuge chi^2
    bool inside_prior = true;
    for(int i = 0; i < nparam; i++){
      if(parameters[i] > prior_high[i]) inside_prior = false;
      if(parameters[i] < prior_low[i]) inside_prior = false;
    }
    if(not inside_prior) return std::numeric_limits<double>::max();

    //=========================================================================================
    //======= Here we set up the cosmology class and solve to get the distance functions ======
    //=========================================================================================
    // Set parameters and compute background
    double param_OmegaB   = 0.05;           // Not important as we just sample and are sensitive to OmegaM = OmegaB+OmegaCDM
    double param_Neff     = 0.0;            // Not relevant at late times
    double param_TCMB     = 2.7255;         // Temperature of the CMB
    double param_h        = parameters[0];
    double param_OmegaCDM = parameters[1] - param_OmegaB; // OmegaCDM = OmegaM - OmegaB
    double param_OmegaK   = parameters[2];
    BackgroundCosmology cosmo(param_h, param_OmegaB, param_OmegaCDM, param_OmegaK, param_Neff, param_TCMB);
    cosmo.solve();
    //=========================================================================================

    // Compute chi^2
    double chi2 = 0.0;
    for (size_t i = 0; i < z_arr.size(); i++){
      double x = -std::log(1.0+z_arr[i]);
      //======= Here we call your distance function ======
      double L = cosmo.get_luminosity_distance_of_x(x) / Gpc; // Luminosity function in units of Gpc
      //==================================================
      chi2 += (L - L_arr[i]) * (L - L_arr[i]) / (dL_arr[i] * dL_arr[i]);
    }
    return chi2;
  };

  // Generic Metropolis MCMC algorithm. This requires no changes
  int steps = 0;
  int nsample = 0;
  double oldchi2 = std::numeric_limits<double>::max();
  std::ofstream out(outfile.c_str());
  out << "#          chi2            h           OmegaM           OmegaK               Acceptrate\n";
  while(nsample < maxsteps){
    steps++;
    
    // Generate new set of parameters
    std::array<double, nparam> new_parameters;
    for(int i = 0; i < nparam; i++)
      new_parameters[i] = parameters[i] + ndist(gen) * stepsize[i];

    // Compute chi^2
    double chi2 = 0.0;
    try {
      chi2 = comp_chi2(new_parameters);
    } catch(...){
      chi2 = std::numeric_limits<double>::max();
    }

    // Always accept sample if lower chi^2
    bool accept_step = chi2 < oldchi2;
    if(not accept_step){
      // Metropolis step - draw random number and accept if low enough
      if(std::exp(-(chi2-oldchi2)/2.0) > udist(gen)){
        accept_step = true;
      }
    }

    // If the step is accepted change parameters to new parameters and update the old chi^2 value
    if(accept_step){
      nsample++;
      oldchi2 = chi2;
      parameters = new_parameters;

      // Write sample to file and screen
      std::cout << "#          chi2            h           OmegaM           OmegaK               Acceptrate\n";
      std::cout << std::setw(15) << chi2 << " ";
      out       << std::setw(15) << chi2 << " ";
      for(int i = 0; i < nparam; i++){
        std::cout << std::setw(15) << parameters[i] << " ";
        out << std::setw(15) << parameters[i] << " ";
      }
      out << "\n";
      std::cout << std::setw(15) << " " << nsample/double(steps)*100.0 << "%\n";
      
      // Record new best-fit 
      if(chi2 < chi2_min){
        chi2_min = chi2;
        best_parameters = parameters;
      }
    }
  }

  // Print best-fit
  std::cout << "Minimum chi^2 found " << chi2_min << " ";
  for(int i = 0; i < nparam; i++)
    std::cout << best_parameters[i] << " ";
  std::cout << "\n";
}
