import numpy as np
from   matplotlib import pyplot as plt
from   Global import const
import Global
import BackgroundCosmology
import RecombinationHistory
import Perturbations
import PowerSpectrum

"""

The whole project runs from this file. The task is split into four milestones
for which you have to implement four classes. Each of them has a solve, plot and info
method and methods to retrieve things you will compute. On top of that add what you want/need.
Do them in order and remove the exit after each milestone when you are ready to move on

Any physical constants and units you need can be found in [const]

"""

# Show plots below
show_plots = True

#============================================
# Set up the constants and units class 
#============================================
const = Global.ConstantsAndUnits("SI")
const.info()

#============================================
# Milestone 1: Solve the background
#============================================
cosmo = BackgroundCosmology.BackgroundCosmology(
    name      = "LCDM",              # Label 
    h         = 0.7,                 # Hubble parameter
    OmegaB    = 0.046,               # Baryon density
    OmegaCDM  = 0.224,               # CDM density
    OmegaK    = 0.0,                 # Curvature density parameter
    TCMB_in_K = 2.7255,              # Temperature of CMB today in Kelvin
    Neff      = 0.0)                 # Effective number of neutrinos

# Solve and plot
cosmo.info()
cosmo.solve()
if show_plots: cosmo.plot() 

# Remove when done with milestone
exit()

#============================================
# Milestone 2: Solve the recombination history
#============================================
rec = RecombinationHistory.RecombinationHistory(
    BackgroundCosmology  = cosmo,
    Yp                   = 0.24,     # Primordial helium fraction
    reionization         = True,     # Include reionization
    z_reion              = 11.0,     # Reionization redshift
    delta_z_reion        = 0.5,      # Reionization width
    helium_reionization  = True,     # Helium double reionization
    z_helium_reion       = 3.5,      # Helium double reionization redshift
    delta_z_helium_reion = 0.5)      # Helium double reionization width

# Solve and plot
rec.info()
rec.solve()
if show_plots: rec.plot()

# Remove when done with milestone
exit()

#============================================
# Milestone 3: Integrate the perturbations
#============================================
pert = Perturbations.Perturbations(
    BackgroundCosmology  = cosmo,
    RecombinationHistory = rec,
    n_ell_theta          = 10,       # Number of ells (0,1,...,n-1) to include in the Boltzmann hierarchy
    keta_max             = 1500.0,   # Set kmax based on keta0. 3000 typically enough for Cell, lower for testing
    npts_k               = 100)      # 100-200 typically enough for Cell, lower for testing

# Solve and plot
pert.info()
pert.solve()
if show_plots: pert.plot(pert.k_min)
if show_plots: pert.plot(pert.k_max)

# Remove when done with milestone
exit()

#============================================
# Milestone 4: PowerSpectrum evaluation
#============================================
power = PowerSpectrum.PowerSpectrum(
    BackgroundCosmology  = cosmo,
    RecombinationHistory = rec,
    Perturbations        = pert,
    kpivot_mpc           = 0.05,     # Pivot scale in 1/Mpc
    n_s                  = 0.96,     # Spectral index
    A_s                  = 2e-9,     # Primordial amplitude
    ell_max              = 1500)     # Maximum ell to compute Cell up to

# Solve and plot
power.info()
power.solve()
if show_plots: power.plot()

