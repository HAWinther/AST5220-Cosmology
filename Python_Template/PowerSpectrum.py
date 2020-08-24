import numpy as np
from   matplotlib import pyplot as plt
from   scipy.special import spherical_jn
from   scipy.interpolate import CubicSpline
from   scipy.integrate import solve_ivp

from   Global import const
import BackgroundCosmology
import RecombinationHistory
import Perturbations

class PowerSpectrum:
  """
  This is a class for solving for power-spectra
  After solving it holds functions for C_ell and P(k)
  
  Input Parameters: 
    cosmo (BackgroundCosmology) : The cosmology we use to integrate perturbations
    rec   (RecombinationHistory): The recombination history we use to integrate perturbations
    pert  (Perturbations)       : The perturbations in the Universe
    kpivot_mpc          (float) : Pivot scale in unit of 1/Mpc
    A_s                 (float) : Primordial amplitude
    n_s                 (float) : Spectral index
    ell_max               (int) : Maximum ell for which we compute Cell for

  Attributes:
    kpivot (float): Pivot scale
  
  Functions:
    cell_TT                 (ell float->float): Temperature power-spectrum l(l+1)/2pi C_ell as function of ell 
    get_matter_power_spectrum (k float->float): Matter power-spectrum P(k) as function of wave-number k
  """
  
  def __init__(self, BackgroundCosmology, RecombinationHistory, Perturbations, kpivot_mpc = 0.05, n_s = 0.96, A_s = 2e-9, ell_max = 1500):
    
    self.cosmo = BackgroundCosmology
    self.rec   = RecombinationHistory
    self.pert  = Perturbations
    self.n_s    = n_s
    self.A_s    = A_s
    self.kpivot = kpivot_mpc / const.Mpc

    # The ells we compute Theta_ell (and then Cell) for with LOS integration
    ell_list = np.array([2, 3, 4, 5, 6, 7, 8, 10, 12, 15,  
        20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 
        120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 
        400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 
        900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,
        1900, 1950, 2000, 2150, 2300, 2450, 2600, 2750, 2900, 3000,
        3200, 3400, 3600, 3800, 4000])
    self.ells = ell_list[ell_list < ell_max]
    self.nells = len(self.ells)
    self.ell_max = ell_max

  #=========================================================================
  # Functions availiable after solving
  #=========================================================================
  
  def cell_TT(self,ell):
    """
    The CMB angular power-spectrum l(l+1)/2pi Cell in units on (muK)^2
    """
    if not hasattr(self, 'cell_TT_spline'):
      raise NameError('The spline [cell_TT_spline] has not been created') 
    return self.cell_TT_spline(ell)
  
  def matter_power_spectrum(self,k,x):
    """
    The matter power-spectrum at wavenumber k at time x = log(a)
    """
    # Compute and return P(k)
    # XXX TODO XXX
    return 1.0
  
  def primordial_power_spectrum(self,k):
    """
    The priordial power-spectrum Delta(k) in P(k) = 2pi^2/k^3 Delta(k)
    """
    return self.A_s * ( k / self.kpivot )**(self.n_s-1.0)

  #=========================================================================
  #=========================================================================
  #=========================================================================
  
  def info(self):
    """
    Print some useful info
    """
    print("")
    print("Powerspectrum:")
    print("A_s: ", self.A_s)
    print("n_s: ", self.n_s)
    print("kpivot (1/Mpc): ", self.kpivot * const.Mpc)
    print("ell_max: ", self.ell_max)
    
  def solve(self): 
    """
    Solve for the CMB power-spectrum
    1) Generate j_ell splines for all ells
    2) Do line of sight integration for all ells
    3) Compute Cell for all ells
    4) Spline it up
    """

    # Create splines of Bessel functions j_ell(.) needed below for all ells in self.ells
    # XXX TODO XXX

    # Set up a k-array to evaluate Theta_ell on
    # XXX TODO XXX

    # Solve for theta_ell(k) for all k in k_array for all ells in self.ells
    # XXX TODO XXX

    # Make splines of theta_ell(k) for all the ells
    # XXX TODO XXX

    # Integrate up to get Cell's for al the ells
    # XXX TODO XXX
    
    # Make spline of Cell
    # XXX TODO XXX
    #self.cell_TT_spline = CubicSpline(self.ells, cell)
  
  def plot(self):
    """
    Make plots of P(k) and Cell
    """
    # Make k-array 
    k_min = self.pert.k_min
    k_max = self.pert.k_max
    npts_k = 100
    k_array = np.exp(np.linspace(np.log(k_min), np.log(k_max), npts_k))

    # Plot matter power-spectrum today with k in 1/Mpc and P(k) in Mpc^3
    pofk = np.array([self.matter_power_spectrum(k,0.0) for k in k_array]).flatten()
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(k_array * const.Mpc, pofk / const.Mpc**3)
    plt.show()
    
    # Plot angular power-spectrum today
    plt.xscale('log')
    plt.yscale('log')
    ells = np.exp(np.linspace(np.log(2.),np.log(self.ell_max), 200))
    cells = self.cell_TT(ells)
    plt.plot(ells,cells)
    plt.show()
    return
  
  #=========================================================================
  #=========================================================================
  #=========================================================================

