import numpy as np
from   matplotlib import pyplot as plt
from   scipy.interpolate import CubicSpline
import scipy.integrate as integrate

from   Global import const

class BackgroundCosmology:
  """
  This is a class for the cosmology at the background level.
  It holds cosmological parameters and functions relevant for the background.
  
  Input Parameters: 
    h           (float): The little Hubble parameter h in H0 = 100h km/s/Mpc
    OmegaB      (float): Baryonic matter density parameter at z = 0
    OmegaCDM    (float): Cold dark matter density parameter at z = 0
    OmegaK      (float,optional): Curvative density parameter at z = 0
    name        (float,optional): A name for describing the cosmology
    TCMB        (float,optional): The temperature of the CMB today in Kelvin. Fiducial value is 2.725K
    Neff        (float,optional): The effective number of relativistic neutrinos

  Attributes:    
    OmegaR      (float): Radiation matter density parameter at z = 0
    OmegaNu     (float): Massless neutrino density parameter at z = 0
    OmegaM      (float): Total matter (CDM+b+mnu) density parameter at z = 0
    OmegaK      (float): Curvature density parameter at z = 0
  
  Functions:
    eta_of_x             (float->float) : Conformal time times c (units of length) as function of x=log(a) 
    H_of_x               (float->float) : Hubble parameter as function of x=log(a) 
    dHdx_of_x            (float->float) : First derivative of hubble parameter as function of x=log(a) 
    Hp_of_x              (float->float) : Conformal hubble parameter H*a as function of x=log(a) 
    dHpdx_of_x           (float->float) : First derivative of conformal hubble parameter as function of x=log(a) 
  """
  
  # Settings for integration and splines of eta
  x_start = np.log(1e-8)
  x_end = np.log(1.0)
  n_pts_splines = 1000

  def __init__(self, h = 0.7, OmegaB = 0.046, OmegaCDM = 0.224, OmegaK = 0.0, 
      name = "FiducialCosmology", TCMB_in_K = 2.725, Neff = 3.046):
    self.OmegaB      = OmegaB
    self.OmegaCDM    = OmegaCDM
    self.OmegaK      = OmegaK
    self.h           = h
    self.H0          = const.H0_over_h * h
    self.name        = name
    self.TCMB        = TCMB_in_K * const.K
    self.Neff        = Neff
 
    # Set the constants
    # XXX TODO XXX
    self.rhoc0       = 1.0 # Critical density today
    self.OmegaR      = 1.0 # Radiation
    self.OmegaNu     = 1.0 # Neutrino radiation
    self.OmegaM      = 1.0 # Total matter
    self.OmegaRtot   = 1.0 # Total radiation
    self.OmegaLambda = 1.0 # Dark energy (from Sum Omega_i = 1)
  
  #=========================================================================
  # Methods availiable after solving
  #=========================================================================
  
  def eta_of_x(self,x):
    if not hasattr(self, 'eta_of_x_spline'):
      raise NameError('The spline eta_of_x_spline has not been created') 
    return self.eta_of_x_spline(x)
  def H_of_x(self,x):
    # XXX TODO XXX
    return 1.0
  def Hp_of_x(self,x):
    # XXX TODO XXX
    return 1.0
  def dHdx_of_x(self,x):
    # XXX TODO XXX
    return 1.0
  def dHpdx_of_x(self,x):
    # XXX TODO XXX
    return 1.0

  #=========================================================================
  #=========================================================================
  #=========================================================================
  
  def info(self):
    """
    Print some useful info about the class
    """
    print("")
    print("Background Cosmology:")
    print("OmegaB:        %8.7f" % self.OmegaB)
    print("OmegaCDM:      %8.7f" % self.OmegaCDM)
    print("OmegaLambda:   %8.7f" % self.OmegaLambda)
    print("OmegaR:        %8.7e" % self.OmegaR)
    print("OmegaNu:       %8.7e" % self.OmegaNu)
    print("OmegaK:        %8.7f" % self.OmegaK)
    print("TCMB (K):      %8.7f" % (self.TCMB / const.K))
    print("h:             %8.7f" % self.h)
    print("H0:            %8.7e" % self.H0)
    print("H0 (km/s/Mpc): %8.7f" % (self.H0 / (const.km / const.s / const.Mpc)))
    print("Neff:          %8.7f" % self.Neff)
    print("OmegaM:        %8.7f" % self.OmegaM)
    print("OmegaRtot:     %8.7e" % self.OmegaRtot)
  
  def solve(self):
    """
    Main driver for all the solving.
    For LCDM we only need to solve for the conformal time eta(x)
    """

    # Make scale factor array (logspaced)
    # XXX TODO XXX
    # x_array  = np.linspace(self.x_start, self.x_end, num = self.n_pts_splines)
    # eta_array = np.zeros(self.n_pts_splines)

    # Compute and spline conformal time eta = Int_0^t dt/a = Int da/(a^2 H(a)) =  Int dx/[ exp(x) * H(exp(x)) ] where x = log a
    # XXX TODO XXX

    # Spline up result
    # XXX TODO XXX
    # self.eta_of_x_spline = CubicSpline(x_array, eta_array)
  
  def plot(self):
    """
    Plot some useful quantities
    """

    # Plot eta(x) * Hp(x) as a test
    npts = 2000
    xarr = np.linspace(self.x_start, self.x_end, num = npts)
    eta  = [self.eta_of_x(xarr[i])*self.H_of_x(xarr[i])*np.exp(xarr[i])/const.c for i in range(npts)]
    fac  = [self.dHpdx_of_x(xarr[i])/self.Hp_of_x(xarr[i]) for i in range(npts)]
    plt.title('EtaHp/c and dlogHpdx')
    plt.plot(xarr, eta)
    plt.plot(xarr, fac)
    plt.show()
  
  #=========================================================================
  #=========================================================================
  #=========================================================================
  

