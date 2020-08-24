import numpy as np
import scipy.integrate as integrate
from   scipy.interpolate import CubicSpline
from   scipy.integrate import solve_ivp
from   matplotlib import pyplot as plt
from   scipy.interpolate import RectBivariateSpline

from   Global import const
import BackgroundCosmology
import RecombinationHistory

class Perturbations:
  """
  This is a class for solving the perturbation history of the Universe
  After solving it holds functions for all perturbations and source functions
  
  Input Parameters: 
    cosmo (BackgroundCosmology) : The cosmology we use to integrate perturbations
    rec   (RecombinationHistory): The recombination history we use to integrate perturbations
    keta_max             (float): The maximum k*eta(x=0) to integrate up to
    npts_k               (int)  : Number of k-values from kmin to kmax
    n_ell_theta          (int)  : Number of theta multipoles to include when solving the ODE

  Functions:
    deltaCDM             (k float,x float->float) : CDM density contrast as function of wavenumber k and x = log(a) 
    deltaB               (k float,x float->float) : Baryon density contrast as function of wavenumber k and x = log(a)  
    vCDM                 (k float,x float->float) : CDM velocity as function of wavenumber k and x = log(a)  
    vB                   (k float,x float->float) : Baryon velocity as function of wavenumber k and x = log(a)  
    Phi                  (k float,x float->float) : gii=(1+2Phi)a^2 metric potential as function of wavenumber k and x = log(a)  
    Psi                  (k float,x float->float) : g00=-(1+2Psi) metric potential as function of wavenumber k and x = log(a)  
    sourceT              (k float,x float->float) : Temperature source function as function of wavenumber k and x = log(a)  
    Theta                (k float,x float, ell int->float) : Photon multipoles Theta_ell (ell=0,1,2) as function of wavenumber k and x = log(a) 
  """

  # Settings x-integration
  x_start = np.log(1e-8)
  x_end   = np.log(1.0)
  npts_x  = 1000
  
  def __init__(self, BackgroundCosmology, RecombinationHistory, keta_max = 3000.0, npts_k = 100, n_ell_theta = 10):
    """
    Intitialize the object
    """
    self.cosmo  = BackgroundCosmology
    self.rec    = RecombinationHistory
    self.k_min  = 1.0 / self.cosmo.eta_of_x(0.0)
    self.k_max  = keta_max / self.cosmo.eta_of_x(0.0)
    self.npts_k = npts_k

    # Number of ells to include and the total number of quantities in the ODE system
    self.n_ell_theta = n_ell_theta
    self.n_tot_tight = 5 + 2
    self.n_tot_full  = 5 + n_ell_theta
  
    # A system to have control over where each 
    # quantity is in the ODE arrays (optional to use)
    # Theta_ell at index index_theta + ell. 
    # These things work in both tight and full ODE
    self.index_deltaCDM = 0
    self.index_vCDM     = 1
    self.index_deltaB   = 2
    self.index_vB       = 3
    self.index_Phi      = 4
    self.index_theta    = 5

    # PhD: you need to add polarization and neutrinos here and in the ODEs below
    # (and don't forget the differences between the tight coupling and full ODE)
    # ...
  
  #=========================================================================
  # Functions availiable after solving
  #=========================================================================
  
  def sourceT(self,k,x):
    if not hasattr(self, 'sourceT_spline'):
      raise NameError('The spline [sourceT_spline] has not been created') 
    return self.sourceT_spline(k,x)
  def deltaCDM(self,k,x):
    if not hasattr(self, 'deltaCDM_spline'):
      raise NameError('The spline [deltaCDM_spline] has not been created') 
    return self.deltaCDM_spline(k,x)
  def deltaB(self,k,x):
    if not hasattr(self, 'deltaB_spline'):
      raise NameError('The spline [deltaB_spline] has not been created') 
    return self.deltaB_spline(k,x)
  def vCDM(self,k,x):
    if not hasattr(self, 'vCDM_spline'):
      raise NameError('The spline [vCDM_spline] has not been created') 
    return self.vCDM_spline(k,x)
  def vB(self,k,x):
    if not hasattr(self, 'vB_spline'):
      raise NameError('The spline [vB_spline] has not been created') 
    return self.vB_spline(k,x)
  def Phi(self,k,x):
    if not hasattr(self, 'Phi_spline'):
      raise NameError('The spline [Phi_spline] has not been created') 
    return self.Phi_spline(k,x)
  def Psi(self,k,x):
    if not hasattr(self, 'Psi_spline'):
      raise NameError('The spline [Psi_spline] has not been created') 
    return self.Psi_spline(k,x)
  def Theta(self,k,x,ell):
    if(ell == 0): 
      if not hasattr(self, 'Theta1_spline'):
        raise NameError('The spline [Theta0_spline] has not been created') 
      return self.Theta0_spline(k,x)
    if(ell == 1): 
      if not hasattr(self, 'Theta1_spline'):
        raise NameError('The spline [Theta0_spline] has not been created') 
      return self.Theta1_spline(k,x)
    if(ell == 2): 
      if not hasattr(self, 'Theta2_spline'):
        raise NameError('The spline [Theta2_spline] has not been created') 
      return self.Theta2_spline(k,x)
    print("Error Theta has not been splined for ell = ", ell)
    return 0.0
  
  #=========================================================================
  #=========================================================================
  #=========================================================================

  def info(self):
    """
    Print some useful info
    """
    print("")
    print("Perturbations:")
    print("kmin (1/Mpc): ", self.k_min * const.Mpc)
    print("kmax (1/Mpc): ", self.k_max * const.Mpc)
    print("Number of ell's: ", self.n_ell_theta)
    print("Number of k-points: ", self.npts_k)

  def solve(self):
    """
    The main driver for doing all the solving 
    """
    self.integrate_perturbations()
  
  def plot(self,kval):
    """
    Plot the perturbations and source function as function of x = log(a) for a single value of k
    """
    x_array = np.linspace(self.x_start, self.x_end, self.npts_x)
    kmpc = '{:.3g}'.format(kval * const.Mpc)

    # Fetch data from splines
    data_deltaCDM = self.deltaB  (kval, x_array)[0,:]
    data_deltaB   = self.deltaB  (kval, x_array)[0,:]  
    data_vCDM     = self.vCDM    (kval, x_array)[0,:]
    data_vB       = self.vB      (kval, x_array)[0,:]
    data_Phi      = self.Phi     (kval, x_array)[0,:]
    data_Psi      = self.Psi     (kval, x_array)[0,:]
    data_Theta0   = self.Theta   (kval, x_array, 0)[0,:]
    data_Theta1   = self.Theta   (kval, x_array, 1)[0,:]
    
    plt.yscale('log')
    plt.title('Density perturbations k = ' + str(kmpc) + '/ Mpc')
    plt.plot(x_array, np.abs(data_deltaB),label='deltaB')
    plt.plot(x_array, data_deltaCDM,label='deltaCDM')
    plt.plot(x_array, np.abs(3.0*data_Theta0),label='deltaR')
    plt.legend()
    plt.show()

    plt.yscale('log')
    plt.title('Velocity perturbations k = ' + str(kmpc) + '/ Mpc')
    plt.plot(x_array, np.abs(data_vB),label='|vB|')
    plt.plot(x_array, np.abs(data_vCDM),label='|vCDM|')
    plt.plot(x_array, np.abs(-3.0*data_Theta1),label='|vR|')
    plt.show()

    plt.title('Potentials k = ' + str(kmpc) + '/ Mpc')
    plt.plot(x_array, data_Phi,label='Phi')
    plt.plot(x_array, data_Psi,label='Psi')
    plt.show()
    
    # Plot the temperature source for the min and max k
    plt.title('Source function k = ' + str(kmpc) + '/ Mpc')
    sourceT_data = self.sourceT_spline(kval, x_array)[0,:]
    plt.plot(x_array,sourceT_data,label='Source')
    plt.show()
  
  #=========================================================================
  #=========================================================================
  #=========================================================================

  def integrate_perturbations(self):
    """
    Integrate the tight coupling and the full system
    Combine the data from both and use this to make splines
    for all the quantities. Use this to make splines of the
    source function(s)
    """

    # Set up k-array
    k_array = np.exp( np.linspace(np.log(self.k_min), np.log(self.k_max), num=self.npts_k) )

    print("")
    print("Start integrating perturbations")

    # Set up large x-array
    x_array = np.linspace(self.x_start, self.x_end, self.npts_x)

    # 2D arrays to store the data
    # XXX TODO XXX
    # e.g. deltaCDM_data = np.zeros((self.npts_k, self.npts_x))

    # Loop over all k-values
    for ik in range(self.npts_k):
      # Current k value
      k = k_array[ik]
      self.k_current = k

      # Compute tight coupling time and set up x-arrays
      x_tc_end = self.get_x_end_tight_coupling(k)
      x_tight  = np.linspace(self.x_start, x_tc_end,   npts_tight)
      x_full   = np.linspace(x_tc_end,     self.x_end, self.npts_x - npts_tight)
      
      print("Current k:", k * const.Mpc, " x: ", x_tc_end)

      # Compute IC for tight coupling
      y_tc_ic = self.set_ic(self.x_start, k)

      # Solve the tight coupling ODE
      sol_tight = solve_ivp(self.rhs_tight_coupling, [self.x_start, x_tc_end], y_tc_ic, t_eval=x_tight, rtol=1e-6, atol=1e-6)
     
      # Get solution after tight coupling
      y_after_tight_coupling = sol_tight.y[:,-1]
      
      # Set initial conditions
      y_full_ic = self.set_ic_after_tight_coupling(x_tc_end, k, y_after_tight_coupling)

      # Solve the full ODE
      sol_full = solve_ivp(self.rhs_full, [x_tc_end, self.x_end], y_full_ic, t_eval=x_full, rtol=1e-6, atol=1e-6)
      
      # Compute quantities not solved for in the tight coupling regime
      # XXX TODO XXX

      # Combine arrays from the two regimes and store the data
      # XXX TODO XXX
      # e.g. cur_deltaCDM = np.concatenate((sol_tight.y[self.index_deltaCDM], sol_full.y[self.index_deltaCDM]))

      # Store the data
      # XXX TODO XXX
      # e.g. deltaCDM_data[ik,:] = cur_deltaCDM

    # Make 2D splines
    # XXX TODO XXX
    # self.deltaCDM_spline = RectBivariateSpline(k_array, x_array, deltaCDM_data)
    # self.vCDM_spline     = RectBivariateSpline(k_array, x_array, vCDM_data    )
    # self.deltaB_spline   = RectBivariateSpline(k_array, x_array, deltaB_data  )
    # self.vB_spline       = RectBivariateSpline(k_array, x_array, vB_data      )
    # self.Phi_spline      = RectBivariateSpline(k_array, x_array, Phi_data     )
    # self.Psi_spline      = RectBivariateSpline(k_array, x_array, Psi_data     )
    # self.Theta0_spline   = RectBivariateSpline(k_array, x_array, Theta0_data  )
    # self.Theta1_spline   = RectBivariateSpline(k_array, x_array, Theta1_data  )
    # self.Theta2_spline   = RectBivariateSpline(k_array, x_array, Theta2_data  )
    # self.Pi_spline       = RectBivariateSpline(k_array, x_array, Pi_data      )
   
    # Compute and spline the source-function (milestone 4)
    # XXX TODO XXX
    # self.sourceT_spline = RectBivariateSpline(k_array, x_array, sourceT_data)

    return

  def set_ic(self,x,k):
    """ 
    Set IC for the tight coupling system at given x = log(a) and wavenumber k 
    """
    y = np.zeros(self.n_tot_tight)

    # Cosmological variables
    Hp        = self.cosmo.Hp_of_x(x)
    ckoverHp  = const.c * k / Hp;
    OmegaNu   = self.cosmo.OmegaNu
    OmegaRtot = self.cosmo.OmegaRtot
    f_nu      = OmegaNu/OmegaRtot

    # Set IC
    # XXX TODO XXX
    Psi                    =  -1.0/(1.5 + 2.0*f_nu/5.0)
    y[self.index_deltaCDM] =  1.0
    y[self.index_vCDM]     =  1.0
    y[self.index_deltaB]   =  1.0
    y[self.index_vB]       =  1.0
    y[self.index_Phi]      =  1.0                 
    y[self.index_theta+0]  =  1.0
    y[self.index_theta+1]  =  1.0

    return y

  def set_ic_after_tight_coupling(self,x,k,y):
    """ 
    Set IC for full system at given x = log(a) and wavenumber k 
    y is the solution we have from the tight coupling system
    """
    yfull = np.zeros(self.n_tot_full)
    
    # Set IC to the full regime
    # XXX TODO XXX

    return yfull

  def rhs_full(self,x,y):
    """ 
    Set the right hand side of the full ODE system dy/dx = RHS
    for a given value of x. The wavenumber k is set in the global variable k_current
    """
    
    # Current value of wavenumber k
    k = self.k_current

    # The array we are to fill and return
    dydx = np.zeros(self.n_tot_full)

    # Set the right hand side
    # XXX TODO XXX
    # e.g. dydx[self.index_deltaCDM] = ...

    return dydx

  def rhs_tight_coupling(self,x,y):
    """ 
    Set the right hand side of the tight coupling ODE system dy/dx = RHS
    for a given value of x. The wavenumber k is set in the global variable k_current
    """
    
    # Current value of wavenumber k
    k = self.k_current
    
    # The array we are to fill and return
    dydx = np.zeros(self.n_tot_tight)

    # Set the right hand side
    # XXX TODO XXX
    dydx[self.index_deltaCDM] = 1.0
    dydx[self.index_vCDM]     = 1.0
    dydx[self.index_deltaB]   = 1.0
    dydx[self.index_vB]       = 1.0
    dydx[self.index_Phi]      = 1.0
    dydx[self.index_theta+0]  = 1.0
    dydx[self.index_theta+1]  = 1.0
   
    return dydx

  def get_x_end_tight_coupling(self,k):
    """
    Compute when (x = log(a)) for when tight coupling ends for a given wavenumber k
    """

    # Return x when tight coupling ends
    # XXX TODO XXX
    return 1.0

