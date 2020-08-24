import numpy as np

class ConstantsAndUnits:
  """
  Very simple units class that holds units and physical constants in those units
  SI:              good old SI units
  ParticlePhysics: the fundamental unit is eV and c=hbar=kb=1
  Planck:          c=hbar=kb=G=1
  ReducedPlanck:   c=hbar=kb=8piG=1
  Cosmology:       we use lengths in Mpc, c=1 and masses in Msun
  User:            supply the size of your units relative to SI units

  Input Parameters: 
    typeOfUnits (string): The type (SI, ParticlePhysics, Planck, ReducedPlanck, Cosmology, User) units to use
    UserLengthinMeter  (float): Optional and irrelevant unless typeOfUnits = User 
    UserTimeinSec      (float): Optional and irrelevant unless typeOfUnits = User
    UserMassinkg       (float): Optional and irrelevant unless typeOfUnits = User
    UserTempinKelvin   (float): Optional and irrelevant unless typeOfUnits = User
    UserChargeinColumb (float): Optional and irrelevant unless typeOfUnits = User

  Attributes:
    Fundamental constants:
    G                  (float): Newtons constant
    c                  (float): Speed of light
    k_b                (float): Boltzmanns constant
    k_e                (float): Coloumbs constant
    hbar               (float): Plancks (reduced) constant
    
    Measured quantities:
    Msun               (float): Mass of the sun
    m_e                (float): Electron mass
    m_H                (float): Hydrogen mass
    sigma_T            (float): Thompson cross section
    epsilon_0          (float): Hydrogen ionization energy (positive)
    xhi0               (float): Helium0 ionization energy
    xhi1               (float): Helium+ ionization energy
    lambda_2s1s        (float): Hydrogen 2s->1s decay rate
    H0_over_h          (float): The Hubble parameter 100km/s/Mpc without little 'h'

    m                  (float): how many user length units there are in one meter
    s                  (float): how many user time units there are in one second
    kg                 (float): how many user mass units there are in one kilo
    K                  (float): how many user temperature units there are in one Kelvin
    Co                 (float): how many user charge units there are in one Columb
    J                  (float): how many user energy units there are in one Joule
    N                  (float): how many user force  units there are in one Newton
    pc                 (float): how many user length units there are in one parsec 
    kpc                (float): how many user length units there are in one kiloparsec 
    Mpc                (float): how many user length units there are in one megaparsec 
    Gpc                (float): how many user length units there are in one gigaparsec 
    eV                 (float): how many user energy units there are in one electronvolt
    km                 (float): how many user length units there are in one kilometer
  
  Functions:
    length_to_SI       (float): Convert user length to SI (m)
    time_to_SI         (float): Convert user time to SI (s)
    mass_to_SI         (float): Convert user mass to SI (kg)
    temperature_to_SI  (float): Convert user temperature to SI (K)
    velocity_to_SI     (float): Convert user unit velocity to SI (m/s)
  """

  # Physical constants in SI units
  _G_SI             = 6.67430e-11
  _c_SI             = 2.99792458e8 
  _k_b_SI           = 1.38064852e-23
  _k_e_SI           = 8.9875517923e9
  _hbar_SI          = 1.054571817e-34
  _eV_SI            = 1.60217653e-19
  _m_e_SI           = 9.10938356e-31
  _m_H_SI           = 1.6735575e-27
  _Msun_SI          = 1.98847e30
  _sigma_T_SI       = 6.6524587158e-29
  _lambda_2s1s_SI   = 8.227
  _epsilon_0_eV     = 13.605693122994
  _xhi0_eV          = 24.587387
  _xhi1_eV          = 4.0 * _epsilon_0_eV
  _pc_SI            = 3.08567758e16
  
  def __init__(self, typeOfUnits = "SI", 
      UserLengthinMeter  = 1.0, 
      UserTimeinSec      = 1.0, 
      UserMassinkg       = 1.0, 
      UserTempinKelvin   = 1.0, 
      UserChargeinColumb = 1.0):
      
    # Set the base units (meter, seconds, kilogram, Kelvin and Coloumb)
    # m,s,kg,K,Co is the value of your base units in SI
    self.name = typeOfUnits
    if (typeOfUnits == "SI"):
      # Normal SI units
      self.m  = 1.0
      self.s  = 1.0
      self.kg = 1.0
      self.K  = 1.0
      self.Co = 1.0
    elif (typeOfUnits == "Planck" or typeOfUnits == "ReducedPlanck"):
      # Planck units where hbar=c=G=kb=1 or 8piG=1 in Reduced Planck units
      factor = (8.0*np.pi) if typeOfUnits == "ReducedPlanck" else 1.0
      self.m  = np.sqrt((self._c_SI**3) / (self._hbar_SI * self._G_SI)) / factor**0.5
      self.s  = self._c_SI * self.m
      self.kg = self._G_SI * self.m**3 / self.s**2 * factor
      self.K  = self._k_b_SI * self.kg * self.m**2 / self.s**2
      self.Co = np.sqrt(self._k_e_SI * self.kg * self.m**3 / self.s**2)
    elif (typeOfUnits == "ParticlePhysics"):
      # Particle physics units c=hbar=kb=1 with eV = 1
      self.m  = self._eV_SI / self._hbar_SI / self._c_SI
      self.s  = self._c_SI * self.m
      self.kg = 1.0 / self._hbar_SI * self.s / self.m**2
      self.K  = self._k_b_SI * self._G_SI * self.m**5 / self.s**2
      self.Co = 1.0
    elif (typeOfUnits == "Cosmology"):
      # Cosmology units with length in Mpc, c = 1 and masses in solar-masses
      self.m  = 1.0/(1e6 * self._pc_SI)
      self.s  = self._c_SI * self.m
      self.kg = 1.0 / self._Msun_SI
      self.K  = 1.0
      self.Co = 1.0
    elif (typeOfUnits == "User"):
      self.m  = 1.0/UserLengthinMeter
      self.s  = 1.0/UserTimeinSec
      self.kg = 1.0/UserMassinkg
      self.K  = 1.0/UserTempinKelvin
      self.Co = 1.0/UserChargeinColumb
    else:
      error  = "The units [" + typeOfUnits +"] is not recognized. "
      error += "Expected: SI, Planck, ReducedPlanck, ParticlePhysics, Cosmology or User"
      raise ValueError(error)

    # Derived units
    self.km          = 1e3 * self.m                   # Kilometers
    self.N           = self.kg*self.m/self.s**2       # Newton
    self.J           = self.N*self.m                  # Joule
    self.Mpc         = 1e6 * self._pc_SI*self.m		    # Megaparsec
    self.eV          = self._eV_SI*self.J		          # Electronvolt
    
    # Physical constants in the desired units   
    self.k_b         = self._k_b_SI * self.J/self.K	            # Bolzmanns constant
    self.m_e         = self._m_e_SI * self.kg	                  # Mass of electron
    self.m_H         = self._m_H_SI * self.kg	                  # Mass of hydrogen atom
    self.c           = self._c_SI * self.m/self.s	              # Speed of light
    self.G           = self._G_SI * self.N*(self.m/self.kg)**2	# Gravitational constant
    self.hbar        = self._hbar_SI * self.J*self.s		        # Reduced Plancks constant
    self.sigma_T     = self._sigma_T_SI * self.m**2    	        # Thomas scattering cross-section
    self.lambda_2s1s = self._lambda_2s1s_SI / self.s            # Transition time is
    self.H0_over_h   = 100. * self.km/self.s/self.Mpc           # H0 / h
    self.epsilon_0   = self._epsilon_0_eV * self.eV	            # Ionization energy for the ground state of hydrogen
    self.xhi0        = self._xhi0_eV * self.eV		              # Ionization energy for neutral Helium
    self.xhi1        = self._xhi1_eV * self.eV		              # Ionization energy for singly ionized Helium
  
  def info(self):
    """
    Print some useful info about the units
    """
    print("")
    print("Units [" + self.name + "]:")
    print("Plancks constant hbar:  ", self.hbar)
    print("Boltzmann constant k_b: ", self.k_b)
    print("Newtons constant G:     ", self.G)
    print("Speed of light c:       ", self.c)
    print("Unit length in m:  ", 1.0/self.m)
    print("Unit time   in s:  ", 1.0/self.s)
    print("Unit mass   in kg: ", 1.0/self.kg)
    print("Unit temp   in K:  ", 1.0/self.K)
    print("Unit charge in Co: ", 1.0/self.Co)

  # Conversion functions
  def velocity_to_SI(self,v):
    return v / (self.m/self.s)
  def length_to_SI(self,L):
    return L / (self.m)
  def time_to_SI(self,T):
    return T / (self.s)
  def mass_to_SI(self,M):
    return M / (self.kg)
  def energy_to_SI(self,E):
    return E / (self.J)
  def temperature_to_SI(self,T):
    return T / (self.K)

# Make a global constants 
const = ConstantsAndUnits("SI")
