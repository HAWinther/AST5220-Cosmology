module params
  use healpix_types
  implicit none

  ! Output options
  integer(i4b), parameter :: lmax = 1200

  ! Units
  real(dp), parameter :: eV  = 1.60217646d-19
  real(dp), parameter :: Mpc = 3.08568025d22

  ! Cosmological parameters
  real(dp), parameter :: Omega_b      = 0.046d0
  real(dp), parameter :: Omega_m      = 0.224d0
  real(dp), parameter :: Omega_r      = 8.3d-5
  real(dp), parameter :: Omega_nu     = 0.d0
  ! If you want neutrinos, use the two below
!  real(dp), parameter :: Omega_r      = 5.04d-5
!  real(dp), parameter :: Omega_nu     = 8.3d-5 - Omega_r
  real(dp), parameter :: Omega_lambda = 1.d0 - Omega_m - Omega_b - Omega_r - Omega_nu
  real(dp), parameter :: T_0          = 2.725d0
  real(dp), parameter :: n_s          = 1.d0
  real(dp), parameter :: A_s          = 1.d0
  real(dp), parameter :: h0           = 0.7d0
  real(dp), parameter :: H_0          = h0 * 100.d0 * 1.d3 / Mpc

  ! General constants
  real(dp), parameter :: c            = 2.99792458d8
  real(dp), parameter :: epsilon_0    = 13.605698d0 * eV
  real(dp), parameter :: m_e          = 9.10938188d-31
  real(dp), parameter :: m_H          = 1.673534d-27
  real(dp), parameter :: sigma_T      = 6.652462d-29
  real(dp), parameter :: G_grav       = 6.67258d-11
  real(dp), parameter :: rho_c        = 3.d0*H_0**2 / (8.d0*pi*G_grav)
  real(dp), parameter :: alpha        = 7.29735308d-3
  real(dp), parameter :: hbar         = 1.05457148d-34
  real(dp), parameter :: k_b          = 1.3806503d-23

end module params
