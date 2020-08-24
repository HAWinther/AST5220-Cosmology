module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop
    logical(lgt) :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H

    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstopo

    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))

    ! Task: Fill in x (rec) grid



    ! Task: Compute X_e and n_e at all grid times
    use_saha = .true.
    do i = 1, n
       if (use_saha) then
          ! Use the Saha equation

          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation

       end if
    end do


    ! Task: Compute splined (log of) electron density function


    ! Task: Compute optical depth at all grid points


    ! Task: Compute splined (log of) optical depth
    ! Task: Compute splined second derivative of (log of) optical depth


    ! Task: Compute splined visibility function
    ! Task: Compute splined second derivative of visibility function


  end subroutine initialize_rec_mod


  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e

  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau

  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

  end function get_ddg


end module rec_mod
