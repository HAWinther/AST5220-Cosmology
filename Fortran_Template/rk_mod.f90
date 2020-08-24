module rk_mod
  use healpix_types
  implicit none

contains

  subroutine rk4(y, dydx, x, h, yout, derivs)
    implicit none

    real(dp), dimension(:), intent(in)  :: y, dydx
    real(dp),               intent(in)  :: x, h
    real(dp), dimension(:), intent(out) :: yout

    interface
       subroutine derivs(x, y, dydx)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    integer(i4b) :: ndum
    real(dp)     :: h6, hh, xh
    real(dp), dimension(size(y)) :: dym, dyt, yt

    ndum = size(y)
    hh   = 0.5d0*h
    h6   = h/6.d0
    xh = x+hh
    yt = y + hh*dydx
    call derivs(xh, yt, dyt)
    yt = y + hh*dyt
    call derivs(xh, yt, dym)
    yt = y + h*dym
    dym = dyt + dym
    call derivs(x+h, yt, dyt)
    yout = y + h6*(dydx+dyt+2.d0*dym)

  end subroutine rk4

  subroutine rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs)
    implicit none

    real(dp), dimension(:), intent(inout) :: y
    real(dp), dimension(:), intent(in)    :: dydx, yscal
    real(dp),               intent(inout) :: x
    real(dp),               intent(in)    :: htry, eps
    real(dp),               intent(out)   :: hdid, hnext

    interface
       subroutine derivs(x, y, dydx)
         use healpix_types
         implicit none

         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    integer(i4b) :: ndum
    real(dp)     :: errmax, h, htemp, xnew
    real(dp), dimension(size(y)) :: yerr, ytemp
    real(dp), parameter :: SAFETY=0.9d0, PGROW=-0.2d0, PSHRNK=-0.25d0, &
         & ERRCON=1.89d-4
    
    ndum = size(y)
    h    = htry
    open(39,file='conv.dat', recl=1000)
    do 
       call rkck(y, dydx, x, h, ytemp, yerr, derivs)
       errmax = maxval(abs(yerr(:)/yscal(:)))/eps
       write(39,*) 
       write(39,*) real(y,sp)
       write(39,*) real(abs(yerr(:)/yscal(:))/eps,sp)
!       write(*,*) real(y,sp)
!       write(*,*) real(abs(yerr(:)/yscal(:))/eps,sp)
       if (errmax <= 1.d0) exit
       htemp = SAFETY*h*(errmax**PSHRNK)
       h = sign(max(abs(htemp),0.1d0*abs(h)),h)
       xnew  = x+h
       if (xnew == x) then
          write(*,*) 'Stepsize underflow in rkqs'
          stop
       end if
    end do
    close(39)
!    write(*,*) h
!    write(*,*) 
!    if (h < 1.d-7) stop
    if (errmax > ERRCON) then
       hnext = SAFETY*h*(errmax**PGROW)
    else
       hnext = 5.d0*h
    end if
    hdid = h
    x    = x+h
    y(:) = ytemp(:)
!    write(*,*) 'hnext = ', hnext

  end subroutine rkqs

  subroutine rkck(y, dydx, x, h, yout, yerr, derivs)
    implicit none

    real(dp), dimension(:), intent(in)  :: y, dydx
    real(dp),               intent(in)  :: x, h
    real(dp), dimension(:), intent(out) :: yout, yerr
    
    interface
       subroutine derivs(x, y, dydx)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    integer(i4b) :: ndum
    real(dp), dimension(size(y)) :: ak2, ak3, ak4, ak5, ak6, ytemp
    real(dp), parameter :: A2=0.2d0, A3=0.3d0, A4=0.6d0, A5=1.d0, &
         & A6=0.875d0, B21=0.2d0, B31=3.d0/40.d0, B32=9.d0/40.d0, &
         & B41=0.3d0, B42=-0.9d0, B43=1.2d0, B51=-11.d0/54.d0, &
         & B52=2.5d0, B53=-70.d0/27.d0, B54=35.d0/27.d0, &
         & B61=1631.d0/55296.d0, B62=175.d0/512.d0, &
         & B63=575.d0/13824.d0, B64=44275.d0/110592.d0, &
         & B65=253.d0/4096.d0, C1=37.d0/378.d0, &
         & C3=250.d0/621.d0, C4=125.d0/594.d0, &
         & C6=512.d0/1771.d0, DC1=C1-2825.d0/27648.d0, &
         & DC3=C3-18575.d0/48384.d0, DC4=C4-13525.d0/55296.d0, &
         & DC5=-277.d0/14336.d0, DC6=C6-0.25d0
    
    ndum = size(y)
    ytemp = y+B21*h*dydx
    call derivs(x+A2*h, ytemp, ak2)
    ytemp = y+h*(B31*dydx+B32*ak2)
    call derivs(x+A3*h, ytemp, ak3)
    ytemp = y+h*(B41*dydx+B42*ak2+B43*ak3)
    call derivs(x+A4*h, ytemp, ak4)
    ytemp = y+h*(B51*dydx+B52*ak2+B53*ak3*B54*ak4)
    call derivs(x+A5*h, ytemp, ak5)
    ytemp = y+h*(B61*dydx+B62*ak2+B63*ak3*B64*ak4*B65*ak5)
    call derivs(x+A6*h, ytemp, ak6)
    
    yout = y + h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr = h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)

  end subroutine rkck


end module rk_mod
