module ode_solver
  use healpix_types
  use rk_mod
  use bs_mod
  implicit none

  integer(i4b)                          :: nok, nbad, kount
  logical(lgt),                 save    :: save_steps = .false.
  real(dp)                              :: dxsav
  real(dp),     dimension(:),   pointer :: xp
  real(dp),     dimension(:,:), pointer :: yp

contains
  
  subroutine odeint(ystart, x1, x2, eps, h1, hmin, derivs, rkqs, output) 
    implicit none

    real(dp), dimension(:), intent(inout) :: ystart
    real(dp),               intent(in)    :: x1, x2, eps, h1, hmin

    interface
       subroutine derivs(x, y, dydx)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs

       subroutine output(x, y)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
       end subroutine output

       subroutine rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs)
         use healpix_types
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
       end subroutine rkqs
    end interface
    
    real(dp),     parameter :: TINY=1.d-30
    integer(i4b), parameter :: MAXSTP=1000000

    integer(i4b) :: nstp, i
    real(dp)     :: h, hdid, hnext, x, xsav
    real(dp), dimension(size(ystart)) :: dydx, y, yscal

    real(dp), dimension(0:6) :: theta
    
    x     = x1
    h     = sign(h1,x2-x1)
    nok   = 0
    nbad  = 0
    kount = 0
    y(:)  = ystart(:)
    
    nullify(xp, yp)
    
    if (save_steps) then
       xsav=x-2.d0*dxsav
       allocate(xp(256))
       allocate(yp(size(ystart),size(xp)))
    end if

    do nstp = 1, MAXSTP
!       write(*,*) 'new step', x
       call derivs(x, y, dydx)

       yscal(:) = abs(y(:)) + abs(h*dydx(:)) + TINY
!       if (size(y) > 7) then
!          yscal(4:size(y)) = 1.d0
!       end if

       if (save_steps .and. (abs(x-xsav) > abs(dxsav))) call save_a_step
       if ((x+h-x2)*(x+h-x1) > 0.d0) h = x2-x
       call rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, derivs)
       if (hdid == h) then
          nok = nok+1
       else
          nbad = nbad+1
       end if
       if ((x-x2)*(x2-x1) >= 0.d0) then
          ystart(:) = y(:)
          if (save_steps) call save_a_step
!          if (size(y)>7) write(71,*) real(x,sp), nstp, h
!          write(*,*) 'ferdig!', real(x,sp), nstp, h
          return
       end if
       if (abs(hnext) < hmin) then
          write(*,*) 'Stepsize smaller than minimum in odeint'
          stop
       end if
       h = hnext
!       if (size(y) == 7) then
!          do i = 1, 6
!             write(57+i,*) real(x,sp), abs(real(y(i),sp))
!          end do
!          theta(0:1) = y(6:7)
!          theta(2:6) = 0.d0
!          write(71,*) real(x,sp), abs(real(theta,sp))
!       end if
!       if (size(y) > 7) then
!          do i = 1, 6
!             write(57+i,*) real(x,sp), abs(real(y(i),sp))
!          end do
!          write(71,*) real(x,sp), abs(real(y(6:12),sp))
!          write(71,*) real(x,sp)
!       end if
!       if (size(y) >= 7) then
!          call output(x, y)
!       end if
    end do
       
!    do i = 1, 7
!       close(57+i)
!    end do

    write(*,*) 'Too many steps in odeint'
    stop

  contains
    
    subroutine save_a_step
      implicit none

      kount = kount+1
      if (kount > size(xp)) then
         xp => reallocate1(xp,2*size(xp))
         yp => reallocate2(yp,size(y,1),size(xp))
      end if
      xp(kount)   = x
      yp(:,kount) = y(:)
      xsav        = x
    end subroutine save_a_step

  end subroutine odeint

  function reallocate1(p, n)
    implicit none

    real(dp),     dimension(:), pointer             :: p, reallocate1
    integer(i4b),                        intent(in) :: n
    
    integer(i4b) :: nold, ierr
    
    allocate(reallocate1(n),stat=ierr)
    if (ierr /= 0) then
       write(*,*) 'Reallocate:Problem in attempt to allocate memory'
       stop
    end if
    if (.not. associated(p)) return
    nold = size(p)
    reallocate1(1:min(nold,n)) = p(1:min(nold,n))
    deallocate(p)

  end function reallocate1

  function reallocate2(p, n, m)
    implicit none

    real(dp),     dimension(:,:), pointer            :: p, reallocate2
    integer(i4b),                         intent(in) :: n, m
    
    integer(i4b) :: nold, mold, ierr
    
    allocate(reallocate2(n,m),stat=ierr)
    if (ierr /= 0) then
       write(*,*) 'Reallocate:Problem in attempt to allocate memory'
       stop
    end if
    if (.not. associated(p)) return
    nold = size(p(:,1))
    mold = size(p(1,:))
    reallocate2(1:min(nold,n),1:min(mold,m)) = p(1:min(nold,n),1:min(mold,m))
    deallocate(p)
    
  end function reallocate2

end module ode_solver
