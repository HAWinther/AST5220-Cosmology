module bs_mod
  use healpix_types
  implicit none

  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

contains

  subroutine mmid(y,dydx,xs,htot, nstep, yout, derivs)
    implicit none
    
    integer(i4b),               intent(in) :: nstep
    real(dp),                   intent(in) :: xs, htot
    real(dp),     dimension(:), intent(in)  :: y, dydx
    real(dp),     dimension(:), intent(out) :: yout

    interface
       subroutine derivs(x,y,dydx)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    integer(i4b) :: n, ndum
    real(dp)     :: h, h2, x
    real(dp), dimension(size(y)) :: ym, yn

    ndum = size(y)
    h    = htot/nstep
    ym   = y
    yn   = y+h*dydx
    x    = xs+h
    call derivs(x,yn,yout)
    h2   = 2.d0*h
    do n = 2, nstep
       call swap(ym,yn)
       yn = yn+h2*yout
       x  = x+h
       call derivs(x,yn,yout)
    end do
    yout = 0.5d0*(ym+yn+h*yout)

  end subroutine mmid

  subroutine bsstep(y, dydx, x, htry, eps, yscal, hdid, hnext, derivs)
    implicit none

    real(dp), dimension(:), intent(inout) :: y
    real(dp), dimension(:), intent(in)    :: dydx, yscal
    real(dp),               intent(inout) :: x
    real(dp),               intent(in)    :: htry, eps
    real(dp),               intent(out)   :: hdid, hnext

    interface
       subroutine derivs(x,y,dydx)
         use healpix_types
         implicit none
         real(dp),               intent(in)  :: x
         real(dp), dimension(:), intent(in)  :: y
         real(dp), dimension(:), intent(out) :: dydx
       end subroutine derivs
    end interface

    integer(i4b), parameter :: IMAX=9, KMAXX=IMAX-1
    real(dp),     parameter :: SAFE1=0.25d0, SAFE2=0.7d0, REDMAX=1.d-5, &
         & REDMIN=0.7d0, TINY=1.d-30, SCALMX=0.1d0

    integer(i4b) :: k, km, ndum
    integer(i4b), dimension(IMAX) :: nseq = (/ 2,4,6,8,10,12,14,16,18 /)
    integer(i4b), save :: kopt, kmax
    real(dp), dimension(KMAXX,KMAXX), save :: alf
    real(dp), dimension(KMAXX) :: err
    real(dp), dimension(IMAX), save :: a
    real(dp), save :: epsold = -1.d0, xnew
    real(dp) :: eps1, errmax, fact, h, red, scale, wrkmin, xest
    real(dp), dimension(size(y)) :: yerr, ysav, yseq
    logical(lgt) :: reduct
    logical(lgt), save :: first=.true.

    ndum = size(y)
    if (eps /= epsold) then
       hnext = -1.d29
       xnew  = -1.d29
       eps1  = SAFE1*eps
       a(:)  = cumsum(nseq,1)
       where (upper_triangle(KMAXX, KMAXX)) alf=eps1** &
            & (outerdiff(a(2:),a(2:))/outerprod(arth( &
            & 3.d0,2.d0,KMAXX),(a(2:)-a(1)+1.d0)))
       epsold=eps
       do kopt=2, KMAXX-1
          if (a(kopt+1) > a(kopt)*alf(kopt-1,kopt)) exit
       end do
       kmax = kopt
    end if
    h=htry
    ysav(:) = y(:)
    if (h /= hnext .or. x /= xnew) then
       first = .true.
       kopt = kmax
    end if
    reduct = .false.
    main_loop: do
       do k = 1, kmax
          xnew = x+h
          if (xnew == x) then
             write(*,*) 'step siez underflow in bsstep'
             stop
          end if
          call mmid(ysav,dydx, x, h, nseq(k), yseq, derivs)
          xest = (h/nseq(k))**2
          call pzextr(k, xest, yseq, y, yerr)
          if (k /= 1) then
             errmax = maxval(abs(yerr(:)/yscal(:)))
             errmax = max(TINY,errmax) / eps
             km = k-1
             err(km) = (errmax/SAFE1)**(1.d0/(2*km+1))
          end if
          if (k /= 1 .and. (k >= kopt-1 .or. first)) then
             if (errmax < 1.d0) exit main_loop
             if (k == kmax .or. k == kopt+1) then
                red = SAFE2/err(km)
                exit
             else if (k == kopt) then
                if (alf(kopt-1,kopt) < err(km)) then
                   red = 1.d0 / err(km)
                   exit
                end if
             else if (kopt == kmax) then
                if (alf(km,kmax-1) < err(km)) then
                   red = alf(km,kmax-1)*SAFE2/err(km)
                   exit
                end if
             else if (alf(km,kopt) < err(km)) then
                red=alf(km,kopt-1)/err(km)
                exit
             end if
          end if
       end do
       red = max(min(red,REDMIN),REDMAX)
       h = h*red
       reduct = .true.
    end do main_loop
    x=xnew
    hdid=h
    first=.false.
    kopt = 1+iminloc(a(2:km+1)*max(err(1:km),SCALMX))
    scale = max(err(kopt-1), SCALMX)
    wrkmin=scale*a(kopt)
    hnext = h/scale
    if (kopt >= k .and. kopt /= kmax .and. .not. reduct) then
       fact = max(scale/alf(kopt-1,kopt),SCALMX)
       if (a(kopt+1)*fact <= wrkmin) then
          hnext = h/fact
          kopt  = kopt+1
       end if
    end if
  end subroutine bsstep

  subroutine pzextr(iest, xest, yest, yz, dy)
    implicit none
    
    integer(i4b), intent(in) :: iest
    real(dp),     intent(in) :: xest
    real(dp), dimension(:), intent(in) :: yest
    real(dp), dimension(:), intent(out) :: yz, dy

    integer(i4b), parameter :: IEST_MAX=16
    integer(i4b) :: j, nv
    integer(i4b), save :: nvold=-1
    real(dp)     :: delta, f1, f2
    real(dp), dimension(size(yz)) :: d, tmp, q
    real(dp), dimension(IEST_MAX), save :: x
    real(dp), dimension(:,:), allocatable, save :: qcol
    
    nv = size(yz)
    if (iest > IEST_MAX) then
       write(*,*) 'pzextr: probable misuse, too much extrapolation.'
       stop
    end if
    if (nv /= nvold) then
       if (allocated(qcol)) deallocate(qcol)
       allocate(qcol(nv,IEST_MAX))
       nvold = nv
    end if
    x(iest) = xest
    dy(:) = yest(:)
    yz(:) = yest(:)
    if (iest == 1) then
       qcol(:,1) = yest(:)
    else
       d(:) = yest(:)
       do j = 1, iest-1
          delta = 1.d0 / (x(iest-j)-xest)
          f1 = xest*delta
          f2 = x(iest-j)*delta
          q(:) = qcol(:,j)
          qcol(:,j) = dy(:)
          tmp(:) = d(:)-q(:)
          dy(:) = f1*tmp(:)
          d(:) = f2*tmp(:)
          yz(:) = yz(:)+dy(:)
       end do
       qcol(:,iest) = dy(:)
    end if
  end subroutine pzextr



  subroutine swap(a, b)
    implicit none

    real(dp), dimension(:), intent(inout) :: a, b

    integer(i4b) :: i
    real(dp) :: c
    do i = 1, size(a)
       c    = a(i)
       a(i) = b(i)
       b(i) = c
    end do

  end subroutine swap

  RECURSIVE FUNCTION cumsum(arr,seed) RESULT(ans)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: seed
    INTEGER(I4B), DIMENSION(size(arr)) :: ans
    INTEGER(I4B) :: n,j,sd
    n=size(arr)
    if (n == 0_i4b) RETURN
    sd=0_i4b
    if (present(seed)) sd=seed
    ans(1)=arr(1)+sd
    if (n < NPAR_CUMSUM) then
       do j=2,n
          ans(j)=ans(j-1)+arr(j)
       end do
    else
       ans(2:n:2)=cumsum(arr(2:n:2)+arr(1:n-1:2),sd)
       ans(3:n:2)=ans(2:n-1:2)+arr(3:n:2)
    end if
  END FUNCTION cumsum

  function upper_triangle(n,m) result(A)
    implicit none

    integer(i4b) :: n, m
    logical(lgt), dimension(n,m) :: A

    integer(i4b) :: i, j

    A = .false.
    do i = 1, n
       do j = i+1, n
          A(i,j) = .true.
       end do
    end do
    
  end function upper_triangle


  FUNCTION outerdiff(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff
    outerdiff = spread(a,dim=2,ncopies=size(b)) - &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff

  FUNCTION outerprod(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod
    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod

  FUNCTION arth(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth(k)=arth(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth(k)=arth(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth

  FUNCTION iminloc(arr)
    REAL(DP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B), DIMENSION(1) :: imin
    INTEGER(I4B) :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc


end module bs_mod
