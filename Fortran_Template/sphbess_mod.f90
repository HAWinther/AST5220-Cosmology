module sphbess_mod
  use healpix_types
  use spline_1D_mod
  implicit none

contains

  subroutine sphbes(n, x, sj, sy, sjp, syp)
    implicit none

    integer(i4b), intent(in)  :: n
    real(dp),     intent(in)  :: x
    real(dp),     intent(out), optional :: sj, sy, sjp, syp

    real(dp), parameter :: RTPIO2 = 1.253314137315500d0
    real(dp) :: factor, order, rj, rjp, ry, ryp

    order = n + 0.5d0
    call bessjy(x, order, rj, ry, rjp, ryp)
    factor = RTPIO2 / sqrt(x)
    if (present(sj))  sj  = factor * rj
    if (present(sy))  sy  = factor * ry
    if (present(sjp)) sjp = factor * rjp - sj/(2.d0*x)
    if (present(syp)) syp = factor * ryp - sy/(2.d0*x)

  end subroutine sphbes

  subroutine bessjy(x, xnu, rj, ry, rjp, ryp)
    implicit none

    real(dp), intent(in)            :: x, xnu
    real(dp), intent(out), optional :: rj, ry, rjp, ryp

    integer(i4b), parameter :: MAXIT = 10000
    real(dp),     parameter :: XMIN = 2.d0, eps = 1.d-10, FPMIN=1.d-30

    integer(i4b) :: i, isign, l, nl
    real(dp)     :: a, b, c, d, del, del1, e, f, fact, fact2, fact3, ff, gam, gam1, gam2
    real(dp)     :: gammi, gampl, h, p, pimu, pimu2, q, r, rjl, rjl1, rjmu, rjp1
    real(dp)     :: rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1, w, x2, xi
    real(dp)     :: xi2, xmu, xmu2
    complex(dpc) :: aa, bb, cc, dd, dl, pq

    nl = merge(int(xnu+0.5d0), max(0,int(xnu-x+1.5d0)), x < XMIN)
    xmu   = xnu-nl
    xmu2  = xmu*xmu
    xi    = 1.d0 / x
    xi2   = 2.d0 * xi
    w     = xi2/pi
    isign = 1
    h     = xnu*xi
    if (h < FPMIN) h = FPMIN
    b     = xi2*xnu
    d     = 0.d0
    c     = h
    do i = 1, MAXIT
       b = b + xi2
       d = b-d
       if (abs(d) < FPMIN) d = FPMIN
       c = b - 1.d0/c
       if (abs(c) < FPMIN) c = FPMIN
       d = 1.d0 / d
       del = c*d
       h = del*h
       if (d < 0.d0) isign = -isign
       if (abs(del-1.d0) < EPS) exit
    end do
    if (i > MAXIT) then
       write(*,*) "bessjy: x too large, try asymptotic expansion, x = ", x
       stop
    end if
    rjl = isign*FPMIN
    rjpl = h*rjl
    rjl1 = rjl
    rjp1 = rjpl
    fact = xnu*xi
    do l = nl, 1, -1
       rjtemp = min(fact*rjl + rjpl,1.d300)
       fact = fact-xi
       rjpl = min(fact*rjtemp-rjl,1.d300)
       rjl = rjtemp
    end do
    if (rjl == 0.d0) rjl = EPS
    f = rjpl/rjl
    if (x < XMIN) then
       x2 = 0.5d0 * x
       pimu = pi*xmu
       if (abs(pimu) < EPS) then
          fact = 1.d0
       else
          fact = pimu / sin(pimu)
       end if
       d = -log(x2)
       e = xmu*d
       if (abs(e) < EPS) then
          fact2 = 1.d0
       else
          fact2 = sinh(e)/e
       end if
       call beschb(xmu, gam1, gam2, gampl, gammi)
       ff = 2.d0/pi * fact * (gam1*cosh(e)+gam2*fact2*d)
       e = exp(e)
       p = e/(gampl*pi)
       q = 1.d0 / (e*pi*gammi)
       pimu2 = 0.5d0 * pimu
       if (abs(pimu2) < EPS) then
          fact3 = 1.d0
       else
          fact3 = sin(pimu2) / pimu2
       end if
       r = pi * pimu2*fact3*fact3
       c = 1.d0
       d = -x2*x2
       sum = ff+r*q
       sum1 = p
       do i = 1, MAXIT
          ff = (i*ff+p+q) / (i*i-xmu2)
          c = c*d/i
          p = p/(i-xmu)
          q = q/(i+xmu)
          del = c*(ff+r*q)
          sum = sum + del
          del1 = c*p-i*del
          sum1 = sum1 + del1
          if (abs(del) < (1.d0+abs(sum))*EPS) exit
       end do
       if (i > MAXIT) then
          write(*,*) 'bessy series failed to converge'
          stop
       end if
       rymu = -sum
       ry1 =-sum1*xi2
       rymup = xmu*xi*rymu-ry1
       rjmu = w/(rymup-f*rymu)
    else
       a = 0.25d0 - xmu2
       pq = cmplx(-0.5d0*xi, 1.d0, kind=dpc)
       aa = cmplx(0.d0, xi*a, kind=dpc)
       bb = cmplx(2.d0*x,2.d0, kind=dpc)
       cc = bb + aa/pq
       dd = 1.d0 / bb
       pq = cc*dd*pq
       do i = 2, MAXIT
          a = a+2*(i-1)
          bb = bb + cmplx(0.d0, 2.d0, kind=dpc)
          dd = a*dd + bb
          if (absc(dd) < FPMIN) dd = FPMIN
          cc = bb + a/cc
          if (absc(cc) < FPMIN) cc = FPMIN
          dd = 1.d0 / dd
          dl = cc*dd
          pq = pq * dl
          if (absc(dl-1.d0) < EPS) exit
       end do
       if (i > MAXIT) then
          write(*,*) 'cf2 failed in bessjy'
          stop
       end if
       p = real(pq)
       q = aimag(pq)
       gam = (p-f)/q
       rjmu = sqrt(w/((p-f)*gam+q))
       rjmu = sign(rjmu,rjl)
       rymu = rjmu * gam
       rymup = rymu * (p+q/gam)
       ry1 = xmu*xi*rymu - rymup
    end if
    fact = rjmu / rjl
    rj =  rjl1 * fact
    rjp = rjp1 * fact
    do i = 1, nl
       rytemp = (xmu+i)*xi2*ry1-rymu
       rymu = ry1
       ry1 = rytemp
    end do
    ry = rymu
    ryp = xnu*xi*rymu - ry1
    
  contains

    function absc(z)
      implicit none
      
      complex(dpc), intent(in) :: z
      real(dp)                 :: absc

      absc = abs(real(z)) + abs(aimag(z))
    end function absc
  end subroutine bessjy

  subroutine beschb(x, gam1, gam2, gampl, gammi)
    implicit none

    real(dp), intent(in)  :: x
    real(dp), intent(out) :: gam1, gam2, gampl, gammi

    integer(i4b), parameter :: NUSE1=5, NUSE2=5

    real(dp) :: xx
    real(dp), dimension(7) :: c1 = (/1.142022680371168d0, &
         & 6.5165112670737d-3, 3.087090173086d-4, -3.4706269649d-6, &
         & 6.9437664d-9, 3.67795d-11, -1.356d-13 /)
    real(dp), dimension(8) :: c2 = (/1.843740587300905d0, &
         & -7.68528408447867d-2, 1.2719271366546d-3, &
         & -4.9717367042d-6, -3.31261198d-8, 2.423096d-10, &
         & -1.702d-13, -1.49d-15 /)
    
    xx = 8.d0*x*x-1.d0
    gam1 = chebev(-1.d0, 1.d0, c1(1:NUSE1), xx)
    gam2 = chebev(-1.d0, 1.d0, c2(1:NUSE2), xx)
    gampl = gam2-x*gam1
    gammi = gam2+x*gam1

  end subroutine beschb

  function chebev(a,b,c,x)
    implicit none

    real(dp),               intent(in) :: a, b, x
    real(dp), dimension(:), intent(in) :: c
    real(dp)                           :: chebev

    integer(i4b) :: j, m
    real(dp)     :: d, dd, sv, y, y2
    if ((x-a)*(x-b) > 0.d0) then
       write(*,*) 'x not in range in chebev'
       stop
    end if
    m = size(c)
    d = 0.d0
    dd = 0.d0
    y = (2.d0*x-a-b) / (b-a)
    y2 = 2.d0 * y
    do j = m, 2, -1
       sv = d
       d = y2*d-dd+c(j)
       dd = sv
    end do
    chebev = y*d-dd+0.5d0*c(1)
  end function chebev

end module sphbess_mod
