module minima_1d

   !--(  To get d1mach, mail netlib send d1mach from core )

   !============================================================
   !  some standard routines for search of the extremum points
   !  of one-dimensional functions
   !------------------------------------------------------------

   implicit none
   private

   public :: funmin

contains

   function funmin(fun,ax,bx,tol) result(xmin)

      real(kind=8), external :: fun
      real(kind=8), intent(in) :: ax, bx, tol
      real(kind=8) :: xmin

      !  An approximation  x  to the point where fun(x) attains a minimum on
      !  the interval (ax,bx) is determined.
      !
      !  Input..
      !
      !  fun   function subprogram which evaluates  f(x)  for any  x
      !        in the interval  (ax,bx)
      !  ax    left endpoint of initial interval
      !  bx    right endpoint of initial interval
      !  tol   desired length of the interval of uncertainty of the final
      !        result (.ge.0.)
      !
      !  Output..
      !
      !  funmin  abcissa approximating the point where  fun  attains a minimum
      !
      !      The method used is a combination of  golden  section  search  and
      !  successive parabolic interpolation.  convergence is never much slower
      !  than  that  for  a  fibonacci search. If fun  has a continuous second
      !  derivative which is positive at the minimum (which is not  at  ax  or
      !  bx),  then  convergence  is  superlinear, and usually of the order of
      !  about  1.324....
      !      The function fun is never evaluated at two points closer together
      !  than  eps*abs(funmin)+(tol/3), where eps is approximately the  square
      !  root  of  the  relative  machine  precision.   if  fun  is a unimodal
      !  function and the computed values of  fun  are  always  unimodal  when
      !  separated  by  at least  eps*abs(x)+(tol/3), then funmin approximates
      !  the abcissa of the global minimum of fun on the interval (ax,bx) with
      !  an error less than  3*eps*abs(funmin)+tol.  If fun is  not  unimodal,
      !  then funmin may approximate a local, but perhaps non-global, minimum
      !  to the same accuracy.
      !      This function subprogram is a slightly modified  version  of  the
      !  algol  60 procedure  localmin  given in Richard Brent, Algorithms for
      !  minimization without derivatives, Prentice-Hall, Inc. (1973).
      !

      real(kind=8) :: a, b, c, d, e, eps, xm, p, q, r
      real(kind=8) :: tol1, t2, u, v, w, fu, fv, fw, fx, x, tol3

      !-- c is the squared inverse of the golden ratio
      c = 0.5d0*(3.0d0 - dsqrt(5.0d0))

      !--  eps is approximately the square root of the relative machine precision.

   10 eps=d1mach(4)
      tol1=eps+1.0d0
      eps=dsqrt(eps)

      a=ax
      b=bx
      v=a+c*(b-a)
      w=v
      x=v
      e=0.0d0
      fx=fun(x)
      fv=fx
      fw=fx
      tol3=tol/3.0d0

      !-- main loop starts here

   20 xm=0.5d0*(a+b)
      tol1=eps*dabs(x)+tol3
      t2=2.0d0*tol1

      !-- check stopping criterion

      if (dabs(x-xm).le.(t2-0.5d0*(b-a))) goto 190
      p=0.0d0
      q=0.0d0
      r=0.0d0
      if (dabs(e).le.tol1) goto 50

      !-- fit parabola

      r=(x-w)*(fx-fv)
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=2.0d0*(q-r)
      if (q.le.0.0d0) goto 30
      p=-p
      goto 40
   30 q=-q
   40 r=e
      e=d
   50 if ((dabs(p).ge.dabs(0.5d0*q*r)).or.(p.le.q*(a-x)).or.(p.ge.q*(b-x))) goto 60

      !-- a parabolic-interpolation step

      d=p/q
      u=x+d

      !-- fun must not be evaluated too close to ax or bx

      if (((u-a).ge.t2).and.((b-u).ge.t2)) goto 90
      d=tol1
      if (x.ge.xm) d=-d
      goto 90

      !-- a golden-section step

   60 if (x.ge.xm) goto 70
      e=b-x
      goto 80
   70 e=a-x
   80 d=c*e

      !-- fun must not be evaluated too close to x

   90 if (dabs(d).lt.tol1) goto 100
      u=x+d
      goto 120
  100 if (d.le.0.0d0) goto 110
      u=x+tol1
      goto 120
  110 u=x-tol1
  120 fu=fun(u)

      !-- update a, b, v, w, and x

      if (fx.gt.fu) goto 140
      if (u.ge.x) goto 130
      a=u
      goto 140
  130 b=u
  140 if (fu.gt.fx) goto 170
      if (u.ge.x) goto 150
      b=x
      goto 160
  150 a=x
  160 v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
      goto 20
  170 if ((fu.gt.fw).and.(w.ne.x)) goto 180
      v=w
      fv=fw
      w=u
      fw=fu
      goto 20
  180 if ((fu.gt.fv).and.(v.ne.x).and.(v.ne.w)) goto 20
      v=u
      fv=fu
      goto 20

      !-- end of main loop

  190 xmin=x

   end function funmin

!------------------------------------------------------------------------------------

   function d1mach(i) result(value)

      implicit none
      integer :: i
      real(kind=8) :: value
      real(kind=8) :: b, x

      !   D1MACH can be used to obtain machine-dependent parameters for the
      !   local machine environment.  It is a function subprogram with one
      !   (input) argument, and can be referenced as follows:
      !
      !        A = D1MACH(I)
      !
      !   where I=1,...,5.  The (output) value of A above is determined by
      !   the (input) value of I.  The results for various values of I are
      !   discussed below.
      !
      !   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
      !   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
      !   D1MACH(3) = B**(-T), the smallest relative spacing.
      !   D1MACH(4) = B**(1-T), the largest relative spacing.
      !   D1MACH(5) = LOG10(B)
      !
      !   Assume single precision numbers are represented in the T-digit,
      !   base-B form
      !
      !              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
      !
      !   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
      !   EMIN .LE. E .LE. EMAX.
      !
      !   The values of B, T, EMIN and EMAX are provided in I1MACH as
      !   follows:
      !   I1MACH(10) = B, the base.
      !   I1MACH(11) = T, the number of base-B digits.
      !   I1MACH(12) = EMIN, the smallest exponent E.
      !   I1MACH(13) = EMAX, the largest exponent E.
      !
      !
      !***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
      !                 a portable library, ACM Transactions on Mathematical
      !                 Software 4, 2 (June 1978), pp. 177-188.

      x = 1.0d0
      b = radix(x)

      select case (i)
        case (1)
          value = b**(minexponent(x)-1) ! the smallest positive magnitude.
        case (2)
          value = huge(x)               ! the largest magnitude.
        case (3)
          value = b**(-digits(x))       ! the smallest relative spacing.
        case (4)
          value = b**(1-digits(x))      ! the largest relative spacing.
        case (5)
          value = log10(b)
        case default
          write (*,'("ERROR 1 in d1mach - i out of bounds")')
          stop
      end select

   end function d1mach

!------------------------------------------------------------------------------------
end module minima_1d
