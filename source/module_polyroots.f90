module polyroots

   !-----------------------------------------------------------------------
   ! Subroutines for solving cubic equation with real coefficients
   !
   ! Solve a cubic equation where a, b, c, and d are real.
   !   a*x**3 + b*x**2 + c*x + d = 0
   !
   ! Variables used:
   !   a, b, c, d  ... coefficients (input)
   !   x(i)        ... three (generally) complex solutions (output)
   !   nroot       ... number of roots
   !
   ! Progamming Note: same as CUBICEQ1.FOR, except with subroutines
   !-----------------------------------------------------------------------
   ! Instructor: Nam Sun Wang
   !-----------------------------------------------------------------------

   implicit none
   private

   real(kind=8), parameter :: pi=3.14159265358979d0

   public :: quad
   public :: cubic

contains

   !----------------------------------------------------------------------
   subroutine quad(a, b, c, x, nroot)
   !-------------------------------------------------------------------
   ! Solve a quadratic equation where a, b, and c are real.
   !   a*x*x + b*x + c = 0
   !
   ! Public Variables
   !   a, b, c     ... coefficients (input)
   !   x(i)        ... two complex solutions (output)
   !   nroot       ... number of roots (output)
   !-------------------------------------------------------------------

      real(kind=8), intent(in) :: a, b, c
      complex(kind=8), intent(out) :: x(2)
      integer, intent(out) :: nroot

      real(kind=8) :: dd

      nroot = 0
      x = cmplx(0.d0,0.d0)

      if (a .eq. 0.d0) then

         if (b .eq. 0.d0) then
            !-- We have a non-equation; therefore, we have no valid solution
            nroot = 0
         else
            !-- We have a linear equation with 1 root.
            nroot = 1
            x(1) = cmplx(-c/b, 0.d0)
         endif

      else

         !-- We have a true quadratic equation.  Apply the quadratic formula to find two roots

         nroot = 2
         dd = b*b - 4.d0*a*c
         if (dd .ge. 0.d0) then
            x(1) = cmplx((-b+sqrt(dd))/2.d0/a, 0.d0)
            x(2) = cmplx((-b-sqrt(dd))/2.d0/a, 0.d0)
         else
            x(1) = cmplx(-b/2.d0/a, +sqrt(-dd)/2.d0/a)
            x(2) = cmplx(-b/2.d0/a, -sqrt(-dd)/2.d0/a)
         endif

      endif

   end subroutine quad


   !----------------------------------------------------------------------
   subroutine cubic(a, b, c, d, x, nroot)
   !----------------------------------------------------------------------
   ! Solve a cubic equation where a, b, c, and d are real.
   !   a*x**3 + b*x**2 + c*x + d = 0
   !
   ! Public Variables
   !   a, b, c, d  ... coefficients (input)
   !   x(i)        ... three (generally) complex solutions (output)
   !   nroot       ... number of roots (output)
   ! Local Variables:
   !   y1, y2, y3  ... three transformed solutions
   !
   ! Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
   !   (McGraw Hill, 1978).
   !   Step 0: If a is 0. use the quadratic formula to avoid dividing by 0.
   !   Step 1: Calculate p and q
   !           p = ( 3*c/a - (b/a)**2 ) / 3
   !           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
   !   Step 2: Calculate discriminant D
   !           D = (p/3)**3 + (q/2)**2
   !   Step 3: Depending on the sign of D, we follow different strategy.
   !           If D<0, three distinct real roots.
   !           If D=0, three real roots of which at least two are equal.
   !           If D>0, one real and two complex roots.
   !   Step 3a: For D>0 and D=0,
   !           Calculate u and v
   !           u = cubic_root(-q/2 + sqrt(D))
   !           v = cubic_root(-q/2 - sqrt(D))
   !           Find the three transformed roots
   !           y1 = u + v
   !           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
   !           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
   !   Step 3b Alternately, for D<0, a trigonometric formulation is more convenient
   !           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
   !           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
   !           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
   !           where phi = acos(-q/2/sqrt(|p|**3/27))
   !                 pi  = 3.141592654...
   !   Step 4  Finally, find the three roots
   !           x = y - b/a/3
   !----------------------------------------------------------------------

      real(kind=8), intent(in) :: a, b, c, d
      complex(kind=8), intent(out) :: x(3)
      integer, intent(out) :: nroot

      real(kind=8) :: dd, p, q, phi, temp1, temp2, y1, y2, y3, u, v, y2r, y2i

      nroot = 0
      x = cmplx(0.d0,0.d0)

      !-- Step 0: If a is 0 use the quadratic formula
      if (a .eq. 0.d0) then
         call quad(b, c, d, x(1:2), nroot)
         return
      endif

      !-- Cubic equation with 3 roots

      nroot = 3

      !-- Step 1: Calculate p and q

      p  = c/a - b*b/a/a/3.d0
      q  = (2.d0*b*b*b/a/a/a - 9.d0*b*c/a/a + 27.d0*d/a) / 27.d0

      !-- Step 2: Calculate DD (discriminant)
      DD = p*p*p/27.d0 + q*q/4.d0

      !-- Step 3: Branch to different algorithms based on DD

      if (DD.lt.0.d0) then

         !-- Step 3b:
         !-- 3 real unequal roots -- use the trigonometric formulation

         phi = acos(-q/2.d0/sqrt(abs(p*p*p)/27.d0))
         temp1=2.d0*sqrt(abs(p)/3.d0)
         y1 =  temp1*cos(phi/3.d0)
         y2 = -temp1*cos((phi+pi)/3.d0)
         y3 = -temp1*cos((phi-pi)/3.d0)

      else

         !-- Step 3a:
         !-- 1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)

         temp1 = -q/2.d0 + sqrt(DD)
         temp2 = -q/2.d0 - sqrt(DD)
         u = abs(temp1)**(1.d0/3.d0)
         v = abs(temp2)**(1.d0/3.d0)
         if(temp1 .lt. 0.d0) u = -u
         if(temp2 .lt. 0.d0) v = -v
         y1  = u + v
         y2r = -(u + v)/2.d0
         y2i =  (u - v)*sqrt(3.d0)/2.d0

      endif

      !-- Step 4: Final transformation -----------------------------------------

      temp1 = b/a/3.d0
      y1  =  y1 - temp1
      if (DD.lt.0.d0) then
         y2  =  y2 - temp1
         y3  =  y3 - temp1
      else
         y2r = y2r - temp1
      endif

      !-- Assign answers -------------------------------------------------------

      if (DD .lt. 0.d0) then
         x(1) = cmplx( y1,  0.d0)
         x(2) = cmplx( y2,  0.d0)
         x(3) = cmplx( y3,  0.d0)
      elseif (DD .eq. 0.d0) then
         x(1) = cmplx( y1,  0.d0)
         x(2) = cmplx(y2r,  0.d0)
         x(3) = cmplx(y2r,  0.d0)
      else
         x(1) = cmplx( y1,  0.d0)
         x(2) = cmplx(y2r, y2i)
         x(3) = cmplx(y2r,-y2i)
      endif

   end subroutine cubic

end module polyroots

