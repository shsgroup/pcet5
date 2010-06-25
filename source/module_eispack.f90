module eispack
!-------------------------------------------------------------------
!  diagonalization routines from EISPACK package
!-------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  module_eispack.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  module_eispack.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 19:38:42  souda
!  Initial PCET-4.0 Release
!
!
!-------------------------------------------------------------------
   implicit none
   private

   public  :: rs, treql, invers, slen
      
   !-------------------------------------------------------------------!
   contains

   !-------------------------------------------------------------------!
   pure subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)

      ! this subroutine calls the recommended sequence of
      ! subroutines from the eigensystem subroutine package (eispack)
      ! to find the eigenvalues and eigenvectors (if desired)
      ! of a real symmetric matrix.
      !
      ! on input
      !
      !    nm  must be set to the row dimension of the two-dimensional
      !    array parameters as declared in the calling program
      !    dimension statement.
      !
      !    n  is the order of the matrix  a.
      !
      !    a  contains the real symmetric matrix.
      !
      !    matz  is an integer variable set equal to zero if
      !    only eigenvalues are desired.  otherwise it is set to
      !    any non-zero integer for both eigenvalues and eigenvectors.
      !
      ! on output
      !
      !    w  contains the eigenvalues in ascending order.
      !
      !    z  contains the eigenvectors if matz is not zero.
      !
      !    ierr  is an integer output variable set equal to an error
      !       completion code described in the documentation for tqlrat
      !       and tql2.  the normal completion code is zero.
      !
      !    fv1  and  fv2  are temporary storage arrays.
      !
      ! questions and comments should be directed to burton s. garbow,
      ! mathematics and computer science div, argonne national laboratory
      !
      ! this version dated august 1983.

      integer, intent(in)    :: n, nm, matz
      real*8,  intent(inout) :: a(nm,n)
      integer, intent(out)   :: ierr
      real*8,  intent(out)   :: w(n), z(nm,n)
      real*8,  intent(inout) :: fv1(n), fv2(n)

      fv1 = 0.d0
      fv2 = 0.d0
      w   = 0.d0
      z   = 0.d0
      ierr = 0

      if (n .gt. nm) then

         ierr = 10*n

      else

         if (matz .eq. 0) then

            ! find eigenvalues only
            call tred1(nm,n,a,w,fv1,fv2)
            call tql1(n,w,fv1,ierr)

         else

            ! find both eigenvalues and eigenvectors
            call tred2(nm,n,a,w,fv1,z)
            call tql2(nm,n,w,fv1,z,ierr)

         endif
         
      endif

      return

      end subroutine rs

   !-------------------------------------------------------------------!
   pure function pythag(a,b) result(hypo)

      ! finds dsqrt(a**2+b**2) without overflow or destructive underflow

      real*8, intent(in) :: a, b
      real*8 :: hypo
      real*8 :: p, r, s, t, u

      p = dmax1(dabs(a),dabs(b))

      if (p .ne. 0.d0) then

         r = (dmin1(dabs(a),dabs(b))/p)**2
         do
            t = 4.0d0 + r
            if (t .eq. 4.0d0) exit
            s = r/t
            u = 1.0d0 + 2.0d0*s
            p = u*p
            r = (s/u)**2 * r
         enddo

      endif

      hypo = p
      return

   end function pythag

   !-------------------------------------------------------------------!
   pure subroutine tql1(n,d,e,ierr)

      ! this subroutine is a translation of the algol procedure tql1,
      ! num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
      ! wilkinson.
      ! handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
      !
      ! this subroutine finds the eigenvalues of a symmetric
      ! tridiagonal matrix by the ql method.
      !
      ! on input
      !
      !    n is the order of the matrix.
      !
      !    d contains the diagonal elements of the input matrix.
      !
      !    e contains the subdiagonal elements of the input matrix
      !      in its last n-1 positions.  e(1) is arbitrary.
      !
      !  on output
      !
      !    d contains the eigenvalues in ascending order.  if an
      !      error exit is made, the eigenvalues are correct and
      !      ordered for indices 1,2,...ierr-1, but may not be
      !      the smallest eigenvalues.
      ! 
      !     e has been destroyed.
      !
      !    ierr is set to
      !      zero       for normal return,
      !      j          if the j-th eigenvalue has not been
      !                 determined after 30 iterations.
      !
      ! calls pythag for  dsqrt(a*a + b*b) .
      !
      ! questions and comments should be directed to burton s. garbow,
      ! mathematics and computer science div, argonne national laboratory
      !
      ! this version dated august 1983.

      integer, intent(in) :: n
      real*8,  intent(inout), dimension(n) :: d, e
      integer, intent(out) :: ierr

      integer :: i, j, l, m, ii, l1, l2, mml
      real*8  :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2

      ierr = 0
      if (n .eq. 1) return

      do i = 2, n
         e(i-1) = e(i)
      enddo

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      do l = 1, n

         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h

         ! look for small sub-diagonal element ..........
         do m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) exit
            ! e(n) is always zero, so there is no exit
            ! through the bottom of the loop ..........
         enddo

         if (m .ne. l) then

            do

               if (j .eq. 30) then
                  ! set error - no convergence to an
                  ! eigenvalue after 30 iterations
                  ierr = l
                  return
               endif
  
               j = j + 1

               ! form shift

               l1 = l + 1
               l2 = l1 + 1
               g = d(l)
               p = (d(l1) - g) / (2.0d0 * e(l))
               r = pythag(p,1.0d0)
               d(l) = e(l) / (p + dsign(r,p))
               d(l1) = e(l) * (p + dsign(r,p))
               dl1 = d(l1)
               h = g - d(l)

               if (l2 .le. n) then
                  do i = l2, n
                     d(i) = d(i) - h
                  enddo
               endif

               f = f + h

               ! ql transformation
               p = d(m)
               c = 1.0d0
               c2 = c
               el1 = e(l1)
               s = 0.0d0
               mml = m - l

               ! for i=m-1 step -1 until l do

               do ii = 1, mml
                  c3 = c2
                  c2 = c
                  s2 = s
                  i = m - ii
                  g = c * e(i)
                  h = c * p
                  r = pythag(p,e(i))
                  e(i+1) = s * r
                  s = e(i) / r
                  c = p / r
                  p = c * d(i) - s * g
                  d(i+1) = h + s * (c * g + s * d(i))
               enddo

               p = -s * s2 * c3 * el1 * e(l) / dl1
               e(l) = s * p
               d(l) = c * p
               tst2 = tst1 + dabs(e(l))

               if (tst2 .le. tst1) exit
            
            enddo

         endif

         p = d(l) + f

         ! order eigenvalues

         if (l .eq. 1) go to 250

         ! for i=l step -1 until 2 do

         do ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
         enddo

  250    i = 1
  270    d(i) = p

      enddo

      return

   end subroutine tql1

   !-------------------------------------------------------------------!
   pure subroutine tql2(nm,n,d,e,z,ierr)

      ! this subroutine is a translation of the algol procedure tql2,
      ! num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
      ! wilkinson.
      ! handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
      !
      ! this subroutine finds the eigenvalues and eigenvectors
      ! of a symmetric tridiagonal matrix by the ql method.
      ! the eigenvectors of a full symmetric matrix can also
      ! be found if  tred2  has been used to reduce this
      ! full matrix to tridiagonal form.
      !
      ! on input
      !
      !    nm must be set to the row dimension of two-dimensional
      !      array parameters as declared in the calling program
      !      dimension statement.
      !
      !    n is the order of the matrix.
      !
      !    d contains the diagonal elements of the input matrix.
      !
      !    e contains the subdiagonal elements of the input matrix
      !      in its last n-1 positions.  e(1) is arbitrary.
      !
      !    z contains the transformation matrix produced in the
      !      reduction by  tred2, if performed.  if the eigenvectors
      !      of the tridiagonal matrix are desired, z must contain
      !      the identity matrix.
      !
      !  on output
      !
      !    d contains the eigenvalues in ascending order.  if an
      !      error exit is made, the eigenvalues are correct but
      !      unordered for indices 1,2,...,ierr-1.
      !
      !    e has been destroyed.
      !
      !    z contains orthonormal eigenvectors of the symmetric
      !      tridiagonal (or full) matrix.  if an error exit is made,
      !      z contains the eigenvectors associated with the stored
      !      eigenvalues.
      !
      !    ierr is set to
      !      zero       for normal return,
      !      j          if the j-th eigenvalue has not been
      !                 determined after 30 iterations.
      !
      ! calls pythag for  dsqrt(a*a + b*b) .
      !
      ! questions and comments should be directed to burton s. garbow,
      ! mathematics and computer science div, argonne national laboratory
      !
      ! this version dated august 1983.

      integer, intent(in) :: nm, n
      real*8,  intent(inout), dimension(n)    :: d, e
      real*8,  intent(inout), dimension(nm,n) :: z
      integer, intent(out) :: ierr

      integer :: i,j,k,l,m,ii,l1,l2,mml
      real*8 :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2

      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
  100 e(i-1) = e(i)

      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0

      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
         !.......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
            !.......... e(n) is always zero, so there is no exit
            !           through the bottom of the loop ..........
  110    continue

  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
         !.......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145

         do 140 i = l2, n
  140    d(i) = d(i) - h

  145    f = f + h
         !.......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         !.......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            !.......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue

  200    continue

         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
      !.......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)

         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue

         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p

         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue

  300 continue

      go to 1001
      !.......... set error -- no convergence to an
      !           eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return

   end subroutine tql2

   !-------------------------------------------------------------------!
   pure subroutine tred1(nm,n,a,d,e,e2)

      ! this subroutine is a translation of the algol procedure tred1,
      ! num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
      ! handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
      !
      ! this subroutine reduces a real symmetric matrix
      ! to a symmetric tridiagonal matrix using
      ! orthogonal similarity transformations.
      !
      ! on input
      !
      !    nm must be set to the row dimension of two-dimensional
      !      array parameters as declared in the calling program
      !      dimension statement.
      !
      !    n is the order of the matrix.
      !
      !    a contains the real symmetric input matrix.  only the
      !      lower triangle of the matrix need be supplied.
      !
      ! on output
      !
      !    a contains information about the orthogonal trans-
      !      formations used in the reduction in its strict lower
      !      triangle.  the full upper triangle of a is unaltered.
      !
      !    d contains the diagonal elements of the tridiagonal matrix.
      !
      !    e contains the subdiagonal elements of the tridiagonal
      !      matrix in its last n-1 positions.  e(1) is set to zero.
      !
      !    e2 contains the squares of the corresponding elements of e.
      !      e2 may coincide with e if the squares are not needed.
      !
      ! questions and comments should be directed to burton s. garbow,
      ! mathematics and computer science div, argonne national laboratory
      !
      ! this version dated august 1983.

      integer, intent(in) :: nm, n
      real*8,  intent(inout), dimension(nm,n) :: a
      real*8,  intent(out),   dimension(n)    :: d, e, e2

      integer :: i, j, k, l, ii, jp1
      real*8  :: f, g, h, scale

      do i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
      enddo

      !.......... for i=n step -1 until 1 do -- ..........
      do ii = 1, n

         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0

         if (l .lt. 1) go to 130

         !.......... scale row (algol tol then not needed) ..........
         do k = 1, l
            scale = scale + dabs(d(k))
         enddo

         if (scale .ne. 0.0d0) go to 140

         do j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
         enddo

  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         cycle

  140    do k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
         enddo

         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g

         if (l .ne. 1) then

            !.......... form a*u ..........
            do j = 1, l
               e(j) = 0.0d0
            enddo

            do j = 1, l

               f = d(j)
               g = e(j) + a(j,j) * f
               jp1 = j + 1

               if (l .ge. jp1) then
                  do k = jp1, l
                     g = g + a(k,j) * d(k)
                     e(k) = e(k) + a(k,j) * f
                  enddo
               endif

               e(j) = g

            enddo

            !.......... form p ..........
            f = 0.0d0

            do j = 1, l
               e(j) = e(j) / h
               f = f + e(j) * d(j)
            enddo

            h = f / (h + h)

            !.......... form q ..........
            do j = 1, l
               e(j) = e(j) - h * d(j)
            enddo

            !.......... form reduced a ..........
            do j = 1, l
               f = d(j)
               g = e(j)
               do k = j, l
                  a(k,j) = a(k,j) - f * e(k) - g * d(k)
               enddo
            enddo

         endif

         do j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
         enddo

      enddo

      return

   end subroutine tred1

   !-------------------------------------------------------------------!
   pure subroutine tred2(nm,n,a,d,e,z)

      ! this subroutine is a translation of the algol procedure tred2,
      ! num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
      ! handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
      !
      ! this subroutine reduces a real symmetric matrix to a
      ! symmetric tridiagonal matrix using and accumulating
      ! orthogonal similarity transformations.
      !
      ! on input
      !
      !    nm must be set to the row dimension of two-dimensional
      !      array parameters as declared in the calling program
      !      dimension statement.
      !
      !    n is the order of the matrix.
      !
      !    a contains the real symmetric input matrix.  only the
      !      lower triangle of the matrix need be supplied.
      !
      ! on output
      !
      !    d contains the diagonal elements of the tridiagonal matrix.
      !
      !    e contains the subdiagonal elements of the tridiagonal
      !      matrix in its last n-1 positions.  e(1) is set to zero.
      !
      !    z contains the orthogonal transformation matrix
      !      produced in the reduction.
      !
      !    a and z may coincide.  if distinct, a is unaltered.
      !
      ! questions and comments should be directed to burton s. garbow,
      ! mathematics and computer science div, argonne national laboratory
      !
      ! this version dated august 1983.

      integer, intent(in) :: nm, n
      real*8,  intent(inout), dimension(nm,n) :: a
      real*8,  intent(out),   dimension(nm,n) :: z
      real*8,  intent(out),   dimension(n)    :: d, e

      integer :: i, j, k, l, ii, jp1
      real*8  :: f, g, h, hh, scale

      do i = 1, n
         do j = i, n
            z(j,i) = a(j,i)
         enddo
         d(i) = a(n,i)
      enddo

      if (n .eq. 1) go to 510

      !.......... for i=n step -1 until 2 do -- ..........
      do ii = 2, n

         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0

         if (l .lt. 2) go to 130
         !.......... scale row (algol tol then not needed) ..........

         do k = 1, l
            scale = scale + dabs(d(k))
         enddo

         if (scale .ne. 0.0d0) go to 140

  130    e(i) = d(l)

         do j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
         enddo

         go to 290

  140    do k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
         enddo

         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g

         !.......... form a*u ..........
         do j = 1, l
            e(j) = 0.0d0
         enddo

         do j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .ge. jp1) then
               do k = jp1, l
                  g = g + z(k,j) * d(k)
                  e(k) = e(k) + z(k,j) * f
               enddo
            endif
            e(j) = g
         enddo

         !.......... form p ..........
         f = 0.0d0

         do j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
         enddo

         hh = f / (h + h)

         !.......... form q ..........
         do j = 1, l
            e(j) = e(j) - hh * d(j)
         enddo

         !.......... form reduced a ..........
         do j = 1, l
            f = d(j)
            g = e(j)
            do k = j, l
               z(k,j) = z(k,j) - f * e(k) - g * d(k)
            enddo
            d(j) = z(l,j)
            z(i,j) = 0.0d0
         enddo

  290    d(i) = h

      enddo

      !.......... accumulation of transformation matrices ..........
      do i = 2, n

         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)

         if (h .ne. 0.0d0) then

            do k = 1, l
               d(k) = z(k,i) / h
            enddo

            do j = 1, l
               g = 0.0d0
               do k = 1, l
                  g = g + z(k,i) * z(k,j)
               enddo
               do k = 1, l
                  z(k,j) = z(k,j) - g * d(k)
               enddo
            enddo

         endif

         do k = 1, l
            z(k,i) = 0.0d0
         enddo

      enddo

  510 do i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
      enddo

      z(n,n) = 1.0d0
      e(1) = 0.0d0

      return

   end subroutine tred2

   !-------------------------------------------------------------------!
   pure subroutine treql(n,z,d,e,ierr)

      integer, intent(in)                    :: n
      real*8,  intent(inout), dimension(n,n) :: z
      real*8,  intent(out),   dimension(n)   :: d, e
      integer, intent(out)                   :: ierr
      
      integer :: n2, ii, i, l, k, j, j1, m, l1, ml, it
      real*8  :: cr1, cr2, f, g, h, p, r, b, p1, c, s

      cr2=1.d-15
      cr1=1.d-15
      n2=n+2
      ierr=0

      do 100 ii=2,n
      i=n2-ii
      l=i-2
      f=z(i,i-1)
      g=0.d0
      if(l) 2,2,123
  123 do 1 k=1,l
    1 g=g+z(i,k)**2
      h=g+f*f
      if(g-cr1)2,3,3
    2 e(i)=f
      h=0.d0
      goto 110
    3 l=l+1
      if(f)4,5,5
    4 g=dsqrt(h)
      goto 6
    5 g=-dsqrt(h)
    6 e(i)=g
      h=h-f*g
      z(i,i-1)=f-g
      f=0.d0
      do 9 j=1,l
      z(j,i)=z(i,j)/h
      g=0.d0
      do 7 k=1,j
    7 g=g+z(j,k)*z(i,k)
      j1=j+1
      if(j1.gt.l) goto 124
      do 8 k=j1,l
    8 g=g+z(k,j)*z(i,k)
  124 e(j)=g/h
      f=f+g*z(j,i)
    9 continue
      p=f/(2.d0*h)
      do 10 j=1,l
      f=z(i,j)
      g=e(j)-p*f
      e(j)=g
      do 10 k=1,j
   10 z(j,k)=z(j,k)-f*e(k)-g*z(i,k)
  110 d(i)=h
  100 continue

      d(1)=0.d0
      e(1)=0.d0
      d(1)=z(1,1)
      z(1,1)=1.d0
      do 200 i=2,n
      l=i-1
      if(d(i))11,15,11
   11 do 14 j=1,l
      g=0.d0
      do 12 k=1,l
   12 g=g+z(i,k)*z(k,j)
      do 13 k=1,l
   13 z(k,j)=z(k,j)-g*z(k,i)
   14 continue
   15 d(i)=z(i,i)
      z(i,i)=1.d0
      do 16 j=1,l
      z(i,j)=0.d0
   16 z(j,i)=0.d0
  200 continue
      do 17 i=2,n
   17 e(i-1)=e(i)
      e(n)=0.d0
      b=0.d0
      f=0.d0
      do 300 l=1,n
      j=0
      h=(dabs(d(l))+dabs(e(l)))*cr2
      if(b.lt.h) b=h
      do 18 m=l,n
      if(dabs(e(m)).le.b) goto 19
   18 continue
   19 if(m.eq.l) goto 111
   20 if(j-30)22,21,21
   21 ierr=-1
      goto 444
   22 j=j+1
      g=d(l)
      p=(d(l+1)-g)/(2.d0*e(l))
      r=dsqrt(p*p+1.d0)
      if(p) 23,24,24
   23 p1=p-r
      goto 25
   24 p1=p+r
   25 continue
      d(l)=e(l)/p1
      h=g-d(l)
      l1=l+1
      if(l1.gt.n)  goto 125
      do 26 i=l1,n
   26 d(i)=d(i)-h
  125 f=f+h
      p=d(m)
      c=1.d0
      s=0.d0
      ml=m-l
      do 30 it=1,ml
      i=m-it
      g=c*e(i)
      h=c*p
      if(dabs(p)-dabs(e(i)))28,27,27
   27 c=e(i)/p
      r=dsqrt(c*c+1.d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.d0/r
      goto 29
   28 c=p/e(i)
      r=dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r
      s=1.d0/r
      c=c/r
   29 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 30 k=1,n
      h=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*h
   30 z(k,i)=c*z(k,i)-s*h
      e(l)=s*p
      d(l)=c*p
      if(dabs(e(l)).gt.b) goto 20
  111 d(l)=d(l)+f
  300 continue

      !n1=n-1
      !do 400 i=1,n1
      !j1=i+1
      !k=i
      !p=d(i)
      !do 41 j=j1,n
      !if(p-d(j))41,41,31
      !  31 p=d(j)
      !k=j
      !  41 continue
      !  32 if(i.eq.k) goto 400
      !d(k)=d(i)
      !d(i)=p
      !do 33 j=1,n
      !p=z(j,i)
      !z(j,i)=z(j,k)
      !  33 z(j,k)=p
      ! 400 continue

  444 return

   end subroutine treql

   !-----------------------------------------------------------------------
   pure subroutine invers(n,a,v,ierr)
 
      integer, intent(in)                    :: n
      real*8,  intent(inout), dimension(n,n) :: a
      real*8,  intent(out),   dimension(n)   :: v
      integer, intent(out)                   :: ierr

      integer :: m, k, i, i1, j
      real*8  :: p, y

      m = n - 1

      do k=1,n

         if (a(1,1).eq.0.d0) then
            ierr = -1
            v = 0.d0
            return
         endif

         p = 1.D0/a(1,1)

         do i=2,n
            v(i-1) = a(1,i)
         enddo

         do i=1,m
            a(i,n) = -v(i)*p
            y = a(i,n)
            i1 = i + 1
            do j=i,m
               a(i,j) = a(i1,j+1) + v(j)*y
            enddo
         enddo

         a(n,n) = -p

      enddo

      do i=1,n
         do j=i,n
            p = -a(i,j)
            a(i,j) = p
            a(j,i) = p
         enddo
      enddo

      ierr = 0
      return

   end subroutine invers

   !-----------------------------------------------------------------------
   pure function slen(n,s) result(sss)

      integer, intent(in) :: n
      real*8,  intent(in), dimension(n) :: s
      real*8 :: sss
      
      real*8 :: ss
      integer :: i

      ss = 0.d0
      do i=1,n
         ss=ss+s(i)*s(i)
      enddo

      sss = dsqrt(ss)

   end function slen

   !-----------------------------------------------------------------------
end module eispack
