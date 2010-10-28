pure subroutine bubbli(n,a,l)
!=======================================================================
!  Sorts real and integer arrays (from somewhere...)
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:35 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!=======================================================================
   implicit none

   integer, intent(in)                :: n
   real(8), intent(in),  dimension(n) :: a
   integer, intent(out), dimension(n) :: l

   integer :: i, m, n1, k, j
   real(8) :: u
   real(8), dimension(n) :: a0

   do i=1,n
      l(i) = i
   enddo
   a0 = a

   n1 = n - 1

   !if (n1 == 0) return

outer: do i=1,n1,1

          k = i

inner:    do while (k.ge.1)
             j = k + 1
             if (a0(k).ge.a0(j)) cycle outer
             u = a0(k)
             m = l(k)
             a0(k) = a0(j)
             l(k) = l(j)
             a0(j) = u
             l(j) = m
             k = k - 1
          enddo inner

       enddo outer

   return

end subroutine bubbli

