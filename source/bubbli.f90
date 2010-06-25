pure subroutine bubbli(n,a,l)
!=======================================================================
!  Sorts real and integer arrays (from somewhere...)
!-------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  bubbli.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  bubbli.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.5  2004/01/14 05:27:24  souda
!  -
!
!  Revision 1.4  2004/01/14 05:23:47  souda
!  -
!
!  Revision 1.3  2004/01/14 05:15:07  souda
!  -
!
!  Revision 1.2  2004/01/14 04:51:32  souda
!  nothing
!
!  Revision 1.1.1.1  2003/12/19 16:48:08  souda
!  Initial PCET-4.0 Release
!
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

