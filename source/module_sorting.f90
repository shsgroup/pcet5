module sorting


   implicit none
   private

   public :: qsort, bubbli
 
contains
 
   !-------------------------------------------------------------
   ! adopted to a real array from:
   ! http://rosettacode.org/wiki/Sorting_algorithms/Quicksort
   !-------------------------------------------------------------
   recursive subroutine qsort(a,order)

      real(8), intent(in out) :: a(:)
      integer, intent(in out) :: order(:)
      integer :: split
 
      if(size(a) > 1) then
         call partition(a, order, split)
         call qsort(a(:split-1),order(:split-1))
         call qsort(a(split:),order(split:))
      endif
 
   end subroutine qsort
 
   subroutine partition(a, order, marker)
 
      real(8), intent(in out) :: a(:)
      integer, intent(in out) :: order(:)
      integer, intent(out) :: marker
      real(8) :: pivot, temp
      integer :: left, right, itemp
 
      pivot = (a(1) + a(size(a))) / 2  ! average of first and last elements to prevent quadratic 
      left = 0                         ! behavior with sorted of reverse sorted data
      right = size(a) + 1
 
      do while (left < right)
         right = right - 1
         do while (a(right) > pivot)
            right = right - 1
         enddo
         left = left + 1
         do while (a(left) < pivot)
            left = left + 1
         enddo
         if (left < right) then 
            temp = a(left)
            itemp = order(left)
            a(left) = a(right)
            order(left) = order(right)
            a(right) = temp
            order(right) = itemp
         endif
      enddo

      if (left == right) then
         marker = left + 1
      else
         marker = left
      endif
 
   end subroutine partition

   !-----------------------------------------------------------------------
   !  sorts real and integer arrays (from somewhere...)
   !-----------------------------------------------------------------------
   subroutine bubbli(n,a,l)

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
         inner: do while (k.ge.1)
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

   end subroutine bubbli

end module sorting
