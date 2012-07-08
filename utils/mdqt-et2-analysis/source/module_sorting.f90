!!c*******************************************************************
!!c
!!c   MODULE:      sorting
!!c   PKG VERSION: DLPROTEIN-2.1w
!!c   ACTION:      contains subroutines and functions implementing
!!c                various sorting algorithms
!!c
!!c------------------------------------------------------------------
!!c
!!c   $Author$
!!c   $Date$
!!c   $Revision$
!!c   $Log: not supported by cvs2svn $
!!c   Revision 3.1  2006/12/14 23:28:16  souda
!!c   cleaning (no major changes)
!!c
!!c   Revision 3.0  2006/07/17 23:07:55  souda
!!c   Revision upgraded to 3.0
!!c
!!c   Revision 1.1.1.1  2006/07/17 22:49:45  souda
!!c   Final release of DLPROTEIN_EVB Version 3.0
!!c
!!c   Revision 3.0  2006/07/15 06:04:54  souda
!!c   Final release of version 3.0 (EVB+PI)
!!c
!!c   Revision 2.0  2006/07/05 20:36:33  souda
!!c   revision upgraded to 2.0
!!c
!!c   Revision 1.1  2006/04/19 18:46:05  souda
!!c   added quicksort recursive routine;
!!c   sorting routines combined into a module sorting_mod
!!c
!!c
!!c*******************************************************************

module sorting

   implicit none

contains

   !!------------------------------------------------------
   subroutine quicksort(numbers, array_size, order)

      integer, intent(in) :: array_size
      integer, intent(inout), dimension(:) :: numbers
      integer, intent(inout), dimension(:) :: order

      integer :: left, right

      left = 1
      right = array_size
      call q_sort(numbers, order, left, right)

   end subroutine quicksort

   !!------------------------------------------------------
   recursive subroutine q_sort(numbers, order, left, right)

      integer :: left, right
      integer, dimension(:) :: numbers
      integer, dimension(:) :: order

      integer :: pivot, pivot_order, l_hold, r_hold

      l_hold = left
      r_hold = right
      pivot = numbers(left)
      pivot_order = order(left)

      do while (left.lt.right)

         do while (numbers(right).ge.pivot.and.left.lt.right)
            right = right - 1
         enddo

         if (left.ne.right) then
            numbers(left) = numbers(right)
            order(left) = order(right)
            left = left + 1
         endif

         do while (numbers(left).le.pivot.and.left.lt.right)
            left = left + 1
         enddo

         if (left.ne.right) then
            numbers(right) = numbers(left)
            order(right) = order(left)
            right = right - 1
         endif

      enddo

      numbers(left) = pivot
      order(left) = pivot_order
      pivot = left
      pivot_order = order(left)
      left = l_hold
      right = r_hold

      if(left.lt.pivot) call q_sort(numbers, order, left, pivot-1)
      if(right.gt.pivot) call q_sort(numbers, order, pivot+1, right)

   end subroutine q_sort

   !!------------------------------------------------------
   !!c
   !!c     dlpoly shell sort routine. 
   !!c     Sorts an array of integers into ascending order
   !!c     
   !!c     copyright daresbury laboratory 1993
   !!c     author - t.forester   november 1993
   !!c
   !!c     adopted to Fortran 90 standard by A. Soudackov
   !!c     Penn State, 2006
   !!c
   !!------------------------------------------------------
   subroutine shellsort(n,list)

      integer, intent(in) :: n
      integer, dimension(:), intent(inout) :: list

      integer :: nl, nn, i, imax, ix, j

      !-- set up sort

      if (n.gt.1) then

         !-- number of lists
         nl = n/2

         !-- iterate shell sort

      10 do nn=1,nl

            !-- begin insertion sort on nnth list
          
            do i=nn+nl,n,nl

               imax = list(i)
               ix = i

               !-- find location for insertion
            
               j = i
           100 j = j - nl

               if(j.lt.1) goto 110

               if (list(j).gt.imax) then
                  ix = j
               else
                  j = 1
               endif
            
               goto 100
           110 continue
            
               !-- insert in index array

               do j=i,ix+nl,-nl
                  list(j) = list(j-nl)
               enddo
            
               list(ix) = imax

            enddo

         enddo
        
         nl = nl/2
         if(nl.gt.0) goto 10
        
      endif

   end subroutine shellsort

end module sorting
