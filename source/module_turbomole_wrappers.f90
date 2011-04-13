module turbomole_wrappers

   !-------------------------------------------------------------------
   !  wrappers for diagonalization routines from TURBOMOLE (rdiag)
   !-------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2011-04-13 23:49:48 $
   !  $Revision: 5.1 $
   !  $Log: not supported by cvs2svn $
   !
   !-------------------------------------------------------------------

   use sorting

   implicit none
   private

   real(8), parameter :: epsln=1.d-14

   public  :: rdiag_wrapper
      
!-------------------------------------------------------------------!
contains

   !-------------------------------------------------------------------!
   subroutine rdiag_wrapper(n,a,w,z,ierr)

      ! this subroutine calls the subroutine RDIAG from TURBOMOLE
      ! to find the eigenvalues and eigenvectors
      ! of a real symmetric matrix.
      !
      ! on input
      !
      !    n  is the order of the matrix  a.
      !
      !    a  contains the real symmetric matrix.
      !
      ! on output
      !
      !    w  contains the eigenvalues in ascending order.
      !
      !    z  contains the eigenvectors
      !
      !    ierr  is an integer error code:
      !         = 0:  successful exit
      !         < 0:  if ierr = -i, the i-th argument had an illegal value
      !         > 0:  if ierr = i, then i eigenvectors failed to converge.

      integer, intent(in)    :: n
      real*8,  intent(in)    :: a(n,n)
      integer, intent(out)   :: ierr
      real*8,  intent(out)   :: w(n), z(n,n)

      !-- temporary arrays for a lower triangle of the input matrix,
      !   eigenvalues, and eigenvectors
      integer, dimension(:),   allocatable :: order
      real(8), dimension(:),   allocatable :: ap
      real(8), dimension(:),   allocatable :: wp
      real(8), dimension(:,:), allocatable :: zp
      real(8), dimension(:,:), allocatable :: zp_sorted

      !-- local variables
      integer :: i, j, ij, jj

      w   = 0.d0
      z   = 0.d0
      ierr = 0

      !-- allocate arrays for RDIAG

      allocate (ap(n*(n+1)/2))
      allocate (wp(n), zp(n,n))
      ap = 0.d0
      zp = 0.d0
      do i=1,n
         zp(i,i) = 1.d0
      enddo

      !-- pack the input matrix
      !   upper triangle: ap(i + (j-1)*j/2) = a(i,j) for 1<=i<=j

      !do j=1,n
      !   do i=1,j
      !      ij = i + (j-1)*j/2
      !      ap(ij) = a(i,j)
      !   enddo
      !enddo

      jj=0
      do j=1,n
         do i=1,j
            ap(jj+i)=a(i,j)
         enddo
         jj=jj+j
      enddo

      !-- call to RDIAG
      call rdiag(ap,zp,wp,n,epsln)

      allocate (order(n))
      do i=1,n
         order(i) = i
      enddo

      !-- sort the eigenvalues
      call qsort(wp,order)

      !-- sort eigenvectors accordingly
      allocate (zp_sorted(n,n))
      do i=1,n
         zp_sorted(:,i) = zp(:,order(i))
      enddo

      !-- delivering the output
      w = wp
      z = zp_sorted
      ierr = 0

      !-- release temporary storage

      deallocate(ap,wp,zp,zp_sorted,order)

   end subroutine rdiag_wrapper

   !-----------------------------------------------------------------------
end module turbomole_wrappers
