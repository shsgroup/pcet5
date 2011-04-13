module lapack_wrappers

   !-------------------------------------------------------------------
   !  wrappers for diagonalization routines from LAPACK package
   !-------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2011-04-13 23:49:48 $
   !  $Revision: 5.1 $
   !  $Log: not supported by cvs2svn $
   !
   !-------------------------------------------------------------------

   implicit none
   private

   public  :: dspevx_wrapper
   public  :: dspevd_wrapper
   public  :: dsyevr_wrapper
   !public  :: evvrsp_wrapper
      
!-------------------------------------------------------------------!
contains

   !-------------------------------------------------------------------!
   subroutine dspevx_wrapper(n,a,w,z,ierr)

      ! this subroutine calls the LAPACK subroutine DSPEVX
      ! to find the eigenvalues and eigenvectors (if desired)
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
      real(8), dimension(:),   allocatable :: ap
      real(8), dimension(:),   allocatable :: wp
      real(8), dimension(:,:), allocatable :: zp

      !-- dspevx arguments
      character(len=1) :: uplo, jobz, range
      integer :: info, il, iu, num_ev_found
      integer, allocatable, dimension(:) :: iwork, ifail    ! work arrays
      real(8), allocatable, dimension(:) :: work            ! work array
      real(8) :: abstol, vl, vu
      
      !-- local variables
      integer :: i, j, ij

      w   = 0.d0
      z   = 0.d0
      ierr = 0

      !-- allocate arrays for DSPEVX

      allocate (ap(n*(n+1)/2))
      allocate (wp(n), zp(n,n))
      allocate (iwork(5*n))
      allocate (work(8*n))
      allocate (ifail(n))
      iwork = 0
      ifail = 0
      work = 0.d0

      !-- set parameters of DSPEVX

      jobz  = "V"    ! eigenvalues and eigenvectors
      range = "A"    ! all eigenvalues and eigenvectors
      uplo  = "U"    ! upper triangle packing

      !-- interval for searching eigenvalues
      !   (not referenced if range="A")
      vl = 0.d0
      vu = 0.d0
      
      !-- The absolute error tolerance for the eigenvalues
      abstol = 0.d0
      
      !-- range of eigenvalue indices to find
      !   (not referenced if range="A")
      il = 1
      iu = n

      !-- pack the input matrix

      if (uplo.eq."U") then

         !-- upper triangle: ap(i + (j-1)*j/2) = a(i,j) for 1<=i<=j

         do j=1,n
            do i=1,j
               ij = i + (j-1)*j/2
               ap(ij) = a(i,j)
            enddo
         enddo

      elseif (uplo.eq."L") then

         !-- lower triangle: ap(i + (j-1)*(2*n-j)/2) = a(i,j) for j<=i<=n

         do j=1,n
            do i=j,n
               ij = i + (j-1)*(2*n-j)/2
               ap(ij) = a(i,j)
            enddo
         enddo

      endif

      !-- call to DSPEVX
      call dspevx(jobz,range,uplo,n,ap,vl,vu,il,iu,abstol,num_ev_found,wp,zp,n,work,iwork,ifail,info)

      !-- delivering the output
      w = wp
      z = zp
      ierr = info

      !-- release temporary storage

      deallocate(ap,wp,zp,work)
      deallocate(iwork,ifail)

   end subroutine dspevx_wrapper


   !-------------------------------------------------------------------!
   subroutine dspevd_wrapper(n,a,w,z,ierr)

      ! this subroutine calls the LAPACK subroutine DSPEVD
      ! (divide and conquer algorithm)
      ! to find the eigenvalues and eigenvectors (if desired)
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
      real(8), dimension(:),   allocatable :: ap
      real(8), dimension(:),   allocatable :: wp
      real(8), dimension(:,:), allocatable :: zp

      !-- dspevx arguments
      character(len=1) :: uplo, jobz
      integer :: lwork, liwork, info
      integer, allocatable, dimension(:) :: iwork    ! work array
      real(8), allocatable, dimension(:) :: work     ! work array
      
      !-- local variables
      integer :: i, j, ij

      w   = 0.d0
      z   = 0.d0
      ierr = 0

      !-- allocate arrays for DSPEVX

      lwork = 1 + 6*n + n*n
      liwork = 3 + 5*n
      allocate (ap(n*(n+1)/2))
      allocate (wp(n), zp(n,n))
      allocate (iwork(liwork))
      allocate (work(lwork))
      iwork = 0
      work = 0.d0

      !-- set parameters of DSPEVD

      jobz  = "V"    ! eigenvalues and eigenvectors
      uplo  = "U"    ! upper triangle packing

      !-- pack the input matrix

      if (uplo.eq."U") then

         !-- upper triangle: ap(i + (j-1)*j/2) = a(i,j) for 1<=i<=j

         do j=1,n
            do i=1,j
               ij = i + (j-1)*j/2
               ap(ij) = a(i,j)
            enddo
         enddo

      elseif (uplo.eq."L") then

         !-- lower triangle: ap(i + (j-1)*(2*n-j)/2) = a(i,j) for j<=i<=n

         do j=1,n
            do i=j,n
               ij = i + (j-1)*(2*n-j)/2
               ap(ij) = a(i,j)
            enddo
         enddo

      endif

      !-- call to DSPEVX
      call dspevd(jobz,uplo,n,ap,wp,zp,n,work,lwork,iwork,liwork,info)

      !-- delivering the output
      w = wp
      z = zp
      ierr = info

      !-- release temporary storage

      deallocate(ap,wp,zp,work)
      deallocate(iwork)

   end subroutine dspevd_wrapper


   !-------------------------------------------------------------------!
   subroutine dsyevr_wrapper(n,a,w,z,ierr)

      ! this subroutine calls the LAPACK subroutine DSYEVR
      ! (relatively robust representations algorithm, RRR)
      ! to find the eigenvalues and eigenvectors (if desired)
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

      !-- temporary arrays for the input matrix,
      !   eigenvalues, and eigenvectors
      real(8), dimension(:,:), allocatable :: ap
      real(8), dimension(:),   allocatable :: wp
      real(8), dimension(:,:), allocatable :: zp

      !-- dspevx arguments
      character(len=1) :: uplo, jobz, range
      integer :: info, il, iu, num_ev_found, lwork, liwork
      integer, allocatable, dimension(:) :: isuppz   ! work array
      integer, allocatable, dimension(:) :: iwork    ! work array
      real(8), allocatable, dimension(:) :: work     ! work array
      real(8) :: abstol, vl, vu
      
      w   = 0.d0
      z   = 0.d0
      ierr = 0

      !-- allocate arrays for DSPEVX

      lwork = 26*n
      liwork = 10*n
      allocate (ap(n,n))
      allocate (wp(n), zp(n,n))
      allocate (iwork(liwork))
      allocate (work(lwork))
      allocate (isuppz(2*n))

      iwork = 0
      work = 0.d0
      isuppz = 0

      !-- set parameters of DSPEVX

      jobz  = "V"    ! eigenvalues and eigenvectors
      range = "A"    ! all eigenvalues and eigenvectors
      uplo  = "U"    ! upper triangle packing

      !-- interval for searching eigenvalues
      !   (not referenced if range="A")
      vl = 0.d0
      vu = 0.d0
      
      !-- The absolute error tolerance for the eigenvalues
      abstol = 0.d0
      
      !-- range of eigenvalue indices to find
      !   (not referenced if range="A")
      il = 1
      iu = n

      !-- copy the input matrix

      ap = a

      !-- call to DSYEVR
      call dsyevr(jobz,range,uplo,n,ap,n,vl,vu,il,iu,abstol,num_ev_found,wp,zp,n,isuppz,work,lwork,iwork,liwork,info)

      !-- delivering the output
      w = wp
      z = zp
      ierr = info

      !-- release temporary storage

      deallocate(ap,wp,zp,work)
      deallocate(iwork,isuppz)

   end subroutine dsyevr_wrapper

   !-----------------------------------------------------------------------
end module lapack_wrappers
