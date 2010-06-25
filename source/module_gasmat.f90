module gasmat
!=======================================================================
!     Gas phase Hamiltonian matrices on the proton grid (kcal/mol)
!-----------------------------------------------------------------------
!     H0    - gas phase Hamiltonian matrix in EVB basis
!     DH0   - first derivative of the gas phase Hamiltonian matrix
!             with respect to the proton coordinate
!     D2H0  - second derivative of the gas phase Hamiltonian matrix
!             with respect to the proton coordinate
!     DGH0  - first derivative of the gas phase Hamiltonian matrix
!             with respect to the gating coordinate
!     DG2H0 - second derivative of the gas phase Hamiltonian matrix
!             with respect to the gating coordinate
!-----------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_gasmat.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_gasmat.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 19:59:49  souda
!  Initial PCET-4.0 Release
!
!
!=======================================================================
   implicit none
   public
   save

   real*8, allocatable, dimension(:,:,:,:) :: h0, dh0, d2h0, dgh0, dg2h0

   !=======================================================================
   ! member subroutines and functions
   !=======================================================================
   contains

   subroutine alloc_gasmat(ndim,ngdim,deriv,derivg)
   
      implicit none
      integer, intent(in) :: ndim, ngdim
      logical, intent(in) :: deriv, derivg

      real*8 :: zero = 0.d0

      if (ndim.le.0.or.ngdim.le.0) then
         write(*,*) "ERROR in alloc_gasmat: Wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(h0) ) then
         allocate (h0(4,4,ndim,ngdim))
         h0 = zero
      else
         write(*,*) "ERROR in alloc_gasmat: h0 is already allocated..."
         stop
      endif

      if (deriv) then

         if (.not.allocated  (dh0).and. &
             .not.allocated (d2h0)        ) then

            allocate ( dh0(4,4,ndim,ngdim))
            allocate (d2h0(4,4,ndim,ngdim))
            dh0  = zero
            d2h0 = zero
            
         else

            write(*,*) "ERROR in alloc_gasmat: either dh0 or d2h0 is already allocated..."
            stop

         endif

      endif

      if (derivg) then

         if (.not.allocated  (dgh0).and. &
             .not.allocated (dg2h0)        ) then

            allocate ( dgh0(4,4,ndim,ngdim))
            allocate (dg2h0(4,4,ndim,ngdim))
            dgh0 = zero
            dg2h0 = zero
            
         else

            write(*,*) "ERROR in alloc_gasmat: either dgh0 or dg2h0 is already allocated..."
            stop

         endif

      endif
      
      return

   end subroutine alloc_gasmat


   subroutine dealloc_gasmat
      implicit none
      if (allocated(   h0)) deallocate (   h0)
      if (allocated(  dh0)) deallocate (  dh0)
      if (allocated( d2h0)) deallocate ( d2h0)
      if (allocated( dgh0)) deallocate ( dgh0)
      if (allocated(dg2h0)) deallocate (dg2h0)
      return
   end subroutine dealloc_gasmat

end module gasmat
