module geogas

!======================================================================
!  Geometry of the reaction complex for gas phase calculations
!----------------------------------------------------------------------
!  NATGAS - number of atoms
!  IPTGAS - the numbers of atoms in PT interface (Dp-H-Ap)
!  IETGAS - the numbers of atoms in ET interface (De-Ae)
!  LABGAS - atomic numbers of atoms in the complex
!  XYZGAS - Cartesian coordinates of atoms in the complex
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!======================================================================

   use pardim

   implicit none
   public
   save

   integer :: natgas
   integer :: iptgas(3),ietgas(2)
   integer, allocatable :: labgas(:)                       ! (numatm)
   real*8,  allocatable :: xyzgas(:,:), chrgas(:,:)

   !======================================================================
   contains

   subroutine alloc_geogas(ndim)

      implicit none
      integer, intent(in) :: ndim

      real*8 :: zero = 0.d0
      integer :: izero = 0

      if (ndim.le.0) then
         write(*,*) "ERROR in alloc_geogas: Wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(labgas) ) then
         allocate (labgas(ndim))
         labgas = izero
      else
         write(*,*) "ERROR in alloc_geogas: labgas is already allocated..."
         stop
      endif

      if ( .not.allocated(xyzgas) ) then
         allocate (xyzgas(3,ndim))
         xyzgas = zero
      else
         write(*,*) "ERROR in alloc_geogas: xyzgas is already allocated..."
         stop
      endif

      if ( .not.allocated(chrgas) ) then
         allocate (chrgas(nelst,ndim))
         chrgas = zero
      else
         write(*,*) "ERROR in alloc_geogas: chrgas is already allocated..."
         stop
      endif

      return

   end subroutine alloc_geogas


   subroutine dealloc_geogas
      implicit none
      if (allocated(labgas)) deallocate (labgas)
      if (allocated(xyzgas)) deallocate (xyzgas)
      if (allocated(chrgas)) deallocate (chrgas)
      return
   end subroutine dealloc_geogas

end module geogas
