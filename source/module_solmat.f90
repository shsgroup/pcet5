module solmat
!====================================================================
!  Reorganization energy matrices and related quantities
!  on the proton grid (kcal/mol)
!--------------------------------------------------------------------
!  T      - inertial reorganization energy matrix
!  TINF   - electronic reorganization energy matrix
!  TR     - reduced inertial reorganization energy matrix
!  TRINF  - reduced electronic reorganization energy matrix
!  T1     - inversed inertial reorganization energy matrix
!  DT     - derivative of the inertial reorganization energy matrix
!           with respect to the proton coordinate
!  DTINF  - derivative of the electronic reorganization energy matrix
!           with respect to the proton coordinate
!  D2TINF - second derivative of the electronic reorganization
!           energy matrix with respect to the proton coordinate
!  ERPT   - reorganization energy for PT process
!  ERET   - reorganization energy for ET process
!  ERX    - ET/PT cross-term matrix element
!--------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_solmat.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_solmat.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:07:09  souda
!  Initial PCET-4.0 Release
!
!
!=======================================================================
   implicit none
   public
   save

   real*8                                  :: erpt, eret, erx
   real*8, allocatable, dimension(:,:,:,:) :: t, tinf, tr, trinf, t1
   real*8, allocatable, dimension(:,:,:)   :: dt, dtr, dtinf, d2tinf

   !=======================================================================
   contains

   subroutine alloc_solmat(ndim,ngdim,deriv)

      implicit none
      integer, intent(in) :: ndim, ngdim
      logical, intent(in) :: deriv

      real*8 :: zero = 0.d0

      if (ndim.le.0.or.ngdim.le.0) then
         write(*,*) "error in alloc_solmat: wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(t) ) then
         allocate (t(4,4,ndim+1,ngdim+1))
         t = zero
      else
         write(*,*) "error in alloc_solmat: t is already allocated..."
         stop
      endif

      if ( .not.allocated(tinf) ) then
         allocate (tinf(4,4,ndim+1,ngdim+1))
         tinf = zero
      else
         write(*,*) "error in alloc_solmat: tinf is already allocated..."
         stop
      endif

      if ( .not.allocated(tr) ) then
         allocate (tr(4,4,ndim+1,ngdim+1))
         tr = zero
      else
         write(*,*) "error in alloc_solmat: tr is already allocated..."
         stop
      endif

      if ( .not.allocated(trinf) ) then
         allocate (trinf(4,4,ndim+1,ngdim+1))
         trinf = zero
      else
         write(*,*) "error in alloc_solmat: trinf is already allocated..."
         stop
      endif

      if ( .not.allocated(t1) ) then
         allocate (t1(2,2,ndim+1,ngdim+1))
         t1 = zero
      else
         write(*,*) "error in alloc_solmat: t1 is already allocated..."
         stop
      endif


      if (deriv) then

         if (.not.allocated (dt)    .and. &
             .not.allocated (dtr)   .and. &
             .not.allocated (dtinf) .and. &
             .not.allocated (d2tinf)        ) then

            allocate (dt    (4,ndim+1,ngdim+1))
            allocate (dtr   (4,ndim+1,ngdim+1))
            allocate (dtinf (4,ndim+1,ngdim+1))
            allocate (d2tinf(4,ndim+1,ngdim+1))
            dt     = zero
            dtr    = zero
            dtinf  = zero
            d2tinf = zero

         else

            write(*,*) "error in alloc_solmat: some of derivative matrices are already allocated..."
            stop

         endif

      endif

   end subroutine alloc_solmat

   subroutine dealloc_solmat
      implicit none
      if (allocated(t    ))  deallocate (t    )
      if (allocated(tinf ))  deallocate (tinf )
      if (allocated(tr   ))  deallocate (tr   )
      if (allocated(trinf))  deallocate (trinf)
      if (allocated(t1   ))  deallocate (t1   )
      if (allocated(dt    )) deallocate (dt    )
      if (allocated(dtr   )) deallocate (dtr   )
      if (allocated(dtinf )) deallocate (dtinf )
      if (allocated(d2tinf)) deallocate (d2tinf)
   end subroutine dealloc_solmat

end module solmat

