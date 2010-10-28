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
!  ERX    - ET/PT cross-reorganization energy
!  ER1    - lowest eigenvalue of the truncated reorganization energy matrix
!  ER2    - highest eigenvalue of the truncated reorganization energy matrix
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!=======================================================================

   implicit none
   public
   save

   real(8) :: erpt, eret, erx
   real(8) :: er1, er2  ! er1<er2
   real(8) :: sq1, sq2  ! scaling factors
   real(8) :: cos_theta, sin_theta, delta_1, delta_2
   real(8), allocatable, dimension(:,:,:,:) :: t, tinf, tr, trinf, t1
   real(8), allocatable, dimension(:,:,:)   :: dt, dtr, dtinf, d2tinf
   real(8), allocatable, dimension(:,:)     :: ztmat, gtmat

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

      if ( .not.allocated(ztmat) ) then
         allocate (ztmat(2,2))
         ztmat = zero
      else
         write(*,*) "error in alloc_solmat: ztmat is already allocated..."
         stop
      endif

      if ( .not.allocated(gtmat) ) then
         allocate (gtmat(2,2))
         gtmat = zero
      else
         write(*,*) "error in alloc_solmat: ztmat is already allocated..."
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
      if (allocated(t     )) deallocate (t    )
      if (allocated(tinf  )) deallocate (tinf )
      if (allocated(tr    )) deallocate (tr   )
      if (allocated(trinf )) deallocate (trinf)
      if (allocated(t1    )) deallocate (t1   )
      if (allocated(dt    )) deallocate (dt    )
      if (allocated(dtr   )) deallocate (dtr   )
      if (allocated(dtinf )) deallocate (dtinf )
      if (allocated(d2tinf)) deallocate (d2tinf)
      if (allocated(ztmat )) deallocate (ztmat)
      if (allocated(gtmat )) deallocate (gtmat)
   end subroutine dealloc_solmat


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (ZP,ZE) to (z1,z2)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine zpze_to_z1z2(zp,ze,z1,z2)
      real(8), intent(in)  :: zp, ze
      real(8), intent(out) :: z1, z2
      z1 = ( cos_theta*zp + sin_theta*ze)/sq1 + delta_1
      z2 = (-sin_theta*zp + cos_theta*ze)/sq2 + delta_2
   end subroutine zpze_to_z1z2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (z1,z2) to (ZP,ZE)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine z1z2_to_zpze(z1,z2,zp,ze)
      real(8), intent(in)  :: z1, z2
      real(8), intent(out) :: zp, ze
      zp = sq1*cos_theta*(z1-delta_1) - sq2*sin_theta*(z2-delta_2)
      ze = sq1*sin_theta*(z1-delta_1) + sq2*cos_theta*(z2-delta_2)
   end subroutine z1z2_to_zpze

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation of the gradient from (ZP,ZE) to (z1,z2)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine gpge_to_g1g2(gp,ge,g1,g2)
      real(8), intent(in)  :: gp, ge
      real(8), intent(out) :: g1, g2
      g1 = sq1*( cos_theta*gp + sin_theta*ge)
      g2 = sq2*(-sin_theta*gp + cos_theta*ge)
   end subroutine gpge_to_g1g2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation of the gradient from (z1,z2) to (ZP,ZE)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine g1g2_to_gpge(g1,g2,gp,ge)
      real(8), intent(in)  :: g1, g2
      real(8), intent(out) :: gp, ge
      gp = cos_theta*g1/sq1 - sin_theta*g2/sq2
      ge = sin_theta*g1/sq1 + cos_theta*g2/sq2
   end subroutine g1g2_to_gpge

end module solmat

