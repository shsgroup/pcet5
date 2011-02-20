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
!  DIPOLE_MOMENT_DIAB - matrix of the dipole moment in the basis of
!                       electronic diabatic states
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-02-20 00:58:11 $
!  $Revision: 5.5 $
!  $Log: not supported by cvs2svn $
!  Revision 5.4  2011/02/09 20:51:41  souda
!  added two subroutines for transformation of velocities
!  (affects only the output of velocities in zp-ze frame)
!
!  Revision 5.3  2010/11/04 22:43:09  souda
!  Next iteration... and two additional Makefiles for building the code with debug options.
!
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
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
   real(8), allocatable, dimension(:,:,:)   :: dipole_moment_diab_x
   real(8), allocatable, dimension(:,:,:)   :: dipole_moment_diab_y
   real(8), allocatable, dimension(:,:,:)   :: dipole_moment_diab_z
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

      allocate (dipole_moment_diab_x(4,4,ndim+1))
      allocate (dipole_moment_diab_y(4,4,ndim+1))
      allocate (dipole_moment_diab_z(4,4,ndim+1))
      dipole_moment_diab_x = zero
      dipole_moment_diab_y = zero
      dipole_moment_diab_z = zero

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
      if (allocated(dipole_moment_diab_x)) deallocate (dipole_moment_diab_x)
      if (allocated(dipole_moment_diab_y)) deallocate (dipole_moment_diab_y)
      if (allocated(dipole_moment_diab_z)) deallocate (dipole_moment_diab_z)
      if (allocated(ztmat )) deallocate (ztmat)
      if (allocated(gtmat )) deallocate (gtmat)
   end subroutine dealloc_solmat

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !  Calculates the self-energy of the inertial polarization
   !  (summation over linearly independent solvent variables)
   !   K     - the proton position (GRID POINT);
   !   KG    - the gating position (GRID POINT);
   !   ZP,ZE - the (PT) and (ET) medium coordinates (kcal/mole);
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function selfen(k,kg,zp,ze) result(sw)

      integer, intent(in)  :: k, kg
      real*8,  intent(in)  :: zp, ze
      real*8               :: sw

      integer :: i, j
      real*8, dimension(2) :: xm

      xm(1) = zp
      xm(2) = ze
      sw = 0.d0
      do i=1,2
         do j=1,2
            sw = sw + (xm(i)+tr(1,i+1,k,kg))*t1(i,j,k,kg)*(xm(j)+tr(1,j+1,k,kg))
         enddo
      enddo

      sw = 0.5d0*sw

   end function selfen

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (ZP,ZE) to (z1,z2) - coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine zpze_to_z1z2(zp,ze,z1,z2)
      real(8), intent(in)  :: zp, ze
      real(8), intent(out) :: z1, z2
      z1 = ( cos_theta*zp + sin_theta*ze)/sq1 + delta_1
      z2 = (-sin_theta*zp + cos_theta*ze)/sq2 + delta_2
   end subroutine zpze_to_z1z2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (VP,VE) to (v1,v2) - velocities
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine vpve_to_v1v2(vp,ve,v1,v2)
      real(8), intent(in)  :: vp, ve
      real(8), intent(out) :: v1, v2
      v1 = ( cos_theta*vp + sin_theta*ve)/sq1
      v2 = (-sin_theta*vp + cos_theta*ve)/sq2
   end subroutine vpve_to_v1v2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (z1,z2) to (ZP,ZE) - coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine z1z2_to_zpze(z1,z2,zp,ze)
      real(8), intent(in)  :: z1, z2
      real(8), intent(out) :: zp, ze
      zp = sq1*cos_theta*(z1-delta_1) - sq2*sin_theta*(z2-delta_2)
      ze = sq1*sin_theta*(z1-delta_1) + sq2*cos_theta*(z2-delta_2)
   end subroutine z1z2_to_zpze

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (v1,v2) to (VP,VE) - velocities
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine v1v2_to_vpve(v1,v2,vp,ve)
      real(8), intent(in)  :: v1, v2
      real(8), intent(out) :: vp, ve
      vp = sq1*cos_theta*v1 - sq2*sin_theta*v2
      ve = sq1*sin_theta*v1 + sq2*cos_theta*v2
   end subroutine v1v2_to_vpve

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

