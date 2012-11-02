module data_et2

!====================================================================
!  Hamiltonian, reorganization energy matrices and related quantities
!  for a two-state ET model (kcal/mol)
!--------------------------------------------------------------------
!
!  HG2    - [2x2] gas-phase Hamiltonian
!  T2     - [2x2] inertial reorganization energy matrix
!  T2INF  - [2x2] electronic reorganization energy matrix
!  T2R    - [2x2] reduced inertial reorganization energy matrix
!  T2RINF - [2x2] reduced electronic reorganization energy matrix
!
!  LAMBDA       - outer-sphere reorganization energy
!  SCALE_FACTOR - scaling factor for the energy gap coordinate
!  DELTA_SHIFT  - shift for the scaled solvent coordinate
!
!  DIPOLE_MOMENT_DIAB - [2x2] matrix of the dipole moment
!                       in the basis of electronic diabatic states
!
!--------------------------------------------------------------------
!
!  $Author$
!  $Id$
!  $Revision$
!
!=======================================================================

   implicit none
   public
   save

   real(8) :: lambda, gsolv_1, gsolv_2, z_1, z_2
   real(8) :: scale_factor
   real(8) :: delta_shift

   real(8), dimension(2,2) :: t2, t2inf, t2r, t2rinf, hg2

!=======================================================================
contains

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !  Calculates the self-energy of the inertial polarization
   !  (summation over linearly independent solvent variables)
   !   ZE - the ET medium coordinate (energy gap) (kcal/mole);
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function selfen_et2(ze) result(sw)
      real(kind=8), intent(in) :: ze
      real(kind=8)             :: sw
      real(kind=8)             :: z_shifted
      z_shifted = ze + t2r(1,2)
      sw = z_shifted*z_shifted/lambda/4.d0
   end function selfen_et2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (Z) to (z) - coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine ze_to_z1(ze,z1)
      real(8), intent(in)  :: ze
      real(8), intent(out) :: z1
      z1 = ze/scale_factor + delta_shift
   end subroutine ze_to_z1

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (VE) to (v1) - velocity
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine ve_to_v1(ve,v1)
      real(8), intent(in)  :: ve
      real(8), intent(out) :: v1
      v1 = ve/scale_factor
   end subroutine ve_to_v1

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (z1) to (ZE) - coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine z1_to_ze(z1,ze)
      real(kind=8), intent(in)  :: z1
      real(kind=8), intent(out) :: ze
      ze = scale_factor*(z1-delta_shift)
   end subroutine z1_to_ze

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation from (v1) to (VE) - velocity
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine v1_to_ve(v1,ve)
      real(kind=8), intent(in)  :: v1
      real(kind=8), intent(out) :: ve
      ve = scale_factor*v1
   end subroutine v1_to_ve

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation of the gradient from (ZE) to (z1)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine ge_to_g1(ge,g1)
      real(kind=8), intent(in)  :: ge
      real(kind=8), intent(out) :: g1
      g1 = scale_factor*ge
   end subroutine ge_to_g1

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Transformation of the gradient from (z1) to (ZE)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine g1_to_ge(g1,ge)
      real(kind=8), intent(in)  :: g1
      real(kind=8), intent(out) :: ge
      ge = g1/scale_factor
   end subroutine g1_to_ge

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Diabatic and adiabatic free energies, interaction gradients,
   ! and nonadiabatic couplings for a two-state ET model.
   ! Note: gradients do not include contribution from self energy
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine fes_et2(mode,energy_gap,free_energy,gradient,eigenvectors,nacoupling)

      character(len=5), intent(in) :: mode
      real(kind=8), intent(in)  :: energy_gap
      real(kind=8), dimension(2),   intent(out) :: free_energy
      real(kind=8), optional, dimension(2),   intent(out) :: gradient
      real(kind=8), optional, dimension(2,2), intent(out) :: eigenvectors
      real(kind=8), optional, dimension(2,2), intent(out) :: nacoupling

      real(kind=8) :: self_energy, h11, h22, h12, dh12, sh12, sqd, vnorm
      real(kind=8) :: v1, v2, v1norm, v2norm

      free_energy = 0.d0
      if (present(gradient))     gradient = 0.d0
      if (present(eigenvectors)) eigenvectors = 0.d0
      if (present(nacoupling))   nacoupling = 0.d0

      self_energy = selfen_et2(energy_gap)

      !-- surrogate Hamiltonian
      h11 = hg2(1,1) + gsolv_1
      h22 = hg2(2,2) + gsolv_1 + energy_gap
      h12 = hg2(1,2)

      if (mode.eq."ADIAB".and.h12.ne.0.d0) then

         sh12 = h11 + h22
         dh12 = h22 - h11
         sqd = sqrt(dh12*dh12 + 4.d0*h12*h12)

         !-- free energies
         free_energy(1) = 0.5d0*(sh12 - sqd) + self_energy
         free_energy(2) = 0.5d0*(sh12 + sqd) + self_energy

         !-- gradients (without self-energy contributions)
         if (present(gradient)) then
            gradient(1) = 0.5d0*(1.d0 - dh12/sqd)
            gradient(2) = 0.5d0*(1.d0 + dh12/sqd)
         endif

         !-- eigenvectors
         if (present(eigenvectors)) then
            v1 = -0.5d0*(dh12 + sqd)/h12
            v2 = -0.5d0*(dh12 - sqd)/h12
            v1norm = sqrt(v1*v1 + 1.d0)
            v2norm = sqrt(v2*v2 + 1.d0)
            eigenvectors(1,1) = v1/v1norm
            eigenvectors(2,1) = 1.d0/v1norm
            eigenvectors(1,2) = v2/v2norm
            eigenvectors(2,2) = 1.d0/v2norm
         endif

         !-- nonadiabatic couplings
         if (present(nacoupling).and.present(eigenvectors)) then
            nacoupling(1,2) = eigenvectors(2,1)*eigenvectors(2,2)/(free_energy(2)-free_energy(1))
            nacoupling(2,1) = -nacoupling(1,2)
         endif

      elseif (mode.eq."DIAB4".or.h12.eq.0.d0) then

         !-- free energies
         free_energy(1) = h11 + self_energy
         free_energy(2) = h22 + self_energy

         !-- gradients (without self-energy contributions)
         if (present(gradient)) then
            gradient(1) = 0.d0
            gradient(2) = 1.d0
         endif

         !-- eigenvectors
         if (present(eigenvectors)) then
            eigenvectors(1,1) = 1.d0
            eigenvectors(2,1) = 0.d0
            eigenvectors(1,2) = 0.d0
            eigenvectors(2,2) = 1.d0
         endif

      endif

   end subroutine fes_et2

end module data_et2

