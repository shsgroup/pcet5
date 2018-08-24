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

   use parsol, only: f0
   use minima_1d

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
   !  Calculates the derivative of self-energy of the inertial
   !  polarization with respect to energy gap coordinate
   !  (summation over linearly independent solvent variables)
   !   ZE - the ET medium coordinate (energy gap) (kcal/mole);
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function selfen_et2_der(ze) result(sw_der)
      real(kind=8), intent(in) :: ze
      real(kind=8)             :: sw_der
      real(kind=8)             :: z_shifted
      z_shifted = ze + t2r(1,2)
      sw_der = z_shifted/lambda/2.d0
   end function selfen_et2_der

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
   subroutine fes_et2(mode,energy_gap,free_energy,gradient,eigenvectors)

      character(len=5), intent(in) :: mode
      real(kind=8), intent(in)  :: energy_gap
      real(kind=8), dimension(2),   intent(out) :: free_energy
      real(kind=8), optional, dimension(2),   intent(out) :: gradient
      real(kind=8), optional, dimension(2,2), intent(out) :: eigenvectors

      real(kind=8) :: self_energy, h11, h22, h12, dh12, sh12, sqd, vnorm
      real(kind=8) :: v1, v2, v1norm, v2norm

      free_energy = 0.d0
      if (present(gradient))     gradient = 0.d0
      if (present(eigenvectors)) eigenvectors = 0.d0

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

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extremum points and its characteristics on the ground state
   ! adiabatic free energy surface for a two-state ET model.
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine adiabatic_fes_et2_parameters(double_well,fr,fb,eb)

      logical, intent(out) :: double_well
      real(kind=8), intent(out) :: fr, fb, eb

      integer :: n_extrema, k, i
      real(kind=8), parameter :: tol=1.d-9
      real(kind=8) :: x, xr, xc, xp, x_right, x_left, dx, x_prev
      real(kind=8) :: xminr, xmaxb, fes_der, fes_der_prev, dg, vet
      real(kind=8), dimension(10) :: extremum

      fr = 0.d0
      fb = 0.d0
      eb = 0.d0
      extremum = 0.d0

      !-- Total reaction free energy
      dg = (hg2(2,2) + 0.5d0*t2inf(2,2) + gsolv_2) - (hg2(1,1) + 0.5d0*t2inf(1,1) + gsolv_1)

      !-- electronic coupling
      vet = hg2(1,2)

      !-- minima and the crossing point for reactant and product diabatic FES
      xr = 0.d0
      xp = -sqrt(2.d0*lambda/f0)
      xc = -(lambda + dg)/sqrt(2.d0*lambda*f0)

      !-- check if we are in the normal region (double-well character)
      if (xc.gt.xr.or.xc.lt.xp) then
         double_well = .false.
         return
      endif

      !-- Start backward scanning of the adiabatic ground state FES
      !   to identify approximate locations of the extrema

      x_right = xr + abs(xp)/2.d0
      x_left = xp - abs(xp)/2.d0
      dx = (x_left - x_right)/100.d0
      x_prev = x_right + abs(dx)

      n_extrema = 0

      do i=1,101

         x = x_right + (i-1)*dx

         fes_der_prev = (2*f0*x_prev + sqrt(2.d0)*sqrt(f0*lambda))/2.d0 &
                 & - (sqrt(f0*lambda)*(dg + lambda + sqrt(2.d0)*x_prev*sqrt(f0*lambda)))/ &
                 &   (sqrt(2.d0)*sqrt(4.d0*vet*vet + (dg + lambda + sqrt(2.d0)*x_prev*sqrt(f0*lambda))**2))

         fes_der = (2*f0*x + sqrt(2.d0)*sqrt(f0*lambda))/2.d0 &
                 & - (sqrt(f0*lambda)*(dg + lambda + sqrt(2.d0)*x*sqrt(f0*lambda)))/ &
                 &   (sqrt(2.d0)*sqrt(4.d0*vet*vet + (dg + lambda + sqrt(2.d0)*x*sqrt(f0*lambda))**2))

         if (fes_der*fes_der_prev.lt.0.d0) then
            n_extrema = n_extrema + 1
            if (fes_der.lt.0) then
               extremum(n_extrema) = funmin(fesx,x,x_prev,tol)
            else
               extremum(n_extrema) = funmin(fesx_neg,x,x_prev,tol)
            endif
         endif

         if (n_extrema.ge.3) exit
         x_prev = x

      enddo

      !-- analyze extremum points

      if (n_extrema.lt.3) then
         double_well = .false.
         return
      else
         double_well = .true.
      endif

      !-- calculate second derivatives at the extremum points

      xminr = extremum(1)
      fr = f0*(1.d0 - (4.d0*lambda*vet*vet)/(4.d0*vet*vet + 2.d0*f0*lambda*xminr*xminr & 
                  & + (dg + lambda)*(dg + lambda + 2.d0*sqrt(2.d0*f0*lambda)*xminr))**1.5d0)

      xmaxb = extremum(2)
      fb = -f0*(1.d0 - (4.d0*lambda*vet*vet)/(4.d0*vet*vet + 2.d0*f0*lambda*xmaxb*xmaxb & 
                  & + (dg + lambda)*(dg + lambda + 2.d0*sqrt(2.d0*f0*lambda)*xmaxb))**1.5d0)

      !-- barrier height
      eb = fesx(xmaxb) - fesx(xminr)

      write(*,*)
      write(*,'(1x,"===================================================================")')
      write(*,'(1x,"         Adiabatic ground state FES characteristics                ")')
      write(*,'(1x,"-------------------------------------------------------------------")')
      write(*,'(1x,"Reactant minimum at x(scaled) = ",f12.6," sqrt(kcal/mol)")') xminr
      write(*,'(1x,"Transition state at x(scaled) = ",f12.6," sqrt(kcal/mol)")') xmaxb
      write(*,'(1x,"Force constant at the reactant minimum: ",f12.6)') fr
      write(*,'(1x,"Force constant at the transition state: ",f12.6)') fb
      write(*,'(1x,"Activation free energy (kcal/mol):      ",f12.6)') eb
      write(*,'(1x,"-------------------------------------------------------------------")')

   end subroutine adiabatic_fes_et2_parameters


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Ground state adiabatic free energy (wrapper)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function fesx(x) result(ground_state_energy)

      implicit none
      real(kind=8) :: x
      real(kind=8) :: ground_state_energy
      real(kind=8), dimension(2) :: fes
      real(kind=8) :: egap

      !-- transform x to energy gap coordinate
      call z1_to_ze(x,egap)

      !-- calculate adiabatic free energy
      call fes_et2("ADIAB",egap,free_energy=fes)

      ground_state_energy = fes(1)

   end function fesx

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Negative ground state adiabatic free energy (wrapper)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function fesx_neg(x) result(ground_state_energy)

      implicit none
      real(kind=8) :: x
      real(kind=8) :: ground_state_energy
      real(kind=8), dimension(2) :: fes
      real(kind=8) :: egap

      !-- transform x to energy gap coordinate
      call z1_to_ze(x,egap)

      !-- calculate adiabatic free energy
      call fes_et2("ADIAB",egap,free_energy=fes)

      ground_state_energy = -fes(1)

   end function fesx_neg

end module data_et2

