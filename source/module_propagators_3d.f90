module propagators_3d

   !---------------------------------------------------------------------
   ! Contains the routines for Langevin and MDQT propagators
   !---------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2012-04-06 22:38:46 $
   !  $Revision: 5.18 $
   !  $Log: not supported by cvs2svn $
   !  Revision 5.17  2012/03/13 22:07:18  souda
   !  Langevin propagator for Onodera-2 model added (4-dimensional)
   !
   !  Revision 5.16  2011/06/03 05:05:45  souda
   !  Kinetic energy components are calculated and added to the trajectory data
   !
   !  Revision 5.15  2011/04/13 23:49:48  souda
   !  Minor restructuring of modules and addition of LAPACK diagonalization wrappers
   !
   !  Revision 5.14  2011/03/28 20:58:46  souda
   !  changes caused by the changes in random_generators module
   !
   !  Revision 5.13  2011/03/01 23:54:03  souda
   !  Variable timestep for quantum propagation implemented (thanks to Sharon) -
   !  that fixes the problems with the conservation of the norm
   !  of the time-dependent wavefunction.
   !
   !  Revision 5.12  2011/02/24 00:56:36  souda
   !  - Additional routines for generating a coherent mixture for initial state in MDQT;
   !  - Franck-Condon factors are not renormalized initially;
   !  - initial wavefunction (coherent mixture) is renormalized.
   !
   !  Revision 5.11  2011/02/23 01:26:51  souda
   !  Another bug in FC part... Hopefully, last one.
   !
   !  Revision 5.10  2011/02/22 22:02:29  souda
   !  Bug fixed (yet another one) in the calculation of FC factors
   !
   !  Revision 5.9  2011/02/21 02:41:41  souda
   !  removed Debug printout
   !
   !  Revision 5.8  2011/02/21 02:40:37  souda
   !  Bug fixed in calculation of Franck-Condon probabilities
   !
   !  Revision 5.7  2011/02/20 00:58:11  souda
   !  Major additions/modifications:
   !  (1) precalculation of the proton vibrational basis functions for METHOD=1
   !  (2) Franck-Condon initial excitation added to DYNAMICS3
   !  (3) addition of the module timers: module_timers.f90 (changes in Makefile!)
   !
   !  Revision 5.6  2011/02/10 23:00:33  souda
   !  improvement (due to Anirban Hazra): avoiding the (redundant) calculation
   !  of vibronic states in the beginning of classical propagation step.
   !
   !  Revision 5.5  2011/02/08 00:46:50  souda
   !  unimportant bug in evaluating the kinetic energy for output - fixed
   !
   !  Revision 5.4  2010/11/10 21:14:21  souda
   !  Last addition/changes. IMPORTANT: Solvation keyword is back to SOLV( for compatibility reasons
   !
   !  Revision 5.3  2010/11/04 22:43:08  souda
   !  Next iteration... and two additional Makefiles for building the code with debug options.
   !
   !  Revision 5.2  2010/10/28 21:29:36  souda
   !  First (working and hopefully bug-free) source of PCET 5.x
   !
   !  Revision 5.1  2010/10/26 21:06:21  souda
   !  new routines/modules
   !
   !
   !---------------------------------------------------------------------

   use cst
   use parsol
   use solmat
   use quantum, only: nprst, npnts, wp_wavefunction
   use feszz_3d
   use random_generators
   use rk_parameters
   use laser

   !---------------------------------------------------------------------
   implicit none
   private

   complex(kind=8), parameter :: ii=(0.d0,1.d0)

   character(len=5) :: mode
   integer :: iset, nstates, nzdim, ielst

   !-- arrays for the energies and wavefunctions
   !------------------------------------------------------
   real(8), allocatable, dimension(:)     :: fe
   real(8), allocatable, dimension(:)     :: fe_prev
   real(8), allocatable, dimension(:,:)   :: z
   real(8), allocatable, dimension(:,:)   :: z_prev
   real(8), allocatable, dimension(:,:,:) :: psiel, psipr
   real(8), allocatable, dimension(:,:)   :: enel, envib

   !-- arrays for absorption probabilities (vibronic spectrum)
   !---------------------------------------------------------
   real(8), allocatable, dimension(:) :: absorption_prob
   real(8), allocatable, dimension(:) :: fc_prob

   !-- array for EVB weights
   !---------------------------------------------
   real(8), allocatable, dimension(:,:), public :: wght

   !-- arrays for nonadaiabatic coupling vectors
   !---------------------------------------------
   real(8), allocatable, dimension(:,:) :: coupz1
   real(8), allocatable, dimension(:,:) :: coupz2
   real(8), allocatable, dimension(:,:) :: coupz1_prev
   real(8), allocatable, dimension(:,:) :: coupz2_prev

   !-- array for the dot product of the classical velosity and coupling
   !-------------------------------------------------------------------
   real(8), allocatable, dimension(:,:) :: v_dot_d
   real(8), allocatable, dimension(:,:) :: v_dot_d_mid
   real(8), allocatable, dimension(:,:) :: v_dot_d_prev

   !-- Array for the renormalized complex amplitudes
   !---------------------------------------------------------
   complex(8), allocatable, dimension(:) :: amplitude
   complex(8), allocatable, dimension(:) :: amplitude_copy

   !-- Array for the density matrix
   !---------------------------------------------------------
   complex(8), allocatable, dimension(:,:) :: density_matrix

   !-- arrays for the MDQT transition probabilities for current state (b_jk)
   !---------------------------------------------------------------------
   real(8), allocatable, dimension(:) :: b_prob

   !-- Array for the switching probabilities from current state
   !-----------------------------------------------------------
   real(8), allocatable, dimension(:) :: switch_prob

   !-- interpolation coefficients for kinetic energy
   !---------------------------------------------------------
   real(8) :: a0_ekin, a1_ekin, a2_ekin

   !-- arrays with interpolation coefficients
   !---------------------------------------------------------
   real(8), allocatable, dimension(:) :: a0_fe, a1_fe, a2_fe
   real(8), allocatable, dimension(:,:) :: a0_vdotd, a1_vdotd, a2_vdotd

   public :: set_mode
   public :: get_free_energy
   public :: allocate_vibronic_states, deallocate_vibronic_states
   public :: allocate_evb_weights, deallocate_evb_weights
   public :: allocate_mdqt_arrays, deallocate_mdqt_arrays
   public :: deallocate_all_arrays
   public :: calculate_vibronic_states
   public :: calculate_absorption_prob
   public :: calculate_fc_prob
   public :: print_vibronic_spectrum
   public :: print_vibronic_spectrum_conv
   public :: calculate_vibronic_couplings
   public :: calculate_v_dot_d
   public :: interpolate_vdotd
   public :: interpolate_energy
   public :: interpolate_kinenergy
   public :: assign_initial_state
   public :: assign_fc_state
   public :: set_initial_amplitudes_pure
   public :: set_initial_amplitudes_laser
   public :: set_initial_amplitudes_fc
   public :: print_initial_amplitudes
   public :: set_initial_density_pure
   public :: set_initial_density_laser
   public :: set_initial_density_fc
   public :: calculate_v_dot_d_mid
   public :: calculate_bprob_amp
   public :: calculate_bprob_den
   public :: calculate_density_matrix
   public :: calculate_population
   public :: tdwf_norm
   public :: store_vibronic_couplings
   public :: store_vibronic_energies
   public :: store_wavefunctions
   public :: get_evb_weights
   public :: get_vibronic_coupling
   public :: langevin_debye_2d
   public :: langevin_onodera_2d
   public :: langevin_onodera2_2d
   public :: print_propagators_3d
   public :: tdse_derivatives
   public :: propagate_amplitudes_rk4
   public :: propagate_amplitudes_phcorr_rk4
   public :: propagate_density_rk4
   public :: switch_state
   public :: adjust_velocities
   public :: reset_switch_prob
   public :: accumulate_switch_prob
   public :: normalize_switch_prob
   public :: save_amplitudes
   public :: restore_amplitudes

contains

   !===DEBUG printout===
   subroutine print_propagators_3d
      integer :: k
      write(*,*)
      write(*,*) "-----Contents of the propagators_3d module"
      write(*,*)
      write(*,*) "mode:    ",mode
      write(*,*) "iset:    ",iset
      write(*,*) "nstates: ",nstates
      write(*,*) "nzdim:   ",nzdim
      write(*,*) "ielst:   ",ielst
      write(*,*)
      write(*,*) "Splittings: ",(fe(k)-fe(k-1),k=2,nstates)
      write(*,*)
      write(*,*) "Vibronic couplings z1: ",coupz1
      write(*,*)
      write(*,*) "Vibronic couplings z2: ",coupz2
      write(*,*)
      write(*,*) "V_dot_d: ",v_dot_d
      write(*,*)
      write(*,*) "a0_fe: ",a0_fe
      write(*,*)
      write(*,*) "a1_fe: ",a1_fe
      write(*,*)
      write(*,*) "a2_fe: ",a2_fe
      write(*,*)
      write(*,*) "a0_vdotd: ",a0_vdotd
      write(*,*)
      write(*,*) "a1_vdotd: ",a1_vdotd
      write(*,*)
      write(*,*) "a2_vdotd: ",a2_vdotd
      write(*,*)
      write(*,*) "TDSE amplitudes: ",amplitude
      write(*,*)
   end subroutine print_propagators_3d

   subroutine set_mode(mode_,iset_,nstates_,nzdim_,ielst_)
      character(len=5), intent(in) :: mode_
      integer, intent(in) :: iset_, nstates_, nzdim_, ielst_
      mode = mode_
      iset = iset_
      nstates = nstates_
      nzdim = nzdim_
      ielst = ielst_
   end subroutine set_mode

   function get_free_energy(i_) result(fe_)
      integer, intent(in) :: i_
      real(8) :: fe_
      integer :: isize
      isize = size(fe)
      if (i_.gt.0.and.i_.le.isize) then
         fe_ = fe(i_)
      else
         fe_ = 0.d0
         write(*,*) "ERROR in get_free_energy: index of state (",i_,") is out of bounds"
         call deallocate_all_arrays
         stop
      endif
   end function get_free_energy

   !---------------------------------------------------------------------
   ! array allocation routines
   !---------------------------------------------------------------------

   subroutine allocate_vibronic_states
      integer :: nz
      nz = ielst*nprst
      allocate (fe(nstates))
      allocate (fe_prev(nstates))
      allocate (absorption_prob(nstates))
      allocate (fc_prob(nstates))
      allocate (z(nz,nz))
      allocate (z_prev(nz,nz))
      allocate (psiel(ielst,npnts,ielst))
      allocate (psipr(ielst,nprst,npnts))
      allocate (enel(ielst,npnts))
      allocate (envib(ielst,nprst))
   end subroutine allocate_vibronic_states

   subroutine allocate_evb_weights
      allocate (wght(ielst,nstates))
   end subroutine allocate_evb_weights

   subroutine allocate_mdqt_arrays
      allocate (coupz1(nstates,nstates))
      allocate (coupz2(nstates,nstates))
      allocate (coupz1_prev(nstates,nstates))
      allocate (coupz2_prev(nstates,nstates))
      allocate (amplitude(nstates),amplitude_copy(nstates))
      allocate (density_matrix(nstates,nstates))
      allocate (switch_prob(nstates))
      allocate (b_prob(nstates))
      allocate (v_dot_d(nstates,nstates))
      allocate (v_dot_d_mid(nstates,nstates))
      allocate (v_dot_d_prev(nstates,nstates))
      allocate (a0_fe(nstates),a1_fe(nstates),a2_fe(nstates))
      allocate (a0_vdotd(nstates,nstates))
      allocate (a1_vdotd(nstates,nstates))
      allocate (a2_vdotd(nstates,nstates))
   end subroutine allocate_mdqt_arrays

   !---------------------------------------------------------------------
   ! array deallocation routines
   !---------------------------------------------------------------------

   subroutine deallocate_vibronic_states
      if (allocated(fe))      deallocate (fe)
      if (allocated(fe_prev)) deallocate (fe_prev)
      if (allocated(z))       deallocate (z)
      if (allocated(z_prev))  deallocate (z_prev)
      if (allocated(psiel))   deallocate (psiel)
      if (allocated(psipr))   deallocate (psipr)
      if (allocated(enel))    deallocate (enel)
      if (allocated(envib))   deallocate (envib)
      if (allocated(absorption_prob)) deallocate (absorption_prob)
      if (allocated(fc_prob)) deallocate (fc_prob)
   end subroutine deallocate_vibronic_states

   subroutine deallocate_evb_weights
      if (allocated(wght)) deallocate (wght)
   end subroutine deallocate_evb_weights

   subroutine deallocate_mdqt_arrays
      if (allocated(coupz1))            deallocate (coupz1)
      if (allocated(coupz2))            deallocate (coupz2)
      if (allocated(coupz1_prev))       deallocate (coupz1_prev)
      if (allocated(coupz2_prev))       deallocate (coupz2_prev)
      if (allocated(amplitude))         deallocate (amplitude)
      if (allocated(amplitude_copy))    deallocate (amplitude_copy)
      if (allocated(density_matrix))    deallocate (density_matrix)
      if (allocated(b_prob))            deallocate (b_prob)
      if (allocated(switch_prob))       deallocate (switch_prob)
      if (allocated(v_dot_d))           deallocate (v_dot_d)
      if (allocated(v_dot_d_mid))       deallocate (v_dot_d_mid)
      if (allocated(v_dot_d_prev))      deallocate (v_dot_d_prev)
      if (allocated(a0_fe))             deallocate (a0_fe)
      if (allocated(a1_fe))             deallocate (a1_fe)
      if (allocated(a2_fe))             deallocate (a2_fe)
      if (allocated(a0_vdotd))          deallocate (a0_vdotd)
      if (allocated(a1_vdotd))          deallocate (a1_vdotd)
      if (allocated(a2_vdotd))          deallocate (a2_vdotd)
   end subroutine deallocate_mdqt_arrays

   subroutine deallocate_all_arrays
      call deallocate_vibronic_states
      call deallocate_evb_weights
      call deallocate_mdqt_arrays
   end subroutine deallocate_all_arrays

   !--------------------------------------------------------------------
   !-- calculate absorption spectrum
   !   (oscillator strengths of vibronic transitions)
   !--------------------------------------------------------------------
   subroutine calculate_absorption_prob(iground,kg,z1,z2)

      integer, intent(in) :: iground, kg
      real(8), intent(in) :: z1, z2

      integer :: ndabf, k, i, j, mu, nu, imu, jnu, kp, i0, j0
      real(8) :: zp, ze, s, dx, dy, dz, davx, davy, davz
      real(8) :: pnorm

      !-- transform to zp,ze frame
      call z1z2_to_zpze(z1,z2,zp,ze)

      !-- calculate vibronic states
      call feszz3(mode,iset,kg,zp,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)

      !-- calculate transition dipole moments between the state iground and the higher states

      absorption_prob = 0.d0
      pnorm = 0.d0

      do k=iground+1,nstates

         dx = 0.d0
         dy = 0.d0
         dz = 0.d0

         imu = 0
         do i=1,ielst
            i0 = i + 2*iset - 2
            do mu=1,nprst

               imu = imu + 1

               jnu = 0
               do j=1,ielst
                  j0 = j + 2*iset - 2
                  do nu=1,nprst

                     jnu = jnu + 1

                     !-- here we assume that the dipole moment might have non-zero off-diagonal elements
                     !   even in the basis of the diabatic electronic states

                     davx = 0.d0
                     davy = 0.d0
                     davz = 0.d0
                     do kp=1,npnts
                        davx = davx + dipole_moment_diab_x(i0,j0,kp)*psipr(i,mu,kp)*psipr(j,nu,kp)
                        davy = davy + dipole_moment_diab_y(i0,j0,kp)*psipr(i,mu,kp)*psipr(j,nu,kp)
                        davz = davz + dipole_moment_diab_z(i0,j0,kp)*psipr(i,mu,kp)*psipr(j,nu,kp)
                     enddo

                     dx = dx + z(imu,iground)*z(jnu,k)*davx
                     dy = dy + z(imu,iground)*z(jnu,k)*davy
                     dz = dz + z(imu,iground)*z(jnu,k)*davz

                  enddo 
               enddo

            enddo
         enddo

         absorption_prob(k) = dx*dx + dy*dy + dz*dz
         pnorm = pnorm + absorption_prob(k)

      enddo

      !-- normalization
      absorption_prob = absorption_prob/pnorm

   end subroutine calculate_absorption_prob

   !--------------------------------------------------------------------
   !-- Output the vibronic spectrum (energies and oscillator strength)
   !--------------------------------------------------------------------
   subroutine print_vibronic_spectrum(ichannel,iground)

      integer, intent(in) :: ichannel, iground

      integer :: i
      real(8) :: excitation_energy, pnorm

      !-- print the header of the table

      write(ichannel,'("#",71("="))')
      write(ichannel,'("#","     Vibronic spectrum: excitations from state ",i3," to ",i3," higher states")') iground, nstates-iground
      write(ichannel,'("#",71("-"))')
      write(ichannel,'("#",t7,"Excitation energy (eV)",t35,"Normalized intensity")')
      write(ichannel,'("#",71("-"))')

      pnorm = 0.d0
      do i=iground+1,nstates
         write(ichannel,'(i5,f20.6,t35,f15.6)') i, (fe(i) - fe(iground))*cal2ev, absorption_prob(i)
         pnorm = pnorm + absorption_prob(i)
      enddo

      write(ichannel,'("#",71("-"))')
      write(ichannel,'("#",t30,"Norm:",f15.6)') pnorm
      write(ichannel,'("#",71("="))')

   end subroutine print_vibronic_spectrum

   !--------------------------------------------------------------------
   !-- Output the convoluted vibronic spectrum
   !--------------------------------------------------------------------
   subroutine print_vibronic_spectrum_conv(ichannel,iground,iconv,width)

      integer, intent(in) :: ichannel, iground, iconv
      real(8), intent(in) :: width

      character(len=10), parameter, dimension(2) :: conv=(/"Lorentzian","Gaussian"/)
      integer, parameter :: npoints = 1000

      integer :: i, k
      real(8) :: energy, exc_en_min, exc_en_max, range, den, p
      real(8), dimension(nstates) :: excitation_energy

      if (iconv.lt.1.or.iconv.gt.2) then
         write(ichannel,'("#  Unknown convoluting function: no output")')
         return
      endif

      !-- print the header of the table

      write(ichannel,'("#",71("="))')
      write(ichannel,'("#","     Vibronic spectrum: ",a," convolution")') trim(conv(iconv))
      write(ichannel,'("#",71("-"))')
      write(ichannel,'("#",t7,"Excitation energy (eV)",t35,"Normalized intensity")')
      write(ichannel,'("#",71("-"))')

      !-- set the excitation energy range

      exc_en_min = cal2ev*(fe(iground+1) - fe(iground))
      exc_en_max = cal2ev*(fe(nstates) - fe(iground))

      range = exc_en_max - exc_en_min
      exc_en_min = exc_en_min - 5.d0*width
      if (exc_en_min.lt.0.d0) exc_en_min = 0.d0
      exc_en_max = exc_en_max + 5.d0*width
      range = exc_en_max - exc_en_min
      den = range/(npoints-1)

      !-- precompute excitation energies (in eV)
      do i=iground+1,nstates
         excitation_energy(i) = (fe(i) - fe(iground))*cal2ev
      enddo

      !-- Loop over energy grid points
      
      do k=1,npoints

         energy = exc_en_min + (k-1)*den
         
         !-- evaluate convolution
         
         p = 0.d0

         select case(iconv)

            case(1)
            do i=iground+1,nstates
               p = p + absorption_prob(i)*lorentzian(energy,excitation_energy(i),width)
            enddo

            case(2)
            do i=iground+1,nstates
               p = p + absorption_prob(i)*gaussian(energy,excitation_energy(i),width)
            enddo

         end select

         write(ichannel,'(f20.6,t35,f15.6)') energy, p

      enddo

      write(ichannel,'("#",71("="))')

   end subroutine print_vibronic_spectrum_conv

   !--------------------------------------------------------------------
   !-- Assign initial state according to the normalized magnitude
   !   of the transition dipole moment matrix element multiplied
   !   by the laser pulse intensity
   !--------------------------------------------------------------------
   function assign_initial_state(iground) result(istate)

      integer, intent(in) :: iground
      integer :: istate

      integer :: i
      real(4) :: r
      real(8) :: s, excitation_energy, pnorm
      real(8), dimension(nstates) :: p

      istate = 0

      p = 0.d0
      pnorm = 0.d0

      do i=iground+1,nstates
         excitation_energy = (fe(i) - fe(iground))*cal2ev
         p(i) = absorption_prob(i)*spectral_shape(excitation_energy)
         pnorm = pnorm + p(i)
      enddo

      !-- normalize p
      p = p/pnorm

      !-- generate a random number between 0 and 1
      !   from a uniform distribution
      r = ran2nr()

      s = 0.d0
      
      do i=iground+1,nstates

         s = s + p(i)

         if (s.gt.r) then
            istate = i
            exit
         endif

      enddo

   end function assign_initial_state

   !--------------------------------------------------------------------
   !-- calculate Franck-Condon transition probabilities
   !--------------------------------------------------------------------
   subroutine calculate_fc_prob(fc_state,kg,z1,z2)

      integer, intent(in) :: fc_state, kg
      real(8), intent(in) :: z1, z2

      integer :: ndabf, n, i, i0, ifc, mu, imu, imu_start, kp
      real(8) :: zp, ze, fc_ovlp, coef, pnorm

      fc_prob = 0

      !-- Check if the specified photoexcited diabatic electronic state
      !   (fc_state) is present in the electronic basis.
      !   If not then fc_prob = 0

      ifc = 0
      do i=1,ielst
         i0 = i + 2*iset - 2
         if (i0.eq.fc_state) then
            ifc = i
            imu_start = (ifc-1)*nprst
            exit
         endif
      enddo

      if (ifc.eq.0) return

      !-- transform to zp,ze frame
      call z1z2_to_zpze(z1,z2,zp,ze)

      !-- calculate vibronic states
      call feszz3(mode,iset,kg,zp,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)

      pnorm = 0.d0

      do n=1,nstates

         coef = 0.d0
         imu = imu_start

         do mu=1,nprst

            imu = imu + 1

            fc_ovlp = 0.d0
            do kp=1,npnts
               fc_ovlp = fc_ovlp + wp_wavefunction(kp)*psipr(ifc,mu,kp)
            enddo

            coef = coef + z(imu,n)*fc_ovlp

         enddo

         fc_prob(n) = coef*coef

         pnorm = pnorm + fc_prob(n)

      enddo

      !-- Print normalized Franck-Condon probabilities

      write(*,'(/1x,130("="))')
      write(*,'( 1x," Franck-Condon factors (probabilities) for the initial state")')
      write(*,'( 1x,130("-"))')
      write(*,'( 1x," Norm of the initial wavefunction before renormalization: ",g20.10)') pnorm
      write(*,'( 1x,130("-"))')

      do n=1,nstates/10
         write(6,'(/1x,10i13)') (i,i=(n-1)*10+1,n*10)
         write(6,'(1x,10f13.6)') (fc_prob(i),i=(n-1)*10+1,n*10)
      enddo

      n = mod(nstates,10)
      if (n.gt.0) then
         write(6,'(/1x,10i13)') (i,i=nstates-n+1,nstates)
         write(6,'(1x,10f13.6)') (fc_prob(i),i=nstates-n+1,nstates)
      endif

      write(*,'(/1x,130("-")/)')

      !-- renormalization
      !   (we don't normalize here: normalize the initial wavefunction later)
      !fc_prob = fc_prob/pnorm

      !-- stop the program if the norm is less than 0.99 or greater than 1 (unphysical situation)
      
      if (pnorm.lt.0.98d0) then
         write(*,'(/1x,"WARNING: Bad normalization (<0.99): increase the number of states NSTATES")')
         call deallocate_all_arrays
         stop
      elseif (pnorm-1.d0.gt.1.d-8) then
         write(*,'(/1x,"WARNING: Bad normalization (>1.0): increase the accuracy of the wavefunctions (NPRST or NPNTS)")')
         call deallocate_all_arrays
         stop
      endif

   end subroutine calculate_fc_prob

   !--------------------------------------------------------------------
   !-- Assign initial state according to the normalized
   !   Franck-Condon probabilities
   !--------------------------------------------------------------------
   function assign_fc_state() result(istate)

      integer :: istate

      integer :: i
      real(4) :: r
      real(8) :: s

      istate = 0

      !-- generate a random number between 0 and 1
      !   from a uniform distribution
      r = ran2nr()

      s = 0.d0
      
      do i=1,nstates

         s = s + fc_prob(i)

         if (s.gt.r) then
            istate = i
            exit
         endif

      enddo

   end function assign_fc_state

   !--------------------------------------------------------------------
   !-- Set initial amplitudes to correspond to a pure adiabatic state
   !--------------------------------------------------------------------
   subroutine set_initial_amplitudes_pure(istate)
      integer, intent(in) :: istate
      amplitude = 0.d0
      amplitude(istate) = cmplx(1.d0,0.d0)
   end subroutine set_initial_amplitudes_pure

   !--------------------------------------------------------------------
   !-- Set initial amplitudes to correspond to a coherent mixture
   !   formed after a laser pulse excitation
   !--------------------------------------------------------------------
   subroutine set_initial_amplitudes_laser
      integer :: i
      amplitude = 0.d0
      do i=1,nstates
         amplitude(i) = cmplx(sqrt(absorption_prob(i)),0.d0)
      enddo
   end subroutine set_initial_amplitudes_laser

   !--------------------------------------------------------------------
   !-- Set initial amplitudes to correspond to a coherent mixture
   !   formed after a Franck-Condon excitation
   !--------------------------------------------------------------------
   subroutine set_initial_amplitudes_fc
      integer :: i
      real(8) :: wfnorm
      amplitude = 0.d0
      wfnorm = 0.d0
      do i=1,nstates
         amplitude(i) = cmplx(sqrt(fc_prob(i)),0.d0)
         wfnorm = wfnorm + fc_prob(i)
      enddo
      amplitude = amplitude/sqrt(wfnorm)
   end subroutine set_initial_amplitudes_fc

   !--------------------------------------------------------------------
   !-- Print initial amplitudes of the time-dependent wavefunction
   !--------------------------------------------------------------------
   subroutine print_initial_amplitudes(ichannel)

      integer, intent(in) :: ichannel
      integer :: i, n
      real(8) :: pnorm, reamp

      write(ichannel,'("#",100("="))')
      write(ichannel,'("# Initial amplitudes of the time-dependent wavefunction")')
      write(ichannel,'("#",100("-"))')

      pnorm = 0.d0
      do i=1,nstates
         reamp = real(amplitude(i))
         pnorm = pnorm + reamp*reamp
      enddo

      do n=1,nstates/10
         write(ichannel,'(/1x,10i13)') (i,i=(n-1)*10+1,n*10)
         write(ichannel,'(1x,10f13.6)') (real(amplitude(i)),i=(n-1)*10+1,n*10)
      enddo

      n = mod(nstates,10)
      if (n.gt.0) then
         write(ichannel,'(/1x,10i13)') (i,i=nstates-n+1,nstates)
         write(ichannel,'(1x,10f13.6)') (real(amplitude(i)),i=nstates-n+1,nstates)
      endif
      write(6,'(/)')

      write(ichannel,'("#",130("-"))')
      write(ichannel,'("#",t10,"Norm:",f15.6)') pnorm
      write(ichannel,'("#",130("="))')

   end subroutine print_initial_amplitudes

   !---------------------------------------------------------------------
   !-- save amplitudes
   !---------------------------------------------------------------------
   subroutine save_amplitudes
      amplitude_copy = amplitude   !array operation
   end subroutine save_amplitudes

   !---------------------------------------------------------------------
   !-- restore amplitudes
   !---------------------------------------------------------------------
   subroutine restore_amplitudes
      amplitude = amplitude_copy   !array operation
   end subroutine restore_amplitudes

   !---------------------------------------------------------------------
   !-- Set initial density matrix to correspond to a pure adiabatic state
   !---------------------------------------------------------------------
   subroutine set_initial_density_pure(istate)
      integer, intent(in) :: istate
      density_matrix = 0.d0
      density_matrix(istate,istate) = cmplx(1.d0,0.d0)
   end subroutine set_initial_density_pure

   !---------------------------------------------------------------------
   !-- Set initial density matrix to correspond to a coherent mixture
   !   formed after a laser pulse excitation
   !---------------------------------------------------------------------
   subroutine set_initial_density_laser
      integer :: i, j
      real(8) :: rhoij
      density_matrix = 0.d0
      do i=1,nstates
         do j=1,nstates
            rhoij = sqrt(absorption_prob(i))*sqrt(absorption_prob(j))
            density_matrix(i,j) = cmplx(rhoij,0.d0)
         enddo
      enddo
   end subroutine set_initial_density_laser

   !---------------------------------------------------------------------
   !-- Set initial density matrix to correspond to a coherent mixture
   !   formed after a laser pulse excitation
   !---------------------------------------------------------------------
   subroutine set_initial_density_fc

      integer :: i, j
      real(8) :: rhoij, wfnorm

      density_matrix = 0.d0

      !-- calculate the norm of the initial wavefunction
      wfnorm = 0.d0
      do i=1,nstates
         wfnorm = wfnorm + fc_prob(i)
      enddo

      do i=1,nstates
         do j=1,nstates
            rhoij = sqrt(fc_prob(i))*sqrt(fc_prob(j))/wfnorm
            density_matrix(i,j) = cmplx(rhoij,0.d0)
         enddo
      enddo

   end subroutine set_initial_density_fc

   !--------------------------------------------------------------------
   !-- Interface to the calculation of the vibronic states
   !--------------------------------------------------------------------
   subroutine calculate_vibronic_states(kg,z1,z2)
      integer, intent(in) :: kg
      real(8), intent(in) :: z1, z2
      integer :: ndabf
      real(8) :: zp, ze
      !-- transform to zp,ze frame
      call z1z2_to_zpze(z1,z2,zp,ze)
      !-- calculate vibronic states
      call feszz3(mode,iset,kg,zp,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
   end subroutine calculate_vibronic_states

   !--------------------------------------------------------------------
   !-- Interface to the calculation of the gradients
   !--------------------------------------------------------------------
   subroutine calculate_gradient(istate,g1,g2)
      integer, intent(in) :: istate
      real(8), intent(out) :: g1, g2
      real(8) :: gp, ge
      !-- calculate gradients in zp,ze frame
      call dvdzz3(mode,iset,istate,nzdim,z,ielst,gp,ge)
      !-- transform gradients
      call gpge_to_g1g2(gp,ge,g1,g2)
   end subroutine calculate_gradient

   !--------------------------------------------------------------------
   !-- Interface to the calculation of the vibronic couplings
   !   for current vibronic wavefunctions
   !--------------------------------------------------------------------
   subroutine calculate_vibronic_couplings
      real(8), dimension(nstates,nstates) :: coupzp, coupze
      real(8) :: dzp, dze, dz1, dz2
      integer :: i, j
      call coupzz3(mode,iset,nstates,fe,nzdim,z,ielst,coupzp,coupze)
      !-- transform the couplings to z1,z2 frame
      do i=1,nstates
         do j=i+1,nstates
            dzp = coupzp(i,j)
            dze = coupze(i,j)
            call gpge_to_g1g2(dzp,dze,dz1,dz2)
            coupz1(i,j) =  dz1
            coupz2(i,j) =  dz2
            coupz1(j,i) = -dz1
            coupz2(j,i) = -dz2
         enddo
      enddo
   end subroutine calculate_vibronic_couplings

   !--------------------------------------------------------------------
   !-- extract vibronic coupling vector
   !--------------------------------------------------------------------
   function get_vibronic_coupling(i_,j_) result(coupz)
      integer, intent(in) :: i_, j_
      real(8), dimension(2) :: coupz
      coupz(1) = coupz1(i_,j_)
      coupz(2) = coupz2(i_,j_)
   end function get_vibronic_coupling

   !--------------------------------------------------------------------
   !-- store vibronic couplings
   !   for current vibronic wavefunctions
   !--------------------------------------------------------------------
   subroutine store_vibronic_couplings
      coupz1_prev = coupz1
      coupz2_prev = coupz2
   end subroutine store_vibronic_couplings

   !--------------------------------------------------------------------
   !-- store energies of vibronic states
   !--------------------------------------------------------------------
   subroutine store_vibronic_energies
      fe_prev = fe
   end subroutine store_vibronic_energies

   !--------------------------------------------------------------------
   !-- calculate population (diagonal element of the density matrix)
   !--------------------------------------------------------------------
   function calculate_population(istate) result(pop)
      integer, intent(in) :: istate
      real(8) :: pop
      pop = amplitude(istate)*conjg(amplitude(istate))
   end function calculate_population

   !--------------------------------------------------------------------
   !-- calculate norm of the time-dependent wavefunction
   !--------------------------------------------------------------------
   function tdwf_norm() result(wfnorm)
      real(8) :: wfnorm
      integer :: i
      wfnorm = 0.d0
      do i=1,nstates
      	wfnorm = wfnorm + amplitude(i)*conjg(amplitude(i))
      enddo
   end function tdwf_norm

   !--------------------------------------------------------------------
   !-- calculate trace of the density matrix
   !--------------------------------------------------------------------
   function density_trace() result(trace)
      complex(8) :: ctrace
      real(8) :: trace
      integer :: i
      ctrace = 0.d0
      do i=1,nstates
      	ctrace = ctrace + density_matrix(i,i)
      enddo
      trace = real(ctrace)
   end function density_trace

   !--------------------------------------------------------------------
   !-- construct the full density matrix
   !--------------------------------------------------------------------
   subroutine calculate_density_matrix
      integer :: i, j
      do i=1,nstates
         do j=1,nstates
            density_matrix(i,j) = amplitude(i)*conjg(amplitude(j))
         enddo
      enddo
   end subroutine calculate_density_matrix

   !--------------------------------------------------------------------
   !-- calculate transition probabilities b_jk
   !   from renormalized amplitudes
   !--------------------------------------------------------------------
   subroutine calculate_bprob_amp(istate,t_)
      integer, intent(in) :: istate
      real(8), intent(in) :: t_
      real(8) :: vdji
      integer :: j
      do j=1,nstates
         vdji = a0_vdotd(j,istate) + a1_vdotd(j,istate)*t_ + a2_vdotd(j,istate)*t_*t_
         b_prob(j) = -2.d0*Real(amplitude(j)*conjg(amplitude(istate))*vdji)
      enddo
   end subroutine calculate_bprob_amp

   !--------------------------------------------------------------------
   !-- calculate transition probabilities b_jk
   !   from density matrix
   !--------------------------------------------------------------------
   subroutine calculate_bprob_den(istate,t_)
      integer, intent(in) :: istate
      real(8), intent(in) :: t_
      real(8) :: vdji
      integer :: j
      do j=1,nstates
         vdji = a0_vdotd(j,istate) + a1_vdotd(j,istate)*t_ + a2_vdotd(j,istate)*t_*t_
         b_prob(j) = -2.d0*Real(density_matrix(j,istate)*vdji)
      enddo
   end subroutine calculate_bprob_den

   !--------------------------------------------------------------------
   !-- zero out switching probabilities
   !--------------------------------------------------------------------
   subroutine reset_switch_prob
      switch_prob = 0.d0
   end subroutine reset_switch_prob

   !--------------------------------------------------------------------
   !-- accumulate switching probabilities
   !--------------------------------------------------------------------
   subroutine accumulate_switch_prob(qtstep_)
      real(8), intent(in) :: qtstep_
      switch_prob = switch_prob + b_prob*qtstep_   ! (array operation)
   end subroutine accumulate_switch_prob

   !--------------------------------------------------------------------
   !-- normalize and clean switching probabilities (if <0, set to zero)
   !--------------------------------------------------------------------
   subroutine normalize_switch_prob(factor)
      real(8), intent(in) :: factor
      switch_prob = switch_prob/(factor+1.d-10)
      where (switch_prob < 0.d0) switch_prob = 0.d0
   end subroutine normalize_switch_prob

   !--------------------------------------------------------------------
   !-- store expansion coefficients of vibronic states
   !--------------------------------------------------------------------
   subroutine store_wavefunctions
      z_prev = z
   end subroutine store_wavefunctions

   !--------------------------------------------------------------------
   !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
   !   at t and t+dt
   !--------------------------------------------------------------------
   subroutine calculate_v_dot_d(vz1,vz1_prev,vz2,vz2_prev)
      real(8), intent(in) :: vz1, vz1_prev, vz2, vz2_prev
      real(8) :: vd, vd_prev
      integer :: i, j
      v_dot_d = 0.d0
      v_dot_d_prev = 0.d0
      do i=1,nstates
         do j=i+1,nstates
            vd = vz1*coupz1(i,j) + vz2*coupz2(i,j)
            vd_prev = vz1_prev*coupz1_prev(i,j) + vz2_prev*coupz2_prev(i,j)
            v_dot_d(i,j) =  vd
            v_dot_d(j,i) = -vd
            v_dot_d_prev(i,j) =  vd_prev
            v_dot_d_prev(j,i) = -vd_prev
         enddo
      enddo
   end subroutine calculate_v_dot_d

   !--------------------------------------------------------------------
   !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
   !   at the half step t+dt/2
   !--------------------------------------------------------------------
   subroutine calculate_v_dot_d_mid(tstep_)
      real(8), intent(in) :: tstep_
      real(8) :: vd
      integer :: i, j, k
      v_dot_d_mid = 0.d0
      do i=1,nstates
         do j=i+1,nstates
            vd = 0.d0
            do k=1,nzdim
               vd  = vd + z_prev(k,i)*z(k,j) - z(k,i)*z_prev(k,j)
            enddo
            vd = vd/(2.d0*tstep_)
            v_dot_d_mid(i,j) =  vd
            v_dot_d_mid(j,i) = -vd
         enddo
      enddo
   end subroutine calculate_v_dot_d_mid

   !--------------------------------------------------------------------
   !-- calculate linear interpolation coefficients
   !   for the adiabatic energies
   !--------------------------------------------------------------------
   subroutine interpolate_energy(t_prev_,t_)
      real(8), intent(in) :: t_prev_, t_
      real(8) :: dt0, f1, f2, t01, t02
      integer :: i
      t01 = t_prev_
      t02 = t_
      dt0 = t_ - t_prev_
      do i=1,nstates
         f2 = fe(i)
         f1 = fe_prev(i)
         a0_fe(i) = (t02*f1 - t01*f2)/dt0
         a1_fe(i) = (f2 - f1)/dt0
         a2_fe(i) = 0.d0
      enddo
   end subroutine interpolate_energy


   !--------------------------------------------------------------------
   !-- calculate linear interpolation coefficients
   !   for the velocities (for phase-corrected scheme)
   !--------------------------------------------------------------------
   subroutine interpolate_kinenergy(interpolation,t_prev_,t_,ekin_,ekin_prev_,ekin_half_)

      character(len=*), intent(in) :: interpolation
      real(8), intent(in) :: t_prev_, t_, ekin_, ekin_prev_, ekin_half_
      real(8) :: x0, x1, x2, x01, x02, x12
      real(8) :: y0, y1, y2, y01, y02, y12
      real(8) :: xdenom, a0, a1, a2

      x0 = t_prev_
      x2 = t_
      x1 = 0.5d0*(x0 + x2)
      
      x01 = x0 - x1
      x02 = x0 - x2
      x12 = x1 - x2
      xdenom = x01*x02*x12

      if (interpolation.eq."QUADRATIC") then

         y0 = ekin_prev_
         y1 = ekin_half_
         y2 = ekin_
         y01 = y0 - y1
         y02 = y0 - y2
         y12 = y1 - y2
         
         a0 = (y0*x1*x2*x12 - x0*y1*x2*x02 + x0*x1*y2*x01)/xdenom
         a1 = (x2*x2*y01 - x1*x1*y02 + x0*x0*y12)/xdenom
         a2 = (-x2*y01 + x1*y02 - x0*y12)/xdenom

         a0_ekin = a0
         a1_ekin = a1
         a2_ekin = a2

      else

         y0 = ekin_prev_
         y2 = ekin_
         y02 = y0 - y2
         
         a0 = (x0*y2 - x2*y0)/x02
         a1 = y02/x02
         a2 = 0.d0

         a0_ekin =  a0
         a1_ekin =  a1
         a2_ekin =  a2

      endif

   end subroutine interpolate_kinenergy



   !--------------------------------------------------------------------
   !-- calculate interpolation coefficients
   !   for the nonadiabatic coupling terms v*d_{kl}
   !--------------------------------------------------------------------
   subroutine interpolate_vdotd(interpolation,t_prev_,t_)

      character(len=*), intent(in) :: interpolation
      real(8), intent(in) :: t_prev_, t_
      integer :: i, j
      real(8) :: x0, x1, x2, x01, x02, x12
      real(8) :: y0, y1, y2, y01, y02, y12
      real(8) :: xdenom, a0, a1, a2
      
      x0 = t_prev_
      x2 = t_
      x1 = 0.5d0*(x0 + x2)
      
      x01 = x0 - x1
      x02 = x0 - x2
      x12 = x1 - x2
      xdenom = x01*x02*x12
      
      if (interpolation.eq."QUADRATIC") then
      
         !-- Quadratic interpolation

         do i=1,nstates
            do j=i+1,nstates

               y0 = v_dot_d_prev(i,j)
               y1 = v_dot_d_mid(i,j)
               y2 = v_dot_d(i,j)
               y01 = y0 - y1
               y02 = y0 - y2
               y12 = y1 - y2
            
               a0 = (y0*x1*x2*x12 - x0*y1*x2*x02 + x0*x1*y2*x01)/xdenom
               a1 = (x2*x2*y01 - x1*x1*y02 + x0*x0*y12)/xdenom
               a2 = (-x2*y01 + x1*y02 - x0*y12)/xdenom

               a0_vdotd(i,j) =  a0
               a0_vdotd(j,i) = -a0
               a1_vdotd(i,j) =  a1
               a1_vdotd(j,i) = -a1
               a2_vdotd(i,j) =  a2
               a2_vdotd(j,i) = -a2

            enddo
         enddo

      else

         !-- Linear interpolation

         do i=1,nstates
            do j=i+1,nstates

               y0 = v_dot_d_prev(i,j)
               y2 = v_dot_d(i,j)
               y02 = y0 - y2
            
               a0 = (x0*y2 - x2*y0)/x02
               a1 = y02/x02
               a2 = 0.d0

               a0_vdotd(i,j) =  a0
               a0_vdotd(j,i) = -a0
               a1_vdotd(i,j) =  a1
               a1_vdotd(j,i) = -a1
               a2_vdotd(i,j) =  a2
               a2_vdotd(j,i) =  a2

            enddo
         enddo


      endif

   end subroutine interpolate_vdotd

   !--------------------------------------------------------------------
   !-- Interface to the calculation of the EVB weights
   !   for current vibronic wavefunctions
   !--------------------------------------------------------------------
   subroutine get_evb_weights
      call evbwei(mode,iset,nstates,nzdim,z,ielst,psiel,psipr,wght)
   end subroutine get_evb_weights

   !--------------------------------------------------------------------
   !-- Runge-Cutta 2-nd order for overdamped Langevin equation
   !   [Honeycutt, Phys. Rev. A, 1992, 45, 600]
   !--------------------------------------------------------------------
   subroutine langevin_debye_2d(istate,kg,z1,z2,vz1,vz2,dt,temp,ekin1,ekin2,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, vz1, vz2
      real(8), intent(out)   :: ekin1, ekin2, efes

      integer :: i, ndabf
      real(8) :: d, psi, psi_factor             !, zp, ze, gp, ge
      real(8), dimension(2) :: x, v, xnew, dx
      real(8), dimension(2) :: fr, f1, f2, g

      d = kb*temp/f0/taul
      psi_factor = sqrt(2.d0*d*dt)

      x(1) = z1
      x(2) = z2
      v(1) = vz1
      v(2) = vz2

      !-- calculate vibronic states and vibronic free energies
      !   REUSE THE VIBRONIC STATES FROM THE END OF THE PREVIOUS TIMESTEP
      !---(AVS)-call calculate_vibronic_states(kg,x(1),x(2))-------------

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      do i=1,2
         psi = gaussdist_boxmuller()
         f1(i) = -x(i)/taul
         !-- add forces from the vibronic surface
         f1(i) = f1(i) - g(i)/(f0*taul)
         fr(i) = psi_factor*psi
         xnew(i) = x(i) + f1(i)*dt + fr(i)
      enddo

      !-- calculate vibronic states and vibronic free energies at new positions
      call calculate_vibronic_states(kg,xnew(1),xnew(2))

      !-- calculate new gradients
      call calculate_gradient(istate,g(1),g(2))

      !-- propagate positions and velocities
      do i=1,2
         f2(i) = -xnew(i)/taul
         !-- add forces from the vibronic surface
         f2(i) = f2(i) - g(i)/(f0*taul)
         dx(i) = 0.5d0*dt*(f1(i) + f2(i)) + fr(i)
         x(i) = x(i) + dx(i)
         v(i) = dx(i)/dt
      enddo

      z1 = x(1)
      z2 = x(2)
      vz1 = v(1)
      vz2 = v(2)

      !-- calculate vibronic states and vibronic free energies at final positions
      call calculate_vibronic_states(kg,x(1),x(2))

      !-- current free energy (PMF)
      efes = fe(istate)

      !-- calculate kinetic energies for both degrees of freedom
      ekin1 = half*f0*tau0*taul*v(1)*v(1)
      ekin2 = half*f0*tau0*taul*v(2)*v(2)

   end subroutine langevin_debye_2d

   !--------------------------------------------------------------------
   !-- Vanden-Eijnden and Ciccotti 2-nd order propagator
   !   for Onodera model
   !   [Chem. Phys. Lett., 2006, 429, 310-316, Eq.23]
   !--------------------------------------------------------------------
   subroutine langevin_onodera_2d(istate,kg,z1,z2,vz1,vz2,dt,temp,ekin1,ekin2,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, vz1, vz2
      real(8), intent(out)   :: ekin1, ekin2, efes

      real(8), parameter :: half=0.5d0
      real(8), parameter :: eighth=1.d0/8.d0
      real(8), parameter :: fourth=0.25d0
      real(8), parameter :: e32=1.5d0

      integer :: i, ndabf
      real(8) :: sq3
      real(8) :: gamma, sigma
      real(8), dimension(2) :: x, v, f, g
      real(8), dimension(2) :: ksi, eta, vhalf
      real(8) :: dt2, sqdt, dt32                       !, zp, ze, gp, ge

      !-- DEBUG variables
      !real(8) :: zinc, tr1b, tr2a, fplus, fminus, gzp, gze, gselfzp, gselfze

      sq3 = sqrt(3.d0)

      dt2 = dt*dt
      sqdt = sqrt(dt)
      dt32 = dt*sqdt

      x(1) = z1
      x(2) = z2
      v(1) = vz1
      v(2) = vz2

      !-- define Ciccotti constants
      
      gamma = (tau0 + tauD)/(tau0*tauD)
      sigma = sqrt(2*kb*temp*(tau0l+taul)/f0)/(tau0*taul)

      !-- generate two independent random variables ksi and eta
      !   from univariate Gaussian distributions

      do i=1,2
         ksi(i) = gaussdist_boxmuller()
         eta(i) = gaussdist_boxmuller()
      enddo

      !-- regular forces from the self-energy
      do i=1,2
         f(i) = -x(i)/(tau0*taul)
      enddo

      !-- calculate vibronic states and vibronic free energies
      !   REUSE THE VIBRONIC STATES FROM THE END OF THE PREVIOUS TIMESTEP
      !---(AVS)-call calculate_vibronic_states(kg,x(1),x(2))-------------

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      !=== DEBUG ===> Check numerical gradients ----------------------------------------------------
      !-zinc = 1.d-4
      !-tr1b = tr(1,2,32,1)
      !-tr2a = tr(1,3,32,1)
      !-call feszz3(mode,iset,kg,zp+zinc/2.d0,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fplus = fe(istate)
      !-call feszz3(mode,iset,kg,zp-zinc/2.d0,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fminus = fe(istate)
      !-gzp = (fplus-fminus)/zinc
      !-gselfzp = 0.5d0*(erx*(tr2a+ze) - eret*(tr1b+zp))/(erx*erx - erpt*eret)
      !-gzp = gzp - gselfzp
      !-call feszz3(mode,iset,kg,zp,ze+zinc/2.d0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fplus = fe(istate)
      !-call feszz3(mode,iset,kg,zp,ze-zinc/2.d0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fminus = fe(istate)
      !-gze = (fplus-fminus)/zinc
      !-gselfze = 0.5d0*(erx*(tr1b+zp) - erpt*(tr2a+ze))/(erx*erx - erpt*eret)
      !-gze = gze - gselfze
      !-write(*,*) "======== DEBUG =====> numerical vs. analytical gradient"
      !-write(*,'(1x,"Analytical (gp,ge) : ",2f20.10)') gp, ge
      !-write(*,'(1x,"Numerical (gzp,gze): ",2f20.10)') gzp, gze
      !=== END DEBUG ===> Check numerical gradients --------------------------------------------------


      !-- add forces from the vibronic surface
      do i=1,2
         f(i) = f(i) - g(i)/(f0*tau0*taul)
      enddo

      !-- propagate velocities for half-step
      
      do i=1,2
         vhalf(i) = v(i) + half*dt*f(i) - half*dt*gamma*v(i) &
         &               + half*sqdt*sigma*ksi(i) &
         &               - eighth*dt2*gamma*(f(i) - gamma*v(i)) &
         &               - fourth*dt32*gamma*sigma*(half*ksi(i) + eta(i)/sq3)

         !-- propagate positions for one step
         x(i) = x(i) + dt*vhalf(i) + half*dt32*sigma*eta(i)/sq3
      enddo

      !-- calculate regular forces at a new position
      do i=1,2
         f(i) = -x(i)/(tau0*taul)
      enddo

      !-- calculate vibronic states and vibronic free energies at new positions
      call calculate_vibronic_states(kg,x(1),x(2))

      !-- current free energy (PMF)
      efes = fe(istate)

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      !-- add regular forces from the vibronic surface
      do i=1,2
         f(i) = f(i) - g(i)/(f0*tau0*taul)
      enddo

      !-- propagate velocities for a remaining half step
      do i=1,2
         v(i) = vhalf(i) + half*dt*f(i) - half*dt*gamma*vhalf(i) &
         &               + half*sqdt*sigma*ksi(i) &
         &               - eighth*dt2*gamma*(f(i) - gamma*vhalf(i)) &
         &               - fourth*dt32*gamma*sigma*(half*ksi(i) + eta(i)/sq3)
      enddo

      z1 = x(1)
      z2 = x(2)
      vz1 = v(1)
      vz2 = v(2)

      !-- calculate kinetic energies for both degrees of freedom
      ekin1 = half*f0*tau0*taul*v(1)*v(1)
      ekin2 = half*f0*tau0*taul*v(2)*v(2)

   end subroutine langevin_onodera_2d
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !-- Vanden-Eijnden and Ciccotti 2-nd order propagator
   !   for Onodera2 model (with two auxiliary variables)
   !   [Chem. Phys. Lett., 2006, 429, 310-316, Eq.23]
   !--------------------------------------------------------------------
   subroutine langevin_onodera2_2d(istate,kg,z1,z2,y1,y2,vz1,vz2,dt,temp,ekin1,ekin2,ekinhalf1,ekinhalf2,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, y1, y2, vz1, vz2
      real(8), intent(out)   :: ekin1, ekin2, ekinhalf1, ekinhalf2, efes

      real(8), parameter :: half=0.5d0
      real(8), parameter :: eighth=1.d0/8.d0
      real(8), parameter :: fourth=0.25d0
      real(8), parameter :: e32=1.5d0

      integer :: i, ndabf
      real(8) :: sq3
      real(8) :: d, dy, psi_factor, f1, f2, fr, ynew
      real(8), dimension(2) :: x, y, v, vy, f, fy, g, mass, gammaz, sigmaz
      real(8), dimension(2) :: ksiz, etaz, psiy, vhalf
      real(8) :: dt2, sqdt, dt32                       !, zp, ze, gp, ge

      !-- DEBUG variables
      !real(8) :: zinc, tr1b, tr2a, fplus, fminus, gzp, gze, gselfzp, gselfze

      sq3 = sqrt(3.d0)

      dt2 = dt*dt
      sqdt = sqrt(dt)
      dt32 = dt*sqdt

      x(1) = z1
      x(2) = z2
      v(1) = vz1
      v(2) = vz2
      y(1) = y1
      y(2) = y2
      mass(1) = effmass1
      mass(2) = effmass2

      !-- define Ciccotti constants
      gammaz(1) = etax/effmass1
      gammaz(2) = etax/effmass2
      sigmaz(1) = sqrt(2*kb*temp*gammaz(1)/effmass1)
      sigmaz(2) = sqrt(2*kb*temp*gammaz(2)/effmass2)

      !-- define Honeycutt constants
      if (abs(etay).gt.1.d-20) then
         d = kb*temp/etay
      else
         d = 0.d0
      endif
      psi_factor = sqrt(2.d0*d*dt)

      !-- generate two independent random variables ksiz and etaz
      !   from univariate Gaussian distributions
      do i=1,2
         ksiz(i) = gaussdist_boxmuller()
         etaz(i) = gaussdist_boxmuller()
      enddo

      !-- generate independent random variables psiy
      !   from univariate Gaussian distributions (for auxiliary variables)
      do i=1,2
         psiy(i) = gaussdist_boxmuller()
      enddo


      !-- regular forces from the self-energy
      do i=1,2
         f(i) = -f0*x(i)/mass(i)
      enddo


      !-- regular forces from the effective potential
      do i=1,2
         f(i) = f(i) - gamma*(x(i) + y(i))/mass(i)
      enddo

      !-- calculate vibronic states and vibronic free energies
      !   REUSE THE VIBRONIC STATES FROM THE END OF THE PREVIOUS TIMESTEP
      !---(AVS)-call calculate_vibronic_states(kg,x(1),x(2))-------------

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      !=== DEBUG ===> Check numerical gradients ----------------------------------------------------
      !-zinc = 1.d-4
      !-tr1b = tr(1,2,32,1)
      !-tr2a = tr(1,3,32,1)
      !-call feszz3(mode,iset,kg,zp+zinc/2.d0,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fplus = fe(istate)
      !-call feszz3(mode,iset,kg,zp-zinc/2.d0,ze,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fminus = fe(istate)
      !-gzp = (fplus-fminus)/zinc
      !-gselfzp = 0.5d0*(erx*(tr2a+ze) - eret*(tr1b+zp))/(erx*erx - erpt*eret)
      !-gzp = gzp - gselfzp
      !-call feszz3(mode,iset,kg,zp,ze+zinc/2.d0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fplus = fe(istate)
      !-call feszz3(mode,iset,kg,zp,ze-zinc/2.d0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
      !-fminus = fe(istate)
      !-gze = (fplus-fminus)/zinc
      !-gselfze = 0.5d0*(erx*(tr1b+zp) - erpt*(tr2a+ze))/(erx*erx - erpt*eret)
      !-gze = gze - gselfze
      !-write(*,*) "======== DEBUG =====> numerical vs. analytical gradient"
      !-write(*,'(1x,"Analytical (gp,ge) : ",2f20.10)') gp, ge
      !-write(*,'(1x,"Numerical (gzp,gze): ",2f20.10)') gzp, gze
      !=== END DEBUG ===> Check numerical gradients --------------------------------------------------


      !-- add forces from the vibronic surface
      do i=1,2
         f(i) = f(i) - g(i)/mass(i)
      enddo

      !-- regular forces for auxiliary variables
      do i=1,2
         fy(i) = -gamma*(x(i) + y(i))
      enddo

      !-- propagate velocities for half-step
      do i=1,2
         vhalf(i) = v(i) + half*dt*f(i) - half*dt*gammaz(i)*v(i) &
         &               + half*sqdt*sigmaz(i)*ksiz(i) &
         &               - eighth*dt2*gammaz(i)*(f(i) - gammaz(i)*v(i)) &
         &               - fourth*dt32*gammaz(i)*sigmaz(i)*(half*ksiz(i) + etaz(i)/sq3)
      enddo

      !-- calculate kinetic energies for both degrees of freedom
      !   at half time step
      ekinhalf1 = half*effmass1*vhalf(1)*vhalf(1)
      ekinhalf2 = half*effmass2*vhalf(2)*vhalf(2)

      !-- Honeycutt RK2 propagation for auxiliary variables
      if (abs(etay).gt.1.d-20) then
         do i=1,2
            f1 = fy(i)/etay
            fr = psi_factor*psiy(i)
            ynew = y(i) + f1*dt + fr
            f2 = -gamma*(x(i) + ynew)/etay
            dy = 0.5d0*dt*(f1 + f2) + fr
            y(i) = y(i) + dy
         enddo
      endif

      !-- propagate positions for one step
      do i=1,2
         x(i) = x(i) + dt*vhalf(i) + half*dt32*sigmaz(i)*etaz(i)/sq3
      enddo

      !-- calculate regular forces at a new position
      do i=1,2
         f(i) = -f0*x(i)/mass(i) - gamma*(x(i) + y(i))/mass(i)
      enddo

      !-- calculate vibronic states and vibronic free energies at new positions
      call calculate_vibronic_states(kg,x(1),x(2))

      !-- current free energy (PMF)
      efes = fe(istate)

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      !-- add regular forces from the vibronic surface
      do i=1,2
         f(i) = f(i) - g(i)/mass(i)
      enddo

      !-- propagate velocities for a remaining half step
      do i=1,2
         v(i) = vhalf(i) + half*dt*f(i) - half*dt*gammaz(i)*vhalf(i) &
         &               + half*sqdt*sigmaz(i)*ksiz(i) &
         &               - eighth*dt2*gammaz(i)*(f(i) - gammaz(i)*vhalf(i)) &
         &               - fourth*dt32*gammaz(i)*sigmaz(i)*(half*ksiz(i) + etaz(i)/sq3)
      enddo

      z1 = x(1)
      z2 = x(2)
      vz1 = v(1)
      vz2 = v(2)
      y1 = y(1)
      y2 = y(2)

      !-- calculate kinetic energies for both degrees of freedom
      ekin1 = half*effmass1*v(1)*v(1)
      ekin2 = half*effmass2*v(2)*v(2)

   end subroutine langevin_onodera2_2d
   !--------------------------------------------------------------------
































   !--------------------------------------------------------------------
   !-- Runge-Cutta 2-nd order for coupled overdamped Langevin equation
   !   with exponential memory kernel
   !   [Basilevsky, Chudinov, Mol. Phys. 1990?]
   !   [Honeycutt, Phys. Rev. A, 1992, 45, 600]
   !--------------------------------------------------------------------
   subroutine langevin_debye2_2d(istate,kg,z1,z2,y1,y2,vz1,vz2,dt,temp,ekin1,ekin2,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, y1, y2, vz1, vz2
      real(8), intent(out)   :: ekin1, ekin2, efes

      real(8), parameter :: half=0.5d0

      integer :: i
      real(8) :: ddx, ddy, dx, dy, psi_factor_x, psi_factor_y
      real(8), dimension(2) :: x, y, v, xtmp, ytmp, g, f1, f2, f1y, f2y
      real(8), dimension(2) :: psix, psiy

      x(1) = z1
      x(2) = z2
      v(1) = vz1
      v(2) = vz2
      y(1) = y1
      y(2) = y2

      !-- define Honeycutt constants for x and y coordinates

      if (abs(etax).gt.1.d-20) then
         ddx = kb*temp/etax
      else
         ddx = 0.d0
      endif
      psi_factor_x = sqrt(2.d0*ddx*dt)

      if (abs(etay).gt.1.d-20) then
         ddy = kb*temp/etay
      else
         ddy = 0.d0
      endif
      psi_factor_y = sqrt(2.d0*ddy*dt)

      !-- generate independent random variables psiy and psix
      !   from univariate Gaussian distributions (for z and auxiliary variables)
      do i=1,2
         psix(i) = gaussdist_boxmuller()
         psiy(i) = gaussdist_boxmuller()
      enddo

      !-- regular forces from the self-energy
      do i=1,2
         f1(i) = -f0*x(i)
      enddo

      !-- regular forces from the effective potential
      do i=1,2
         f1(i) = f1(i) - gamma*(x(i) + y(i))
      enddo

      !-- calculate vibronic states and vibronic free energies
      !   REUSE THE VIBRONIC STATES FROM THE END OF THE PREVIOUS TIMESTEP
      !---(AVS)-call calculate_vibronic_states(kg,x(1),x(2))-------------

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      !-- add forces from the vibronic surface
      do i=1,2
         f1(i) = f1(i) - g(i)
      enddo

      !-- divide forces by friction coefficient
      do i=1,2
         f1(i) = f1(i)/etax
      enddo

      !-- regular forces for auxiliary variables
      do i=1,2
         f1y(i) = -gamma*(x(i) + y(i))
         f1y(i) = f1y(i)/etay
      enddo

      !-- propagate coordinates to obtain an intermediate point for force evaluation

      do i=1,2
         xtmp(i) = x(i) + dt*f1(i)  + psi_factor_x*psix(i)
         ytmp(i) = y(i) + dt*f1y(i) + psi_factor_y*psiy(i)
      enddo

      !-- calculate forces at the intermidiate point (F_2 in Honeycutt notations)
      do i=1,2
         f2(i)  = - f0*xtmp(i) - gamma*(xtmp(i) + ytmp(i))
         f2y(i) = - gamma*(xtmp(i) + ytmp(i))
         f2y(i) = f2y(i)/etay
      enddo

      call calculate_vibronic_states(kg,xtmp(1),xtmp(2))
      call calculate_gradient(istate,g(1),g(2))

      !-- add forces from the vibronic surface
      do i=1,2
         f2(i) = f2(i) - g(i)
         f2(i) = f2(i)/etax
      enddo

      !-- Honeycutt RK2 propagation for x and y variables

      do i=1,2
         dx = half*dt*(f1(i)+f2(i))   + psi_factor_x*psix(i)
         dy = half*dt*(f1y(i)+f2y(i)) + psi_factor_y*psiy(i)
         x(i) = x(i) + dx
         y(i) = y(i) + dy
         v(i) = dx/dt
      enddo

      !-- calculate vibronic states and vibronic free energies at new positions
      call calculate_vibronic_states(kg,x(1),x(2))

      !-- current free energy (PMF)
      efes = fe(istate)

      z1 = x(1)
      z2 = x(2)
      vz1 = v(1)
      vz2 = v(2)
      y1 = y(1)
      y2 = y(2)

      !-- kinetic energies are zero for overdamped mption
      ekin1 = 0.d0
      ekin2 = 0.d0

   end subroutine langevin_debye2_2d
   !--------------------------------------------------------------------













































   !--------------------------------------------------------------------
   !-- Right-hand side (time derivative) of von-Neumann equation
   !   for density matrix
   !--------------------------------------------------------------------
   subroutine tdvn_derivatives(t_,ro,dro)

      real(8), intent(in) :: t_
      complex(8), intent(in),  dimension(nstates,nstates) :: ro
      complex(8), intent(out), dimension(nstates,nstates) :: dro

      integer :: i, j, k
      real(8) :: ei, ej, vdik, vdkj
      complex(8) :: droij

      do i=1,nstates
         do j=i,nstates

            droij = (0.d0, 0.d0)
            
            if (i.eq.j) then

               do k=1,nstates
                  !-- interpolate the nonadiabatic coupling term
                  vdik = a0_vdotd(i,k) + a1_vdotd(i,k)*t_ + a2_vdotd(i,k)*t_*t_
                  droij = droij - 2.d0*real(ro(k,i))*vdik
               enddo

            else

               !-- interpolate adiabatic energies
               ei = a0_fe(i) + a1_fe(i)*t_ + a2_fe(i)*t_*t_
               ej = a0_fe(j) + a1_fe(j)*t_ + a2_fe(j)*t_*t_
               droij = droij - ro(i,j)*(ei - ej)/(ii*hbarps)

               do k=1,nstates
                  !-- interpolate the nonadiabatic coupling terms
                  vdik = a0_vdotd(i,k) + a1_vdotd(i,k)*t_ + a2_vdotd(i,k)*t_*t_
                  vdkj = a0_vdotd(k,j) + a1_vdotd(k,j)*t_ + a2_vdotd(k,j)*t_*t_
                  droij = droij - (ro(k,j)*vdik - ro(i,k)*vdkj)
               enddo

               dro(i,j) = droij
               dro(j,i) = conjg(droij)

            endif

         enddo
      enddo

   end subroutine tdvn_derivatives

   !--------------------------------------------------------------------
   !-- Quantum propagator for the density matrix (von Neumann equation)
   !   I. Based on the classic 4-th order Runge-Kutta algorithm.
   !--------------------------------------------------------------------
   subroutine propagate_density_rk4(t_,qtstep_)

      real(8), intent(in) :: t_, qtstep_

      complex(8), dimension(nstates,nstates) :: y, dy1, dy2, dy3, dy4

      !-- initial amplitudes
      y = density_matrix

      !-- Runge-Kutta steps
      call tdvn_derivatives(t_,y,dy1)
      call tdvn_derivatives(t_+0.5d0*qtstep_,y+0.5d0*qtstep_*dy1,dy2)
      call tdvn_derivatives(t_+0.5d0*qtstep_,y+0.5d0*qtstep_*dy2,dy3)
      call tdvn_derivatives(t_+qtstep_,y+qtstep_*dy3,dy4)

      !-- final amplitude (array operation)
      density_matrix = y + qtstep_*(dy1 + 2.d0*dy2 + 2.d0*dy3 + dy4)/6.d0

   end subroutine propagate_density_rk4

   !--------------------------------------------------------------------
   !-- Right-hand side (time derivative) of TDSE
   !--------------------------------------------------------------------
   subroutine tdse_derivatives(t_,c,dc)

      real(8), intent(in) :: t_
      complex(8), intent(in),  dimension(nstates) :: c
      complex(8), intent(out), dimension(nstates) :: dc

      integer :: i, j
      real(8) :: e0, de, vdij
      complex(8) :: dci

      !-- interpolated value of the ground state energy
      e0 = a0_fe(1) + a1_fe(1)*t_ + a2_fe(1)*t_*t_

      do i=1,nstates

         !-- calculate the interpolated value of energy
         !   relative to the ground state energy
         de = a0_fe(i) + a1_fe(i)*t_ + a2_fe(i)*t_*t_ - e0

         dci = -ii*c(i)*de/hbarps
         
         do j=1,nstates
            !-- calculate the interpolated value of the
            !   nonadiabatic coupling term
            vdij = a0_vdotd(i,j) + a1_vdotd(i,j)*t_ + a2_vdotd(i,j)*t_*t_
            dci = dci - c(j)*vdij
         enddo

         dc(i) = dci

      enddo

   end subroutine tdse_derivatives





   !--------------------------------------------------------------------
   !-- Right-hand side (time derivative) of TDSE with phase correction
   !--------------------------------------------------------------------
   subroutine tdse_derivatives_phcorr(istate_,t_,c,dc)

      integer, intent(in) :: istate_
      real(8), intent(in) :: t_
      complex(8), intent(in),  dimension(nstates) :: c
      complex(8), intent(out), dimension(nstates) :: dc

      integer :: i, j
      real(8) :: ekinocc, efeocc, efei, hocc, hii, conservation, vdij
      complex(8) :: dci

      !-- interpolated value of the kinetic energy in the occupied state
      ekinocc = a0_ekin + a1_ekin*t_ + a2_ekin*t_*t_

      hocc = -2.d0*ekinocc

      !-- interpolated value of the potential (free) energy of the occupied state
      efeocc = a0_fe(istate_) + a1_fe(istate_)*t_ + a2_fe(istate_)*t_*t_

      do i=1,nstates

         if (i.eq.istate_) then

            hii = hocc

         else

            !-- calculate the interpolated value of energy
            !   of the current state
            efei = a0_fe(i) + a1_fe(i)*t_ + a2_fe(i)*t_*t_
            
            conservation = ekinocc + efeocc - efei
            
            if (conservation.gt.0.d0) then
               hii = -2.d0*sqrt(ekinocc*conservation)
            else
               hii = 0.d0
            endif

         endif

         dci = -ii*c(i)*hii/hbarps

         do j=1,nstates
            !-- calculate the interpolated value of the
            !   nonadiabatic coupling term
            vdij = a0_vdotd(i,j) + a1_vdotd(i,j)*t_ + a2_vdotd(i,j)*t_*t_
            dci = dci - c(j)*vdij
         enddo

         dc(i) = dci

      enddo

   end subroutine tdse_derivatives_phcorr



   !--------------------------------------------------------------------
   !-- Quantum propagator for renormalized amplitudes
   !   I. Based on the classic 4-th order Runge-Kutta algorithm.
   !--------------------------------------------------------------------
   subroutine propagate_amplitudes_rk4(t_,qtstep_)

      real(8), intent(in) :: t_, qtstep_

      complex(8), dimension(nstates) :: y, dy1, dy2, dy3, dy4
      complex(8), dimension(nstates) :: y2, y3, y4

      !-- initial amplitudes
      y = amplitude

      !-- Runge-Kutta steps
      call tdse_derivatives(t_,y,dy1)
      y2 = y + 0.5d0*qtstep_*dy1
      call tdse_derivatives(t_+0.5d0*qtstep_,y2,dy2)
      y3 = y + 0.5d0*qtstep_*dy2
      call tdse_derivatives(t_+0.5d0*qtstep_,y3,dy3)
      y4 = y + qtstep_*dy3
      call tdse_derivatives(t_+qtstep_,y4,dy4)

      !-- final amplitude (array operation)
      amplitude = y + qtstep_*(dy1 + 2.d0*dy2 + 2.d0*dy3 + dy4)/6.d0

   end subroutine propagate_amplitudes_rk4


   !--------------------------------------------------------------------
   !-- Quantum propagator for renormalized amplitudes
   !   (for phase-corrected algorithm)
   !   I. Based on the classic 4-th order Runge-Kutta algorithm.
   !--------------------------------------------------------------------
   subroutine propagate_amplitudes_phcorr_rk4(istate_,t_,qtstep_)

      integer, intent(in) :: istate_
      real(8), intent(in) :: t_, qtstep_

      complex(8), dimension(nstates) :: y, dy1, dy2, dy3, dy4
      complex(8), dimension(nstates) :: y2, y3, y4

      !-- initial amplitudes
      y = amplitude

      !-- Runge-Kutta steps
      call tdse_derivatives_phcorr(istate_,t_,y,dy1)
      y2 = y + 0.5d0*qtstep_*dy1
      call tdse_derivatives_phcorr(istate_,t_+0.5d0*qtstep_,y2,dy2)
      y3 = y + 0.5d0*qtstep_*dy2
      call tdse_derivatives_phcorr(istate_,t_+0.5d0*qtstep_,y3,dy3)
      y4 = y + qtstep_*dy3
      call tdse_derivatives_phcorr(istate_,t_+qtstep_,y4,dy4)

      !-- final amplitude (array operation)
      amplitude = y + qtstep_*(dy1 + 2.d0*dy2 + 2.d0*dy3 + dy4)/6.d0

   end subroutine propagate_amplitudes_phcorr_rk4


   !--------------------------------------------------------------------
   !-- Check if surface hoppping should take place
   !--------------------------------------------------------------------
   function switch_state(istate) result(new_state)
      integer, intent(in) :: istate
      integer :: new_state
      integer :: i
      real(8) :: total_prob
      real(4) :: s
      new_state = istate
      s = ran2nr()
      total_prob = 0.d0
      do i=1,nstates
         if (i.ne.istate) total_prob = total_prob + switch_prob(i)
         if (total_prob.gt.s) then
            new_state = i
            exit
         endif
      enddo
   end function switch_state

   !--------------------------------------------------------------------
   !-- Check if surface hoppping should take place
   !--------------------------------------------------------------------
   subroutine adjust_velocities(istate,new_state,vz1,vz2,success)

      integer, intent(in)    :: istate, new_state
      real(8), intent(inout) :: vz1, vz2
      logical, intent(out)   :: success

      real(8) :: dkl_z1, dkl_z2
      real(8) :: akl, bkl, ekl, discr
      real(8) :: gamma

      dkl_z1 = coupz1(istate,new_state)
      dkl_z2 = coupz2(istate,new_state)

      akl = dkl_z1*dkl_z1/effmass1 + dkl_z2*dkl_z2/effmass2

      bkl  = v_dot_d(istate,new_state)
      ekl = fe(istate) - fe(new_state)

      discr = bkl*bkl + 2.d0*akl*ekl

      if (discr.lt.0.d0) then

         !-- not enough energy for a switch: do nothing
         success = .false.

         !=== DEBUG ===================================================================
         write(*,*)
         write(*,'(137("-"))')
         write(*,*) "Rejected switch:",istate," --->",new_state
         write(*,*) "Energy gap:", -ekl, " kcal/mol"
         write(*,'(137("-"))')
         !=== end DEBUG ===============================================================

      else

         !-- calculate adjustment coefficient
         gamma = bkl/akl
         if (bkl.lt.0) then
            gamma = gamma + sqrt(discr)/akl
         else
            gamma = gamma - sqrt(discr)/akl
         endif

         !-- adjust velocities
         vz1 = vz1 - gamma*dkl_z1/effmass1
         vz2 = vz2 - gamma*dkl_z2/effmass2

         success = .true.

         !=== DEBUG ===================================================================
         write(*,*)
         write(*,'(137("-"))')
         write(*,*) "Switch: ",istate," --->", new_state
         write(*,*) "Nonadiabatic coupling: d_z1 = ", dkl_z1
         write(*,*) "                       d_z2 = ", dkl_z2
         write(*,*) "                       |d|  = ", sqrt(dkl_z1*dkl_z1 + dkl_z2*dkl_z2)
         write(*,*) "Switch probability: ", switch_prob(new_state)
         write(*,*) "Kinetic energy gain (loss if negative): ", ekl, " kcal/mol"
         write(*,*) "Velocity adjustments: vz1 = vz1 + ", -gamma*dkl_z1/effmass1
         write(*,*) "                      vz2 = vz2 + ", -gamma*dkl_z2/effmass2
         write(*,'(137("-"))')
         !=== end DEBUG ===============================================================

      endif

   end subroutine adjust_velocities

!===============================================================================
end module propagators_3d
