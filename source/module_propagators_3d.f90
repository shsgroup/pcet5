module propagators_3d

   !---------------------------------------------------------------------
   ! Contains the routines for langevin propagators
   !---------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2010-10-26 21:06:21 $
   !  $Revision: 5.1 $
   !  $Log: not supported by cvs2svn $
   !
   !---------------------------------------------------------------------

   use cst
   use parsol
   use solmat
   use feszz_3d
   use random_generators
   use rk_parameters

   !---------------------------------------------------------------------
   implicit none
   private

   complex(kind=8), parameter :: ii=(zero,one)

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

   !-- Array for the density matrix
   !---------------------------------------------------------
   complex(8), allocatable, dimension(:,:) :: density_matrix

   !-- arrays for the transition probabilities for current state (b_jk)
   !---------------------------------------------------------------------
   real(8), allocatable, dimension(:) :: b_prob

   !-- Array for the switching probabilities from current state
   !-----------------------------------------------------------
   real(8), allocatable, dimension(:) :: switch_prob

   !-- arrays with interpolation coefficients
   !---------------------------------------------------------
   real(8), allocatable, dimension(:) :: a0_fe, a1_fe, a2_fe
   real(8), allocatable, dimension(:,:) :: a0_vdotd, a1_vdotd, a2_vdotd

   public :: set_mode
   public :: allocate_vibronic_states, deallocate_vibronic_states
   public :: allocate_evb_weights, deallocate_evb_weights
   public :: allocate_mdqt_arrays, deallocate_mdqt_arrays
   public :: calculate_vibronic_states
   public :: calculate_vibronic_couplings
   public :: calculate_v_dot_d
   public :: interpolate_vdotd
   public :: interpolate_energy
   public :: set_initial_amplitudes
   public :: set_initial_density
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
   public :: langevin_debye_2d
   public :: langevin_onodera_2d
   public :: print_propagators_3d
   public :: tdse_derivatives
   public :: propagate_amplitudes_rk4
   public :: propagate_density_rk4
   public :: switch_state
   public :: adjust_velocities
   public :: reset_switch_prob
   public :: accumulate_switch_prob
   public :: normalize_switch_prob

contains

   !===DEBUG printout===
   subroutine print_propagators_3d
      write(*,*)
      write(*,*) "-----Contents of the propagators_3d module"
      write(*,*) "mode:    ",mode
      write(*,*) "iset:    ",iset
      write(*,*) "nstates: ",nstates
      write(*,*) "nzdim:   ",nzdim
      write(*,*) "ielst:   ",ielst
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

   !---------------------------------------------------------------------
   ! array allocation routines
   !---------------------------------------------------------------------

   subroutine allocate_vibronic_states(nprst_,npnts_)
      integer, intent(in) :: nprst_, npnts_
      integer :: nz
      nz = ielst*nprst_
      allocate (fe(nstates))
      allocate (fe_prev(nstates))
      allocate (z(nz,nz))
      allocate (z_prev(nz,nz))
      allocate (psiel(ielst,npnts_,ielst))
      allocate (psipr(ielst,nprst_,npnts_))
      allocate (enel(ielst,npnts_))
      allocate (envib(ielst,nprst_))
   end subroutine allocate_vibronic_states

   subroutine allocate_evb_weights
      allocate (wght(ielst,nstates))
   end subroutine allocate_evb_weights

   subroutine allocate_mdqt_arrays
      allocate (coupz1(nstates,nstates))
      allocate (coupz2(nstates,nstates))
      allocate (coupz1_prev(nstates,nstates))
      allocate (coupz2_prev(nstates,nstates))
      allocate (amplitude(nstates))
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
      deallocate (fe)
      deallocate (fe_prev)
      deallocate (z)
      deallocate (z_prev)
      deallocate (psiel)
      deallocate (psipr)
      deallocate (enel)
      deallocate (envib)
   end subroutine deallocate_vibronic_states

   subroutine deallocate_evb_weights
      deallocate (wght)
   end subroutine deallocate_evb_weights

   subroutine deallocate_mdqt_arrays
      deallocate (coupz1)
      deallocate (coupz2)
      deallocate (coupz1_prev)
      deallocate (coupz2_prev)
      deallocate (amplitude)
      deallocate (density_matrix)
      deallocate (b_prob)
      deallocate (switch_prob)
      deallocate (v_dot_d)
      deallocate (v_dot_d_mid)
      deallocate (v_dot_d_prev)
      deallocate (a0_fe,a1_fe,a2_fe)
      deallocate (a0_vdotd)
      deallocate (a1_vdotd)
      deallocate (a2_vdotd)
   end subroutine deallocate_mdqt_arrays

   !--------------------------------------------------------------------
   !-- Set initial amplitudes to correspond to a pure adiabatic state
   !--------------------------------------------------------------------
   subroutine set_initial_amplitudes(istate)
      integer, intent(in) :: istate
      amplitude = 0.d0
      amplitude(istate) = (1.d0,0.d0)
   end subroutine set_initial_amplitudes

   !---------------------------------------------------------------------
   !-- Set initial density matrix to correspond to a pure adiabatic state
   !---------------------------------------------------------------------
   subroutine set_initial_density(istate)
      integer, intent(in) :: istate
      density_matrix = 0.d0
      density_matrix(istate,istate) = (1.d0,0.d0)
   end subroutine set_initial_density

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
   subroutine calculate_bprob_amp(istate,t)
      integer, intent(in) :: istate
      real(8), intent(in) :: t
      real(8) :: vdji
      integer :: j
      do j=1,nstates
         vdji = a0_vdotd(j,istate) + a1_vdotd(j,istate)*t + a2_vdotd(j,istate)*t*t
         b_prob(j) = -2.d0*Real(amplitude(j)*conjg(amplitude(istate))*vdji)
      enddo
   end subroutine calculate_bprob_amp

   !--------------------------------------------------------------------
   !-- calculate transition probabilities b_jk
   !   from density matrix
   !--------------------------------------------------------------------
   subroutine calculate_bprob_den(istate,t)
      integer, intent(in) :: istate
      real(8), intent(in) :: t
      real(8) :: vdji
      integer :: j
      do j=1,nstates
         vdji = a0_vdotd(j,istate) + a1_vdotd(j,istate)*t + a2_vdotd(j,istate)*t*t
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
   subroutine accumulate_switch_prob(qtstep)
      real(8), intent(in) :: qtstep
      switch_prob = switch_prob + b_prob*qtstep   ! (array operation)
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
   subroutine interpolate_energy(t_prev,t)
      real(8), intent(in) :: t_prev, t
      real(8) :: dt, f1, f2, t1, t2
      integer :: i
      t1 = t_prev
      t2 = t
      dt = t - t_prev
      do i=1,nstates
         f2 = fe(i)
         f1 = fe_prev(i)
         a0_fe(i) = (t2*f1 - t1*f2)/dt
         a1_fe(i) = (f2 - f1)/dt
         a2_fe(i) = 0.d0
      enddo
   end subroutine interpolate_energy

   !--------------------------------------------------------------------
   !-- calculate interpolation coefficients
   !   for the nonadiabatic coupling terms v*d_{kl}
   !--------------------------------------------------------------------
   subroutine interpolate_vdotd(interpolation,t_prev,t)

      character(len=*), intent(in) :: interpolation
      real(8), intent(in) :: t_prev, t
      integer :: i, j
      real(8) :: x0, x1, x2, x01, x02, x12
      real(8) :: y0, y1, y2, y01, y02, y12
      real(8) :: xdenom, a0, a1, a2
      
      x0 = t_prev
      x2 = t
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
   subroutine langevin_debye_2d(istate,kg,z1,z2,vz1,vz2,dt,temp,ekin,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, vz1, vz2
      real(8), intent(out)   :: ekin, efes

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
      call calculate_vibronic_states(kg,x(1),x(2))

      !-- calculate gradients
      call calculate_gradient(istate,g(1),g(2))

      do i=1,2
         call gaussdist_boxmuller(psi,iseed)
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

      !-- calculate kinetic energy
      do i=1,2
         ekin = ekin + v(i)*v(i)
      enddo
      ekin = half*f0*tau0*taul*ekin

   end subroutine langevin_debye_2d

   !--------------------------------------------------------------------
   !-- Vanden-Eijnden and Ciccotti 2-nd order propagator
   !   for Onodera model
   !   [Chem. Phys. Lett., 2006, 429, 310-316, Eq.23]
   !--------------------------------------------------------------------
   subroutine langevin_onodera_2d(istate,kg,z1,z2,vz1,vz2,dt,temp,ekin,efes)

      implicit none
      integer, intent(in)    :: istate, kg
      real(8), intent(in)    :: dt, temp
      real(8), intent(inout) :: z1, z2, vz1, vz2
      real(8), intent(out)   :: ekin, efes

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
         call gaussdist_boxmuller(ksi(i),iseed)
         call gaussdist_boxmuller(eta(i),iseed)
      enddo

      !-- regular forces from the self-energy
      do i=1,2
         f(i) = -x(i)/(tau0*taul)
      enddo

      !-- calculate vibronic states and vibronic free energies
      call calculate_vibronic_states(kg,x(1),x(2))

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

      !-- calculate kinetic energy
      do i=1,2
         ekin = ekin + v(i)*v(i)
      enddo
      ekin = half*f0*tau0*taul*ekin

   end subroutine langevin_onodera_2d
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !-- Right-hand side (time derivative) of von-Neumann equation
   !   for density matrix
   !--------------------------------------------------------------------
   subroutine tdvn_derivatives(t,ro,dro)

      real(8), intent(in) :: t
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
                  vdik = a0_vdotd(i,k) + a1_vdotd(i,k)*t + a2_vdotd(i,k)*t*t
                  droij = droij - 2.d0*real(ro(k,i))*vdik
               enddo

            else

               !-- interpolate adiabatic energies
               ei = a0_fe(i) + a1_fe(i)*t + a2_fe(i)*t*t
               ej = a0_fe(j) + a1_fe(j)*t + a2_fe(j)*t*t
               droij = droij - ro(i,j)*(ei - ej)/(ii*hbarps)

               do k=1,nstates
                  !-- interpolate the nonadiabatic coupling terms
                  vdik = a0_vdotd(i,k) + a1_vdotd(i,k)*t + a2_vdotd(i,k)*t*t
                  vdkj = a0_vdotd(k,j) + a1_vdotd(k,j)*t + a2_vdotd(k,j)*t*t
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
   subroutine propagate_density_rk4(t,qtstep)

      real(8), intent(in) :: t, qtstep

      complex(8), dimension(nstates,nstates) :: y, dy1, dy2, dy3, dy4

      !-- initial amplitudes
      y = density_matrix

      !-- Runge-Kutta steps
      call tdvn_derivatives(t,y,dy1)
      call tdvn_derivatives(t+0.5d0*qtstep,y+0.5d0*qtstep*dy1,dy2)
      call tdvn_derivatives(t+0.5d0*qtstep,y+0.5d0*qtstep*dy2,dy3)
      call tdvn_derivatives(t+qtstep,y+qtstep*dy3,dy4)

      !-- final amplitude (array operation)
      density_matrix = y + qtstep*(dy1 + 2.d0*dy2 + 2.d0*dy3 + dy4)/6.d0

   end subroutine propagate_density_rk4

   !--------------------------------------------------------------------
   !-- Right-hand side (time derivative) of TDSE
   !--------------------------------------------------------------------
   subroutine tdse_derivatives(t,c,dc)

      real(8), intent(in) :: t
      complex(8), intent(in),  dimension(nstates) :: c
      complex(8), intent(out), dimension(nstates) :: dc

      integer :: i, j
      real(8) :: e0, de, vdij
      complex(8) :: dci

      !-- interpolated value of the ground state energy
      e0 = a0_fe(1) + a1_fe(1)*t + a2_fe(1)*t*t

      do i=1,nstates

         !-- calculate the interpolated value of energy
         !   relative to the ground state energy
         de = a0_fe(i) + a1_fe(i)*t + a2_fe(i)*t*t - e0

         dci = c(i)*de/(ii*hbarps)
         
         do j=1,nstates
            !-- calculate the interpolated value of the
            !   nonadiabatic coupling term
            vdij = a0_vdotd(i,j) + a1_vdotd(i,j)*t + a2_vdotd(i,j)*t*t
            dci = dci - c(j)*vdij
         enddo

         dc(i) = dci

      enddo

   end subroutine tdse_derivatives

   !--------------------------------------------------------------------
   !-- Quantum propagator for renormalized amplitudes
   !   I. Based on the classic 4-th order Runge-Kutta algorithm.
   !--------------------------------------------------------------------
   subroutine propagate_amplitudes_rk4(t,qtstep)

      real(8), intent(in) :: t, qtstep

      complex(8), dimension(nstates) :: y, dy1, dy2, dy3, dy4

      !-- initial amplitudes
      y = amplitude

      !-- Runge-Kutta steps
      call tdse_derivatives(t,y,dy1)
      call tdse_derivatives(t+0.5d0*qtstep,y+0.5d0*qtstep*dy1,dy2)
      call tdse_derivatives(t+0.5d0*qtstep,y+0.5d0*qtstep*dy2,dy3)
      call tdse_derivatives(t+qtstep,y+qtstep*dy3,dy4)

      !-- final amplitude (array operation)
      amplitude = y + qtstep*(dy1 + 2.d0*dy2 + 2.d0*dy3 + dy4)/6.d0

   end subroutine propagate_amplitudes_rk4

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
      call ran2nr(s,iseed)
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

      real(8) :: dkl_z1, dkl_z2, dkl_sq
      real(8) :: vdkl, ekl, discr
      real(8) :: gamma

      dkl_z1 = coupz1(istate,new_state)
      dkl_z2 = coupz2(istate,new_state)
      dkl_sq = dkl_z1*dkl_z1 + dkl_z2*dkl_z2

      vdkl  = v_dot_d(istate,new_state)
      ekl = fe(istate) - fe(new_state)

      discr = vdkl*vdkl - 2.d0*dkl_sq*ekl/effmass

      if (discr.lt.0.d0) then

         !-- not enough energy for a switch: do nothing
         success = .false.

         !=== DEBUG ===================================================================
         write(*,*) "-----------------------------------------------------------------"
         write(*,*) "Rejected switch:",istate," --->",new_state
         write(*,*) "Energy gap (positive is switch down):",fe(istate)-fe(new_state)
         write(*,*) "-----------------------------------------------------------------"
         !=== end DEBUG ===============================================================

      else

         !-- calculate adjustment coefficient
         gamma = vdkl/dkl_sq
         if (vdkl.lt.0) then
            gamma = gamma + sqrt(discr)/dkl_sq
         else
            gamma = gamma - sqrt(discr)/dkl_sq
         endif

         !-- adjust velocities
         vz1 = vz1 - gamma*dkl_z1
         vz2 = vz2 - gamma*dkl_z2

         success = .true.

         !=== DEBUG ===================================================================
         write(*,*) "-----------------------------------------------------------------"
         write(*,*) "Switch:",istate," --->",new_state
         write(*,*) "Switch probability:",switch_prob(new_state)
         write(*,*) "Energy gap (positive is switch down):",fe(istate)-fe(new_state)
         write(*,*) "Velocity adjustments: vz1 = vz1 - ",gamma*dkl_z1
         write(*,*) "                      vz2 = vz2 - ",gamma*dkl_z2
         write(*,*) "-----------------------------------------------------------------"
         !=== end DEBUG ===============================================================

      endif

   end subroutine adjust_velocities

!===============================================================================
end module propagators_3d
