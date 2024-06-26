subroutine dynamicset2
!===================================================================C
!
!  Driver for solvent dynamics on two one-dimensional
!  electronic free energy surfaces (Marcus ET model).
!
!  OPTIONS:
!
!  MDQT - employ Tully's surface hopping algorithm to incorporate
!         non-adiabatic transitions between electronic states
!         (otherwise classical dynamics on a single electronic
!         free energy surface is simulated)
!
!  ADIAB - move on adiabatic electronic free energy surfaces
!
!  DIABATIC - move on the ET diabatic electronic free energy surfaces
!
!  ISTATE=<N> - index of the occupied adiabatic electronic state
!               (1 for ground state) at t=0.
!
!  PURESTATEDENSITY - the initial density matrix corresponds to a pure state
!
!  MIXEDSTATEDENSITY - the initial density matrix corresponds to a coherent mixture of states
!
!  NOWEIGHTS - do not calculate the evb weights along the trajectory (default is YES)
!
!  SCALED - initial value of the solvent coordinate is given in the scaled frame (z).
!
!  Z0=<float> - center of the initial distribution along Z=X (or z=x) coordinate
!
!  FIXEDINIT - use center of the initial distribution as an initial coordinate value
!              (no sampling from Gaussian distribution)
!
!  VZE0=<float> - initial velocity along Z (energy gap) coordinate
!
!  VZ10=<float> - initial velocity along z (scaled energy gap) coordinate
!
!  TRAJOUT=<string> - name of the output files with trajectory data
!                     (default filename "trajectory_n.dat" where n is the
!                     trajectory index starting from 1)
!
!  SEED=<int> - random seed for RAN2NR random number generator.
!               (if SEED=0 then clock will be used to generate seed)
!               Default value is generated using the current time.
!
!  DUNISEEDS=i/j/k/l - random seeds for DUNI random number generator.
!               (if DUNISEEDS=CLOCK then clock will be used to generate seeds)
!
!  RESET_RANDOM - reset the ran2nr() sequence when starting each trajectory
!
!  NTRAJ=<int> - number of trajectories to generate (default is 1)
!
!  EPSMODEL=<KEY> - dielectric relaxation model (no defaults!).
!
!            DEBYE - Simple Debye relaxation model
!                    (overdamped Langevin dynamics).
!                    Not compatible with MDQT option!
!
!            ONODERA - Onodera model with short time correction
!                      (ordinary Langevin dynamics).
!
!            ONODERA2 - Onodera model with short time correction and two relaxation timescales
!                       (generalized Langevin dynamics).
!
!  EPSPARS - specify parameters of the dielectric function for ONODERA2 model:
!            TAU0/EFFMASS, EPS1, TAU1, TAU2
!
!  GLEPARS - specify parameters of the GLE for ONODERA2 model:
!            EFFMASS, GAMMA, TAUALPHA, ETA
!
!            Note that the EPS0 and EPS8 parameters of the corresponding dielectric
!            function can be specified within the SOLV keyword.
!
!  TSTEP=<float> - timestep for solvent dynamics in picoseconds (default=0.0005)
!
!  NSTEPS=<int> - number of steps (length of the trajectory)
!
!  NQSTEPS=<int> - number of TDSE steps per classical step in MDQT
!
!  MAXNQSTEPS=<int> - maximum number of TDSE steps per classical step in MDQT
!
!  INTERPOLATION=<KEY> - interpolation scheme for TDSE
!
!            LINEAR - linear interpolation of the nonadiabatic coupling term
!                     between t and t+dt (default)
!
!            QUADRATIC - quadratic interpolation of the nonadiabatic coupling term
!                        using values at t, t+dt/2, t+dt
!
!  NDUMP=<int> - trajectory output frequency (every NDUMP steps)
!
!  NDUMP6=<int> - trajectory screen output frequency (every NDUMP6 steps)
!
!  NDUMP777=<int> - populations and coherences output frequency
!                   (every NDUMP777 steps; no output if NDUMP777=0, default)
!
!  NDUMP888=<int> - couplings output frequency
!                   (every NDUMP888 steps; no output if NDUMP888=0, default)
!
!  NDUMP999=<int> - moments output frequency (AFSSH)
!                   (every NDUMP999 steps; no output if NDUMP999=0, default)
!
!  T=<float> - temperature in K
!
!  RESTART[=<file>] - restart trajectories (and random number sequences)
!                     File (default name is dynamics_checkpoint) must be
!                     present and not empty.
!
!  PHASE - Phase correction algorithm
!          [N. Shenvi, J. E. Subotnik, and W. Yang, JCP 135, 024101 (2011)]
!
!  AFSSH - Decoherence algorithm, Augmented Fewest Switches Surface Hopping (AFSSH)
!                [N. Shenvi, J. E. Subotnik, 2012, original algorithm]
!
!  COLLAPSE_REGION_COUPLING - A simple decoherence algorithm: wavefunction is collapsed
!                             to occupied state each time the trajectory leaves the interaction
!                             region defined by the cutoff value for the magnitude of the largest
!                             nonadiabatic coupling vector
!
!  COUPLING_CUTOFF=<float> - cutoff for the magnitude of nonadiabatic coupling vector
!
!  !!!COUPLE - couple TDSE to EOM for moments in decoherence algorithm (Eq.18),
!  !!!         Augmented Fewest Switches Surface Hopping (AFSSH) [N. Shenvi, J. E. Subotnik, 2012]
!  !!!         (this option is questionable because the density matrix evolving according
!  !!!          to this equation is not positive definite)
!  !!!         !!! ABANDONED !!!
!
!  ADJUST_ALONG_MOMENTS - adjust velocities along the difference of vectors of moments in decoherence algorithm,
!           Augmented Fewest Switches Surface Hopping (AFSSH-0) [N. Shenvi, J. E. Subotnik, 2011]
!
!  IDS - "Instantaneous Decoherence fror Succesful Hops" decoherence algorithm
!
!  IDA - "Instantaneous Decoherence fror Any Hops" decoherence algorithm
!
!  GEDC - Granucci's Energy-based Decoherence Correction (with C = 1, E0 = 0.1 a.u.)
!
!  REVVEL - Use Subotnik (Truhlar’s Nabla-V) approach for velocity reversal in case of frustrated hops
!
!  REVVEL0 - Velocity reversal for all frustrated hops
!
!--------------------------------------------------------------------------------
!  REACTIVE_FLUX - Starts trajectory at dividing surface.  Propagates
!  "backwards" in time until reactant or product is reached.  Follow stored trajectory
!  forward and update time dependent Schroedinger equation coefficients until
!  the dividing surface is reached.  Start the trajectory at the dividing
!  surface if reactants were reached in the reverse trajectory
!  and propagate normally until reactants or products are reached, at which time the
!  trajecotory is stopped.
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Id: dynamicset2.f90 176 2013-03-19 21:30:29Z souda $
!  $Revision: 176 $
!
!===================================================================C

   use pardim
   use keys
   use cst
   use timers
   use strings
   use control
   use control_dynamics
   use random_generators
   use data_et2
   use rates_et2
   use geogas, only: iptgas, xyzgas
   use parsol
   use propagators_et2
   use potential
   use frcm
   use elcm

   implicit none

   character(len=1024) :: options, solvoptions
   character(len=  80) :: fname
   character(len=  15) :: zgapdim="kcal/mol      "
   character(len=  15) :: zscadim="(kcal/mol)^1/2"
   character(len=   5) :: mode_dyn
   character(len=   6) :: traj_suffix
   character(len=   1) :: state_index
   character(len=  20) :: str
   character(len=1), dimension(4) :: iset_char_diab4=(/"1","2","3","4"/)

   logical :: adiab, diab4, weights
   logical :: switch=.false.
   logical :: success=.false.
   logical :: scaled=.false.
   logical :: initial_state_pure=.true.
   logical :: initial_state_diab=.false.
   logical :: purestate=.true.
   logical :: revvel_all=.false.
   logical :: revvel=.false.
   logical :: revvel_cond1=.false.
   logical :: revvel_cond2=.false.

   logical :: interaction_region_prev=.false.
   logical :: interaction_region=.false.
   logical :: double_well=.false.
   logical :: fixedinit=.false.

   logical :: reset_random=.false.

   integer :: nstates_dyn, nzdim_dyn, ielst_dyn, iseed_inp, iset_dyn
   integer :: initial_state=-1, iground=1, idecoherence=0
   integer :: istate, new_state, ndump6, ndump777, ndump888, ndump999
   integer :: number_of_skipped_trajectories=0
   integer :: number_of_failed_trajectories=0
   integer :: itraj_start=1
   integer :: ntraj_valid
   integer :: nqsteps_var

   integer :: ikey, ikeya, ikeyb, ioption, ioption2, ioption3, ioptionve, ioptionv1, kg0, ifrom, ito
   integer :: islash1, islash2, islash3, idum1, idum2, idum3, idum4
   integer :: ioutput, lenf, ispa, idash, icount, ireac, iprod
   integer :: itraj, istep, iqstep, k, itmp, itmp1
   integer :: number_of_switches=0, number_of_rejected=0, number_of_reversals=0
   integer :: itraj_channel=11
   integer :: ifes_channel=12

   real(kind=8) :: sigma, sigma1, sample, population_current
   real(kind=8) :: wf_norm, zmom1_norm, pzmom1_norm
   real(kind=8) :: zeit_start, zeit_end, zeit_total, traj_time_start, traj_time_end
   real(kind=8) :: zeit, zeit_prev
   real(kind=8) :: zeitq, zeitq_prev
   real(kind=8) :: z1, ze, vz10, vz1, vze0, vze, z10, ze0, y1, force1, force2, coup12
   real(kind=8) :: ekin, ekin1, ekin_prev, ekinhalf1, efes
   real(kind=8) :: vz1_prev, z1_prev, ze_prev
   real(kind=8) :: qtstep_var

   real(kind=8) :: hreac, hprod, vet, zi, zi_scaled
   real(kind=8) :: dg_reaction, dg_activation, et_marcus_rate
   real(kind=8) :: kappa_ad, prefactor1, prefactor2
   real(kind=8) :: et_rips_jortner_rate, adiabatic_rate
   real(kind=8) :: fr, fb, eb
   real(kind=8), dimension(4,4) :: h0k
   real(kind=8), dimension(4,4) :: tk, tinfk, trk, trinfk
   real(kind=8), dimension(2)   :: fe_diab, fe_adiab

   !-- Variables for the reactive_flux method
   !   Note: reactive_flux works for the two state system
   !   (the calculate_w_mu routine is coded to only work for two states)

   real(kind=8) :: N_tot, prob, rand          ! For calculating the initial adiabat
   integer :: i, final_step_reverse           ! keep track of time steps in reverse time loop

   real(kind=8), allocatable, dimension(:) :: z1_storage, vz1_storage        ! Storage of energy gap/vel coordinate
   real(kind=8), allocatable, dimension(:) :: vz1_storage_beforehop

   real(kind=8), allocatable, dimension(:,:,:) :: fnj_storage                ! Storage of the "fake" hopping probability
   integer, allocatable, dimension(:) :: switch_attempt                      ! has value 1 if switch attempted, 0 if not attempted
   integer, allocatable, dimension(:) :: occupied_adiabat                    ! holds the occupied adiabat at end of reverse time step istep
   integer, allocatable, dimension(:) :: occupied_adiabat_beforehop          ! holds occupied adiabat at beginning of reverse timestep
   integer, allocatable, dimension(:) :: switch_attempt_state                ! holds the state to which a switch was attempted

   real(kind=8) :: W, alpha, dividing_surface_egap, Fn, Fd, alpha_limit      ! quantities for reactive_flux

   logical :: successfulreverse     ! true if trajectory did not run out of time steps in the reverse trajectory before reaching reactants
   logical :: normalforward         ! true is should propagate from dividing surface foward

   integer :: Fnfile, Fdfile        ! files to hold values of Fn and Fd for all trajectories

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Start of executable statements
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   adiab   = .true.
   diab4   = .false.
   weights = .true.

   !~~~~~~~~~~~~~~
   ! Print banner
   !~~~~~~~~~~~~~~
   write(6,*)
   write(6,'("================================================================")')
   write(6,'("             ET SOLVENT DYNAMICS MODULE (optional MDQT)         ")')
   write(6,'("         (single solvent coordinate: Marcus 2-state model)      ")')
   write(6,'("================================================================"/)')

   !~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' DYNAMICSET2(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify options for DYNAMICSET2 keyword ***"/)')
      call clean_exit
   else
      call getopt(keywrd,ikey+13,options)
   endif

   !~~~~~~~~~~~~~
   ! Temperature
   !~~~~~~~~~~~~~

   ioption = index(options," T=")
   if (ioption.ne.0) then
      temp = reada(options,ioption+3)
   else
      temp = 300.d0
   endif

   write(6,'(1x,"Temperature: ",f10.3," K"/)') temp

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Dielectric relaxation model
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options,' EPSMODEL=')

   if (ioption.ne.0) then

      !-- extract the keyword for the model
      ispa = index(options(ioption+10:),space)
      solvent_model = options(ioption+10:ioption+ispa+8)
      write(6,'(/1x,"Dielectric model for dynamical calculations: ",a)') trim(solvent_model)

      if (solvent_model.eq."DEBYE") then
         write(6,'(1x,"(Overdamped solvent dynamics with a single relaxation period)")')
      elseif (solvent_model.eq."DEBYE2") then
         write(6,'(1x,"(Overdamped solvent dynamics with two relaxation periods)")')
      elseif (solvent_model.eq."ONODERA") then
         write(6,'(1x,"(Solvent dynamics with a single relaxation period and effective solvent mass)")')
      elseif (solvent_model.eq."ONODERA2") then
         write(6,'(1x,"(Solvent dynamics with two relaxation periods and effective solvent mass)")')
      elseif (solvent_model.eq."NEWTON") then
         write(6,'(1x,"(Ballistic Solvent dynamics with effective solvent mass and no random forces)")')
      else
         write(6,'(1x,"(UNKNOWN model: abort calculation)")')
         call clean_exit
      endif

   else

      solvent_model = "DEBYE"
      write(6,'(/1x,"Dielectric model for dynamical calculations (default): ",a)') trim(solvent_model)
      write(6,'(1x,"(Overdamped solvent dynamics with a single relaxation period)")')

   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Parameters of the dielectric relaxation model
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(6,'(/1x,"Parameters of the dielectric function ",a/)') trim(solvent_model)

   if (solvent_model.eq."NEWTON") then

      ioption = index(options,' EFFMASS=')
      if (ioption.ne.0) then
         effmass1 = reada(options,ioption+9)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify EFFMASS= option for NEWTON model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' EPS0=')
      if (ioption.ne.0) then
         eps0_dyn = reada(options,ioption+6)
      else
         eps0_dyn = eps0
      endif

      ioption = index(options,' EPS8=')
      if (ioption.ne.0) then
         eps8_dyn = reada(options,ioption+6)
      else
         eps8_dyn = eps8
      endif

      call set_newton_model_parameters()
      write(6,'(1x,"Static dielectric constant EPS0:           ",f15.6)') eps0_dyn
      write(6,'(1x,"Optical dielectric constant EPS_inf:       ",f15.6)') eps8_dyn
      write(6,'(1x,"Inverse Pekar factor f_0 (force constant): ",f15.6)') f0
      write(6,'(1x,"Effective mass of the solvent (ps^2):      ",f15.6)') effmass1

   elseif (solvent_model.eq."DEBYE") then

      ioption = index(options,' TAUD=')
      if (ioption.ne.0) then
         taud = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAUD= option for DEBYE model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' EPS0=')
      if (ioption.ne.0) then
         eps0_dyn = reada(options,ioption+6)
      else
         eps0_dyn = eps0
      endif

      ioption = index(options,' EPS8=')
      if (ioption.ne.0) then
         eps8_dyn = reada(options,ioption+6)
      else
         eps8_dyn = eps8
      endif

      call set_debye_model_parameters()
      write(6,'(1x,"Static dielectric constant EPS0         ",f15.6)') eps0_dyn
      write(6,'(1x,"Optical dielectric constant EPS_inf     ",f15.6)') eps8_dyn
      write(6,'(1x,"Inverse Pekar factor f_0                ",f15.6)') f0
      write(6,'(1x,"Debye relaxation time TAUD (ps):        ",f15.6)') taud
      write(6,'(1x,"Longitudianl relaxation time TAUL (ps): ",f15.6)') taul
      write(6,'(1x,"Effective mass of the solvent (ps^2):   ",f15.6)') effmass1


   elseif (solvent_model.eq."DEBYE2") then

      ioption = index(options,' TAU1=')
      if (ioption.ne.0) then
         tau1 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAU1= option for DEBYE2 model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' TAU2=')
      if (ioption.ne.0) then
         tau2 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAU2= option for DEBYE2 model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' EPS1=')
      if (ioption.ne.0) then
         eps1_dyn = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify EPS1= option for DEBYE2 model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' EPS0=')
      if (ioption.ne.0) then
         eps0_dyn = reada(options,ioption+6)
      else
         eps0_dyn = eps0
      endif

      ioption = index(options,' EPS8=')
      if (ioption.ne.0) then
         eps8_dyn = reada(options,ioption+6)
      else
         eps8_dyn = eps8
      endif

      call set_debye2_model_parameters()
      write(6,'(1x,"Static dielectric constant EPS0         ",f15.6)') eps0_dyn
      write(6,'(1x,"Optical dielectric constant EPS_inf     ",f15.6)') eps8_dyn
      write(6,'(1x,"Additional dielectric constant EPS1     ",f15.6)') eps1_dyn
      write(6,'(1x,"Inverse Pekar factor f_0                ",f15.6)') f0
      write(6,'(1x,"First  relaxation time TAU1 (ps):       ",f15.6)') tau1
      write(6,'(1x,"Second relaxation time TAU2 (ps):       ",f15.6)') tau2
      write(6,'(1x,"Longitudianl relaxation time TAUL2(ps): ",f15.6)') taul
      write(6,'(1x,"Effective mass of the solvent (ps^2):   ",f15.6)') effmass1

   elseif (solvent_model.eq."ONODERA") then

      ioption = index(options,' TAUD=')
      if (ioption.ne.0) then
         taud = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAUD= option for ONODERA model ***"/)')
         call clean_exit
      endif

      ioption = index(options,' TAU0=')
      ioption2 = index(options,' EFFMASS=')

      if (ioption.ne.0) then

         tau0 = reada(options,ioption+6)
         if (tau0.eq.0.d0) then
            write(*,'(/1x,"*** (in DYNAMICSET2): Use EFFMASS option instead of setting TAU0 to zero ***"/)')
            call clean_exit
         endif

      elseif (ioption2.ne.0) then

         tau0 = 0.d0
         ikey = ioption2 + 9
         effmass1 = reada(options,ikey)

      else

         write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify either TAU0 or EFFMASS option for ONODERA model ***"/)')
         call clean_exit

      endif

      !-- dielectric constants

      ioption = index(options,' EPS0=')
      if (ioption.ne.0) then
         eps0_dyn = reada(options,ioption+6)
      else
         eps0_dyn = eps0
      endif

      ioption = index(options,' EPS8=')
      if (ioption.ne.0) then
         eps8_dyn = reada(options,ioption+6)
      else
         eps8_dyn = eps8
      endif

      call set_onodera_model_parameters()
      write(6,'(/1x,"Dielectric function parameters")')
      write(6,'( 1x,"------------------------------------------------------------")')
      write(6,'( 1x,"Static dielectric constant EPS0:         ",f15.6)') eps0_dyn
      write(6,'( 1x,"Optical dielectric constant EPS_inf:     ",f15.6)') eps8_dyn
      write(6,'( 1x,"Inverse Pekar factor f_0:                ",f15.6)') f0
      write(6,'( 1x,"Debye   relaxation time TAUD (ps):       ",f15.6)') taud
      write(6,'( 1x,"Onodera inertial   time TAU0 (ps):       ",f15.6)') tau0
      write(6,'( 1x,"Longitudinal relaxation time TAUL  (ps): ",f15.6)') taul
      write(6,'( 1x,"Longitudinal relaxation time TAU0L (ps): ",f15.6)') tau0l
      write(6,'( 1x,"Longitudinal relaxation time TAUL~ (ps): ",f15.6)') taul_tilde

      write(6,'(/1x,"Langevin equation parameters")')
      write(6,'(/1x,"mu*dv/dt + eta*v = F_x(t) + R(t)")')
      write(6,'( 1x,"------------------------------------------------------------")')
      write(6,'( 1x,"Effective mass of the solvent (ps^2):    ",f15.6)') effmass1
      write(6,'( 1x,"Friction parameter eta (ps):             ",f15.6)') f0*taul_tilde
      write(6,'( 1x,"Friction coefficient g=eta/mass (1/ps):  ",f15.6)') (tau0+taud)/tau0/taud

   !-- Onodera model with two relaxation periods
   elseif (solvent_model.eq."ONODERA2") then

      ioption = index(options,' EPS0=')
      if (ioption.ne.0) then
         eps0_dyn = reada(options,ioption+6)
      else
         eps0_dyn = eps0
      endif

      ioption = index(options,' EPS8=')
      if (ioption.ne.0) then
         eps8_dyn = reada(options,ioption+6)
      else
         eps8_dyn = eps8
      endif

      if (index(options,' EPSPARS').ne.0) then

         !-- parameters of the dielectric function

         ioption = index(options,' TAU0=')
         ioption2 = index(options,' EFFMASS=')

         if (ioption.ne.0) then

            tau0 = reada(options,ioption+6)
            if (tau0.eq.0.d0) then
               write(*,'(/1x,"*** (in DYNAMICSET2): Use EFFMASS option instead of setting TAU0 to zero ***"/)')
               call clean_exit
            endif

         elseif (ioption2.ne.0) then

            tau0 = 0.d0
            ikey = ioption2 + 9
            effmass1 = reada(options,ikey)

         else

            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify either TAU0 or EFFMASS option for ONODERA-2 model ***"/)')
            call clean_exit

         endif

         ioption = index(options,' TAU1=')
         if (ioption.ne.0) then
            tau1 = reada(options,ioption+6)
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAU1= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         ioption = index(options,' TAU2=')
         if (ioption.ne.0) then
            tau2 = reada(options,ioption+6)
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAU2= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         ioption = index(options,' EPS1=')
         if (ioption.ne.0) then
            eps1_dyn = reada(options,ioption+6)
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify EPS1= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         call set_onodera2_model_parameters()
         write(6,'(/1x,"Dielectric function parameters")')
         write(6,'( 1x,"------------------------------------------------------------")')
         write(6,'( 1x,"Static dielectric constant EPS0        ",f15.6)') eps0_dyn
         write(6,'( 1x,"Optical dielectric constant EPS_inf    ",f15.6)') eps8_dyn
         write(6,'( 1x,"Inverse Pekar factor f_0               ",f15.6)') f0
         write(6,'( 1x,"Additional dielectric constant EPS1    ",f15.6)') eps1_dyn
         write(6,'( 1x,"Onodera relaxation time TAU0 (ps)      ",f15.6)') tau0
         write(6,'( 1x,"First   relaxation time TAU1 (ps)      ",f15.6)') tau1
         write(6,'( 1x,"Second  relaxation time TAU2 (ps)      ",f15.6)') tau2
         write(6,'( 1x,"Longitudinal relaxation time TAUL2 (ps)",f15.6)') eps8_dyn*tau2/eps0_dyn

         write(6,'(/1x,"Generalized Langevin equation parameters")')
         write(6,'(/1x,"mu*dv/dt + eta*v + int[gamma*exp[-(t-tau)/tau_alpha)]*v(tau)*dtau = F_x(t) + R(t)")')
         write(6,'( 1x,"----------------------------------------------------------------------------------")')
         write(6,'( 1x,"Effective mass of the solvent (ps^2)    ",f15.6)') effmass1
         write(6,'( 1x,"Friction parameter eta (ps):            ",f15.6)') etax
         write(6,'( 1x,"Friction coefficient eta/mass (1/ps):   ",f15.6)') etax/effmass1
         write(6,'( 1x,"Friction parameter eta_y (ps):          ",f15.6)') etay
         write(6,'( 1x,"Memory function parameter gamma:        ",f15.6)') gamma
         write(6,'( 1x,"Memory function timescale tau_a (ps):   ",f15.6)') taualpha


      elseif (index(options,' GLEPARS').ne.0) then

         !-- parameters of the GLE

         ioption = index(options,' EFFMASS=')
         if (ioption.ne.0) then
            effmass1 = reada(options,ioption+9)
            effmass2 = effmass1
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify EFFMASS= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         ioption = index(options,' GAMMA=')
         if (ioption.ne.0) then
            gamma = reada(options,ioption+7)
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify GAMMA= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         ioption = index(options,' TAUALPHA=')
         if (ioption.ne.0) then
            taualpha = reada(options,ioption+10)
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify TAUALPHA= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         ioption = index(options,' ETA=')
         if (ioption.ne.0) then
            etax = reada(options,ioption+5)
            etay = gamma*taualpha
         else
            write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify ETA= option for ONODERA-2 model ***"/)')
            call clean_exit
         endif

         write(6,'(/1x,"Dielectric function parameters")')
         write(6,'( 1x,"------------------------------------------------------------")')
         write(6,'( 1x,"Static dielectric constant EPS0        ",f15.6)') eps0_dyn
         write(6,'( 1x,"Optical dielectric constant EPS_inf    ",f15.6)') eps8_dyn
         write(6,'( 1x,"Inverse Pekar factor f_0               ",f15.6)') f0

         write(6,'(/1x,"Generalized Langevin equation parameters")')
         write(6,'(/1x,"mu*dv/dt + eta*v + int[gamma*exp[-(t-tau)/tau_alpha)]*v(tau)*dtau = F_x(t) + R(t)")')
         write(6,'( 1x,"----------------------------------------------------------------------------------")')
         write(6,'( 1x,"Effective mass of the solvent (ps^2)    ",f15.6)') effmass1
         write(6,'( 1x,"Friction parameter eta (ps):            ",f15.6)') etax
         write(6,'( 1x,"Friction coefficient eta/mass (1/ps):   ",f15.6)') etax/effmass1
         write(6,'( 1x,"Friction parameter eta_y (ps):          ",f15.6)') etay
         write(6,'( 1x,"Memory function parameter gamma:        ",f15.6)') gamma
         write(6,'( 1x,"Memory function timescale tau_a (ps):   ",f15.6)') taualpha

      else

         write(*,'(/1x,"*** (in DYNAMICSET2): You must specify either EPSPARS or GLEPARS, check your input ***"/)')
         call clean_exit

      endif


      if (effmass1.eq.0.d0) then
         write(*,'(/1x,"*** (in DYNAMICSET2): The effective solvent mass MUST NOT BE ZERO, check your input ***"/)')
         call clean_exit
      endif

   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Type of dynamics (classical or MDQT)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' MDQT').ne.0) then

      if (solvent_model.eq."DEBYE".or.solvent_model.eq."DEBYE2") then
         write(6,'(/1x,"MDQT is not compatible with overdamped diffusive dynamics (DEBYE and DEBYE2).")')
         write(6,'( 1x,"Please choose ONODERA or ONODERA2 model for solvent dynamics.")')
         write(6,'( 1x,"Check your input and make up your mind!!! Aborting...")')
         call clean_exit
      endif

      mdqt = .true.
      write(6,'(/1x,"Mixed quantum-classical dynamics on two electronic free energy surfaces.",/,&
                &1x,"Tullys fewest switches surface hopping algorithm (MDQT) will be utilized."/)')

      if (index(options,' PHASE').ne.0) then
         phase_corr = .true.
         write(6,'(/1x,"Phase correction algorithm will be used.",/,&
                   &1x,"[N. Shenvi, J. E. Subotnik, and W. Yang, J. Chem. Phys. 135, 024101 (2011) ]"/)')
      endif

      !-- velocity reversal option

      if (index(options,' REVVEL').ne.0) then
         revvel = .true.
         if (index(options,' REVVEL0').ne.0) revvel_all = .true.
         write(6,'(/1x,"Velocity reversal will be performed in case of frustrated hops.")')
         if (revvel_all) then
            write(6,'(1x,"Original algorithm: Velocity reversal will be ALWAYS performed in case of frustrated hops.",/)')
         else
            write(6,'(1x,"Truhlar criterion:",/,&
                   &1x,"[A. Jain and J. E. Subotnik. Surface hopping, transition state theory and decoherence",/,&
                   &1x," 2: Thermal rate constants and detailed balance. J. Chem. Phys., Preprint (2015)]"/)')
         endif
      endif

      !-- decoherence options

      idecoherence = 0
      if (index(options,' AFSSH').ne.0) idecoherence = idecoherence + 1
      if (index(options,' COLLAPSE_REGION_COUPLING').ne.0) idecoherence = idecoherence + 1
      if (index(options,' IDS').ne.0) idecoherence = idecoherence + 1
      if (index(options,' IDA').ne.0) idecoherence = idecoherence + 1
      if (index(options,' GEDC').ne.0) idecoherence = idecoherence + 1


      if (idecoherence.gt.0) then

         !-- make sure that only one decoherence option has been chosen
         if (idecoherence.gt.1) then
            write(6,'(/1x,"Only one decoherence option can be specified.")')
            write(6,'( 1x,"Your input contains more then one decoherence option in DYNAMICSET2 keyword.")')
            write(6,'( 1x,"Choose one of: AFSSH, COLLAPSE_REGION_COUPLING, IDS, IDA, or GEDC. ")')
            write(6,'( 1x,"Check your input and make up your mind!!! Aborting...")')
            call clean_exit
         endif

         decoherence = .true.

      endif


      if (index(options,' AFSSH').ne.0) then

         if (.not.phase_corr) then

            afssh = .true.
            collapse_region_coupling = .false.
            ids = .false.
            ida = .false.
            gedc = .false.

            ioption = index(options," DZETA=")
            if (ioption.ne.0) then
               dzeta = reada(options,ioption+7)
            else
               dzeta = 1.d0
            endif
            write(6,'(/1x,"Decoherence algorithm (AFSSH) with dzeta =",f8.3," will be used.",/,&
                      &1x,"[B. R. Landry, N. Shenvi, J. E. Subotnik, 2012]"/)') dzeta

            if (index(options," COUPLE").ne.0) then
               decouple = .false.
               write(6,'(/1x,"The TDSE in decoherence algorithm (AFSSH) will be coupled to EOM for the moments (Eq.18)")')
            else
               decouple = .true.
               write(6,'(/1x,"The TDSE in decoherence algorithm (AFSSH) will be decoupled from EOM for the moments (original algorithm)")')
            endif

            if (index(options," ADJUST_ALONG_MOMENTS").ne.0) then
               along_moments = .true.
               write(6,'(/1x,"The adjustment of velocities will be performed along the direction of the difference")')
               write(6,'( 1x,"between the moments of momenta in decoherence algorithm (AFSSH-0)")')
            else
               along_moments = .false.
            endif

         else

            write(6,'(/1x,"In current version decoherence (AFSSH) and phase-correction algorithms are incompatible")')
            write(6,'( 1x,"Your input contains both PHASE and DECOHERENCE options in DYNAMICS3 keyword.")')
            write(6,'( 1x,"Check your input and make up your mind!!! Aborting...")')
            call clean_exit

         endif


      elseif (index(options,' COLLAPSE_REGION_COUPLING').ne.0) then

         collapse_region_coupling = .true.
         afssh = .false.
         ids = .false.
         ida = .false.
         gedc = .false.

         ioption = index(options," COUPLING_CUTOFF=")
         if (ioption.ne.0) then
            coupling_cutoff = reada(options,ioption+17)
         else
            coupling_cutoff = 1.d-5
         endif
         write(6,'(/1x,"Simple decoherence algorithm with collapsing events occuring upon leaving the interaction region")')
         write(6,'( 1x,"The interaction region is defined as region where the largest nonadiabatic coupling")')
         write(6,'( 1x,"is smaller in magnitude than the cutoff value of ",g15.6," (kcal/mol)^{-1/2}"/)') coupling_cutoff


      elseif (index(options,' IDS').ne.0) then

         ids = .true.
         collapse_region_coupling = .false.
         afssh = .false.
         ida = .false.
         gedc = .false.
         write(6,'(/1x,"Instantaneous decoherence algorithm with collapsing events occuring upon succesful hops (IDS)")')


      elseif (index(options,' IDA').ne.0) then

         ida = .true.
         collapse_region_coupling = .false.
         afssh = .false.
         ids = .false.
         gedc = .false.
         write(6,'(/1x,"Instantaneous decoherence algorithm with collapsing events occuring upon any hop (IDA)")')

      elseif (index(options,' GEDC').ne.0) then

         gedc = .true.
         ida = .false.
         collapse_region_coupling = .false.
         afssh = .false.
         ids = .false.
         write(6,'(/1x,"Granucci-s Energy-based Decoherence correction (GEDC) with C=1 and E0=0.1au")')

      endif


   else

      mdqt = .false.
      write(6,'(/1x,"Simulation of classical dynamics on a single free energy surface (default)."/)')

   endif

   if (index(options,' REACTIVE_FLUX').ne.0) then
      reactive_flux = .true.
      Fnfile = 15
      Fdfile = 16
      open(unit=Fnfile,file='Fn',action='write')
      open(unit=Fdfile,file='Fd',action='write')
      write(6,'(1X,"Begin the calculation at the dividing surface. Propagate trajectories backward and forward in time.")')
      if (.not.mdqt) then
         write(6,'(1x,"The reactive_flux option can only be used when MDQT is turned on."/)')
         write(6,'(1x,"MDQT is not turned on.  Please remove REACTIVE_FLUX or include MDQT in the DYNAMICSET2 group."/)')
         Call clean_exit
      end if
   else
      reactive_flux = .false.
      normalforward = .false.
      successfulreverse = .false.
   end if


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! flag for EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (index(options,' NOWEIGHTS').ne.0) then
      weights = .false.
      write(6,'(/1x,"EVB weights WILL NOT be calculated along the trajectory."/)')
   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Surface type (ADIAB, DIAB)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then

      write(6,'(/1x,"Solvent dynamics on adiabatic free energy surfaces."/)')

      mode_dyn  = 'ADIAB'
      adiab = .true.
      diab4 = .false.
      ielst_dyn = 2
      iset_dyn = 1

      ioption = index(options,' ADIAB=')

      if (ioption.eq.0.) then

         ireac = 1
         iprod = 3

      else

         ikey = ioption + 7
         ispa = index(options(ikey:),space)
         islash1 = index(options(ikey:ikey+ispa-1),'/')
         if (islash1.ne.0) then
            ireac = reada(options(ikey:ikey+islash1-2),1)
            iprod = reada(options,ikey+islash1)
         else
            write(*,'(//1x,"*** (in DYNAMICSET2): BOTH reactant and product states should be specified in ADIAB= keyword."/)')
            write(*,'(//1x,"*** (in DYNAMICSET2): Example: ADIAB=1/3"/)')
            call clean_exit
         endif

      endif

   elseif (index(options,' DIAB=').ne.0) then

      ioption = index(options,' DIAB=')
      mode_dyn  = 'DIAB4'
      adiab = .false.
      diab4 = .true.
      ielst_dyn = 1
      iset_dyn = reada(options,ioption+6)

      if (iset_dyn.eq.1) then
         write(6,'(/1x,"Classical solvent dynamics on the 1-st (1a) diabatic surface")')
      elseif (iset_dyn.eq.2) then
         write(6,'(/1x,"Classical solvent dynamics on the 2-nd (1b) diabatic surface")')
      elseif (iset_dyn.eq.3) then
         write(6,'(/1x,"Classical solvent dynamics on the 3-rd (2a) diabatic surface")')
      elseif (iset_dyn.eq.4) then
         write(6,'(/1x,"Classical solvent dynamics on the 4-th (2b) diabatic surface")')
      else
         write(*,'(/1x,"*** (in DYNAMICSET2): subset in DIAB keyword must be 1 or 2 or 3 or 4 ***"/)')
         call clean_exit
      endif

      if (mdqt) then
         mdqt = .false.
         write(6,'(/1x,"****** Only classical solvent dynamics on the diabatic free energy surface is available")')
         write(6,'( 1x,"****** WARNING: MDQT option has been turned off")')
      endif

   else

      mode_dyn  = 'ADIAB'
      adiab = .true.
      diab4 = .false.
      ielst_dyn = 2
      iset_dyn = 1
      ireac = 1
      iprod = 3
      write(6,'(/1x,"Solvent dynamics on ET adiabatic free energy surfaces"/)')
      write(6,'( 1x,"(by default, the diabatic states 1 (1a) and 3 (2a) form the basis)."/)')

   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initialize 2x2 gas-phase Hamiltonian
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (igas.eq.0) then

      ! constant potential for a fixed solute
      call h0mat_constant(h0k)

   else

      write(*,'(/1x,"****** Only constant gas-phase potential can be used in DYNAMICSET2")')
      call clean_exit

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Form the gas phase 2x2 Hamiltonian matrix for the ET
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hreac = h0k(ireac,ireac)
   hprod = h0k(iprod,iprod)
   vet   = h0k(ireac,iprod)
   hg2(1,1) = hreac
   hg2(2,2) = hprod
   hg2(1,2) = vet
   hg2(2,1) = vet
   call primat(hg2,2,2,2,5,6,0,'[2x2] gas phase ET Hamiltonian (kcal/mol)')


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initialize 2x2 reorganization energy matrices
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   if (isolv.eq.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !  Read partial reorganization energies and reconstruct
      !  full reorganization energy matrices (TSET option).
      !
      !  There are two alternative sets of values which must be
      !  specified to fully define the reorganization energy matrices.
      !
      !  Set A:
      !  - outer-sphere reorganization energy (lambda)
      !  - solvation energy of the reactant state (gsolv_1)
      !  - the value of the equilibrium interaction energy gap
      !    in the reactant state (z_1)
      !
      !  Set B:
      !  - outer-sphere reorganization energy (lambda)
      !  - solvation energy of the reactant state (gsolv_1)
      !  - solvation energy of the product state (gsolv_2)
      !
      !  One of these sets must be specified in the input.
      !  Note that only
      !
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ikey = index(keywrd,' SOLV(')
      call getopt(keywrd,ikey+6,solvoptions)

      write(6,'(/1x,''[2x2] Reorganization energy matrix will be reconstructed from the following input.'')')
      write(6,'( 1x,''(It is assumed that electronic reorganization energies are zero...)'')')

      !-- Read reorganization energy

      ikey = index(solvoptions,' ER=')
      if (ikey.ne.0) then
         lambda = reada(solvoptions,ikey+4)
      else
         write(*,'(/1x,"*** (INPUT ERROR in DYNAMICSET2): you must specify the ER reorganization energy ***"/)')
         call clean_exit
      endif

      !-- Read solvation free energy of reactant diabatic state
      !   (default value is zero)

      ikey = index(solvoptions,' GSOLV1=')
      if (ikey.ne.0) then
         gsolv_1 = reada(solvoptions,ikey+8)
      else
         gsolv_1 = 0.d0
      endif

      write(6,'(/1x,"Solvation free energy of the reactant diabatic state (from input): ",f15.3," kcal/mol")') gsolv_1

      !-- Read energy gap value at the minimum of the reactant diabatic free energy surface
      !   or equilibrium solvation energy of the product state
      !   (default values are zero)

      ikeya = index(solvoptions,' MIN1=')
      ikeyb = index(solvoptions,' GSOLV2=')

      if (ikeya.ne.0) then

         ikey = ikeya + 6
         z_1 = reada(solvoptions,ikey)
         z_2 = z_1 - 2.d0*lambda

         t2r(1,1) = -2.d0*gsolv_1
         t2r(1,2) = -z_1
         t2r(2,1) = -z_1
         t2r(2,2) = 2.d0*lambda

         gsolv_2 = z_1 + gsolv_1 - lambda

         write(6,'(1x,"Reorganization free energy (from input):                                ",f12.3," kcal/mol")') lambda
         write(6,'(1x,"Solvation free energy of the product diabatic state (calculated):       ",f12.3," kcal/mol")') gsolv_2
         write(6,'(1x,"Difference in solvation energies of diabatic states (calculated):       ",f12.3," kcal/mol")') gsolv_2-gsolv_1
         write(6,'(1x,"Value of the equilibrium energy gap in the reactant state (from input): ",f12.3," kcal/mol")') z_1
         write(6,'(1x,"Value of the equilibrium energy gap in the product  state (calculated): ",f12.3," kcal/mol")') z_2
         write(6,'(1x,"(GSOLV2= option is ignored if present)")')

      elseif (ikeyb.ne.0) then

         ikey = ikeyb + 8
         gsolv_2 = reada(solvoptions,ikey)

         z_1 = gsolv_2 - gsolv_1 + lambda
         z_2 = z_1 - 2.d0*lambda

         t2r(1,1) = -2.d0*gsolv_1
         t2r(1,2) = -z_1
         t2r(2,1) = -z_1
         t2r(2,2) = 2.d0*lambda

         write(6,'(1x,"Reorganization free energy (from input):                                ",f12.3," kcal/mol")') lambda
         write(6,'(1x,"Solvation free energy of the product diabatic state (from input):       ",f12.3," kcal/mol")') gsolv_2
         write(6,'(1x,"Difference in solvation energies of diabatic states (calculated):       ",f12.3," kcal/mol")') gsolv_2-gsolv_1
         write(6,'(1x,"Value of the equilibrium energy gap in the reactant state (calculated): ",f12.3," kcal/mol")') z_1
         write(6,'(1x,"Value of the equilibrium energy gap in the product  state (calculated): ",f12.3," kcal/mol")') z_2
         write(6,'(1x,"(MIN1= option is ignored if present)")')

      else

         write(*,'(//1x,"*** (INPUT ERROR in DYNAMICSET2): you must specify either MIN1= or GSOLV2= (not both!) ***"/)')
         call clean_exit

      endif

      !-- What about electronic reorganization energy matrices?

      !-- (for now we set them to zero assuming that the electronic solvation energies
      !    are included in the diagonal elements of the gas phase Hamiltonian)

      t2rinf(1,1) = 0.d0
      t2rinf(1,2) = 0.d0
      t2rinf(2,1) = 0.d0
      t2rinf(2,2) = 0.d0

      !-- initialize the reorganization energy matrices

      t2(1,1) = t2r(1,1)
      t2(1,2) = t2r(1,2) + t2r(1,1)
      t2(2,1) = t2(1,2)
      t2(2,2) = t2r(2,2) + t2r(1,2) + t2r(1,2) + t2r(1,1)

      t2inf(1,1) = t2rinf(1,1)
      t2inf(1,2) = t2rinf(1,2) + t2rinf(1,1)
      t2inf(2,1) = t2inf(1,2)
      t2inf(2,2) = t2rinf(2,2) + t2rinf(1,2) + t2rinf(1,2) + t2rinf(1,1)


   elseif (isolv.eq.1.or.isolv.eq.2) then

      if (isolv.eq.1) then

         !-- Ellipsoidal model
         call tmat(tk,tinfk,trk,trinfk)

      elseif (isolv.eq.2) then

         !-- FRCM model
         call solint(tk,tinfk,trk,trinfk)

      endif

      t2r(1,1) = trk(ireac,ireac)
      t2r(1,2) = trk(ireac,iprod)
      t2r(2,1) = trk(ireac,iprod)
      t2r(2,2) = trk(iprod,iprod)

      t2rinf(1,1) = trinfk(ireac,ireac)
      t2rinf(1,2) = trinfk(ireac,iprod)
      t2rinf(2,1) = trinfk(iprod,ireac)
      t2rinf(2,2) = trinfk(iprod,iprod)

      t2(1,1) = tk(ireac,ireac)
      t2(1,2) = tk(ireac,iprod)
      t2(2,1) = tk(iprod,ireac)
      t2(2,2) = tk(iprod,iprod)

      t2inf(1,1) = tinfk(ireac,ireac)
      t2inf(1,2) = tinfk(ireac,iprod)
      t2inf(2,1) = tinfk(iprod,ireac)
      t2inf(2,2) = tinfk(iprod,iprod)

      gsolv_1 = -0.5d0*t2(1,1)
      gsolv_2 = -0.5d0*t2(2,2)
      lambda  = 0.5d0*t2r(2,2)

      z_1 = gsolv_2 - gsolv_1 + lambda
      z_2 = z_1 - 2.d0*lambda

      write(6,'(/1x,"Reorganization free energy:                                ",f12.3," kcal/mol")') lambda
      write(6,'( 1x,"Solvation free energy of the product diabatic state:       ",f12.3," kcal/mol")') gsolv_2
      write(6,'( 1x,"Difference in solvation energies of diabatic states:       ",f12.3," kcal/mol")') gsolv_2-gsolv_1
      write(6,'( 1x,"Value of the equilibrium energy gap in the reactant state: ",f12.3," kcal/mol")') z_1
      write(6,'( 1x,"Value of the equilibrium energy gap in the product  state: ",f12.3," kcal/mol")') z_2

   else

      write(*,'(//1x,"ERROR in DYNAMICSET2: unrecognized solvation model ISOLV =",i2)') isolv
      write(*,'(  1x,"ISOLV should be 0(TSET), 1(ELLIPSE), or 2(FRCM)")')
      call clean_exit

   endif

   call primat(t2r,2,2,2,5,6,0,'[2x2] inertial reorganization energy matrix (kcal/mol)')
   call primat(t2r,2,2,2,5,6,0,'[2x2] reduced inertial reorganization energy matrix (kcal/mol)')

   !-- calculate parameters for transformation between the energy gap coordinate and
   !   dynamical solvent coordinate

   lambda = 0.5d0*t2r(2,2)
   scale_factor = sqrt(2.d0*f0*lambda)
   delta_shift  = t2r(1,2)/scale_factor

   !-- Calculate Marcus nonadiabatic rate constant

   dg_reaction    = (hg2(2,2) + 0.5*t2inf(2,2) + gsolv_2) - (hg2(1,1) + 0.5*t2inf(1,1) + gsolv_1)
   dg_activation  = (dg_reaction + lambda)**2/(4.d0*lambda)
   et_marcus_rate = 1.d12*(sqrt(pi/(kb*temp))/hbarps)*vet*vet*exp(-dg_activation/(kb*temp))/dsqrt(lambda)

   write(*,'(/1x,"========================================================")')
   write(*,'( 1x,"Marcus model parameters and nonadiabatic rates")')
   write(*,'( 1x,"--------------------------------------------------------")')
   write(*,'( 1x,"Reorganization free energy: ",f12.3," kcal/mol")') lambda
   write(*,'( 1x,"ET reaction free energy:    ",f12.3," kcal/mol")') dg_reaction
   write(*,'( 1x,"ET activation free energy:  ",f12.3," kcal/mol")') dg_activation
   write(*,'( 1x,"ET Marcus rate constant:    ",e16.9," sec^-1")')   et_marcus_rate


   !-- Calculate Rips-Jortner nonadiabatic rate constant (solvent control)

   if (solvent_model.eq."ONODERA2".and.index(options,' GLEPARS').ne.0) then

      et_rips_jortner_rate = 0.d0
      kappa_ad = 0.d0
      write(*,'( 1x,"--------------------------------------------------------")')
      write(*,'( 1x,"Longitudinal relaxation time is not defined")')

   else

      kappa_ad = 4.d0*pi*(vet*vet*cal2au*cal2au)*taul*ps2au/(lambda*cal2au)
      prefactor1 = (2.d0*pi/au2ps)*(vet*vet*cal2au*cal2au)/sqrt(4.d0*pi*cal2au*cal2au*lambda*kb*temp)
      prefactor2 = 1.d0 + 4.d0*pi*(vet*vet*cal2au*cal2au)*taul*ps2au/(lambda*cal2au)
      et_rips_jortner_rate = 1.d12*(prefactor1/prefactor2)*exp(-dg_activation/(kb*temp))
      write(*,'( 1x,"--------------------------------------------------------")')
      write(*,'( 1x,"Rips-Jortner adiabaticity:  ",f12.3," kcal/mol")') kappa_ad
      write(*,'( 1x,"Rips-Jortner rate constant: ",e16.9," sec^-1")')   et_rips_jortner_rate

   endif

   write(*,'(1x,"========================================================")')

   !-- Calculate Kramers (For Debye-1, Onodera-1) or
   !   Kramers-Grote-Hynes (For Debye-2, Onodera-2)
   !   adiabatic rate constant

   call adiabatic_fes_et2_parameters(double_well,fr,fb,eb)

   if (double_well) then

      if (solvent_model.eq."NEWTON") then
         adiabatic_rate = tst_rate(temp,fr,eb)
      elseif (solvent_model.eq."DEBYE") then
         adiabatic_rate = kgh_rate_debye1(temp,fr,fb,eb)
      elseif (solvent_model.eq."ONODERA") then
         adiabatic_rate = kgh_rate_onodera1(temp,fr,fb,eb)
      elseif (solvent_model.eq."DEBYE2") then
         adiabatic_rate = kgh_rate_debye2(temp,fr,fb,eb)
      elseif (solvent_model.eq."ONODERA2") then
         adiabatic_rate = kgh_rate_onodera2(temp,fr,fb,eb)
      endif

      write(*,'( 1x,"Adiabatic rate constant:    ",e16.9," sec^-1")')   adiabatic_rate

   else

      write(*,'( 1x,"Single-well adiabatic FES: adiabatic rate is not defined")')

   endif

   write(*,'(1x,"========================================================")')

   !-- print out the free energy profiles to the external file

   fname = job(1:ljob)//"/"
   open(ifes_channel,file=trim(fname)//"ET_free_energy_profiles.dat")

   write(ifes_channel,'("#",74("="))')
   write(ifes_channel,'("#   ET free energy profiles (units are based on kcal/mol)")')
   write(ifes_channel,'("#   Reorganization energy:      ",f12.3," kcal/mol")') lambda
   write(ifes_channel,'("#   Scaling factor:             ",f12.3," sqrt(kcal/mol)")') scale_factor
   write(ifes_channel,'("#   Shift:                      ",f12.3," sqrt(kcal/mol)")') delta_shift
   write(ifes_channel,'("#   ET reaction free energy:    ",f12.3," kcal/mol")') dg_reaction
   write(ifes_channel,'("#   ET activation free energy:  ",f12.3," kcal/mol")') dg_activation
   write(ifes_channel,'("#   Gap at the reactant minimum:",f12.3," kcal/mol")') z_1
   write(ifes_channel,'("#   Gap at the product  minimum:",f12.3," kcal/mol")') z_2
   write(ifes_channel,'("#   ET Marcus rate constant:    ",e20.9," 1/sec")') et_marcus_rate
   write(ifes_channel,'("#",74("-"))')

   if (solvent_model.ne."ONODERA2".or.index(options,' GLEPARS').eq.0) then
      write(ifes_channel,'("#   Rips-Jortner adiabaticity:  ",e20.9," 1/sec")') kappa_ad
      write(ifes_channel,'("#   Rips-Jortner rate constant: ",e20.9," 1/sec")') et_rips_jortner_rate
      write(ifes_channel,'("#",74("-"))')
   endif

   if (double_well) then
      write(ifes_channel,'("#   Adiabatic activation energy:  ",e20.9," kcal/mol")') eb
      write(ifes_channel,'("#   Adiabatic rate constant:      ",e20.9," 1/sec")') adiabatic_rate
      write(ifes_channel,'("#",74("-"))')
   endif

   write(ifes_channel,'("#",t7,"Z(gap)",t19,"z(scaled)",t30,"U1(diab)",t42,"U2(diab)",t54,"U1(adiab)",t66,"U2(adiab)")')
   write(ifes_channel,'("#",74("-"))')

   do i=1,101
      zi = 3.d0*lambda - (i-1)*6.d0*lambda/100.d0
      call ze_to_z1(zi,zi_scaled)
      call fes_et2(mode="DIAB4",energy_gap=zi,free_energy=fe_diab)
      call fes_et2(mode="ADIAB",energy_gap=zi,free_energy=fe_adiab)
      write(ifes_channel,'(f13.6,5f12.5)') zi, zi_scaled, fe_diab(1), fe_diab(2), fe_adiab(1), fe_adiab(2)
   enddo

   write(ifes_channel,'("#",74("="))')
   close(ifes_channel)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of states to include in dynamics
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   nstates_dyn = 2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of trajectories
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NTRAJ=")

   if (ioption.ne.0) then
      ntraj = reada(options,ioption+7)
   else
      ntraj = 1
   endif
   write(6,'(1x,"Number of trajectories to generate: ",i4/)') ntraj

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Timesteps
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," TSTEP=")
   if (ioption.ne.0) then
      tstep = reada(options,ioption+7)
   else
      tstep = 0.0005d0
   endif

   write(6,'(1x,"Timestep for solvent dynamics: ",g15.6," ps"/)') tstep

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of steps in each trajectory
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NSTEPS=")

   if (ioption.ne.0) then
      nsteps = reada(options,ioption+8)
      write(6,'(1x,"Number of steps in each trajectory: ",i10/)') nsteps
   else
      nsteps = 100
      write(6,'(1x,"Number of steps in each trajectory (default value): ",i10/)') nsteps
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of steps in TDSE (for MDQT)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (mdqt) then

      ioption = index(options," NQSTEPS=")
      if (ioption.ne.0) then
         itmp = reada(options,ioption+9)
         write(6,'(1x,"Number of TDSE steps per classical step in MDQT: ",i10/)') itmp
      else
         itmp = 100
         write(6,'(1x,"Number of TDSE steps per classical step in MDQT (default value): ",i10/)') itmp
      endif
      write(6,'(1x,"Timestep for TDSE: ",g15.6," ps"/)') tstep/real(itmp,kind=8)

      ioption = index(options," MAXNQSTEPS=")
      if (ioption.ne.0) then
         itmp1 = reada(options,ioption+12)
         write(6,'(1x,"Maximum number of TDSE steps per classical step in MDQT: ",i10/)') itmp1
      else
         itmp1 = 10000
         write(6,'(1x,"Maximum number of TDSE steps per classical step in MDQT (default value): ",i10/)') itmp1
      endif
      write(6,'(1x,"Minimum timestep for TDSE: ",g15.6," ps"/)') tstep/real(itmp1,kind=8)

      call set_tdse_timestep(itmp,itmp1,tstep)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Interpolation scheme for nonadiabatic coupling term in TDSE
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (mdqt) then

      ioption = index(options,' INTERPOLATION=')

      if (ioption.ne.0) then

         !-- extract the keyword for the interpolation scheme
         ispa = index(options(ioption+15:),space)
         interpolation = options(ioption+15:ioption+ispa+13)
         write(6,'(/1x,"Interpolation scheme for the nonadiabatic coupling term in TDSE: ",a)') trim(interpolation)

         if (interpolation.eq."LINEAR") then
            write(6,'(1x,"(linear interpolation using the values at t and t+dt)")')
         elseif (interpolation.eq."QUADRATIC") then
            write(6,'(1x,"(quadratic interpolation using the values at t, t+dt/2, and t+dt)")')
         else
            write(6,'(1x,"(ERROR in DYNAMICS3: UNKNOWN interpolation scheme. Check your input file.)")')
            call clean_exit
         endif

      else

         interpolation = "LINEAR"
         write(6,'(/1x,"Interpolation scheme for the nonadiabatic coupling term in TDSE (default): ",a)') trim(interpolation)
         write(6,'(1x,"(linear interpolation using the values at t and t+dt)")')

      endif

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Trajectory output frequency (every NDUMP steps)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NDUMP=")

   if (ioption.ne.0) then
      ndump = reada(options,ioption+7)
      write(6,'(1x,"Dump trajectory data every ",i10," steps"/)') ndump
   else
      ndump = 1
      write(6,'(1x,"Dump trajectory data every single step (default)"/)')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Trajectory screen output frequency (every NDUMP6 steps)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NDUMP6=")

   if (ioption.ne.0) then
      ndump6 = reada(options,ioption+8)
      write(6,'(1x,"Dump trajectory data to screen every ",i10," steps"/)') ndump6
   else
      ndump6 = 0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Populations and coherences output frequency (every NDUMP777 steps)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NDUMP777=")

   if (ioption.ne.0) then
      ndump777 = reada(options,ioption+10)
      write(6,'(1x,"Dump populations and coherences every ",i10," steps"/)') ndump777
   else
      ndump777 = 0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Couplings output frequency (every NDUMP888 steps)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NDUMP888=")

   if (ioption.ne.0) then
      ndump888 = reada(options,ioption+10)
      write(6,'(1x,"Dump couplings every ",i10," steps"/)') ndump888
   else
      ndump888 = 0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! AFSSH moments output frequency (every NDUMP999 steps)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NDUMP999=")

   if (ioption.ne.0) then
      ndump999 = reada(options,ioption+10)
      write(6,'(1x,"Dump expectation values of moments for all states every ",i10," steps"/)') ndump999
   else
      ndump999 = 0
   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Random seed for RAN2NR (negative to reinitialize)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," SEED=")

   if (ioption.ne.0) then

      if (options(ioption+6:ioption+10).eq."PBSID") then

         call set_random_seed("PBS_JOBID")
         write(6,'(1x,"Random seed for RAN2NR (from PBS variable PBS_JOBID): ",i6/)') iseed

      elseif (options(ioption+6:ioption+10).eq."SGEID") then

         call set_random_seed("JOB_ID")
         write(6,'(1x,"Random seed for RAN2NR (from SGE variable JOB_ID): ",i6/)') iseed

      elseif (options(ioption+6:ioption+10).eq."CLOCK") then

         !-- use clock to generate random seed
         call set_random_seed()
         write(6,'(1x,"Random seed for RAN2NR (from current clock value): ",i6/)') iseed

      else

         iseed_inp = reada(options,ioption+6)
         call set_random_seed(iseed_inp)
         write(6,'(1x,"Random seed for RAN2NR (from input): ",i6/)') iseed

      endif

   else

      !-- use clock to generate random seed (default)
      call set_random_seed()
      write(6,'(1x,"Random seed for RAN2NR was generated based on the clock (by default): ",i6/)') iseed

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! reset ran2nr() random sequence when starting each trajectory
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' RESET_RANDOM').ne.0) then
      reset_random = .true.
      write(6,'(/1x,"Reset ran2nr() random sequence when starting each new trajectory."/)')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Random seeds for DUNI
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," DUNISEEDS=")
   ioption2 = ioption + 11

   if (ioption.ne.0) then

      if (options(ioption2:ioption2+4).eq."CLOCK") then

         !-- set seeds using current time
         call set_duni_random_seeds()
         write(6,'(/1x,"Random seeds for DUNI generated from the clock:")')

      else

         !-- set seeds manually (from input)

         islash1 = index(options(ioption2:),'/')
         islash2 = index(options(ioption2+islash1:),'/')
         islash3 = index(options(ioption2+islash1+islash2:),'/')
         idum1 = reada(options(ioption2:ioption2+islash1-2),1)
         idum2 = reada(options(ioption2+islash1:ioption2+islash1+islash2-1),1)
         idum3 = reada(options(ioption2+islash1+islash2:ioption2+islash1+islash2+islash3-1),1)
         idum4 = reada(options,ioption2+islash1+islash2+islash3)
         call set_duni_random_seeds(idum1,idum2,idum3,idum4)
         write(6,'(/1x,"Random seeds for DUNI read from input:")')

      endif

   else

      write(6,'(/1x,"Random seeds for DUNI have fixed default values:")')

   endif

   write(6,'(/1x,"Random seeds for DUNI random number generator:")')
   write(6,'( 1x,"i_seed = ",i6)')  i_seed
   write(6,'( 1x,"j_seed = ",i6)')  j_seed
   write(6,'( 1x,"k_seed = ",i6)')  k_seed
   write(6,'( 1x,"l_seed = ",i6/)') l_seed


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! The solvent coordinate frame used for initial values
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ioption = index(options,"SCALED")
   if (ioption.ne.0) then
      scaled = .true.
   else
      scaled = .false.
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initial solvent coordinate
   ! (center of the initial distribution)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption  = index(options,' ZE0=')
   ioptionve = index(options,' VZE0=')
   ioptionv1 = index(options,' VZ10=')
   ioption2 = index(options,' ZE0R')
   ioption3 = index(options,' ZE0P')

   if (ioptionve.ne.0.and.ioptionv1.ne.0) then
      write(*,'(/1x,"*** (in DYNAMICSET2): You can specify either of the VZE0= or VZ10=, not both ***"/)')
      call clean_exit
   endif

   fixedinit = index(options,' FIXEDINIT').ne.0

   if (ioption.ne.0) then

      ioption = ioption + 5
      ze0 = reada(options,ioption)
      if (scaled) then
         z10 = ze0
         call z1_to_ze(z10,ze0)
      endif

   elseif (ioption2.ne.0) then

      !-- reactant minimum
      ze0 = z_1

   elseif (ioption3.ne.0) then

      !-- product minimum
      ze0 = z_1 - 2.d0*lambda

   else

      write(*,'(/1x,"*** (in DYNAMICSET2): You MUST specify one of the ZE0=, ZE0R, ZE0P options for DYNAMICSET2 keyword ***"/)')
      call clean_exit

   endif


   if (ioptionve.ne.0) then
      ioptionve = ioptionve + 6
      vze0 = reada(options,ioptionve)
      call ve_to_v1(vze0,vz10)
   endif

   if (ioptionv1.ne.0) then
      ioptionv1 = ioptionv1 + 6
      vz10 = reada(options,ioptionve)
      call v1_to_ve(vz10,vze0)
   endif

   if (reactive_flux) then
      if (scaled) then 
         Call z1_to_ze(z10,ze0)
         dividing_surface_egap = ze0
      else
         dividing_surface_egap = ze0
      end if
   end if

   call ze_to_z1(ze0,z10)

   if (.not.reactive_flux) then 
      write(6,'(/1x,"Center of the initial distribution of the solvent coordinate:",/,&
      &          1x,"Ze(0) = ",F10.3,2X,A,/,&
      &          1x,"z1(0) = ",F10.3,2X,A)') ze0, zgapdim, z10, zscadim
   else
      write(6,'(/1x,"Initial solvent coordinate for each reactive flux trajectory:",/,&
      &          1x,"Ze(0) = ",F10.3,2X,A,/,&
      &          1x,"z1(0) = ",F10.3,2X,A)') ze0, zgapdim, z10, zscadim
   end if

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate dynamics (not implemented yet)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !-- set control variables in the propagators module
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call set_mode_et2(mode_dyn,nstates_dyn,ielst_dyn,iset_dyn)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !-- Allocate arrays in propagators module
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call allocate_electronic_states
   if (mdqt) then
      call allocate_mdqt_arrays
      if (afssh) call allocate_afssh_arrays
   endif
   if (weights) call allocate_evb_weights

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !-- nature of the initial density matrix
   !   (meaningful only for MDQT dynamics)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (index(options,' MIXEDSTATEDENSITY').ne.0) purestate = .false.
   if (index(options,' PURESTATEDENSITY').ne.0) purestate = .true.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Type of the initial condition:
   !
   ! ISTATE - initial state is a pure electronic state
   !          specified by the ISTATE keyword
   !
   ! DSTATE - initial state is a diabatic electronic state
   !          specified by the DSTATE keyword
   !
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   initial_state_pure = index(options,' ISTATE=') .ne. 0
   initial_state_diab = index(options,' DSTATE=') .ne. 0

   !-- These two options should be mutually exclusive.

   icount = 0
   if (initial_state_pure)  icount = icount + 1
   if (initial_state_diab)  icount = icount + 1

   if (icount.eq.0) then

      !-- none of the relevant keywords has been specified: default mode "pure" and ISTATE=1

      initial_state_pure = .true.
      initial_state = 1
      write(6,'(1x,"(Default initial condition) At t=0: Pure ground adiabatic state.")')

   elseif (icount.gt.1) then

      !-- more than one keyword has been specified: ambiguity => stop the program

      write(6,'(1x,"Ambiguious input in DYNAMICSET2: only ONE of (ISTATE, DSTATE) keywords must be specified.)")')
      call clean_exit

   endif

   if (reactive_flux) then
      write(6,'(1x,"Ignoring any ISTATE and DSTATE values specified for the reactive_flux calculation.)")')
      initial_state_diab = .false.
      initial_state_pure = .false.
   end if

   if (initial_state_pure) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !  Initial state is a pure adiabatic state
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (initial_state.lt.0) then
         ioption2 = index(options,' ISTATE=')
         write(6,'(/1x,"(Initial state is a pure adiabatic electronic state)")')
      endif

      initial_state = reada(options,ioption2+8)
      if (initial_state.eq.0) initial_state = 1

      if (mode_dyn.eq.'ADIAB') then

         if (initial_state.le.0) initial_state = 1
         write(6,'(1x,"At t=0: initial adiabatic electronic state: ",i6)') initial_state

      elseif (mode_dyn.eq.'DIAB4') then

         write(6,'(1x,"At t=0: initial diabatic electronic state: ",i6)') initial_state

      endif

   elseif (initial_state_diab) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !  Initial state is a diabatic state (coherent mixture of adiabatic states)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (initial_state.lt.0) then
         ioption2 = index(options,' DSTATE=')
         write(6,'(/1x,"(Initial state is a diabatic electronic state - coherent mixture of adiabatic states)")')
      endif

      initial_state = reada(options,ioption2+8)
      if (initial_state.eq.0) initial_state = 1
      write(6,'(1x,"At t=0: initial diabatic electronic state: ",i6)') initial_state

   elseif (reactive_flux) then 

      write(6,'(1x,"Initial state will be chosen at the beginning of each trajectory for reactive_flux")')

   else

      write(6,'(1x,"No initial state was chosen.  WARNING: You might need to specify ISTATE or DSTATE.")')

   endif


   !-- Allocate arrays for reactive_flux if that is the calculation type

   if (reactive_flux) then

      allocate(z1_storage(1:nsteps))
      allocate(vz1_storage(1:nsteps))
      allocate(vz1_storage_beforehop(1:nsteps))
      allocate(fnj_storage(1:nsteps,1:nstates_dyn,1:nstates_dyn))
      allocate(switch_attempt(1:nsteps))
      allocate(occupied_adiabat(1:nsteps))
      allocate(occupied_adiabat_beforehop(1:nsteps))
      allocate(switch_attempt_state(1:nsteps))

      if (.not.mdqt) then
         write(6,*) 'ERROR: REACTIVE_FLUX can only be turned on when using MDQT.'
         Call clean_exit
      endif

      if (afssh.or.collapse_region_coupling.or.ids.or.ida.or.gedc) then
         write(6,*) 'ERROR: Decoherence (either AFSSH, collapse_coupling_region, IDS, IDA, GEDC) cannot be used with reactive flux algorithm.'
         Call clean_exit
      endif

      !-- Make sure the solvent model is not ONODERA2 for reactive flux (memory of trajectory is unphysical for reverse trajectory)
      if (solvent_model.ne."ONODERA") then
         write(6,*) 'ERROR: Reactive flux can only be performed with the Onodera-1 model.'
         write(6,*) 'Please select this solvent model.'
         call clean_exit
      endif

   endif

   !===DEBUG===
   !call print_propagators_et2


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output files for trajectories
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' TRAJOUT=')

   if (ioutput.eq.0) then

      fname = job(1:ljob)//'/trajectory'
      lenf = ljob + 11

   else

      ispa = index(options(ioutput+9:),space)
      fname = options(ioutput+9:ioutput+ispa+7)
      lenf = ispa - 1
      call locase(fname,lenf)
      fname = job(1:ljob)//'/'//fname(1:lenf)
      lenf = lenf + ljob + 1

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Append "adiab" or "diab" to the trajectory file names
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (adiab) then
      fname = trim(fname)//"_adiab"
   elseif (diab4) then
      fname = trim(fname)//"_diaba_"//iset_char_diab4(iset_dyn)
   endif

   write(6,'(/1x,"Trajectory data are written to the file(s) <",a,">"/)') trim(fname)//"_<xxxxxx>.dat"

   !======================================!
   !      MAIN LOOP OVER TRAJECTORIES     !
   !======================================!

   zeit_start = second()

   number_of_skipped_trajectories = 0
   number_of_failed_trajectories = 0
   itraj_start = 1
   ntraj_valid = ntraj

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Restart trajectories
   ! (restore all random seeds from the checkpoint file)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," RESTART")

   if (ioption.ne.0) then

      if (options(ioption+8:ioption+8).eq."=") then
         ispa = index(options(ioption+9:),space)
         fname = options(ioption+9:ioption+ispa+7)
         call locase(fname,ispa-1)
      else
         fname = job(1:ljob)//".dchk"
      endif

      !-- open checkpoint file

      open(unit=1,file=trim(fname),action="read",form="unformatted",status="old")

      !-- read last trajectory number etc.
      read(1) itraj_start, number_of_skipped_trajectories, number_of_failed_trajectories
      itraj_start = itraj_start + 1
      write(6,'(/1x,"Restarting from trajectory: ",i6)') itraj_start

      !-- restore random seeds
      call restore_random_seeds(1)

      !-- close checkpoint file
      close(1)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !-- initialize DUNI() random number generator
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call initialize_duni()

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !-- start (restart) loop over trajectories
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   loop_over_trajectories: do itraj=itraj_start,ntraj

      !-- initialize the suffix of the output file

      write(traj_suffix,'(i6.6)') itraj

      write(6,'(/1x,60("#"))')
      write(6,'( 1x,"Trajectory #",i10)') itraj
      write(6,'( 1x,60("#")/)')

      !-- pick the initial value of the solvent coordinate
      !   from a gaussian distribution centered at (z10)
      !   or use a fixed value from the input (z10)

      !-- CAS adding this to set initial energy gap to value chosen in input file (no distribution around center)
      if (reactive_flux.or.fixedinit) then
         z1 = z10
      else
         sigma = sqrt(kb*temp/f0)
         sample = gaussdist_boxmuller()
         z1 = z10 + sigma*sample
      endif

      !-- Reactive flux cannot be done for Onodera2 so the y1
      !   (auxiliary solvent coordinate) value will not be modified.

      if (solvent_model.eq."ONODERA2") then
         !-- pick the initial value of auxiliary solvent coordinate
         !   from a gaussian distribution centered at (-z10)
         sigma = sqrt(kb*temp/gamma)
         sample = gaussdist_boxmuller()
         if (fixedinit) then
            y1 = -z10
         else
            y1 = -z10 + sigma*sample
         endif
      endif

      !-- zero out the moments (A-FSSH)
      if (mdqt.and.afssh) then
         call reset_zmoments
         call reset_pzmoments
      endif

      !-- initialize initial velocity

      if (solvent_model.eq."DEBYE".or.solvent_model.eq."DEBYE2") then

         !-- overdamped dynamics
         vz1 = 0.d0

      elseif (solvent_model.eq."ONODERA".or.solvent_model.eq."ONODERA2".or.solvent_model.eq."NEWTON") then

         !-- initial velocity from Maxwell distribution
         sigma1 = sqrt(kb*temp/effmass1)
         sample = gaussdist_boxmuller()
         vz1 = sigma1*sample

      else

         !-- Other models are not implemented yet...
         write(6,'(1x,"(From DYNAMICSET2: solvent model ",a10," is not implemented: abort calculation)")') solvent_model
         call clean_exit

      endif

      if (fixedinit) vz1 = vz10

      !-- For reactive flux, check if velocity is negative.  If it is, reverse it.  
      if (reactive_flux) then
         if (vz1.lt.0.d0) then
            vz1 = -1.d0*vz1
         endif
         write(*,*) 'Initial (scaled) velocity for reactive flux trajectory:', vz1
      endif

      !-- calculate electronic state for the initial value of solvent coordinate
      !   at t=0 (very first time for this trajectory)
      call calculate_electronic_states(z1)

      if (collapse_region_coupling) then
         interaction_region = interaction_region_check()
         interaction_region_prev = interaction_region
      endif


      !-- open the trajectory output file (channel 1)
      open(itraj_channel,file=trim(fname)//"_"//traj_suffix//".dat")

      !-- write the header of the trajectory file

      if (weights) then
         write(itraj_channel,'("#",191("="))')
      else
         write(itraj_channel,'("#",131("="))')
      endif

      if (reactive_flux) then
         write(itraj_channel,'("#   Data for the reactive flux trajectory #",i6.6)') itraj
      else
         write(itraj_channel,'("#   Data for the trajectory #",i6.6)') itraj
      endif

      if (weights) then
         write(itraj_channel,'("#",191("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"VZe",t68,"Ekin",t80,"Efe",t89,"occ.",t100,"EVB weights (1,2)",t127,"diabatic populations")')
         write(itraj_channel,'("#",191("-"))')
      else
         write(itraj_channel,'("#",131("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"vze",t68,"Ekin",t80,"Efe",t89,"occ.")')
         write(itraj_channel,'("#",131("-"))')
      endif



      !-- Assign the initial occupied state

      if (initial_state_pure) then

         !-- always the same initial state
         istate = initial_state

         if (istate.gt.nstates_dyn) then
            write(6,'(/1x,"From DYNAMICS3: index of the initial pure state specified (",i2,") ",/,&
            &          1x,"is greater than the number of states included in dynamics (",i2,").",/,&
            &          1x,"Check the ISTATE option!")') istate, nstates_dyn
            call clean_exit
         endif

         if (mdqt) then
            !-- set initial amplitudes (and/or density matrix) at t=0
            !   (always a pure state)
            call set_initial_amplitudes_pure(istate)
            call set_initial_density_pure(istate)
         endif

      elseif (initial_state_diab) then

         !-- The initial state is sampled according to its amplitude in the coherent mixture
         istate = assign_initial_state(initial_state)

         if (istate.gt.nstates_dyn) then
            write(6,'(/1x,"From DYNAMICS3: index of the sampled initial state specified (",i2,") ",/,&
            &          1x,"is greater than the number of states included in dynamics (",i2,").",/,&
            &          1x,"Check the input option!")') istate, nstates_dyn
            call clean_exit
         endif

         if (mdqt) then
            !-- set initial amplitudes (and/or density matrix) at t=0
            !   (coherent mixture of states or pure state with the largest amplitude)
            if (purestate) then
               call set_initial_amplitudes_pure(istate)
               call set_initial_density_pure(istate)
            else
               call set_initial_amplitudes_mixture(initial_state)
               call set_initial_density_mixture(initial_state)
            endif
         endif

      elseif (reactive_flux) then 

         !-- For reactive_flux, we need to do the backwards trajectory first before we can set the initial amplitudes.
         !   Choose the initial adiabat using the Maxwell-Boltzmann distribution
         call reactiveflux_choose_initial_state(nstates_dyn,temp,istate,initial_state)

         !**********************************for single adiabat propagation only
         !istate = 1
         !initial_state = 1

         !-- Set initial values of W, alpha, Fn, Fd for each trajectory.  Also zero all arrays for each trajectory.

         call initialize_reactiveflux_variables(successfulreverse,normalforward,alpha_limit,W,alpha,Fn,Fd,fnj_storage,&
              &occupied_adiabat,occupied_adiabat_beforehop,switch_attempt,switch_attempt_state,vz1_storage,z1_storage,&
              &vz1_storage_beforehop,nsteps,nstates_dyn)
         call set_initial_amplitudes_pure(istate)
         call set_initial_density_pure(istate)

         call reverse_time_propagation(z1_storage,vz1_storage,fnj_storage,&
               &dg_reaction,switch_attempt,nstates_dyn,occupied_adiabat,switch_attempt_state,final_step_reverse,&
               &vz1,z1,ekin,istate,successfulreverse,vz1_storage_beforehop,occupied_adiabat_beforehop,alpha_limit,&
               &normalforward,dividing_surface_egap)

         if (.not.successfulreverse) then
           !The trajectory did not make it to the reactant or product region.  Cycle for the next trajectory.
           cycle
         endif

         istate = occupied_adiabat_beforehop(final_step_reverse)

         !-- This z1 and vz1 will be stored in z1_prev and vz1_prev upon entering forward_time_propgagtion_quantum_only routine.
         !   This allows for interpolation and integration of the TDSE.

         z1 = z1_storage(final_step_reverse)
         vz1 =-1.d0*vz1_storage_beforehop(final_step_reverse)

         ekin = 0.5d0*f0*tau0*taul*vz1*vz1

         call set_initial_amplitudes_pure(istate)
         call set_initial_density_pure(istate)

         if (weights) then
            call get_evb_weights
            call get_diabatic_populations
            call get_diabatic_populations_lfs(istate)
         endif

         call forward_time_propagation_quantum_only(z1_storage,vz1_storage,switch_attempt,occupied_adiabat,switch_attempt_state,&
                &z1,vz1,itraj,fnj_storage,nstates_dyn,final_step_reverse,alpha,dividing_surface_egap,ekin,&
                &W,vz1_storage_beforehop, occupied_adiabat_beforehop,istate,itraj_channel,ielst_dyn,weights)
         call print_amplitudes(6)

         if (.not.normalforward) then
            write(6,*) 'This trajectory went to products (in the reverse trajectory).'
            write(6,*) 'Cycle to next trajectory.'
            call v1_to_ve(vz1_storage_beforehop(1),vze)
            write(Fdfile,*) "Fd ", vze*W
            write(Fnfile,*) "Fn ", 0.d0
            cycle  !did not make it to reactant in reverse traj, don't do normal traj
         endif

         if (normalforward.and.(W.eq.0.d0)) then
            write(6,*) 'Cycling for next trajectory because Fn and Fd will be 0 (W=0)'
            write(Fdfile,*) 'Fd ', W
            write(Fnfile,*) 'Fn ', W
            cycle
         endif

         if (.not.normalforward) then
            write(6,*) 'This trajectory went to products (in the reverse trajectory).  Stop here.'
            call v1_to_ve(vz1_storage_beforehop(1),vze)
            write(6,*) 'Value of W', W
            write(6,*) "Fd ", vze*W
            write(6,*) "Fn ", 0.d0
            cycle  ! did not make it to reactant in reverse traj, don't do normal traj
         endif

         if (normalforward.and.(W.eq.0.d0)) then
            write(6,*) 'normalforward is true but value of W is', W
            write(6,*) 'Cycling for next trajectory because Fn and Fd will be 0'
            write(6,*) 'Fd ', W
            write(6,*) 'Fn ', W
            cycle
         endif

      else

         write(6,'(/1x,"From DYNAMICS3: For unknown reason (VERY SERIOUS BUG?) no initial condition was chosen. Abort.")')
         call clean_exit

      endif


      if (mdqt) then

         !-- calculate electronic states at t=0 (very first time for this trajectory)
         call calculate_electronic_states(z1)

         !-- print out the initial amplitudes of the time-dependent wavefunction
         if (.not.reactive_flux) then
            call print_initial_amplitudes(6)
         endif

         !-- calculate force matrices (A-FSSH specific)
         if (afssh) call calculate_force_matrices(z1)

      endif


      if (weights) then
         call get_evb_weights
         call get_diabatic_populations
         call get_diabatic_populations_lfs(istate)
      endif


      write(6,'(/1x,"===> Trajectory ",i5," starts on the electronic state ",i3)') itraj, istate
      write(6,'( 1x,"===> Initial solvent coordinate (z1), (kcal/mol)^(1/2):    ",f13.6)') z1
      write(6,'( 1x,"===> Initial solvent velocity  (vz1), (kcal/mol)^(1/2)/ps: ",f13.6)') vz1
      call z1_to_ze(z1,ze)

      write(6,*)
      write(6,'(191("-"))')
      write(6,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"vze",t68,"Ekin",t80,"Efe",t89,"state",t100,"diab populations (evb,wavefun,LFSdensity)")')
      write(6,'(191("-"))')

      !write(6,'(137x,$)')

      !--(DEBUG)--start
      if (ndump777.ne.0) then
         open(777,file=job(1:ljob)//"/populations_"//traj_suffix//".dat")
         open(778,file=job(1:ljob)//"/coherences_"//traj_suffix//".dat")
      endif

      if (ndump888.ne.0) then
         open(888,file=job(1:ljob)//"/couplings_"//traj_suffix//".dat")
         write(888,'("#",132("-"))')
         write(888,'("#",t6,"t(ps)",t22,"|d(1,2)|",t42,"d12",t62,"v*d12",t82,"U(ground)",t102,"U(excited)",t120,"Adiabatic gap")')
         write(888,'("#",132("-"))')
      endif

      if (ndump999.ne.0) then
         do i=1,nstates_dyn
            write(state_index,'(i1)') i
            open(999+i+1,file=job(1:ljob)//"/afssh-moments_state"//state_index//"_traj"//traj_suffix//".dat")
            write(999+i+1,'("#-State ",i1,124("-"))') i
            write(999+i+1,'("#",t6,"t(ps)",t22,"z-moment(Re,Im)",t62,"pz-moment(Re,Im)")')
            write(999+i+1,'("#",132("-"))')
         enddo
      endif
      !--(DEBUG)--end

      !-- transform initial values at time t=0
      call z1_to_ze(z1,ze)
      call v1_to_ve(vz1,vze)

      !-- initial free energy (PMF)
      efes = get_free_energy(istate)

      !-- initial kinetic energy
      ekin1 = half*f0*tau0*taul*vz1*vz1
      ekin = ekin1

      if (weights) then
          write(itraj_channel,'(f13.6,6f12.5,i5,12f15.9)') &
          & 0.d0, z1, vz1, ze, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn), (diabatic_population(k),k=1,ielst_dyn), (diabatic_population_lfs(k),k=1,ielst_dyn)
      else
          write(itraj_channel,'(f13.6,6f12.5,i5)') &
          & 0.d0, z1, vz1, ze, vze, ekin, efes, istate
      endif

      number_of_switches = 0
      number_of_rejected = 0
      number_of_reversals = 0

      traj_time_start = second()

      !===============================!
      !   MAIN LOOP OVER TIME STEPS   !
      !===============================!
      loop_over_time: do istep=1,nsteps

         switch = .false.

         zeit_prev = real(istep-1,kind=8)*tstep
         zeit = real(istep,kind=8)*tstep

         !-- MDQT: store couplings, electronic energies, and velocities
         !         from the previous step (for interpolation)

         z1_prev = z1
         if (mdqt) then
            vz1_prev = vz1
            ekin_prev = ekin
            call store_nonadiabatic_couplings             !  coupz(:,:) -> coupz_prev(:,:)
            if (collapse_region_coupling) interaction_region_prev = interaction_region
            call store_electronic_energies                !  fe(:)      -> fe_prev(:)
            call store_wavefunctions                      !  z(:,:)     -> z_prev(:,:)
            if (afssh) call store_force_matrices    !  fmatz(:,:) -> fmatz_prev(:,:)
         endif

         !-- Propagate solvent coordinates and velocities

         if (solvent_model.eq."NEWTON") then

            !-- Newton equation (ballistic model)
            call velocity_verlet_1d(istate,z1,vz1,tstep,ekin1,efes)
            ekin = ekin1

         elseif (solvent_model.eq."DEBYE") then

            !-- overdamped Langevin equation (pure Debye model)
            call langevin_debye_1d(istate,z1,vz1,tstep,temp,ekin1,efes)
            ekin = ekin1

         elseif (solvent_model.eq."DEBYE2") then

            !-- overdamped Langevin equation with memory friction
            !   (Debye model with two relaxation periods)
            call langevin_debye2_1d(istate,z1,y1,vz1,tstep,temp,ekin1,efes)
            ekin = ekin1
            !write(*,'(/1x,"DYNAMICS3: Debye2 propagator is not coded yet...")')
            !call clean_exit

         elseif (solvent_model.eq."ONODERA") then

            !-- ordinary Langevin equation (Onodera model)
            call langevin_onodera_1d(istate,z1,vz1,tstep,temp,ekin1,efes,reactive_flux,normalforward,istep)
            ekin = ekin1

         elseif (solvent_model.eq."ONODERA2") then

            !-- ordinary Langevin equation with memory friction
            !   (Onodera model with two relaxation periods)
            call langevin_onodera2_1d(istate,z1,y1,vz1,tstep,temp,ekin1,ekinhalf1,efes)
            ekin = ekin1

         endif

         !----------------!
         !-- MDQT stage --!
         !----------------!

         if (mdqt) then

            !-- calculate electronic states and nonadiabatic couplings at t+dt
            !-----------------------------------------------------------------
            call calculate_electronic_states(z1)

            if (collapse_region_coupling) then
               !-- are we still in the interaction region?
               interaction_region = interaction_region_check()
            endif

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at t and t+dt
            !-------------------------------------------------------
            call calculate_v_dot_d(vz1,vz1_prev)

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at half timestep for quadratic interpolation scheme
            !-------------------------------------------------------
            if (interpolation.eq."QUADRATIC") then
               call calculate_v_dot_d_mid(tstep)
            endif

            if (afssh) then

               !-- calculate force matrices at t+dt (A-FSSH)
               !------------------------------------------------------
               call calculate_force_matrices(z1)

               !-- calculate interpolation coefficients for the force matrices
               !------------------------------------------------------------------
               call interpolate_force_matrices(zeit_prev,zeit)

            endif

            !-- calculate interpolation coefficients for the kinetic energy
            !   for phase-corrected surface hopping scheme
            !------------------------------------------------------------------
            if (phase_corr) then
               call interpolate_kinenergy(interpolation,zeit_prev,zeit,ekin,ekin_prev,ekinhalf1)
            endif

            !-- calculate interpolation coefficients for the adiabatic energies
            !------------------------------------------------------------------
            call interpolate_energy(zeit_prev,zeit)

            !-- calculate interpolation coefficients
            !   for the nonadiabatic coupling terms v*d_{kl}
            !-----------------------------------------------
            call interpolate_vdotd(interpolation,zeit_prev,zeit)

            !-- calculate the population of the current state at time t_prev
            !---------------------------------------------------------------
            population_current = calculate_population(istate)

            call reset_switch_prob

            !--(DEBUG)--start
            !if (istep.eq.19419) call print_propagators_et2
            !--(DEBUG)--end

            !-- propagate the amplitudes (or density matrix)
            !   and switching probabilities from t_prev to t
            !------------------------------------------------------------------

            nqsteps_var = nqsteps
            qtstep_var = qtstep
            call save_amplitudes
            if (afssh) then
               call save_density_matrix
               call save_moments
            endif

            24 continue
            call reset_switch_prob

            do iqstep=1,nqsteps_var

               zeitq_prev = (iqstep-1)*qtstep_var + zeit_prev
               zeitq = iqstep*qtstep_var + zeit_prev

               if (afssh) then

                  !-- propagate moments and amplitudes forward in time
                  !-------------------------------------------------------
                  if (phase_corr) then
                     call propagate_moments_and_amplitudes_phcorr(istate,zeitq_prev,qtstep_var)
                  else
                     call propagate_moments_and_amplitudes(istate,zeitq_prev,qtstep_var)
                  endif

               else

                  !-- propagate amplitudes forward in time (AFSSH)
                  !-------------------------------------------------------
                  if (phase_corr) then
                     call propagate_amplitudes_phcorr_rk4(istate,zeitq_prev,qtstep_var)
                  else
                     call propagate_amplitudes_rk4(zeitq_prev,qtstep_var)
                  endif

               endif

               call calculate_density_matrix

               !-- calculate transition probabilities from current state
               !--------------------------------------------------------
               call calculate_bprob_amp(istate,zeitq)

               !-- accumulate swithing probabilities (array operation)
               !------------------------------------------------------
               call accumulate_switch_prob(qtstep_var)

            enddo

            !-- check the norm of the time-dependent wavefunction
            !-----------------------------------------------------
            wf_norm = tdwf_norm()

            if (abs(wf_norm-1.d0).gt.1.d-3) then

               write(*,'(/1x,"-------------------------------------------------------------------------------")')
               write(*,'( 1x,"DYNAMICS3: Amplitudes are not normalized after timestep ",i6)') istep
               write(*,'( 1x,"           Norm of the time-dependent wavefunction:     ",g20.10)') wf_norm
               write(*,'(/1x,"-------------------------------------------------------------------------------")')

               if (nqsteps_var.lt.maxnqsteps) then

                  !-- reduce the TDSE timestep by ten times and repeat the quantum propagation

                  nqsteps_var = nqsteps_var*10
                  qtstep_var = qtstep_var/10.d0

                  write(*,'( 1x,"Number of quantum timesteps is increased ten times to ",i6)') nqsteps_var
                  write(*,'( 1x,"and the quantum propagation will be repeated with a 10 times smaller timestep.")')
                  write(*,'(/1x,"-------------------------------------------------------------------------------")')

                  call restore_amplitudes
                  if (afssh) then
                     call restore_density_matrix
                     call restore_moments
                  endif
                  goto 24

               else

                  !-- discard the failed trajectory

                  write(*,'( 1x,"--- The trajectory ",i6," will be discarded-------------------------------"/)') itraj

                  number_of_failed_trajectories = number_of_failed_trajectories + 1

                  if (weights) then
                     write(itraj_channel,'("#",191("-"))')
                  else
                     write(itraj_channel,'("#",131("-"))')
                  endif

                  write(itraj_channel,'("# Amplitudes are not normalized after timestep ",i6)') istep
                  write(itraj_channel,'("# Norm of the time-dependent wavefunction:     ",g20.10)') wf_norm
                  write(itraj_channel,'("# This trajectory has failed... Even after several tries with smaller TDSE timesteps.")')

                  if (weights) then
                     write(itraj_channel,'("#",191("-"))')
                  else
                     write(itraj_channel,'("#",131("-"))')
                  endif

                  close(itraj_channel)
                  call system("mv "//trim(fname)//"_"//traj_suffix//".dat"//" "//trim(fname)//"_"//traj_suffix//".dat_failed")

                  cycle loop_over_trajectories

               endif

            endif

            !-----------------------------------------------------------------
            !-- check the traces of the matrices of moments (A-FSSH algorithm)
            !   (should be zero?)
            !-----------------------------------------------------------------
            !if (afssh) then
            !   zmom1_norm = zmom1_trace()
            !   pzmom1_norm = pzmom1_trace()
            !   if (abs(zmom1_norm).gt.1.d-6) then
            !      write(*,'(/1x,"-------------------------------------------------------------------------------")')
            !      write(*,'( 1x,"DYNAMICSET2(A-FSSH): Trace(zmom1) is not zero at timestep ",i6)') istep
            !      write(*,'( 1x,"                     The value is: ",g20.10)') zmom1_norm
            !      write(*,'(/1x,"-------------------------------------------------------------------------------")')
            !   endif
            !   if (abs(pzmom1_norm).gt.1.d-6) then
            !      write(*,'(/1x,"-------------------------------------------------------------------------------")')
            !      write(*,'( 1x,"DYNAMICSET2(A-FSSH): Trace(pzmom1) is not zero at timestep ",i6)') istep
            !      write(*,'( 1x,"                     The value is: ",g20.10)') pzmom1_norm
            !      write(*,'(/1x,"-------------------------------------------------------------------------------")')
            !   endif
            !endif
            !-----------------------------------------------------------------

            !--(DEBUG)--start

            !-- print out the populations and coherences (channels 777 and 778)
            if (ndump777.ne.0.and.mod(istep,ndump777).eq.0) then
               if (afssh) then
                  call print_populations_den(777,zeit)
                  call print_coherences_den(778,zeit,istate)
               else
                  call print_populations_amp(777,zeit)
                  call print_coherences_amp(778,zeit,istate)
               endif
            endif

            !-- print out couplings (channel 888)
            if (ndump888.ne.0.and.mod(istep,ndump888).eq.0) then
               call print_couplings_and_splittings(888,zeit)
            endif

            !-- print out afssh moments (channel 999)
            if (afssh) then
               if (ndump999.ne.0.and.mod(istep,ndump999).eq.0) then
                  do i=1,nstates_dyn
                     call print_afssh_moments(999+i+1,zeit,i)
                  enddo
               endif
            endif

            !-- hook for the debugger (for setting up a breakpoint)
            !if (istep.eq.16777216) then
            !   write(*,*)
            !   write(*,*) "===========BREAKPOINT: Step number: ",istep
            !   write(*,*)
            !endif

            !--(DEBUG)--end


            !-- Normalize swithing probabilities by the current state population
            !   and zero out the negative ones
            !-------------------------------------------------------------------
            call normalize_switch_prob(population_current)

            !-- construct the full density matrix (do we really need it?)
            !------------------------------------------------------------
            call calculate_density_matrix

            !-- decision time: should we make a hop?
            !-------------------------------------------
            new_state = switch_state(istate)
            switch = new_state.ne.istate

            if (switch) then

               !-- attempt adjusting velocities (and possibly moments for A-FSSH)

               if (afssh) then
                  if (.not.along_moments) then
                     call adjust_velocities_and_moments(istate,new_state,vz1,success)
                  else
                     call adjust_velocities_and_moments_0(istate,new_state,vz1,success)
                  endif
               else
                  call adjust_velocities(istate,new_state,vz1,success)
               endif

               if (success) then

                  write(itraj_channel,'("#--------------------------------------------------------------------")')
                  write(itraj_channel,'("#  t  = ",f13.6," ps ==> SWITCH ",i3,"  -->",i3)') zeit,istate,new_state
                  write(itraj_channel,'("#  d  = ",f20.6)') get_nonadiabatic_coupling(istate,new_state)

                  istate = new_state
                  number_of_switches = number_of_switches + 1

                  !-- Instantaneous decoherence algorithms (IDS and IDA): collapse the wavefunction after a successful hop
                  if (ids.or.ida) then
                     call collapse_wavefunction(istate)
                     call calculate_density_matrix
                     write(*,'("*** (ID): wavefunction collapsed to pure state ",i2)') istate
                     write(itraj_channel,'("# (ID) Wavefunction collapse to state ",i2," occurred")') istate
                  endif
                  write(itraj_channel,'("#--------------------------------------------------------------------")')

               else

                  write(itraj_channel,'("#--------------------------------------------------------------------")')
                  write(itraj_channel,'("#  t  = ",f13.6," ps ==> REJECTED SWITCH ",i3,"  -->",i3)') zeit,istate,new_state
                  write(itraj_channel,'("#  d  = ",f20.6)') get_nonadiabatic_coupling(istate,new_state)

                  number_of_rejected = number_of_rejected + 1

                  !-- Instantaneous decoherence algorithm (IDA): collapse the wavefunction after a rejected hop
                  if (ida) then
                     call collapse_wavefunction(istate)
                     call calculate_density_matrix
                     write(*,'("*** (ID): wavefunction collapsed to pure state ",i2)') istate
                     write(itraj_channel,'("# (ID) Wavefunction collapse to state ",i2," occurred")') istate
                  endif

                  if (revvel) then

                     if (.not.revvel_all) then

                        !-- calculate both criteria for Truhlar velocity reversal algorithm

                        !--(1) (F1*d12)(F2*d12) < 0
                        !--(2) (P1*d12)(F2*d12) < 0

                        force1 = -get_gradient(istate) - f0*z1
                        force2 = -get_gradient(new_state) - f0*z1
                        coup12 =  get_nonadiabatic_coupling(istate,new_state)

                        revvel_cond1 = (force1*coup12)*(force2*coup12) .lt. 0
                        revvel_cond2 = (vz1*coup12)*(force2*coup12) .lt. 0

                        if (revvel_cond1.and.revvel_cond2) then
                           vz1 = -vz1
                           number_of_reversals = number_of_reversals + 1
                           write(*,'("*** (REVVEL): velocity has been reversed")')
                           write(itraj_channel,'("#  velocity has been reversed")')
                        endif

                     else

                        vz1 = -vz1
                        number_of_reversals = number_of_reversals + 1
                        write(*,'("*** (REVVEL0): velocity has been reversed")')
                        write(itraj_channel,'("#  velocity has been reversed")')

                     endif

                  endif

                  write(itraj_channel,'("#--------------------------------------------------------------------")')

               endif

            endif

            !-- A-FSSH specific part: collapsing events and resetting the moments

            if (afssh) then
               call collapse_and_reset_afssh_erratum(istate,tstep,dzeta)
               call calculate_density_matrix
            endif

            !-- Collapse the wavefunction if leaving the interaction region: "poor man's" decoherence
            if (collapse_region_coupling) then
               if (interaction_region_prev.and.(.not.interaction_region)) then
                  call collapse_wavefunction(istate)
                  call calculate_density_matrix
                  write(*,'("*** Leaving interaction region: wavefunction collapsed to pure state ",i2)') istate
               endif
            endif

            !-- Energy-based decoherence correction: damp the amplitudes using an approximate Granucci's prescription (GEDC)
            if (gedc) then
               call damp_amplitudes_gedc(istate,tstep,ekin,1.d0,0.1d0*au2cal)
               call calculate_density_matrix
            endif


         endif  !mdqt

         !---------------------------!
         !-- end of the MDQT stage --!
         !---------------------------!

         !-- calculate EVB weights
         if (weights) then
            call get_evb_weights
            call get_diabatic_populations
            call get_diabatic_populations_lfs(istate)
         endif

         !-- write the current data to the trajectory file

         call z1_to_ze(z1,ze)
         call v1_to_ve(vz1,vze)
         call z1_to_ze(z1_prev,ze_prev)

         if (reactive_flux) then

            !-- Check if passed through dividing surface toward products
            if ((ze_prev - dividing_surface_egap).ge.0.d0.and.(ze - dividing_surface_egap).le.0.d0) then
               alpha = alpha + 1.d0
            end if

            !-- Check if we are in product region
            !if ((ze.lt.(dg_reaction - lambda)).and.(istate.eq.1))  then
            if ((ze.lt.z_2).and.(istate.eq.1)) then    !z_2 has product energy gap only for solvation terms.  If gas phase bias is zero z_2 = dg_reaction - lambda
               write(6,*) 'Reached products on adiabat 1 in normal forward propagation.'
               write(6,*) 'value of alpha', alpha
               write(6,*) 'value of W', W
               call v1_to_ve(vz1_storage_beforehop(1),vze)
               Fn = vze*W/alpha
               Fd = vze*W
               write(Fnfile,*) 'Fn', Fn
               write(Fdfile,*) 'Fd', Fd
               exit
            end if
 
            !-- Check if we are in the reactant region
            !if ((ze.gt.(lambda + dg_reaction)).and.(istate.eq.1)) then
            if ((ze.gt.z_1).and.(istate.eq.1)) then  !z_1 has the reactant energy for only solvation terms.  If gas phase bias is zero z_1 = lambda + dg_reaction
               write(6,*) 'Reached reactants on adiabat 1 in normal forward propagation.'
               call v1_to_ve(vz1_storage_beforehop(1),vze)
               Fd = vze*W
               Fn = 0.d0
               write(Fnfile,*) 'Fn', Fn
               write(Fdfile,*) 'Fd ', Fd
               exit
            end if

            if (alpha.gt.alpha_limit) then
               !-- W is set to 1, and update only Fd
               call v1_to_ve(vz1_storage_beforehop(1),vze)
               write(6,*) 'Alpha exceeded the upper limit in the normal forward trajectory.  Do not count this trajectory.'
               write(6,*) 'Value of W and initial velocity', W, vze
               !write(*,*) 'Fn ', 0.d0
               !write(*,*) 'Fd ', vze*1.d0
               exit
            end if


            if (istep.eq.nsteps) then
               write(6,*) 'Reached the last time step without reaching products, reactants or exceeding alpha. Disregard this trajectory. Arrays too short'
               call v1_to_ve(vz1_storage_beforehop(1),vze)
               write(6,*) 'Value of W and initial velocity', W, vze
               exit
            end if

         end if  !reactive_flux


         if (mod(istep,ndump).eq.0) then
            if (weights) then
               write(itraj_channel,'(f13.6,6f12.5,i5,12f15.9)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn), (diabatic_population(k),k=1,ielst_dyn), (diabatic_population_lfs(k),k=1,ielst_dyn)
            else
               write(itraj_channel,'(f13.6,6f12.5,i5)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate
            endif
         endif

         !if (ndump6.gt.0.and.mod(istep,ndump6).eq.0) &
         !& write(6,'(137("\b"),f13.6,6f12.5,i5,$)') zeit, z1, vz1, ze, vze, ekin, efes, istate

         if (ndump6.gt.0.and.mod(istep,ndump6).eq.0) &
         & write(6,'(f13.6,6f12.5,i5,12f12.8)') zeit, z1, vz1, ze, vze, ekin, efes, istate, &
         & (wght(k,istate),k=1,ielst_dyn), (diabatic_population(k),k=1,ielst_dyn), (diabatic_population_lfs(k),k=1,ielst_dyn)

      enddo loop_over_time

      traj_time_end = second()


      if (weights) then
         write(itraj_channel,'("#",191("-"))')
      else
         write(itraj_channel,'("#",131("-"))')
      endif
      write(itraj_channel,'("# Number of allowed  switches:  ",i5)') number_of_switches
      write(itraj_channel,'("# Number of rejected switches:  ",i5)') number_of_rejected
      write(itraj_channel,'("# Number of velocity reversals: ",i5)') number_of_reversals
      if (weights) then
         write(itraj_channel,'("#",191("-"))')
      else
         write(itraj_channel,'("#",131("-"))')
      endif
      close(itraj_channel)

      !--(DEBUG)--start
      !-- close files with populations and coherences
      if (ndump777.ne.0) then
         close(777)
         close(778)
      endif
      if (ndump888.ne.0) close(888)
      if (ndump999.ne.0) close(999)
      !--(DEBUG)--end

      write(6,*)
      write(6,'(131("-"))')
      write(6,'("# Total number of allowed  switches: ",i5)') number_of_switches
      write(6,'("# Total number of rejected switches: ",i5)') number_of_rejected
      write(6,'("# CPU time elapsed (sec):            ",f12.3)') traj_time_end - traj_time_start
      write(6,'(131("-"))')
      write(6,*)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !-- save restart data to the binary checkpoint file
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      open(unit=1,file=job(1:ljob)//".dchk",action="write",form="unformatted",status="replace")

      !-- write last trajectory number etc.
      write(1) itraj, number_of_skipped_trajectories, number_of_failed_trajectories

      !-- save random seeds
      call save_random_seeds(1)

      !-- close checkpoint file
      close(1)

      if (reset_random) call reset_random_seed

   enddo loop_over_trajectories

   ntraj_valid = ntraj - number_of_skipped_trajectories - number_of_failed_trajectories

   if (ntraj_valid.le.0) ntraj_valid = 1
   zeit_end = second()
   zeit_total = zeit_end - zeit_start
   write(6,*)
   write(6,'(1x,"================================================================================")')
   write(6,'(1x,"Done. Number of trajectories generated: ",i6)') ntraj_valid
   write(6,'(1x,"      Number of discarded trajectories: ",i6)') number_of_skipped_trajectories
   write(6,'(1x,"      Number of failed trajectories:    ",i6)') number_of_failed_trajectories
   write(6,'(1x,"================================================================================")')
   write(6,'(1x,"Done. Time elapsed         (sec): ",f20.3)') zeit_total
   write(6,'(1x,"      Time per trajectory  (sec): ",f20.3)') zeit_total/ntraj_valid
   write(6,'(1x,"      Time per timestep    (sec): ",f20.3)') zeit_total/(ntraj_valid*nsteps)
   write(6,'(1x,"      Productivity rate (ps/day): ",f20.3)') 3600.d0*24.d0*tstep*ntraj_valid*nsteps/zeit_total
   write(6,'(1x,"================================================================================")')
   write(6,*)

   !-- Deallocate arrays in propagators module

   call deallocate_all_arrays
   if (reactive_flux) then
      close(Fnfile)
      close(Fdfile)
   end if

contains

   subroutine clean_exit
      call deallocate_all_arrays
      stop
   end subroutine clean_exit

end subroutine dynamicset2



!(CAS)&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


   subroutine reverse_time_propagation(z1_storage, vz1_storage, fnj_storage,&
              &dg_reaction,switch_attempt,nstates_dyn,occupied_adiabat,switch_attempt_state,final_step_reverse,vz1,z1,&
              &ekin,istate,successfulreverse,vz1_storage_beforehop,occupied_adiabat_beforehop,alpha_limit,normalforward,dividing_surface_egap)

      use control_dynamics
      use random_generators
      use propagators_et2
      use data_et2
      Implicit NONE

      integer :: nstates_dyn, final_step_reverse
      real(kind=8) :: vz1,z1,ekin
      real(kind=8), dimension(:) :: z1_storage(nsteps), vz1_storage(nsteps), vz1_storage_beforehop(nsteps)
      real(kind=8), dimension(:,:,:) :: fnj_storage(nsteps,nstates_dyn,nstates_dyn)
      integer, dimension(:) :: occupied_adiabat(1:nsteps), occupied_adiabat_beforehop(1:nsteps)
      integer, dimension(:) :: switch_attempt_state(1:nsteps)
      integer, dimension(:) :: switch_attempt(1:nsteps)
      integer :: istep, istate, new_state, number_of_switches, number_of_rejected
      logical :: switch, success, successfulreverse, normalforward
      real(kind=8) :: zeit, zeit_prev
      real(kind=8) :: zeitq, zeitq_prev
      real(kind=8) :: ze, vze, z10, ze0
      real(kind=8) :: vz1_prev, ekin_prev, ekin1, efes, ekinhalf1, dg_reaction, ze_prev, z1_prev
      real(kind=8) :: alpha_limit, dividing_surface_egap
      real(kind=8) :: alpha_temp  !local

      loop_over_time: do istep=1,nsteps

         switch = .false.
         !-- Store energy gap, velocity and current adiabat before checking for a hop
         z1_storage(istep) = z1
         vz1_storage_beforehop(istep) = vz1
         occupied_adiabat_beforehop(istep) = istate

         if (mdqt) then

            !-- calculate electronic states and nonadiabatic couplings at t+dt
            !-----------------------------------------------------------------
            call calculate_electronic_states(z1)

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at t and t+dt
            !-------------------------------------------------------
            call calculate_v_dot_d(vz1,0.d0)
            call reset_switch_prob

            !-- Use f_nj function to estimate switch_prob().
            !   The f_nj values are stored in switch_prob so we can use the switch_state subroutine.
            !---------------------------------------------------------------------------------------
            Call calculate_switch_prob_function(istate,istep,fnj_storage,nsteps,nstates_dyn,tstep)

            !-- decision time: should we make a hop?
            !-------------------------------------------

            new_state = switch_state(istate)
            switch = new_state.ne.istate
            if (switch) then
               !-- record the attempted switch, which is needed to compute w_mu
               !   The following two lines should be uncommented if you want to record each time a switch is attempted (even if it is not successful)
               !switch_attempt(istep) = 1  !1 b/c a switch was attempted, 0 if not
               !switch_attempt_state(istep) = new_state

               !-- attempt adjusting velocities

               call adjust_velocities(istate,new_state,vz1,success)
               if (success) then
                  !-- Uncomment the two following lines if you want to record in switch_attempt only when a hop was actually successful.
                  switch_attempt(istep) = 1
                  switch_attempt_state(istep) = new_state
                  istate = new_state
               endif

            endif

         endif  !mdqt

         !---------------------------!
         !-- end of the MDQT stage --!
         !---------------------------!

         !-- Collect the velocity and occupied adiabat at the end of the hop
         vz1_storage(istep) = vz1
         occupied_adiabat(istep)  = istate

         !-- Propagate solvent coordinates and velocities

         if (solvent_model.eq."DEBYE") then

            !-- overdamped Langevin equation (pure Debye model)
            call langevin_debye_1d(istate,z1,vz1,tstep,temp,ekin1,efes)
            ekin = ekin1

       ! elseif (solvent_model.eq."DEBYE2") then

       !    !-- overdamped Langevin equation with memory friction
       !    !   (Debye model with two relaxation periods)
       !    call langevin_debye2_1d(istate,z1,y1,vz1,tstep,temp,ekin1,efes)
       !    ekin = ekin1

         elseif (solvent_model.eq."ONODERA") then

            !-- ordinary Langevin equation (Onodera model)
            call langevin_onodera_1d(istate,z1,vz1,tstep,temp,ekin1,efes,reactive_flux,normalforward,istep)
            ekin = ekin1

         endif

         !-- Convert scaled energy gaps to unscaled values before determining if
         !   the energy gap is in the reactant or product state

         call z1_to_ze(z1,ze)

         if (mod(istep,10000).eq.0) then
           write(6,*) 'Step and Energy gap:', istep, ze
         endif


         !if ((ze.gt.(lambda + dg_reaction))) then   
         if ((ze.gt.z_1).and.(istate.eq.1)) then !z_1 should have value (lambda + dg_reaction) if the gas phase bias is 0
           !-- z_1 should have the reactant energy gap for only solvation terms
           !   In reactant regime!  We can stop the reverse trajectory.
           write(6,*) 'Final energy gap in reverse trajectory:', ze

           !-- Store the last time step

           final_step_reverse = istep 
           successfulreverse = .true.
           normalforward = .true.
           exit

         end if

         !if ((ze.lt.(dg_reaction - lambda))) then
         if ((ze.lt.z_2).and.(istate.eq.1)) then  !z_2 should have the value of (dg_reaction-lambda) if the gas phase bias is 0
            !-- z_2 should have the product energy gap for only solvation terms
            !   In product regime!  We can stop the reverse trajectory. 
            write(6,*) 'Final energy gap in reverse trajectory:', ze
            final_step_reverse = istep
            normalforward = .false.
            successfulreverse = .true.
            exit
         end if


         if (istep.eq.nsteps) then 
            write(6,'("In the reverse trajectory, the maximum number of time steps was exceeded before reaching the reactant or product.")')
            write(6,'("Disregard this trajectory.  Arrays were not long enough.")')
            !-- Switch to false in order to prevent the forward trajectory from propagating
            successfulreverse = .false.
            normalforward = .false.
         end if

         !if ((ze.lt.(dg_reaction - lambda))) then
         if ((ze.lt.z_2).and.(istate.eq.1)) then  !z_2 should have the value of (dg_reaction-lambda) if the gas phase bias is 0
            !-- z_2 should have the product energy gap for only solvation terms
            !   In product regime!  We can stop the reverse trajectory. 
            write(6,*) 'Final energy gap in reverse trajectory:', ze
            write(6,*) 'istate is', istate
            final_step_reverse = istep
            normalforward = .false.
            successfulreverse = .true.
            exit
         end if

         if (istep.eq.nsteps) then 
            write(6,'("In the reverse trajectory, the maximum number of time steps was exceeded before reaching the reactant or product.")')
            write(6,'("Disregard this trajectory.  Arrays were not long enough.")')
            !-- Switch to false in order to prevent the forward trajectory from propagating
            successfulreverse = .false.
            normalforward = .false.
         end if


      enddo loop_over_time



   end subroutine reverse_time_propagation

!(CAS)&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

   subroutine forward_time_propagation_quantum_only(z1_storage,vz1_storage,switch_attempt,occupied_adiabat,&
              &switch_attempt_state,z1,vz1,itraj,fnj_storage,nstates_dyn,final_step_reverse,alpha,&
              &dividing_surface_egap,ekin,W,vz1_storage_beforehop,occupied_adiabat_beforehop,istate,itraj_channel,ielst_dyn,weights)

      use parsol, only: f0, tau0, taul
      use propagators_et2
      use control_dynamics
      use random_generators
      use data_et2

      Implicit NONE

      integer :: itraj, nstates_dyn, final_step_reverse
      real(kind=8), dimension(:) :: z1_storage(nsteps),  vz1_storage(nsteps), vz1_storage_beforehop(nsteps)
      real(kind=8), dimension(:,:,:) :: fnj_storage(nsteps,nstates_dyn,nstates_dyn)
      integer, dimension(:) :: switch_attempt(nsteps), occupied_adiabat(nsteps), occupied_adiabat_beforehop(nsteps)
      integer, dimension(:) ::  switch_attempt_state(nsteps)
      real(kind=8) :: vz1_prev, ekin_prev, zeit, zeit_prev, qtstep_var, z1_prev, ze_prev
      integer :: istep, istate, nqsteps_var, iqstep, p, newstate, k
      integer, intent(in) :: itraj_channel, ielst_dyn
      logical, intent(in) :: weights
      logical :: switch
      real(kind=8) :: z1,vz1, ekin, ekinhalf1, efes, population_current
      real(kind=8) :: zeitq, zeitq_prev, wf_norm, w_mu, fsum,gsum, W,ze, vze
      integer :: new_state
      real(kind=8) :: dividing_surface_egap, alpha

      !-- Initialize alpha for forward trajectory
      alpha = 0.d0

      !-- LOOP OVER TIME STEPS
      loop_over_time: do istep=1,final_step_reverse-1

         switch = .false.

         zeit_prev = real(istep-1,kind=8)*tstep
         zeit = real(istep,kind=8)*tstep

         !-- MDQT: store couplings, electronic energies, and velocities
         !         from the previous step (for interpolation)

         z1_prev = z1

         if (mdqt) then
            vz1_prev = vz1       ! same as -1.d0*vz1_storage_beforehop(istep+1): want to use the velocity from 
                                 ! before a hop (previous step) in the reverse when interpolating for 
                                 ! the time dependent Schrodinger equation
            ekin_prev = ekin
            call store_nonadiabatic_couplings     !  coupz(:,:) -> coupz_prev(:,:)
            call store_electronic_energies        !  fe(:)      -> fe_prev(:)
            call store_wavefunctions              !  z(:,:)     -> z_prev(:,:)
         endif

         !-- Instead of calling the classical propagation routine, get the position
         !   (z1) and velocity (vz1) coordinates from data stored for the
         !   reverse trajectory

         z1 = z1_storage(final_step_reverse - istep)
         vz1 = -1.d0*vz1_storage(final_step_reverse - istep)
         istate = occupied_adiabat(final_step_reverse - istep)

         ekin_prev = ekin
         ekin = 0.5d0*f0*tau0*taul*vz1*vz1


         !-- Check to see if the trajectory crossed over the dividing surface
         !   (where the energy gap is 0)
         !   Only check this if we are not on the first timestep. If statement prob
         !   not relevant b/c starting the loop one lower than final step for
         !   interpolation purposes anyway.

         call z1_to_ze(z1,ze)
         call z1_to_ze(z1_prev,ze_prev)

         if (((ze_prev-dividing_surface_egap).gt.0.d0).and.((ze-dividing_surface_egap).lt.0.d0)) then
            alpha = alpha + 1.d0
         end if

         !----------------!
         !-- MDQT stage --!
         !----------------!

         if (mdqt) then

            !-- calculate electronic states and nonadiabatic couplings at t+dt
            !-----------------------------------------------------------------
            call calculate_electronic_states(z1)

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at t and t+dt
            !-------------------------------------------------------
            call calculate_v_dot_d(vz1,vz1_prev)

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at half timestep for quadratic interpolation scheme
            !-------------------------------------------------------
            if (interpolation.eq."QUADRATIC") then
               call calculate_v_dot_d_mid(tstep)
            endif

            !-- calculate interpolation coefficients for the kinetic energy
            !   for phase-corrected surface hopping scheme
            !------------------------------------------------------------------
            if (phase_corr) then
               call interpolate_kinenergy(interpolation,zeit_prev,zeit,ekin,ekin_prev,ekinhalf1)
            endif

            !-- calculate interpolation coefficients for the adiabatic energies
            !------------------------------------------------------------------
            call interpolate_energy(zeit_prev,zeit)

            !-- calculate interpolation coefficients
            !   for the nonadiabatic coupling terms v*d_{kl}
            !-----------------------------------------------
            call interpolate_vdotd(interpolation,zeit_prev,zeit)

            !-- calculate the population of the current state at time t_prev
            !---------------------------------------------------------------
            population_current = calculate_population(istate)

            call reset_switch_prob

            nqsteps_var = nqsteps
            qtstep_var = qtstep
            call save_amplitudes

            24 continue
            call reset_switch_prob

            do iqstep=1,nqsteps_var

               zeitq_prev = (iqstep-1)*qtstep_var + zeit_prev
               zeitq = iqstep*qtstep_var + zeit_prev

               !-- propagate amplitudes forward in time 
               !-------------------------------------------------------
               if (phase_corr) then
                  call propagate_amplitudes_phcorr_rk4(istate,zeitq_prev,qtstep_var)
               else
                  call propagate_amplitudes_rk4(zeitq_prev,qtstep_var)
               endif

               call calculate_density_matrix

               !-- calculate transition probabilities from current state
               !--------------------------------------------------------
               call calculate_bprob_amp(istate,zeitq)

               !-- accumulate swithing probabilities (array operation)
               !------------------------------------------------------
               call accumulate_switch_prob(qtstep_var)

            enddo

            !-- check the norm of the time-dependent wavefunction
            !-----------------------------------------------------
            wf_norm = tdwf_norm()

            if (abs(wf_norm-1.d0).gt.1.d-3) then

               write(*,'(/1x,"-------------------------------------------------------------------------------")')
               write(*,'( 1x,"DYNAMICS3: Amplitudes are not normalized after timestep ",i6)') istep
               write(*,'( 1x,"           Norm of the time-dependent wavefunction:     ",g20.10)') wf_norm
               write(*,'(/1x,"-------------------------------------------------------------------------------")')

               if (nqsteps_var.lt.maxnqsteps) then

                  !-- reduce the TDSE timestep by ten times and repeat the quantum propagation

                  nqsteps_var = nqsteps_var*10
                  qtstep_var = qtstep_var/10.d0

                  write(*,'( 1x,"Number of quantum timesteps is increased ten times to ",i6)') nqsteps_var
                  write(*,'( 1x,"and the quantum propagation will be repeated with a 10 times smaller timestep.")')
                  write(*,'(/1x,"-------------------------------------------------------------------------------")')

                  call restore_amplitudes
                  goto 24

               else

                  !-- discard the failed trajectory

                  write(*,'( 1x,"--- The trajectory ",i6," will be discarded-------------------------------"/)') itraj

               endif

            endif


            !-- Normalize swithing probabilities by the current state population
            !   and zero out the negative ones
            !-------------------------------------------------------------------
            call normalize_switch_prob(population_current)

            !-- In the forward_time_propagation we don't decide if we will hop
            !   We just take it from the reverse time trajectory
            !   After quantum transition (or attempt) update the velocity and
            !   state.  If there was not hop, these values will be the same as
            !   before the hop attempt (the "beforehop" array holds the same value as the corresponding array)

            ! istate = occupied_adiabat_beforehop(istep)
            ! vz1 = -1.d0*vz1_storage_beforehop(istep)

            istate = occupied_adiabat_beforehop(final_step_reverse - istep)
            vz1 = -1.d0*vz1_storage_beforehop(final_step_reverse - istep)


            !-- Calculate the w_mu value for this time step.  Need g_nj and f_nj
            !   and knowledge of whether or not there is a switch
            !   Calculate w_mu

            call calculate_w_mu(switch_attempt,w_mu,fnj_storage,nsteps,nstates_dyn,(final_step_reverse - istep),switch_attempt_state,&
                 &istate,occupied_adiabat,occupied_adiabat_beforehop)
            W = W*w_mu

            !**************************** Single adiabat propagation only
            ! W = 1.d0
            !****************************

         endif  !mdqt

         !---------------------------!
         !-- end of the MDQT stage --!
         !---------------------------!

         call z1_to_ze(z1,ze)
         call v1_to_ve(vz1,vze)

         call calculate_electronic_states(z1)
         efes = get_free_energy(istate)

         !-- calculate EVB weights
         if (weights) then
            call get_evb_weights
            call get_diabatic_populations
            call get_diabatic_populations_lfs(istate)
         endif

         !-- write the current data to the trajectory file

         if (mod(istep,ndump).eq.0) then
            if (weights) then
               write(itraj_channel,'(f13.6,6f12.5,i5,12f15.9)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn), (diabatic_population(k),k=1,ielst_dyn), (diabatic_population_lfs(k),k=1,ielst_dyn)
            else
               write(itraj_channel,'(f13.6,6f12.5,i5)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate
            endif
         endif

      enddo loop_over_time

   end subroutine forward_time_propagation_quantum_only

