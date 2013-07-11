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
!-------------------------------------------------------------------
!
!  $Author$
!  $Id$
!  $Revision$
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
   character(len=  20) :: str
   character(len=1), dimension(4) :: iset_char_diab4=(/"1","2","3","4"/)

   logical :: adiab, diab4, weights
   logical :: switch=.false.
   logical :: success=.false.
   logical :: scaled=.false.
   logical :: initial_state_pure=.true.
   logical :: initial_state_diab=.false.
   logical :: purestate=.true.

   logical :: interaction_region_prev=.false.
   logical :: interaction_region=.false.
   logical :: double_well=.false.

   integer :: nstates_dyn, nzdim_dyn, ielst_dyn, iseed_inp, iset_dyn
   integer :: initial_state=-1, iground=1
   integer :: istate, new_state, ndump6, ndump777, ndump888
   integer :: number_of_skipped_trajectories=0
   integer :: number_of_failed_trajectories=0
   integer :: itraj_start=1
   integer :: ntraj_valid
   integer :: nqsteps_var

   integer :: ikey, ikeya, ikeyb, ioption, ioption2, ioption3, kg0, ifrom, ito
   integer :: islash1, islash2, islash3, idum1, idum2, idum3, idum4
   integer :: ioutput, lenf, ispa, idash, icount, ireac, iprod
   integer :: itraj, istep, iqstep, k, itmp, itmp1
   integer :: number_of_switches=0, number_of_rejected=0
   integer :: itraj_channel=11
   integer :: ifes_channel=12

   real(kind=8) :: sigma, sigma1, sample, population_current
   real(kind=8) :: wf_norm, zmom1_norm, pzmom1_norm
   real(kind=8) :: zeit_start, zeit_end, zeit_total, traj_time_start, traj_time_end
   real(kind=8) :: zeit, zeit_prev
   real(kind=8) :: zeitq, zeitq_prev
   real(kind=8) :: z1, ze, vz1, vze, z10, ze0, y1
   real(kind=8) :: ekin, ekin1, ekin_prev, ekinhalf1, efes
   real(kind=8) :: vz1_prev
   real(kind=8) :: qtstep_var

   real(kind=8) :: hreac, hprod, vet, zi, zi_scaled
   real(kind=8) :: dg_reaction, dg_activation, et_marcus_rate
   real(kind=8) :: kappa_ad, prefactor1, prefactor2
   real(kind=8) :: et_rips_jortner_rate, adiabatic_rate
   real(kind=8) :: fr, fb, eb
   real(kind=8), dimension(4,4) :: h0k
   real(kind=8), dimension(4,4) :: tk, tinfk, trk, trinfk
   real(kind=8), dimension(2)   :: fe_diab, fe_adiab


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

   if (solvent_model.eq."DEBYE") then

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

      mdqt = .true.
      write(6,'(/1x,"Mixed quantum-classical dynamics on two electronic free energy surfaces.",/,&
                &1x,"Tullys fewest switches surface hopping algorithm (MDQT) will be utilized."/)')

      if (index(options,' PHASE').ne.0) then
         phase_corr = .true.
         write(6,'(/1x,"Phase correction algorithm will be used.",/,&
                   &1x,"[N. Shenvi, J. E. Subotnik, and W. Yang, J. Chem. Phys. 135, 024101 (2011) ]"/)')
      endif


      !-- decoherence options

      if (index(options,' AFSSH').ne.0.or.&
         &index(options,' COLLAPSE_REGION_COUPLING').ne.0) then

         !-- make sure that only one decoherence option has been chosen
         if (index(options,' AFSSH').ne.0.and.&
            &index(options,' COLLAPSE_REGION_COUPLING').ne.0) then

            write(6,'(/1x,"Only one decoherence option can be specified.")')
            write(6,'( 1x,"Your input contains both AFSSH and COLLAPSE_REGION_COUPLING options in DYNAMICS2 keyword.")')
            write(6,'( 1x,"Check your input and make up your mind!!! Aborting...")')
            call clean_exit

         endif

         decoherence = .true.

      endif


      if (index(options,' AFSSH').ne.0) then

         if (.not.phase_corr) then

            afssh = .true.
            collapse_region_coupling = .false.

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

         ioption = index(options," COUPLING_CUTOFF=")
         if (ioption.ne.0) then
            coupling_cutoff = reada(options,ioption+17)
         else
            coupling_cutoff = 1.d-5
         endif
         write(6,'(/1x,"Simple decoherence algorithm with collapsing events occuring upon leaving the interaction region")')
         write(6,'( 1x,"The interaction region is defined as region where the largest nonadiabatic coupling")')
         write(6,'( 1x,"is smaller in magnitude than the cutoff value of ",g15.6," (kcal/mol)^{-1/2}"/)') coupling_cutoff

      endif

   else

      mdqt = .false.
      write(6,'(/1x,"Simulation of classical dynamics on a single free energy surface (default)."/)')

   endif


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

      if (solvent_model.eq."DEBYE") then
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
      write(ifes_channel,'("#   Adiabatic activation energy:  ",e20.9,"kcal/mol")') eb
      write(ifes_channel,'("#   Adiabatic rate constant:      ",e20.9," 1/sec")') adiabatic_rate
      write(ifes_channel,'("#",74("-"))')
   endif

   write(ifes_channel,'("#",t7,"Z(gap)",t19,"z(scaled)",t30,"U1(diab)",t42,"U2(diab)",t54,"U1(adiab)",t66,"U2(adiab)")')
   write(ifes_channel,'("#",74("-"))')

   do zi=3.d0*lambda,-3.d0*lambda,-6.d0*lambda/101.d0
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
   ioption2 = index(options,' ZE0R')
   ioption3 = index(options,' ZE0P')

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

   call ze_to_z1(ze0,z10)

   write(6,'(/1x,"Center of the initial distribution of the solvent coordinate:",/,&
   &          1x,"Ze(0) = ",F10.3,2X,A,/,&
   &          1x,"z1(0) = ",F10.3,2X,A)') ze0,zgapdim,z10,zscadim

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

      sigma = sqrt(kb*temp/f0)
      sample = gaussdist_boxmuller()
      z1 = z10 + sigma*sample

      if (solvent_model.eq."ONODERA2") then
         !-- pick the initial value of auxiliary solvent coordinate
         !   from a gaussian distribution centered at (-z10)
         sigma = sqrt(kb*temp/gamma)
         sample = gaussdist_boxmuller()
         y1 = -z10 + sigma*sample
      endif

      !--(DEBUG)--start
      !-- ignoring sampling
      !z1 = z10
      !if (solvent_model.eq."ONODERA2") then
      !   y1 = -z10
      !endif
      !--(DEBUG)--end

      !-- zero out the moments (A-FSSH)
      if (mdqt.and.afssh) then
         call reset_zmoments
         call reset_pzmoments
      endif

      !-- initialize initial velocity

      if (solvent_model.eq."DEBYE".or.solvent_model.eq."DEBYE2") then

         !-- overdamped dynamics
         vz1 = 0.d0

      elseif (solvent_model.eq."ONODERA") then

         !-- initial velocity from Maxwell distribution
         sigma1 = sqrt(kb*temp/effmass1)
         sample = gaussdist_boxmuller()
         vz1 = sigma1*sample

      elseif (solvent_model.eq."ONODERA2") then

         !-- initial velocity from Maxwell distribution
         sigma1 = sqrt(kb*temp/effmass1)
         sample = gaussdist_boxmuller()
         vz1 = sigma1*sample

      else

         !-- Other models are not implemented yet...
         write(6,'(1x,"(From DYNAMICS3: solvent model ",a10," is not implemented: abort calculation)")') solvent_model
         call clean_exit

      endif

      !-- calculate electronic state for the initial value of solvent coordinate
      !   at t=0 (very first time for this trajectory)
      call calculate_electronic_states(z1)

      if (collapse_region_coupling) then
         interaction_region = interaction_region_check()
         interaction_region_prev = interaction_region
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

      else

         write(6,'(/1x,"From DYNAMICS3: For unknown reason (VERY SERIOUS BUG?) no initial condition was chosen. Abort.")')
         call clean_exit

      endif


      if (mdqt) then

         !-- calculate electronic states at t=0 (very first time for this trajectory)
         !call calculate_electronic_states(z1)

         !-- print out the initial amplitudes of the time-dependent wavefunction
         call print_initial_amplitudes(6)

         !-- calculate force matrices (A-FSSH specific)
         if (afssh) call calculate_force_matrices(z1)

      endif

      !-- open the trajectory output file (channel 1)

      open(itraj_channel,file=trim(fname)//"_"//traj_suffix//".dat")

      !-- write the header of the trajectory file

      if (weights) then
         call get_evb_weights
         write(itraj_channel,'("#",131("="))')
      else
         write(itraj_channel,'("#",131("="))')
      endif
      write(itraj_channel,'("#   Data for the trajectory #",i6.6)') itraj
      if (weights) then
         write(itraj_channel,'("#",131("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"VZe",t68,"Ekin",t80,"Efe",t89,"occ.",t100,"EVB weights (1,2)")')
         write(itraj_channel,'("#",131("-"))')
      else
         write(itraj_channel,'("#",131("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"vze",t68,"Ekin",t80,"Efe",t89,"occ.")')
         write(itraj_channel,'("#",131("-"))')
      endif

      write(6,'(/1x,"===> Trajectory ",i5," starts on the electronic state ",i3)') itraj, istate
      write(6,'( 1x,"===> Initial solvent coordinate (z1), (kcal/mol)^(1/2): ",f13.6)') z1

      write(6,*)
      write(6,'(131("-"))')
      write(6,'("#",t6,"t(ps)",t20,"z1",t32,"vz1",t44,"Ze",t56,"vze",t68,"Ekin",t80,"Efe",t89,"occ.")')
      write(6,'(131("-"))')

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
          write(itraj_channel,'(f13.6,6f12.5,i5,4f15.9)') &
          & 0.d0, z1, vz1, ze, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn)
      else
          write(itraj_channel,'(f13.6,6f12.5,i5)') &
          & 0.d0, z1, vz1, ze, vze, ekin, efes, istate
      endif

      number_of_switches = 0
      number_of_rejected = 0

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

         if (solvent_model.eq."DEBYE") then

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
            call langevin_onodera_1d(istate,z1,vz1,tstep,temp,ekin1,efes)
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
                     write(itraj_channel,'("#",131("-"))')
                  else
                     write(itraj_channel,'("#",131("-"))')
                  endif

                  write(itraj_channel,'("# Amplitudes are not normalized after timestep ",i6)') istep
                  write(itraj_channel,'("# Norm of the time-dependent wavefunction:     ",g20.10)') wf_norm
                  write(itraj_channel,'("# This trajectory has failed... Even after several tries with smaller TDSE timesteps.")')

                  if (weights) then
                     write(itraj_channel,'("#",131("-"))')
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
                  write(itraj_channel,'("#--------------------------------------------------------------------")')

                  istate = new_state
                  number_of_switches = number_of_switches + 1

               else

                  write(itraj_channel,'("#--------------------------------------------------------------------")')
                  write(itraj_channel,'("#  t  = ",f13.6," ps ==> REJECTED SWITCH ",i3,"  -->",i3)') zeit,istate,new_state
                  write(itraj_channel,'("#  d  = ",f20.6)') get_nonadiabatic_coupling(istate,new_state)
                  write(itraj_channel,'("#--------------------------------------------------------------------")')

                  number_of_rejected = number_of_rejected + 1

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
                  write(*,'("*** Leaving interaction region: wavefunction collapsed to pure state ",i2)') istate
               endif
            endif

         endif  !mdqt

         !---------------------------!
         !-- end of the MDQT stage --!
         !---------------------------!


         !-- calculate EVB weights
         if (weights) call get_evb_weights

         !-- write the current data to the trajectory file

         call z1_to_ze(z1,ze)
         call v1_to_ve(vz1,vze)

         if (mod(istep,ndump).eq.0) then
            if (weights) then
               write(itraj_channel,'(f13.6,6f12.5,i5,4f15.9)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn)
            else
               write(itraj_channel,'(f13.6,6f12.5,i5)') &
               & zeit, z1, vz1, ze, vze, ekin, efes, istate
            endif
         endif

         !if (ndump6.gt.0.and.mod(istep,ndump6).eq.0) &
         !& write(6,'(137("\b"),f13.6,6f12.5,i5,$)') zeit, z1, vz1, ze, vze, ekin, efes, istate

         if (ndump6.gt.0.and.mod(istep,ndump6).eq.0) &
         & write(6,'(f13.6,6f12.5,i5)') zeit, z1, vz1, ze, vze, ekin, efes, istate

      enddo loop_over_time

      traj_time_end = second()

      if (weights) then
         write(itraj_channel,'("#",131("-"))')
      else
         write(itraj_channel,'("#",131("-"))')
      endif
      write(itraj_channel,'("# Number of allowed  switches: ",i5)') number_of_switches
      write(itraj_channel,'("# Number of rejected switches: ",i5)') number_of_rejected
      if (weights) then
         write(itraj_channel,'("#",131("-"))')
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

contains

   subroutine clean_exit
      call deallocate_all_arrays
      stop
   end subroutine clean_exit

end subroutine dynamicset2

