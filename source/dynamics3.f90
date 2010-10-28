subroutine dynamics3
!===================================================================C
!
!  Driver for solvent dynamics on (multiple) two-dimensional
!  vibronic free energy surfaces. Gating coordinate
!  is assumed to be fixed (will be alleviated in later
!  versions).
!
!  OPTIONS:
!
!  MDQT - employ Tully's surface hopping algorithm to incorporate
!         non-adiabatic transitions between vibronic states
!         (otherwise classical dynamics on a single vibronic
!         free energy surface is simulated)
!
!  ADIAB - move on adiabatic electron/proton vibrational free energy surfaces
!
!  DIAB  - move on diabatic electron/proton vibrational free energy surfaces
!
!  KG=<int> - index of the grid point along the gating coordinate
!             (it is assumed that the gating coordinate is frozen
!             during the dynamics). The default is KG=1 which implies
!             that the number of grid points along the gating coordinate
!             is also NPNTSG=1 (gating distance is fixed at the values
!             from the input geometry). Otherwise KG must be smaller
!             than NPNTSG.
!
!  ISTATE[/SET] - index of the occupied vibronic state (1 for ground state) at t=0.
!                 If METHOD=1 then SET specifies one of the diabatic electronic
!                 states ("1" for 1a, "2" for 1b, "3" for 2a, "4" for 2b).
!                 If METHOD=2 then SET specifies one of the adiabatic electronic
!                 states ("1" for ground, "2" for first excited, etc.).
!
!  NOWEIGHTS - do not calculate the evb weights along the trajectory (default is YES)
!
!  TRANSFORM - initial values of solvent coordinates are given in the
!              rotated and scaled frame (z1,z2).
!
!  ZP0=<float> - center of the initial distribution along ZP (or z1) coordinate
!
!  ZE0=<float> - center of the initial distribution along ZE (or z2) coordinate
!
!  NSTATES=<int> - number of states to include in dynamics (for MDQT only)
!
!  TRAJOUT=<string> - name of the output files with trajectory data
!                     (default filename "trajectory_n.dat" where n is the
!                     trajectory index starting from 1)
!
!  SEED=<int> - random seed for random number generator.
!               (if SEED=0 then clock will be used to generate seed)
!               Default value is generated using the current time.
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
!                   Note that the parameters of the corresponding dielectric
!                   function must be specified within the SOLVENT keyword.
!
!  TSTEP=<float> - timestep for solvent dynamics in picoseconds (default=0.0005)
!
!  NSTEPS=<int> - number of steps (length of the trajectory)
!
!  NQSTEPS=<int> - number of TDSE steps per classical step in MDQT
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
!  T=<float> - temperature in K
!
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:35 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!  Revision 5.1  2010/10/26 21:06:20  souda
!  new routines/modules
!
!===================================================================C

   use pardim
   use keys
   use cst
   use strings
   use solmat
   use control
   use control_dynamics
   use random_generators
   use quantum
   use geogas, only: iptgas, xyzgas
   use parsol
   use propagators_3d
   use feszz_3d, only: reset_feszz3_counter

   implicit none

   character(len=1024) :: options
   character(len=  40) :: fname
   character(len=  15) :: zdim ="kcal/mol      "
   character(len=  15) :: z1dim="(kcal/mol)^1/2"
   character(len=   5) :: mode_dyn
   character(len=   5) :: traj_suffix

   character(len=3), dimension(2) :: iset_char_diab2=(/"1ab","2ab"/)
   character(len=2), dimension(4) :: iset_char_diab4=(/"1a","1b","2a","2b"/)

   logical :: adiab, diab2, diab4, weights
   logical :: switch=.false.
   logical :: transform=.false.

   integer :: nstates_dyn, nzdim_dyn, ielst_dyn, iseed_inp
   integer :: istate, new_state, ndump6

   integer :: ikey, ioption, islash, istart, kg0
   integer :: ize1, nze, ioutput, lenf, ispa
   integer :: itraj, istep, iqstep, k, itmp
   integer :: itraj_channel=11

   real(8) :: sigma, sample, population_current, wf_norm
   real(8) :: zeit_start, zeit_end, zeit_total, second
   real(8) :: zeit, zeit_prev
   real(8) :: zeitq, zeitq_prev
   real(8) :: z1, z2, zp, ze, vz1, vz2, vzp, vze, z10, z20, zp0, ze0, ekin, efes
   real(8) :: vz1_prev, vz2_prev

   adiab   = .true.
   diab2   = .false.
   diab4   = .false.
   weights = .true.

   !~~~~~~~~~~~~~~
   ! Print banner
   !~~~~~~~~~~~~~~
   write(6,*)
   write(6,'("================================================================")')
   write(6,'("              SOLVENT DYNAMICS MODULE (optional MDQT)           ")')
   write(6,'("         (two solvent coordinates + fixed gating distance)      ")')
   write(6,'("================================================================"/)')

   !~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' DYNAMICS3(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in DYNAMICS3): You MUST specify options for DYNAMICS3 keyword ***"/)')
      stop
   else
      call getopt(keywrd,ikey+11,options)
   endif

   !~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate
   !~~~~~~~~~~~~~~~~~~~

   ioption = index(options," KG=")
   
   if (ioption.ne.0) then
      kg0 = reada(options,ioption+4)
      if (kg0.gt.npntsg) kg0 = npntsg/2
   else
      kg0 = 1
   endif

   if (npntsg.eq.1) then
      kg0 = 1
      write(6,'(1x,"The gating distance is fixed at the value from the input geometry: ",f8.3," A"/)') abs(xyzgas(1,iptgas(3)) - xyzgas(1,iptgas(1)))
   else
      if (kg0.ge.1.and.kg0.le.npntsg) then
         write(6,'(1x," The gating distance is fixed at the value: ",f8.3," A (grid point #",i3,")"/)') glist(kg0),kg0
      else
         write(6,'(1x," The specified gating grid point #",i3," is outside the allowed range (1,",i3,")"/)') kg0,npntsg
         stop
      endif
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

   write(6,'(1x,"Temperature: ",f8.3," K"/)') temp

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
         stop
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
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify TAUD= option for DEBYE model ***"/)')
         stop
      endif
      call set_debye_model_parameters()
      write(6,'(1x,"Debye relaxation time TAUD (ps):        ",f15.6)') taud
      write(6,'(1x,"Longitudianl relaxation time TAUL (ps): ",f15.6)') taul
      write(6,'(1x,"Effective mass of the solvent (ps^2):   ",f15.6)') effmass

   elseif (solvent_model.eq."DEBYE2") then

      ioption = index(options,' TAU1=')
      if (ioption.ne.0) then
         tau1 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify TAU1= option for DEBYE2 model ***"/)')
         stop
      endif

      ioption = index(options,' TAU2=')
      if (ioption.ne.0) then
         tau2 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify TAU2= option for DEBYE2 model ***"/)')
         stop
      endif

      ioption = index(options,' EPS1=')
      if (ioption.ne.0) then
         eps1 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify EPS1= option for DEBYE2 model ***"/)')
         stop
      endif

      call set_debye2_model_parameters()
      write(6,'(1x,"First  relaxation time TAU1 (ps):       ",f15.6)') tau1
      write(6,'(1x,"Second relaxation time TAU2 (ps):       ",f15.6)') tau2
      write(6,'(1x,"Longitudianl relaxation time TAUL (ps): ",f15.6)') taul
      write(6,'(1x,"Effective mass of the solvent (ps^2):   ",f15.6)') effmass

   elseif (solvent_model.eq."ONODERA") then

      ioption = index(options,' TAUD=')
      if (ioption.ne.0) then
         taud = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify TAUD= option for ONODERA model ***"/)')
         stop
      endif

      ioption = index(options,' TAU0=')
      if (ioption.ne.0) then
         tau0 = reada(options,ioption+6)
      else
         write(*,'(/1x,"*** (in DYNAMICS): You MUST specify TAU0= option for ONODERA model ***"/)')
         stop
      endif

      call set_onodera_model_parameters()
      write(6,'(1x,"Inverse Pekar factor f_0         :       ",f15.6)') f0
      write(6,'(1x,"Debye   relaxation time TAUD (ps):       ",f15.6)') taud
      write(6,'(1x,"Onodera relaxation time TAU0 (ps):       ",f15.6)') tau0
      write(6,'(1x,"Longitudinal relaxation time TAUL  (ps): ",f15.6)') taul
      write(6,'(1x,"Longitudinal relaxation time TAU0L (ps): ",f15.6)') tau0l
      write(6,'(1x,"Effective mass of the solvent (ps^2):    ",f15.6)') effmass

   elseif (solvent_model.eq."ONODERA2") then

      !-- not implemented yet....
      call set_onodera2_model_parameters()

   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Type of dynamics (classical or MDQT)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' MDQT').ne.0) then
      mdqt = .true.
      write(6,'(/1x,"Mixed quantum-classical dynamics on multiple vibronic free energy surfaces.",/,&
                &1x,"Tullys fewest switches surface hopping algorithm (MDQT) will be utilized."/)')
   else
      mdqt = .false.
      write(6,'(/1x,"Classical dynamics on a single vibronic free energy surface (default)."/)')
   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Surface type (ADIAB, DIAB)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then
      mode_dyn  = 'ADIAB'
      adiab = .true.
      diab2 = .false.
      diab4 = .false.
      ielst_dyn = nelst
      write(6,'(/1x,"Solvent dynamics on adiabatic free energy surface(s)."/)')
   elseif (index(options,' DIAB2').ne.0) then
      mode_dyn  = 'DIAB2'
      adiab = .false.
      diab2 = .true.
      diab4 = .false.
      ielst_dyn = 2
      write(6,'(/1x,"Solvent dynamics on ET adiabatic free energy surface(s)."/)')
   elseif (index(options,' DIAB4').ne.0) then
      mode_dyn  = 'DIAB4'
      adiab = .false.
      diab2 = .false.
      diab4 = .true.
      ielst_dyn = 1
      write(6,'(/1x,"Solvent dynamics on a single diabatic free energy surface."/)')
      if (mdqt) then
         write(*,'(/1x,"MDQT in the diabatic representation is not implemented in the current version."/)')
         stop
      endif
   else
      mode_dyn  = 'ADIAB'
      adiab = .true.
      diab2 = .false.
      diab4 = .false.
      ielst_dyn = nelst
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! flag EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (index(options,' NOWEIGHTS').ne.0) then
      weights = .false.
      write(6,'(/1x,"EVB weights WILL NOT be calculated along the trajectory."/)')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initial occupied state at time t=0
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options,' ISTATE=')

   if (ioption.ne.0) then

      istart = ioption + 8
      ispa = index(options(ioption+1:),space)

      if (mode_dyn.eq.'ADIAB') then

         initial_state = reada(options,istart)
         if (initial_state.le.0) initial_state = 1
         initial_set = 1
         write(6,'(1x,"At t=0: initial adiabatic state: ",i6)') initial_state

      elseif (mode_dyn.eq.'DIAB2') then

         islash = index(options(istart:istart+ispa+1),'/')
         if (islash.eq.0) then
            write(*,'(/1x,"For diabatic representation DIAB2 you must specify BOTH the initial state",/,&
                      &1x," and the initial ET subset: two integers separated by a slash."/)')
            stop
         endif
         initial_state = reada(options(istart:istart+islash-2),1)
         initial_set = reada(options,istart+islash)
         if (initial_set.eq.1) then
            write(6,'(1x,"At t=0: initial ET diabatic state: ",i6," within the first (1a,1b) ET subset")') initial_state
         elseif (initial_set.eq.2) then
            write(6,'(1x,"At t=0: initial ET diabatic state: ",i6," within the second (2a,2b) ET subset")') initial_state
         else
            write(*,'(/1x,"*** (in DYNAMICS): subset in ISTATE keyword must be 1 or 2 ***"/)')
            stop
         endif

      elseif (mode_dyn.eq.'DIAB4') then

         islash = index(options(istart:istart+ispa+1),'/')
         if (islash.eq.0) then
            write(*,'(/1x,"For diabatic representation DIAB4 you must specify BOTH the initial state",/,&
                      &1x," and the initial diabatic electronic state: two integers separated by a slash."/)')
            stop
         endif
         initial_state = reada(options(istart:istart+islash-2),1)
         initial_set = reada(options,istart+islash)
         if (initial_set.eq.1) then
            write(6,'(1x,"At t=0: initial diabatic state: ",i6," within the 1a electronic set")') initial_state
         elseif (initial_set.eq.2) then
            write(6,'(1x,"At t=0: initial diabatic state: ",i6," within the 1b electronic set")') initial_state
         elseif (initial_set.eq.3) then
            write(6,'(1x,"At t=0: initial diabatic state: ",i6," within the 2a electronic set")') initial_state
         elseif (initial_set.eq.4) then
            write(6,'(1x,"At t=0: initial diabatic state: ",i6," within the 2b electronic set")') initial_state
         else
            write(*,'(/1x,"*** (in DYNAMICS): subset in ISTATE keyword must be 1 or 2 or 3 or 4 ***"/)')
            stop
         endif

      endif

   else

      initial_state = 1
      initial_set = 1
      write(6,'(1x,"(Default) At t=0: initial: ",i6," within the first electronic set")') initial_state

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of states in MDQT
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," NSTATES=")
   
   if (ioption.ne.0) then
      nstates_dyn = reada(options,ioption+9)
      write(6,'(/1x,"Number of vibronic states to include in MDQT dynamics: ",i4/)') nstates_dyn
   else
      nstates_dyn = nelst*nprst
      write(6,'(/1x,"Number of vibronic states to include in MDQT dynamics (default): ",i4/)') nstates_dyn
   endif

   if (.not.mdqt) then
      nstates_dyn = 1
      write(6,'(/1x,"However, no MDQT keyword was found: number of states is reset to 1 (one)."/)')
   endif

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

      call set_tdse_timestep(itmp,tstep)
      write(6,'(1x,"Timestep for TDSE: ",g15.6," ps"/)') tstep/real(itmp)

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
         elseif (interpolation.eq."DEBYE2") then
            write(6,'(1x,"(quadratic interpolation using the values at t, t+dt/2, and t+dt)")')
         else
            write(6,'(1x,"(ERROR in DYNAMICS3: UNKNOWN interpolation scheme. Check your input file.)")')
            stop
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
      write(6,'(1x,"Dump trajectory to screen data every ",i10," steps"/)') ndump
   else
      ndump6 = 0
   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Random seed
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options," SEED=")
   if (ioption.ne.0) then
      iseed_inp = reada(options,ioption+6)
   else
      !-- use clock to generate random seed (to be done...)
      iseed_inp = -11
   endif
   call set_random_seed(iseed_inp)

   write(6,'(1x,"Random seed: ",i6/)') iseed_inp

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! The solvent coordinate frame used for initial values
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ioption = index(options,"TRANSFORM")
   if (ioption.ne.0) then
      transform = .true.
   else
      transform = .false.
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initial solvent coordinates
   ! (center of the initial distribution)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioption = index(options,' ZP0=')
   if (ioption.ne.0) then
      ioption = ioption + 5
      zp0 = reada(options,ioption)
   else
      write(*,'(/1x,"*** (in DYNAMICS): You MUST specify ZP0= option for DYNAMICS keyword ***"/)')
      stop
   endif

   ioption = index(options,' ZE0=')
   if (ioption.ne.0) then
      ioption = ioption + 5
      ze0 = reada(options,ioption)
   else
      write(*,'(/1x,"*** (in DYNAMICS): You MUST specify ZE0= option for DYNAMICS keyword ***"/)')
      stop
   endif

   if (.not.transform) then
      call zpze_to_z1z2(zp0,ze0,z10,z20)
   else
      z10 = zp0
      z20 = ze0
      call z1z2_to_zpze(z10,z20,zp0,ze0)
   endif

   write(6,'(/1x,"Center of the initial distribution of solvent coordinates:",/,&
   &" ZP(0) = ",F7.3,2X,A," and ZE(0) = ",F7.3,2X,A,/,&
   &" z1(0) = ",F7.3,2X,A," and z2(0) = ",F7.3,2X,A)')&
   &  zp0,zdim,ze0,zdim,z10,z1dim,z20,z1dim

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate dynamics (not implemented yet)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output files for trajectories
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' TRAJOUT=')

   if (ioutput.eq.0) then

      fname = job(1:ljob)//'/trajectory'
      lenf = ljob + 11

   else

      ispa = index(options(ioutput+9:),' ')
      fname = options(ioutput+9:ioutput+ispa+9)
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
   elseif (diab2) then
      fname = trim(fname)//"_diab2_"//iset_char_diab2(initial_set)
   elseif (diab4) then
      fname = trim(fname)//"_diab4_"//iset_char_diab4(initial_set)
   endif

   write(6,'(/1x,"Trajectory data are written to the file(s) <",a,">"/)') trim(fname)//"_<NNNN>.dat"

   !-- set some variables in propagators module
   nzdim_dyn = nprst*ielst_dyn
   call set_mode(mode_dyn,initial_set,nstates_dyn,nzdim_dyn,ielst_dyn)
   istate = initial_state

   !-- Allocate arrays in propagators module
   call allocate_vibronic_states(nprst,npnts)
   if (mdqt) call allocate_mdqt_arrays
   if (weights) call allocate_evb_weights

   !===DEBUG===
   !call print_propagators_3d

   zeit_start = second()

   !======================================!
   !      MAIN LOOP OVER TRAJECTORIES     !
   !======================================!

   loop_over_trajectories: do itraj=1,ntraj

      !-- initialize the suffix of the output file
      write(traj_suffix,'(i5.5)') itraj

      !-- open the trajectory output file (channel 1)
      open(itraj_channel,file=trim(fname)//"_"//traj_suffix//".dat")

      !-- pick the initial values of solvent coordinates
      !   from gaussian distribution centered at ZP0,ZE0

      sigma = sqrt(kb*temp/f0)
      call gaussdist_boxmuller(sample,iseed)
      z1 = z10 + sigma*sample
      call gaussdist_boxmuller(sample,iseed)
      z2 = z20 + sigma*sample

      !-- initialize initial velocities

      if (solvent_model.eq."DEBYE".or.solvent_model.eq."DEBYE2") then

         !-- overdamped dynamics
         vz1 = 0.d0
         vz2 = 0.d0

      elseif (solvent_model.eq."ONODERA") then

         !-- initial velocities from Maxwell distribution
         sigma = sqrt(kb*temp/effmass)
         call gaussdist_boxmuller(sample,iseed)
         vz1 = sigma*sample
         call gaussdist_boxmuller(sample,iseed)
         vz2 = sigma*sample

      else

         !-- Other models are not implemented yet...
         write(6,'(1x,"(From DYNAMICS3: solvent model ",a10," is not implemented yet: abort calculation)")') solvent_model
         stop

      endif

      !-- in case of MDQT trajectory assign the initial state based
      !   on the quantum amplitudes of the initial wavefunction

      if (mdqt) then

         !-- calculate the adiabatic states and vibronic couplings at t=0

         call calculate_vibronic_states(kg0,z1,z2)
         call calculate_vibronic_couplings

         !-- set initial amplitudes (and/or density matrix) at t=0
         !--> to be properly coded (now it is a pure state)
         call set_initial_amplitudes(istate)
         call set_initial_density(istate)

         !-- choose the initial occupied state randomly
         !   according to the initial populations
         !--> to be coded, for now it is just istate <--

      endif

      !-- reset the counter of calls to feszz3
      call reset_feszz3_counter

      !-- write the header of the trajectory file

      if (weights) then
         write(itraj_channel,'("#",168("="))')
      else
         write(itraj_channel,'("#",141("="))')
      endif
      write(itraj_channel,'("#   Data for the trajectory ",i5)') itraj
      if (weights) then
         write(itraj_channel,'("#",168("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"z2",t44,"vz1",t56,"vz2",t68,"zp",t80,"ze",t92,"vzp",t103,"vze",t115,"Ekin",t126,"Efe",t135,"occ.",t142,"EVB weights (1a,1b,2a,2b)")')
         write(itraj_channel,'("#",168("-"))')
      else
         write(itraj_channel,'("#",141("-"))')
         write(itraj_channel,'("#",t6,"t(ps)",t20,"z1",t32,"z2",t44,"vz1",t56,"vz2",t68,"zp",t80,"ze",t92,"vzp",t103,"vze",t115,"Ekin",t126,"Efe",t135,"occ.")')
         write(itraj_channel,'("#",141("-"))')
      endif

      write(6,'(/1x,"===> Trajectory ",i5)') itraj
      write(6,*)
      write(6,'(137x,$)')

      !===============================!
      !   MAIN LOOP OVER TIME STEPS   !
      !===============================!
      loop_over_time: do istep=1,nsteps

         switch = .false.

         zeit_prev = real(istep-1)*tstep
         zeit = real(istep)*tstep

         !-- MDQT: store couplings, vibronic energies, and velocities
         !         from the previous step (for iterpolation)
         
         if (mdqt) then
            vz1_prev = vz1
            vz2_prev = vz2
            call store_vibronic_couplings                                 !  coupz(:,:) -> coupz_prev(:,:)
            call store_vibronic_energies                                  !  fe(:)      -> fe_prev(:)
            if (interpolation.eq."QUADRATIC") call store_wavefunctions    !  z(:,:)     -> z_prev(:,:)
         endif

         !-- Propagate solvent coordinates and velocities
         
         if (solvent_model.eq."DEBYE") then
         
            !-- overdamped Langevin equation (pure Debye model)
            call langevin_debye_2d(istate,kg0,z1,z2,vz1,vz2,tstep,temp,ekin,efes)
            !write(*,'(/1x,"DYNAMICS3: Debye propagator is not coded yet...")')
            !stop

         elseif (solvent_model.eq."DEBYE2") then

            !-- overdamped Langevin equation with memory friction
            !   (Debye model with two relaxation periods)
            !call langevin_debye2_2d(istate,kg0,z1,z2,vz1,vz2,tstep,temp,ekin,efes)
            write(*,'(/1x,"DYNAMICS3: Debye2 propagator is not coded yet...")')
            stop

         elseif (solvent_model.eq."ONODERA") then

            !-- ordinary Langevin equation (Onodera model)
            call langevin_onodera_2d(istate,kg0,z1,z2,vz1,vz2,tstep,temp,ekin,efes)

         elseif (solvent_model.eq."ONODERA2") then

            !-- ordinary Langevin equation with memory friction
            !   (Onodera model with two relaxation periods)
            !call langevin_onodera2_2d(istate,kg0,z1,z2,vz1,vz2,tstep,temp,ekin,efes)
            write(*,'(/1x,"DYNAMICS3: Onodera2 propagator is not coded yet...")')
            stop

         endif

         !----------------!
         !-- MDQT stage --!
         !----------------!

         if (mdqt) then

            !---------------------------------------------------------------------------------------------------
            !write(*,'(/1x,"DYNAMICS3: MDQT is not coded yet... Be patient!")')
            !stop
            !---------------------------------------------------------------------------------------------------

            !-- calculate vibronic couplings at t+dt
            !---------------------------------------
            call calculate_vibronic_couplings

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at t and t+dt
            !-------------------------------------------------------
            call calculate_v_dot_d(vz1,vz1_prev,vz2,vz2_prev)

            !-- Calculate the nonadiabatic coupling terms (v*d_{kl})
            !   at half timestep for quadratic interpolation scheme
            !-------------------------------------------------------
            if (interpolation.eq."QUADRATIC") then
               call calculate_v_dot_d_mid(tstep)
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

            !-- propagate the amplitudes and switch probabilities
            !   from t_prev to t
            !------------------------------------------------------------------
            do iqstep=1,nqsteps

               zeitq_prev = (iqstep-1)*qtstep + zeit_prev
               zeitq = iqstep*qtstep + zeit_prev

               !-- propagate amplitudes forward in time
               !-------------------------------------------------------
               call propagate_amplitudes_rk4(zeitq_prev,qtstep)
               !call propagate_density_rk4(zeitq_prev,qtstep)

               !-- calculate transition probabilities from current state
               !--------------------------------------------------------
               call calculate_bprob_amp(istate,zeitq)
               !call calculate_bprob_den(istate,zeitq)

               !-- accumulate swithing probabilities (array operation)
               !------------------------------------------------------
               call accumulate_switch_prob(qtstep)
               
            enddo

            !-- check the norm of the time-dependent wavefunction
            !-----------------------------------------------------
            wf_norm = tdwf_norm()
            if (abs(wf_norm-1.d0).gt.1.d-6) then
               write(*,*) "DYNAMICS3: Amplitudes are not normalized after timestep", istep
               call deallocate_all_arrays
               stop
            endif

            !-- Normalize swithing probabilities by the current state population
            !   and zero out the negative ones
            !-------------------------------------------------------------------
            call normalize_switch_prob(population_current)

            !-- construct the full density matrix (do we really need it?)
            !------------------------------------------------------------
            !call calculate_density_matrix

            !-- decision time: should we make a hop?
            !-------------------------------------------
            new_state = switch_state(istate)

            if (new_state.ne.istate) then
               !-- attempt adjusting velocities
               call adjust_velocities(istate,new_state,vz1,vz2,switch)
               if (switch) then
                  write(itraj_channel,'("#--------------------------------------------------------------------")')
                  write(itraj_channel,'("#  t  = ",f12.6," ps ==> switch ",i3,"  -->",i3)') zeit,istate,new_state
                  write(itraj_channel,'("#  d  = (",f20.6,",",f20.6,")")') &
                  & get_vibronic_coupling(istate,new_state)
                  write(itraj_channel,'("# |d| = ",f20.6)') &
                  & sqrt(dot_product(get_vibronic_coupling(istate,new_state),get_vibronic_coupling(istate,new_state)))
                  write(itraj_channel,'("#--------------------------------------------------------------------")')
                  istate = new_state
               endif
            endif

         endif  !mdqt

         !---------------------------!
         !-- end of the MDQT stage --!
         !---------------------------!

         !-- calculate EVB weights
         if (weights) call get_evb_weights

         !-- write the current data to the trajectory file

         call z1z2_to_zpze(z1,z2,zp,ze)
         call z1z2_to_zpze(vz1,vz2,vzp,vze)

         if (mod(istep,ndump).eq.0) then
            if (weights) then
               write(itraj_channel,'(f12.6,10f12.5,i5,4f8.3)') &
               & zeit, z1, z2, vz1, vz2, zp, ze, vzp, vze, ekin, efes, istate, (wght(k,istate),k=1,ielst_dyn)
            else
               write(itraj_channel,'(f12.6,10f12.5,i5)') &
               & zeit, z1, z2, vz1, vz2, zp, ze, vzp, vze, ekin, efes, istate
            endif
         endif

         if (ndump6.gt.0.and.mod(istep,ndump6).eq.0) &
         & write(6,'(137("\b"),f12.6,10f12.5,i5,$)') zeit, z1, z2, vz1, vz2, zp, ze, vzp, vze, ekin, efes, istate

      enddo loop_over_time

      write(itraj_channel,'("#",141("-"))')
      close(itraj_channel)
      
      write(6,*)
      write(6,*)

      !-- reset the initial state
      istate = initial_state

   enddo loop_over_trajectories

   zeit_end = second()
   zeit_total = zeit_end - zeit_start
   write(6,*)
   write(6,'(1x,"================================================================================")')
   write(6,'(1x,"Done. Time elapsed         (sec): ",f20.3)') zeit_total
   write(6,'(1x,"      Time per trajectory  (sec): ",f20.3)') zeit_total/ntraj
   write(6,'(1x,"      Time per timestep    (sec): ",f20.3)') zeit_total/(ntraj*nsteps)
   write(6,'(1x,"      Productivity rate (ps/day): ",f20.3)') 3600.d0*24.d0*tstep*ntraj*nsteps/zeit_total
   write(6,'(1x,"================================================================================")')

   !-- Deallocate arrays in propagators module

   call deallocate_all_arrays

contains

   subroutine deallocate_all_arrays
      call deallocate_vibronic_states
      if (mdqt) call deallocate_mdqt_arrays
      if (weights) call deallocate_evb_weights
   end subroutine deallocate_all_arrays

end subroutine dynamics3

