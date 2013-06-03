!======================================================================================
!
!  Analysis of the dynamical MDQT trajectories (two-state ET model)
!  Alexander V. Soudackov, Penn State University
!  July 7, 2012
!
!  $Author$
!  $Id$
!  $Revision$
!
!======================================================================================

program analyze_et2_trajectories

   use marcus
   use string_utilities
   use sorting
   use timers
   implicit none

   integer, parameter :: number_of_bins_z1=100
   integer, parameter :: number_of_bins_ze=100
   integer, parameter :: ndt=27
   integer, parameter :: ndt2 = (ndt-1)/2

   logical, parameter :: calculate_reactive_flux_rate = .false.

   character(len=60)  :: filename
   character(len=200) :: record
   character(len=3)   :: state_suffix

   logical :: found

   integer :: number_of_traj, number_of_timesteps, number_of_states
   integer :: number_of_occ_states
   integer :: itraj, istep, idstep, eof=0, ierr1, ierr2
   integer :: i, ilargest, ii, iocc, k, i1, iargc, nsteps
   integer :: ibin_1, ibin_e
   integer :: ichannel_z1, ichannel_ze

   real(kind=8) :: z1_min, z1_max
   real(kind=8) :: ze_min, ze_max, z_divide, z_crossing
   real(kind=8) :: z1_curr, ze_curr, ze_curr1, ze_curr2
   real(kind=8) :: vz1_curr, vze_curr
   real(kind=8) :: z1_0, ze_0, efe_curr, ekin_curr
   real(kind=8) :: n_aver_r, n_aver_p, wh1_mean_0, w1_mean_0
   real(kind=8) :: time_start, time_end, total_time_start, total_time_end

   real(kind=8) :: bin_width_z1, bin_width_ze

   real(kind=8) :: z1tav
   real(kind=8) :: zetav

   real(kind=8) :: w1tav, w2tav

   real(kind=8) :: ze_cut, ze_crossing, ze_product, dze, time_mfpt
   integer :: number_of_reactive_events

   real(kind=8), dimension(number_of_bins_z1) :: bin_center_z1
   real(kind=8), dimension(number_of_bins_ze) :: bin_center_ze

   real(kind=8), dimension(:,:),   allocatable :: time
   real(kind=8), dimension(:),     allocatable :: time1
   real(kind=8), dimension(:,:),   allocatable :: z1, vz1
   real(kind=8), dimension(:,:),   allocatable :: ze, vze
   real(kind=8), dimension(:,:),   allocatable :: ekin, efe
   real(kind=8), dimension(:,:),   allocatable :: w1, w2
   real(kind=8), dimension(:),     allocatable :: histogram_z1
   real(kind=8), dimension(:),     allocatable :: histogram_ze
   real(kind=8), dimension(:,:),   allocatable :: state_histogram_z1
   real(kind=8), dimension(:,:),   allocatable :: state_histogram_ze
   real(kind=8), dimension(:),     allocatable :: z1_mean, ze_mean
   real(kind=8), dimension(:),     allocatable :: vz1_mean, vze_mean
   real(kind=8), dimension(:),     allocatable :: z1_var, ze_var
   real(kind=8), dimension(:),     allocatable :: z11_tcf
   real(kind=8), dimension(:),     allocatable :: zee_tcf
   real(kind=8), dimension(:),     allocatable :: theta_corr
   real(kind=8), dimension(:),     allocatable :: efe_mean, ekin_mean
   real(kind=8), dimension(:),     allocatable :: w1_mean, w2_mean
   real(kind=8), dimension(:),     allocatable :: wh1_mean, wh2_mean
   real(kind=8), dimension(:,:),   allocatable :: pop_ad

   integer, dimension(:,:), allocatable :: istate
   integer, dimension(:),   allocatable :: all_states
   integer, dimension(:),   allocatable :: istate_occ

   character(len=20), dimension(100) :: carr
   integer,           dimension(100) :: iarr
   real(kind=8),      dimension(100) :: rarr
   real(kind=8),      dimension(2)   :: tmparray

   !-- Marcus parameters

   logical :: marcus_flag
   real(kind=8) :: electronic_coupling
   real(kind=8) :: reorganization_energy
   real(kind=8) :: reaction_free_energy
   real(kind=8) :: eps_0, eps_inf, tau_2
   real(kind=8) :: temperature, k_fit_log, k_fit_odr, rcorr, shift, defect
   namelist /marcus_parameters/ electronic_coupling, reorganization_energy, reaction_free_energy, &
                              & eps_0, eps_inf, tau_2, temperature

   real(kind=8), dimension(2) :: p_marcus, p_rips_jortner, p_zusman, p_fit0, p_fit

   !--------------------------------------------------------------------------------------

   total_time_start = secondi()

   !-- read Marcus parameters from the input file

   call getarg(1,filename)

   open(1,file=trim(filename),status="old",iostat=ierr1)
   read(1,nml=marcus_parameters,iostat=ierr2)

   if (ierr1.ne.0.or.ierr2.ne.0) then
   
      marcus_flag = .false.
      write(*,*) "*** Error reading the file with Marcus parameters."
      write(*,*) "*** Marcus rate constant and populations will not be calculated."
      write(*,*)

   else

      marcus_flag = .true.
      call set_marcus_parameters(electronic_coupling, reorganization_energy, reaction_free_energy, eps_0, eps_inf, tau_2, temperature)

      !-- calculate equilibrium quantities
      call calculate_equilibrium_quantities()

      write(*,'(/1x,"======================================================")')
      write(*,'( 1x,"Equilibrium quantities")')
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Temperature:                     ",f12.3," K")')        temperature
      write(*,'( 1x,"ET reaction free energy:         ",f12.3," kcal/mol")') reaction_free_energy
      write(*,'( 1x,"Equilibrium constant:            ",g20.10)') k_equil
      write(*,'( 1x,"Equilibrium reactant population: ",g12.6)')  n_equil_r
      write(*,'( 1x,"Equilibrium product  population: ",g12.6)')  n_equil_p
      write(*,'( 1x,"======================================================")')

      !-- calculate Marcus nonadiabatic rate constant in ps^(-1)
      call calculate_marcus_rate_constant()

      write(*,'(/1x,"======================================================")')
      write(*,'( 1x,"Marcus model parameters")')
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Temperature:                ",f12.3," K")')        temperature
      write(*,'( 1x,"Electronic coupling:        ",f12.3," kcal/mol")') electronic_coupling
      write(*,'( 1x,"Reorganization free energy: ",f12.3," kcal/mol")') reorganization_energy
      write(*,'( 1x,"ET reaction free energy:    ",f12.3," kcal/mol")') reaction_free_energy
      write(*,'( 1x,"ET activation free energy:  ",f12.3," kcal/mol")') (reorganization_energy+reaction_free_energy)**2.d0/(4.d0*reorganization_energy)
      write(*,'( 1x,"ET Marcus rate constant:    ",e16.9," ps^(-1)")')  k_marcus
      write(*,'( 1x,"R-J Adiabaticity parameter: ",e16.9)')             kappa_ad
      write(*,'( 1x,"======================================================")')


      !-- calculate Rips-Jortner and Zusman nonadiabatic rate constants in ps^(-1)

      call calculate_rips_jortner_rate_constant()
      call calculate_zusman_rate_constant()

      write(*,'(/1x,"======================================================")')
      write(*,'( 1x,"Dielectric relaxation parameters")')
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Static dielectric constant:  ",f12.3)') eps_0
      write(*,'( 1x,"Optical dielectric constant: ",f12.3)') eps_inf
      write(*,'( 1x,"Longest longitudinal relaxation period: ",f12.3," ps")') eps_inf*tau_2/eps_0
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"R-J Adiabaticity parameter: ",e16.9)')             kappa_ad
      write(*,'( 1x,"Rips-Jortner rate constant: ",e16.9," ps^(-1)")')  k_rips_jortner
      write(*,'( 1x,"Zusman rate constant:       ",e16.9," ps^(-1)")')  k_zusman
      write(*,'( 1x,"======================================================")')

   endif

   number_of_traj = iargc() - 1
   if (number_of_traj == 0) then
      write(*,*) "No trajectory files specified... Come again..."
      stop
   endif

   !-- scan the first trajectory file
   !   to get the number of time steps

   call getarg(2,filename)
   filename = trim(filename)
   open(1,file=filename)

   number_of_timesteps = 0
   do
      read(1,'(a)',iostat=eof) record
      !-- end-of-file condition
      if (eof.ne.0) exit
      !-- skip comments and empty lines
      if (record(1:1).eq."#".or.record.eq."") cycle
      number_of_timesteps = number_of_timesteps + 1
   enddo
   close(1)

   write(*,*)
   write(*,'(1x,"Number of timesteps:    ",i6)') number_of_timesteps
   write(*,'(1x,"Number of trajectories: ",i6)') number_of_traj
   write(*,*)

   if (number_of_timesteps.eq.0) then
      write(*,*)
      write(*,*) "The first trajectory file is empty, Something is wrong..."
      write(*,*)
      stop
   endif

   !-- allocate the data arrays

   allocate(time(number_of_traj,number_of_timesteps))
   allocate(z1(number_of_traj,number_of_timesteps))
   allocate(ze(number_of_traj,number_of_timesteps))
   allocate(vz1(number_of_traj,number_of_timesteps))
   allocate(vze(number_of_traj,number_of_timesteps))
   allocate(ekin(number_of_traj,number_of_timesteps))
   allocate(efe(number_of_traj,number_of_timesteps))
   allocate(istate(number_of_traj,number_of_timesteps))
   allocate(w1(number_of_traj,number_of_timesteps))
   allocate(w2(number_of_traj,number_of_timesteps))

   call allocated_memory

   write(*,'(1x,"All data arrays allocated, now reading data files... ",t64,"--> ",$)')
   time_start = secondi()

   !=============================================================
   !-- reading all trajectory files and fill the data arrays
   !=============================================================

   loop_over_trajectories: do itraj=1,number_of_traj

      call getarg(itraj+1,filename)
      open(1,file=trim(filename))

      nsteps = 0

      loop_over_timesteps: do

         read(1,'(a)',iostat=eof) record

         !-- end-of-file condition
         if (eof.ne.0) exit

         !-- skip comments and empty lines
         if (record(1:1).eq."#".or.record.eq."") cycle

         call readrec(record,"rrrrrrrirr",carr,iarr,rarr)

         nsteps = nsteps + 1

         if (nsteps.gt.number_of_timesteps) then
            call deallocate_all_arrays
            write(*,*)
            write(*,*) "Trajectory file ",trim(filename)," has inconsistent number of timesteps."
            write(*,*) "Skipping this trajectory..."
            write(*,*)
            cycle loop_over_trajectories
         endif

         time(itraj,nsteps) = rarr(1)
         z1(itraj,nsteps)   = rarr(2)
         vz1(itraj,nsteps)  = rarr(3)
         ze(itraj,nsteps)   = rarr(4)
         vze(itraj,nsteps)  = rarr(5)
         ekin(itraj,nsteps) = rarr(6)
         efe(itraj,nsteps)  = rarr(7)

         istate(itraj,nsteps) = iarr(1)
 
         w1(itraj,nsteps) = rarr(8)
         w2(itraj,nsteps) = rarr(9)

      enddo loop_over_timesteps

   enddo loop_over_trajectories

   allocate (time1(nsteps))
   time1(:) = time(1,:)
   deallocate(time)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !=============================================================
   !-- Write the time values (ps) to the external file
   !=============================================================

   open(1,file="timesteps.dat",form="formatted")
   do i=1,number_of_timesteps
      write(1,'(i10,2x,f20.6)') i, time1(i)
   enddo
   close(1)

   !=============================================================
   !-- Scan the occupied state array and determine the highest
   !   occupied state for the whole set of trajectories.
   !=============================================================

   call allocated_memory

   write(*,*)
   write(*,'(1x,"Scanning MDQT trajectories to identify unique occupied states... ")')
   time_start = secondi()

   !=============================================================
   !-- Scan the occupied state array and create a list of states
   !   occupied at least once for the whole set of trajectories.
   !=============================================================

   number_of_states = 0
   number_of_occ_states = 0

   !--------------------------------------------------------------------------------
   !--(ET-2 specific) there are only two states in ET, so we already know the answer
   !                  assuming that both states get occupied at least once
   !
   !number_of_states = 2
   !number_of_occ_states = 2
   !--------------------------------------------------------------------------------
 
   !-- unwrap 2D array and store in the temporary array

   allocate (all_states(number_of_traj*number_of_timesteps))

   k = 0
   do itraj=1,number_of_traj
      do istep=1,number_of_timesteps
         k = k + 1
         iocc = istate(itraj,istep)
         all_states(k) = istate(itraj,istep)
         if (iocc.gt.number_of_states) number_of_states = iocc
      enddo
   enddo

   !-- sort the array all_states
   call shellsort(number_of_traj*number_of_timesteps,all_states)

   allocate(istate_occ(number_of_states))
   istate_occ = 0

   !-- find the number of unique occupied states
   !-- and fill the array of unique occupied states

   ii = 1
   istate_occ(1) = all_states(1)
   do i=2,number_of_traj*number_of_timesteps
      if (all_states(i).ne.all_states(i-1)) then
         ii = ii + 1
         istate_occ(ii) = all_states(i)
      endif
   enddo
   number_of_occ_states = ii

   !-- release temporary array
   deallocate(all_states)

   !-- sort the array istate_occ
   !call shellsort(number_of_occ_states,istate_occ)

   !-------------------
   !--(ET-2 specific)
   !-------------------
   !istate_occ(1) = 1
   !istate_occ(2) = 2
   !-------------------

   write(*,'(/1x,"Highest occupied adiabatic state: ",i3)') number_of_states
   write(*,'( 1x,"Number of unique occupied states: ",i3)') number_of_occ_states
   write(*,'( 1x,"Unique occupied states:",/,(10(1x,i3)))') (istate_occ(k),k=1,number_of_occ_states)

   time_end = secondi()
   write(*,'(/1x,"Done in ",f10.3," sec"/)') time_end-time_start

   !======================================
   !-- Start calculating the observables
   !======================================

   !---------------------------------------------------------------
   !--(1)-- Time-dependent distributions of solvent coordinates
   !---------------------------------------------------------------

   !-- find the global ranges for solvent coordinates

   z1_min =  999999.d0
   z1_max = -999999.d0
   ze_min =  999999.d0
   ze_max = -999999.d0

   do itraj=1,number_of_traj
      do istep=1,number_of_timesteps

         z1_curr = z1(itraj,istep)
         ze_curr = ze(itraj,istep)

         if (z1_curr.lt.z1_min) z1_min = z1_curr
         if (z1_curr.gt.z1_max) z1_max = z1_curr

         if (ze_curr.lt.ze_min) ze_min = ze_curr
         if (ze_curr.gt.ze_max) ze_max = ze_curr

      enddo
   enddo

   !-- shift min and max by very small amounts to ensure that
   !   binning will always put the values within bounds

   z1_min = z1_min - 1.d-10
   ze_min = ze_min - 1.d-10

   z1_max = z1_max + 1.d-10
   ze_max = ze_max + 1.d-10

   !-- define the histogram parameters

   bin_width_z1 = (z1_max - z1_min)/number_of_bins_z1
   bin_width_ze = (ze_max - ze_min)/number_of_bins_ze

   do i=1,number_of_bins_z1
      bin_center_z1(i) = z1_min + (i-1)*bin_width_z1 + bin_width_z1/2.d0
   enddo

   do i=1,number_of_bins_ze
      bin_center_ze(i) = ze_min + (i-1)*bin_width_ze + bin_width_ze/2.d0
   enddo

   !---------------------------------------------------------------
   !-- build histograms for occupied adiabatic states
   !---------------------------------------------------------------

   write(*,'(1x,"Building state-resolved histograms for solvent coordinates... ",t64,"--> ",$)')
   time_start = secondi()

   allocate(state_histogram_z1(number_of_states,number_of_bins_z1))
   allocate(state_histogram_ze(number_of_states,number_of_bins_ze))

   !-- open a separate file for each occupied state

   do i=1,number_of_occ_states
      ichannel_z1 = 10 + 2*istate_occ(i) - 1
      ichannel_ze = 10 + 2*istate_occ(i)
      write(state_suffix,'(i3.3)') istate_occ(i)
      open(ichannel_z1,file="z1_distribution_"//state_suffix//".dat")    !,form="unformatted")
      open(ichannel_ze,file="ze_distribution_"//state_suffix//".dat")    !,form="unformatted")
   enddo

   do istep=1,number_of_timesteps

      state_histogram_z1 = 0.d0
      state_histogram_ze = 0.d0

      do itraj=1,number_of_traj

         z1_curr = z1(itraj,istep)
         ze_curr = ze(itraj,istep)

         iocc = istate(itraj,istep)

         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         state_histogram_z1(iocc,ibin_1) = state_histogram_z1(iocc,ibin_1) + 1.d0

         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)
         state_histogram_ze(iocc,ibin_e) = state_histogram_ze(iocc,ibin_e) + 1.d0

      enddo

      !-- normalize distributions
   
      do i1=1,number_of_bins_z1
         do k=1,number_of_occ_states
            state_histogram_z1(istate_occ(k),i1) = state_histogram_z1(istate_occ(k),i1)/number_of_traj
         enddo
      enddo

      do i1=1,number_of_bins_ze
         do k=1,number_of_occ_states
            state_histogram_ze(istate_occ(k),i1) = state_histogram_ze(istate_occ(k),i1)/number_of_traj
         enddo
      enddo

      !-- output to the external files for visualization

      !-- state-resolved distributions

      do i=1,number_of_occ_states

         ichannel_z1 = 10 + 2*istate_occ(i) - 1
         ichannel_ze = 10 + 2*istate_occ(i)

         do i1=1,number_of_bins_z1
            write(ichannel_z1,'(2f15.6,1x,g15.6)') time1(istep), bin_center_z1(i1), state_histogram_z1(istate_occ(i),i1)
            !write(ichannel_z1) bin_center_z1(i1), state_histogram_z1(istate_occ(i),i1)
         enddo
         write(ichannel_z1,*)

         do i1=1,number_of_bins_ze
            write(ichannel_ze,'(2f15.6,1x,g15.6)') time1(istep), bin_center_ze(i1), state_histogram_ze(istate_occ(i),i1)
            !write(ichannel_ze) bin_center_ze(i1), state_histogram_ze(istate_occ(i),i1)
         enddo
         write(ichannel_ze,*)

      enddo

   enddo

   do i=1,number_of_occ_states
      ichannel_z1 = 10 + 2*istate_occ(i) - 1
      ichannel_ze = 10 + 2*istate_occ(i)
      close(ichannel_z1)
      close(ichannel_ze)
   enddo

   deallocate (state_histogram_z1,state_histogram_ze)
   deallocate (istate_occ)
   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !-- build global histograms
   !---------------------------------------------------------------

   write(*,'(1x,"Building global histograms for solvent coordinates... ",t64,"--> ",$)')
   time_start = secondi()

   !-- global distributions
   allocate(histogram_z1(number_of_bins_z1))
   allocate(histogram_ze(number_of_bins_ze))

   open(21,file="z1_distribution_global.dat")    !,form="unformatted")
   open(22,file="ze_distribution_global.dat")    !,form="unformatted")

   open(215,file="z1_initial_distribution.dat")
   open(225,file="ze_initial_distribution.dat")

   do istep=1,number_of_timesteps
   
      histogram_z1 = 0.d0
      histogram_ze = 0.d0
   
      do itraj=1,number_of_traj
   
         z1_curr = z1(itraj,istep)
         ze_curr = ze(itraj,istep)
         
         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         histogram_z1(ibin_1) = histogram_z1(ibin_1) + 1.d0
   
         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)
         histogram_ze(ibin_e) = histogram_ze(ibin_e) + 1.d0
   
      enddo
   
      !-- normalize distributions
   
      do i1=1,number_of_bins_z1
         histogram_z1(i1) = histogram_z1(i1)/number_of_traj
      enddo
   
      do i1=1,number_of_bins_ze
         histogram_ze(i1) = histogram_ze(i1)/number_of_traj
      enddo
   
      !-- output to the external files for visualization
   
      do i1=1,number_of_bins_z1
         write(21,'(3g15.6)') time1(istep), bin_center_z1(i1), histogram_z1(i1)
         !write(21) bin_center_z1(i1), histogram_z1(i1)
      enddo
      write(21,*)
   
      do i1=1,number_of_bins_ze
         write(22,'(3g15.6)') time1(istep), bin_center_ze(i1), histogram_ze(i1)
         !write(22) bin_center_ze(i1), histogram_ze(i1)
      enddo
      write(22,*)
   
      if (istep.eq.1) then
         do i1=1,number_of_bins_z1
            write(215,'(2g15.6)') bin_center_z1(i1), histogram_z1(i1)
         enddo
         close(215)
         do i1=1,number_of_bins_ze
            write(225,'(2g15.6)') bin_center_ze(i1), histogram_ze(i1)
         enddo
         close(225)
      endif
      
   enddo
   
   close(21)
   close(22)
   
   deallocate (histogram_z1,histogram_ze)
   
   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !--(2)-- Time-dependent averages of solvent coordinates
   !        and corresponding velocities
   !---------------------------------------------------------------

   write(*,'(1x,"Building time-dependent averages... ",t64,"--> ",$)')
   time_start = secondi()

   allocate (z1_mean(number_of_timesteps))
   allocate (ze_mean(number_of_timesteps))
   allocate (vz1_mean(number_of_timesteps))
   allocate (vze_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      z1_mean(istep) = 0.d0
      ze_mean(istep) = 0.d0
      vz1_mean(istep) = 0.d0
      vze_mean(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_curr = z1(itraj,istep)
         ze_curr = ze(itraj,istep)
         vz1_curr = vz1(itraj,istep)
         vze_curr = vze(itraj,istep)

         z1_mean(istep) = z1_mean(istep) + z1_curr
         ze_mean(istep) = ze_mean(istep) + ze_curr
         vz1_mean(istep) = vz1_mean(istep) + vz1_curr
         vze_mean(istep) = vze_mean(istep) + vze_curr

      enddo

      z1_mean(istep) = z1_mean(istep)/number_of_traj
      ze_mean(istep) = ze_mean(istep)/number_of_traj
      vz1_mean(istep) = vz1_mean(istep)/number_of_traj
      vze_mean(istep) = vze_mean(istep)/number_of_traj

   enddo

   !-- output to the external files for visualization

   open(2,file="z_mean.dat")
   write(2,'("#",t10,"time(ps)",t30,"<z1>",t50,"<Ze>")')
   do istep=1,number_of_timesteps
       write(2,'(5g20.10)') time1(istep), z1_mean(istep), ze_mean(istep)
   enddo
   close(2)

   open(2,file="z1_phase_mean.dat")
   write(2,'("#",t10,"time(ps)",t30,"<z1>",t50,"<vz1>")')
   do istep=1,number_of_timesteps
       write(2,'(5g20.10)') time1(istep), z1_mean(istep), vz1_mean(istep)
   enddo
   close(2)

   open(2,file="ze_phase_mean.dat")
   write(2,'("#",t10,"time(ps)",t30,"<ze>",t50,"<vze>")')
   do istep=1,number_of_timesteps
       write(2,'(5g20.10)') time1(istep), ze_mean(istep), vze_mean(istep)
   enddo
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !--(3)-- Time-dependent second moments of solvent distribution
   !---------------------------------------------------------------

   write(*,'(1x,"Building time-dependent variances... ",t64,"--> ",$)')
   time_start = secondi()

   allocate (z1_var(number_of_timesteps))
   allocate (ze_var(number_of_timesteps))

   do istep=1,number_of_timesteps

      z1_var(istep) = 0.d0
      ze_var(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_curr = z1(itraj,istep) - z1_mean(istep)
         ze_curr = ze(itraj,istep) - ze_mean(istep)

         z1_var(istep) = z1_var(istep) + z1_curr*z1_curr
         ze_var(istep) = ze_var(istep) + ze_curr*ze_curr

      enddo

      z1_var(istep) = z1_var(istep)/number_of_traj
      ze_var(istep) = ze_var(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_var.dat")
   write(2,'("#",t10,"time(ps)",t30,"sigma(z1)",t50,"sigma(Ze)")')
   do istep=1,number_of_timesteps
       write(2,'(3g20.10)') time1(istep), sqrt(z1_var(istep)), sqrt(ze_var(istep))
   enddo
   close(2)

   deallocate (z1_var, ze_var)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !--------------------------------------------------------------------
   !--(4)-- Time-averages of solvent coordinates
   !--------------------------------------------------------------------

   write(*,'(1x,"Building time-averages for solvent coordinates... ",t64,"--> ",$)')
   time_start = secondi()

   z1tav = z1_mean(1)
   zetav = ze_mean(1)

   open(2,file="z_corr_tav.dat")
   write(2,'("#",t10,"time(ps)",t30,"<z1tav>",t50,"<zetav>")')
   write(2,'(7g20.10)') time1(1), z1tav, zetav

   do istep=2,number_of_timesteps
      z1tav = z1tav + z1_mean(istep)
      zetav = zetav + ze_mean(istep)
      write(2,'(7g20.10)') time1(istep), z1tav/istep, zetav/istep
   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !--------------------------------------------------------------------
   !--(5)-- Running (local) time-averages of solvent coordinates
   !--------------------------------------------------------------------

   write(*,'(1x,"Building running time-averages for solvent coordinates... ",t64,"--> ",$)')
   time_start = secondi()

   !-- time interval for averaging (number of timesteps),
   !   must be an odd number

   open(2,file="z_corr_rtav.dat")
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time1(2)-time1(1))
   write(2,'("#",t10,"time(ps)",t30,"<z1tav>",t50,"<zetav>")')

   do istep=1+ndt2,number_of_timesteps-ndt2

      z1tav  = 0.d0
      zetav  = 0.d0

      do i=istep-ndt2,istep+ndt2
         z1tav = z1tav + z1_mean(i)
         zetav = zetav + ze_mean(i)
      enddo

      z1tav  = z1tav/ndt
      zetav  = zetav/ndt
      write(2,'(7g20.10)') time1(istep), z1tav, zetav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !--(6)-- Non-equilibrium time-correlation functions
   !        of solvent coordinates
   !        (not clear what it means...)
   !---------------------------------------------------------------

   write(*,'(1x,"Building time-correlation functions... ",t64,"--> ",$)')
   time_start = secondi()

   allocate (z11_tcf(number_of_timesteps))
   allocate (zee_tcf(number_of_timesteps))

   do istep=1,number_of_timesteps

      z11_tcf(istep) = 0.d0
      zee_tcf(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_0 = z1(itraj,1) - z1_mean(1)
         ze_0 = ze(itraj,1) - ze_mean(1)

         z1_curr = z1(itraj,istep) - z1_mean(istep)
         ze_curr = ze(itraj,istep) - ze_mean(istep)

         z11_tcf(istep) = z11_tcf(istep) + z1_0*z1_curr
         zee_tcf(istep) = zee_tcf(istep) + ze_0*ze_curr

      enddo

      z11_tcf(istep) = z11_tcf(istep)/number_of_traj
      zee_tcf(istep) = zee_tcf(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_tcf.dat")
   write(2,'("#",t10,"time(ps)",t25,"<z1(0)z1(t)>",t45,"<ze(0)ze(t)>")')
   do istep=1,number_of_timesteps
       write(2,'(9g20.10)') time1(istep), z11_tcf(istep), zee_tcf(istep)
   enddo
   close(2)

   deallocate (z11_tcf)
   deallocate (zee_tcf)

   deallocate (z1, z1_mean, ze_mean)
   deallocate (vz1, vz1_mean, vze_mean)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !--(7)-- Time-dependent average free energy and kinetic energy
   !---------------------------------------------------------------

   write(*,'(1x,"Building average energies... ",t64,"--> ",$)')
   time_start = secondi()

   allocate (efe_mean(number_of_timesteps))
   allocate (ekin_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      efe_mean(istep) = 0.d0
      ekin_mean(istep) = 0.d0

      do itraj=1,number_of_traj
         efe_curr = efe(itraj,istep)
         ekin_curr = ekin(itraj,istep)
         efe_mean(istep) = efe_mean(istep) + efe_curr
         ekin_mean(istep) = ekin_mean(istep) + ekin_curr
      enddo

      efe_mean(istep) = efe_mean(istep)/number_of_traj
      ekin_mean(istep) = ekin_mean(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="ene_mean.dat")
   write(2,'("#",79("-"))')
   write(2,'("#",t10,"time",t25,"Free energy",t45,"Kinetic energy")')
   write(2,'("#",t10,"(ps)",t25,"(kcal/mol) ",t45,"  (kcal/mol)  ")')
   write(2,'("#",79("-"))')
   do istep=1,number_of_timesteps
      write(2,'(5g20.10)') time1(istep), efe_mean(istep), ekin_mean(istep)
   enddo
   close(2)

   deallocate (efe, efe_mean)
   deallocate (ekin, ekin_mean)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !---------------------------------------------------------------
   !--(8)-- Time-dependent EVB weights
   !---------------------------------------------------------------

   write(*,'(1x,"Building average EVB weights... ",t64,"--> ",$)')
   time_start = secondi()

   allocate (w1_mean(number_of_timesteps))
   allocate (w2_mean(number_of_timesteps))

   allocate (wh1_mean(number_of_timesteps))
   allocate (wh2_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      w1_mean(istep) = 0.d0
      w2_mean(istep) = 0.d0

      wh1_mean(istep) = 0.d0
      wh2_mean(istep) = 0.d0

      do itraj=1,number_of_traj

         w1_mean(istep) = w1_mean(istep) + w1(itraj,istep)
         w2_mean(istep) = w2_mean(istep) + w2(itraj,istep)

         !-- assign 1 for the largest weight (whXX arrays)
         tmparray(1) = w1(itraj,istep)
         tmparray(2) = w2(itraj,istep)

         ilargest = 1
         if (tmparray(2).gt.tmparray(1)) ilargest = 2

         select case(ilargest)
            case(1)
               wh1_mean(istep) = wh1_mean(istep) + 1.d0
            case(2)
               wh2_mean(istep) = wh2_mean(istep) + 1.d0
         end select

      enddo

      w1_mean(istep) = w1_mean(istep)/number_of_traj
      w2_mean(istep) = w2_mean(istep)/number_of_traj

      wh1_mean(istep) = wh1_mean(istep)/number_of_traj
      wh2_mean(istep) = wh2_mean(istep)/number_of_traj

   enddo

   !-- output to the external files for visualization

   open(2,file="weights_mean.dat")
   write(2,'("#",80("-"))')
   write(2,'("#   Average diabatic populations (expectation values)")')
   write(2,'("#",80("-"))')
   write(2,'("#",t10,"time",t30,"<1>",t50,"<2>")')
   write(2,'("#",80("-"))')
   do istep=1,number_of_timesteps
      write(2,'(3g20.10)') time1(istep), w1_mean(istep), w2_mean(istep)
   enddo
   close(2)

   open(2,file="weights_assigned_mean.dat")
   write(2,'("#",80("-"))')
   write(2,'("#   Average diabatic populations (values by assignment)")')
   write(2,'("#",80("-"))')
   write(2,'("#",t10,"time",t30,"<1>",t50,"<2>")')
   write(2,'("#",80("-"))')
   do istep=1,number_of_timesteps
      write(2,'(3g20.10)') time1(istep), wh1_mean(istep), wh2_mean(istep)
   enddo
   close(2)

   wh1_mean_0 = wh1_mean(1)
   w1_mean_0  = w1_mean(1)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   deallocate (w1)
   deallocate (w2)

   !-----------------------------------------------------------------------
   !--(9)-- Time-dependent Marcus/Rips-Jortner/Zusman diabatic populations
   !-----------------------------------------------------------------------

   if (marcus_flag) then

      write(*,'(1x,"Writing out Marcus populations... ",t64,"--> ",$)')
      time_start = secondi()

      open(2,file="marcus_diab_pop.dat")
      write(2,'("#",80("-"))')
      write(2,'("#   Marcus diabatic populations")')
      write(2,'("#   Nonadiabatic rate constant: ",g20.10," ps^(-1)")') k_marcus
      write(2,'("#",80("-"))')
      write(2,'("#",t10,"time",t30,"P(1)",t50,"P(2)")')
      write(2,'("#",80("-"))')
      do istep=1,number_of_timesteps
         p_marcus = exp_diabatic_populations(time1(istep),k_marcus,1.d0,0.d0)
         write(2,'(3g20.10)') time1(istep), p_marcus(1), p_marcus(2)
      enddo
      close(2)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start


      write(*,'(1x,"Writing out Rips-Jortner populations... ",t64,"--> ",$)')
      time_start = secondi()

      open(2,file="rips_jortner_diab_pop.dat")
      write(2,'("#",80("-"))')
      write(2,'("#   Rips-Jortner diabatic populations")')
      write(2,'("#   Rips-Jortner nonadiabatic rate constant: ",g20.10," ps^(-1)")') k_rips_jortner
      write(2,'("#",80("-"))')
      write(2,'("#",t10,"time",t30,"P(1)",t50,"P(2)")')
      write(2,'("#",80("-"))')
      do istep=1,number_of_timesteps
         p_rips_jortner = exp_diabatic_populations(time1(istep),k_rips_jortner,1.d0,0.d0)
         write(2,'(3g20.10)') time1(istep), p_rips_jortner(1), p_rips_jortner(2)
      enddo
      close(2)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start

      write(*,'(1x,"Writing out Zusman populations... ",t64,"--> ",$)')
      time_start = secondi()

      open(2,file="zusman_diab_pop.dat")
      write(2,'("#",80("-"))')
      write(2,'("#   Zusman diabatic populations")')
      write(2,'("#   Zusman nonadiabatic rate constant: ",g20.10," ps^(-1)")') k_zusman
      write(2,'("#",80("-"))')
      write(2,'("#",t10,"time",t30,"P(1)",t50,"P(2)")')
      write(2,'("#",80("-"))')
      do istep=1,number_of_timesteps
         p_zusman = exp_diabatic_populations(time1(istep),k_zusman,1.d0,0.d0)
         write(2,'(3g20.10)') time1(istep), p_zusman(1), p_zusman(2)
      enddo
      close(2)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start

      call allocated_memory

      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Electronic coupling:      ",e16.9," kcal/mol")') electronic_coupling
      write(*,'( 1x,"ET reaction free energy:  ",e16.9," kcal/mol")') reaction_free_energy
      write(*,'( 1x,"Reorganization energy:    ",e16.9," kcal/mol")') reorganization_energy
      write(*,'( 1x,"Activation energy:        ",e16.9," kcal/mol")') (reorganization_energy+reaction_free_energy)**2.d0/(4.d0*reorganization_energy)
      write(*,'( 1x,"ET Marcus rate constant:  ",e16.9," ps^(-1) ")') k_marcus
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Equilibrium constant:            ",g20.10)') k_equil
      write(*,'( 1x,"Equilibrium reactant population: ",g12.6)')  n_equil_r
      write(*,'( 1x,"Equilibrium product  population: ",g12.6)')  n_equil_p
      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Adiabaticity parameter:   ",e16.9)')             kappa_ad
      write(*,'( 1x,"ET R-J rate constant:     ",e16.9," ps^(-1) ")') k_rips_jortner
      write(*,'( 1x,"ET Zusman rate constant:  ",e16.9," ps^(-1) ")') k_zusman

      write(*,'(1x,"Fitting the rate constant (log regression)... ",t64,"--> ",$)')
      time_start = secondi()

      call fit_rate_constant_log(number_of_timesteps,time1(:),wh1_mean,k_fit_log,rcorr,shift,defect)

      open(2,file="fitlog_diab_pop.dat")
      write(2,'("#",80("-"))')
      write(2,'("#   Fitted diabatic populations (LOG linear regression)")')
      write(2,'("#   Fitted rate constant:        ",g20.10," ps^(-1)")') k_fit_log
      write(2,'("#   Correlation coefficient:     ",g20.10)') rcorr
      write(2,'("#   Initial condition shift:     ",g20.10)') shift
      write(2,'("#   Vertical shift (fit defect): ",g20.10)') defect
      write(2,'("#",80("-"))')
      write(2,'("#",t10,"time",t30,"Pfit(1)",t50,"Pfit(2)",t70,"Pfit0(1)",t90,"Pfit0(2)")')
      write(2,'("#",80("-"))')
      do istep=1,number_of_timesteps
         p_fit = exp_diabatic_populations(time1(istep),k_fit_log,wh1_mean_0,defect)
         p_fit0 = exp_diabatic_populations(time1(istep),k_fit_log,1.d0,0.d0)
         write(2,'(5g20.10)') time1(istep), p_fit(1), p_fit(2), p_fit0(1), p_fit0(2)
      enddo
      close(2)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start

      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"Log fit -  rate constant:  ",e16.9," ps^(-1) ")') k_fit_log
      write(*,'( 1x,"*correlation coefficient: ",g20.9)') rcorr
      write(*,'( 1x,"*Initial condition shift: ",g20.9)') shift
      write(*,'( 1x,"*vertical shift (defect): ",g20.9)') defect
      write(*,'( 1x,"------------------------------------------------------")')


      write(*,'(1x,"Fitting the rate constant (exponential ODR fit)... ",t64,"--> "/)')
      time_start = secondi()

      call fit_rate_constant_odr(number_of_timesteps,time1(:),wh1_mean,k_fit_odr,rcorr,defect)

      open(2,file="fitexp_diab_pop.dat")
      write(2,'("#",80("-"))')
      write(2,'("#   Fitted diabatic populations (ODR fit)")')
      write(2,'("#   Fitted rate constant:        ",g20.10," ps^(-1)")') k_fit_odr
      write(2,'("#   Error:                       ",g20.10)') rcorr
      write(2,'("#   Vertical shift (fit defect): ",g20.10)') defect
      write(2,'("#",80("-"))')
      write(2,'("#",t10,"time",t30,"Pfit(1)",t50,"Pfit(2)",t70,"Pfit0(1)",t90,"Pfit0(2)")')
      write(2,'("#",80("-"))')
      do istep=1,number_of_timesteps
         p_fit = exp_diabatic_populations(time1(istep),k_fit_odr,wh1_mean_0,defect)
         p_fit0 = exp_diabatic_populations(time1(istep),k_fit_odr,1.d0,0.d0)
         write(2,'(5g20.10)') time1(istep), p_fit(1), p_fit(2), p_fit0(1), p_fit0(2)
      enddo
      close(2)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start

      write(*,'( 1x,"------------------------------------------------------")')
      write(*,'( 1x,"ODR fit -  rate constant:  ",e16.9," ps^(-1) ")') k_fit_odr
      write(*,'( 1x,"*error:                   ",g20.9)') rcorr
      write(*,'( 1x,"*vertical shift (defect): ",g20.9)') defect
      write(*,'( 1x,"------------------------------------------------------")')


      open(2,file="marcus_and_fitted_rates.dat")
      write(2,'("#",t5,"V, kcal/mol",t26,"dG,kcal/mol",t45,"k_fit_log, 1/ps",t65,"k_fit_odr, 1/ps",t85,"k_marcus, 1/ps",t105,"k_Rips-Jortner",t125,"k_Zusman")')
      write(2,'(10g20.10)') electronic_coupling, reaction_free_energy, k_fit_log, k_fit_odr, k_marcus, k_rips_jortner, k_zusman
      close(2)

   endif

   call allocated_memory

   !------------------------------------------------------------------------------
   !--(10)-- Time-averaged EVB weights
   !------------------------------------------------------------------------------

   write(*,'(1x,"Building time-averaged EVB weights... ",t64,"--> ",$)')
   time_start = secondi()

   w1tav = w1_mean(1)
   w2tav = w2_mean(1)

   open(2,file="weights_tav.dat")
   write(2,'("#",80("-"))')
   write(2,'("#   Time averages of diabatic populations (expectation values)")')
   write(2,'("#",80("-"))')
   write(2,'("#",t10,"time(ps)",t30,"<w1tav>",t50,"<w2tav>")')
   write(2,'("#",80("-"))')
   write(2,'(11g20.10)') time1(1), w1tav, w2tav

   do istep=2,number_of_timesteps
      w1tav = w1tav + w1_mean(istep)
      w2tav = w2tav + w2_mean(istep)
      write(2,'(3g20.10)') time1(istep), w1tav/istep, w2tav/istep
   enddo

   close(2)

   !-- for mean weights calculated by assignment

   w1tav = wh1_mean(1)
   w2tav = wh2_mean(1)

   open(2,file="weights_assigned_corr_tav.dat")
   write(2,'("#",t10,"time(ps)",t30,"<w1tav>",t50,"<w2tav>")')
   write(2,'(3g20.10)') time1(1), w1tav, w2tav

   do istep=2,number_of_timesteps
      w1tav = w1tav + wh1_mean(istep)
      w2tav = w2tav + wh2_mean(istep)
      write(2,'(3g20.10)') time1(istep), w1tav/istep, w2tav/istep
   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   !--------------------------------------------------------------------
   !--(11)-- Running (local) time-averages of EVB weights
   !--------------------------------------------------------------------

   write(*,'(1x,"Building running time-averages for EVB weights... ",t64,"--> ",$)')
   time_start = secondi()

   !-- time interval for averaging (number of timesteps),
   !   must be an odd number

   open(2,file="weights_corr_rtav.dat")
   write(2,'("#",80("-"))')
   write(2,'("#    Running time averages of diabatic populations (expectation values)")')
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time1(2)-time1(1))
   write(2,'("#",80("-"))')
   write(2,'("#",t10,"time(ps)",t30,"<w1tav>",t50,"<w2tav>")')
   write(2,'("#",80("-"))')

   do istep=1+ndt2,number_of_timesteps-ndt2

      w1tav  = 0.d0
      w2tav  = 0.d0

      do i=istep-ndt2,istep+ndt2
         w1tav = w1tav + w1_mean(i)
         w2tav = w2tav + w2_mean(i)
      enddo

      w1tav = w1tav/ndt
      w2tav = w2tav/ndt

      write(2,'(3g20.10)') time1(istep), w1tav, w2tav

   enddo

   close(2)


   !-- for mean weights calculated by assignment

   open(2,file="weights_assigned_corr_rtav.dat")
   write(2,'("#",80("-"))')
   write(2,'("#    Running time averages of diabatic populations (by assignment)")')
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time1(2)-time1(1))
   write(2,'("#",80("-"))')
   write(2,'("#",t10,"time(ps)",t30,"<w1tav>",t50,"<w2tav>")')
   write(2,'("#",80("-"))')

   do istep=1+ndt2,number_of_timesteps-ndt2

      w1tav  = 0.d0
      w2tav  = 0.d0

      do i=istep-ndt2,istep+ndt2
         w1tav = w1tav + wh1_mean(i)
         w2tav = w2tav + wh2_mean(i)
      enddo

      w1tav = w1tav/ndt
      w2tav = w2tav/ndt

      write(2,'(3g20.10)') time1(istep), w1tav, w2tav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   deallocate (w1_mean,wh1_mean)
   deallocate (w2_mean,wh2_mean)

   !---------------------------------------------------------------
   !--(10)-- Time-dependent adiabatic populations
   !---------------------------------------------------------------

   write(*,'(1x,"Building time-dependent adiabatic populations... ",t64,"--> ",$)')
   time_start = secondi()

   allocate(pop_ad(number_of_states,number_of_timesteps))

   pop_ad = 0.d0

   do istep=1,number_of_timesteps

      do itraj=1,number_of_traj
         iocc = istate(itraj,istep)
         pop_ad(iocc,istep) = pop_ad(iocc,istep) + 1.d0
      enddo

      do i=1,number_of_states
         pop_ad(i,istep) = pop_ad(i,istep)/number_of_traj
      enddo

   enddo

   !-- output to the external file for visualization

   open(2,file="pop_ad.dat")
   write(2,'("#",40("-"))')
   write(2,'("#   Adiabatic populations")')
   write(2,'("#",40("-"))')
   write(2,'("#",t5,"time(ps)",t19,"<1>",t31,"<2>")')
   write(2,'("#",40("-"))')
   do istep=1,number_of_timesteps
      write(2,'(101f12.6)') time1(istep), (pop_ad(k,istep),k=1,number_of_states)
   enddo
   close(2)

   deallocate (pop_ad)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start

   call allocated_memory

   !---------------------------------------------------------------
   !--(11)-- Calculate Mean First Passage Time
   !---------------------------------------------------------------

   write(*,'(1x,"Calculating Mean First Passage Time... ",t64,"--> ",$)')
   time_start = secondi()

   ze_crossing = 0.d0
   ze_product = -reorganization_energy + reaction_free_energy
   dze = (ze_product - ze_crossing)/10.d0

   open(2,file="mean_first_passage_time.dat")
   write(2,'("#",60("-"))')
   write(2,'("#   Mean First Passage Time (MFPT) and MFPT rate constant")')
   write(2,'("#",60("-"))')
   write(2,'("#",t5,"ze_cut(kcal/mol)",t23,"MFPT (ps)",t38,"k_MFPT (ps^-1)")')
   write(2,'("#",40("-"))')

   do ze_cut=ze_crossing,ze_product,dze

      number_of_reactive_events = 0
      time_mfpt = 0.d0

      do itraj=1,number_of_traj
         do istep=1,number_of_timesteps
            if (ze(itraj,istep).lt.ze_cut.and.istate(itraj,istep).eq.1) then
               time_mfpt = time_mfpt + time1(istep)
               number_of_reactive_events = number_of_reactive_events + 1
               exit
            endif
         enddo
      enddo

      if (number_of_reactive_events.gt.0) then
         time_mfpt = time_mfpt/number_of_reactive_events
         write(2,'(3f15.6)') ze_cut, time_mfpt, 1.d0/time_mfpt
      else
         write(2,'("#",f14.6,"  -- no reactive events --")') ze_cut
      endif

   enddo

   close(2)

   deallocate (istate)

   time_end = secondi()
   write(*,'("Done in ",f10.3," sec")') time_end-time_start


   !---------------------------------------------------------------
   !--(12)-- Calculate equilibrium correlation function of the
   !         reactant Heaviside function with a dividing surface
   !         ze* defined as ze* = ze(crossing) = 0
   !---------------------------------------------------------------

   if (calculate_reactive_flux_rate) then

      call allocated_memory

      write(*,'(/1x,"Calculating reactive flux correlation function... ")')
      time_start = secondi()

      z_divide = z_crossing

      allocate (theta_corr(number_of_timesteps))

      !-- calculate global averages (timesteps AND trajectories)
      !   and compare it to equilibrium averages

      n_aver_r = 0.d0
      n_aver_p = 0.d0
      do itraj=1,number_of_traj
         do istep=1,number_of_timesteps
            ze_curr = ze(itraj,istep)
            n_aver_r = n_aver_r + reactant_heaviside(ze_curr,z_divide)
            n_aver_p = n_aver_p + 1.d0 - reactant_heaviside(ze_curr,z_divide)
         enddo
      enddo
      n_aver_r = n_aver_r/(number_of_timesteps*number_of_traj)
      n_aver_p = n_aver_p/(number_of_timesteps*number_of_traj)

      write(*,'( 1x,"Simulated equilibrium populations: ",2g12.6)')  n_aver_r, n_aver_p
      write(*,'( 1x,"Exact     equilibrium populations: ",2g12.6)')  n_equil_r, n_equil_p

      !-- calculate correlation function

      do idstep=1,number_of_timesteps/2
         theta_corr(idstep) = 0.d0
         do itraj=1,number_of_traj
            do istep=1,number_of_timesteps+1-idstep
               ze_curr1 = ze(itraj,idstep)
               ze_curr2 = ze(itraj,idstep+istep-1)
               theta_corr(idstep) = theta_corr(idstep) + reactant_heaviside(ze_curr1,z_divide)*reactant_heaviside(ze_curr2,z_divide)
            enddo
         enddo
      enddo

      open(2,file="theta_theta_tcf.dat")
      write(2,'("#",60("-"))')
      write(2,'("#   Equilibrium TCF of the reactant population")')
      write(2,'("#   (its time derivative is the time-dependent rate constant)")')
      write(2,'("#",60("-"))')
      write(2,'("#",t10,"time(ps)",t30,"<d_theta[z(0)]*d_theta[z(t)]>/<theta>")')
      write(2,'("#",60("-"))')

      !-- normalise the correlation function by <theta>=n_equil_r and output the result
      do idstep=1,number_of_timesteps
         theta_corr(idstep) = theta_corr(idstep)/real(number_of_timesteps+1-idstep)/real(number_of_traj)
         write(2,'(2g20.10)') time1(idstep), (theta_corr(idstep) - n_aver_r*n_aver_r)/n_aver_r
      enddo

      deallocate (theta_corr)

      time_end = secondi()
      write(*,'("Done in ",f10.3," sec")') time_end-time_start

   endif

   call deallocate_all_arrays
   call allocated_memory

   total_time_end = secondi()
   write(*,'(/"===> All Done in ",f10.3," sec"/)') total_time_end-total_time_start

contains

   subroutine deallocate_all_arrays
      if (allocated(time)) deallocate(time)
      if (allocated(time1)) deallocate(time1)
      if (allocated(z1)) deallocate(z1)
      if (allocated(ze)) deallocate(ze)
      if (allocated(ekin)) deallocate(ekin)
      if (allocated(efe)) deallocate(efe)
      if (allocated(istate)) deallocate(istate)
      if (allocated(istate_occ)) deallocate(istate_occ)
      if (allocated(all_states)) deallocate(all_states)
      if (allocated(w1)) deallocate(w1)
      if (allocated(w2)) deallocate(w2)
      if (allocated(histogram_z1)) deallocate(histogram_z1)
      if (allocated(histogram_ze)) deallocate(histogram_ze)
      if (allocated(state_histogram_z1)) deallocate(state_histogram_z1)
      if (allocated(state_histogram_ze)) deallocate(state_histogram_ze)
      if (allocated(z1_mean)) deallocate(z1_mean)
      if (allocated(ze_mean)) deallocate(ze_mean)
      if (allocated(z1_var)) deallocate(z1_var)
      if (allocated(ze_var)) deallocate(ze_var)
      if (allocated(z11_tcf)) deallocate(z11_tcf)
      if (allocated(zee_tcf)) deallocate(zee_tcf)
      if (allocated(theta_corr)) deallocate(theta_corr)
      if (allocated(ekin_mean)) deallocate(ekin_mean)
      if (allocated(efe_mean)) deallocate(efe_mean)
      if (allocated(w1_mean)) deallocate(w1_mean)
      if (allocated(w2_mean)) deallocate(w2_mean)
      if (allocated(wh1_mean)) deallocate(wh1_mean)
      if (allocated(wh2_mean)) deallocate(wh2_mean)
      if (allocated(pop_ad)) deallocate(pop_ad)
   end subroutine deallocate_all_arrays


   subroutine allocated_memory

      real(kind=8) :: memory, memory_mb, memory_gb

      memory = 0.d0

      if (allocated(istate))             memory = memory + number_of_traj*number_of_timesteps*4.d0
      if (allocated(all_states))         memory = memory + number_of_traj*number_of_timesteps*4.d0
      if (allocated(istate_occ))         memory = memory + number_of_states*4.d0

      if (allocated(time))               memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(z1))                 memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(vz1))                memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(ze))                 memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(vze))                memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(ekin))               memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(efe))                memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(w1))                 memory = memory + number_of_traj*number_of_timesteps*8.d0
      if (allocated(w2))                 memory = memory + number_of_traj*number_of_timesteps*8.d0

      if (allocated(pop_ad))             memory = memory + number_of_states*number_of_timesteps*8.d0
      if (allocated(state_histogram_z1)) memory = memory + number_of_bins_z1*number_of_timesteps*8.d0
      if (allocated(state_histogram_ze)) memory = memory + number_of_bins_ze*number_of_timesteps*8.d0

      if (allocated(histogram_z1))       memory = memory + number_of_bins_z1*8.d0
      if (allocated(histogram_ze))       memory = memory + number_of_bins_ze*8.d0

      if (allocated(z1_mean))            memory = memory + number_of_timesteps*8.d0
      if (allocated(ze_mean))            memory = memory + number_of_timesteps*8.d0
      if (allocated(vz1_mean))           memory = memory + number_of_timesteps*8.d0
      if (allocated(vze_mean))           memory = memory + number_of_timesteps*8.d0
      if (allocated(z1_var))             memory = memory + number_of_timesteps*8.d0
      if (allocated(ze_var))             memory = memory + number_of_timesteps*8.d0
      if (allocated(z11_tcf))            memory = memory + number_of_timesteps*8.d0
      if (allocated(zee_tcf))            memory = memory + number_of_timesteps*8.d0
      if (allocated(efe_mean))           memory = memory + number_of_timesteps*8.d0
      if (allocated(ekin_mean))          memory = memory + number_of_timesteps*8.d0
      if (allocated(w1_mean))            memory = memory + number_of_timesteps*8.d0
      if (allocated(wh1_mean))           memory = memory + number_of_timesteps*8.d0
      if (allocated(w2_mean))            memory = memory + number_of_timesteps*8.d0
      if (allocated(wh2_mean))           memory = memory + number_of_timesteps*8.d0
      if (allocated(theta_corr))         memory = memory + number_of_timesteps*8.d0

      memory_mb = memory/1024.d0/1024.d0
      memory_gb = memory/1024.d0/1024.d0/1024.d0

      write(*,'(/1x,"Current memory allocation: ",f20.3," MB = ",f8.3," GB"/)') memory_mb, memory_gb

   end subroutine allocated_memory


   function reactant_heaviside(z_, z_divide_) result(out)
      implicit none
      real(kind=8), intent(in) :: z_, z_divide_
      real(kind=8) :: out
      real(kind=8) :: deltaz
      deltaz = z_ - z_divide_
      if (deltaz.lt.0.d0) then
         out = 0.d0
      else
         out = 1.d0
      endif
   end function reactant_heaviside


end program analyze_et2_trajectories
