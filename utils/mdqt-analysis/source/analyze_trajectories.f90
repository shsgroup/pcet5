!======================================================================================
!  Analysis of the dynamical trajectories
!  Alexander V. Soudackov, Penn State
!  March 14, 2011
!======================================================================================

program analyze_trajectories

   use string_utilities
   use sorting
   use timers
   implicit none

   integer, parameter :: number_of_bins_z12=50
   integer, parameter :: number_of_bins_zpe=50

   character(len=60)  :: filename
   character(len=200) :: record
   character(len=3)   :: state_suffix

   logical :: found

   integer :: number_of_traj, number_of_timesteps, number_of_states
   integer ::  number_of_occ_states
   integer :: itraj, istep, eof=0
   integer :: i, ii, iocc, k, i1, i2, iargc, nsteps
   integer :: ibin_1, ibin_2, ibin_p, ibin_e
   integer :: ichannel_z12, ichannel_zpe

   real(kind=8) :: z1_min, z1_max, z2_min, z2_max
   real(kind=8) :: zp_min, zp_max, ze_min, ze_max
   real(kind=8) :: z1_curr, z2_curr, zp_curr, ze_curr
   real(kind=8) :: z1_0, z2_0, zp_0, ze_0, efe_curr, ekin_curr
   real(kind=8) :: time_start, time_end, total_time_start, total_time_end

   real(kind=8) :: bin_width_z1, bin_width_z2, bin_width_zp, bin_width_ze

   real(kind=8) :: cd1a1a, cd1b1b, cd2a2a, cd2b2b
   real(kind=8) :: cc1a1b, cc1a2a, cc1a2b, cc1b2a, cc1b2b, cc2a2b
   real(kind=8) :: sq1a, sq1b, sq2a, sq2b

   real(kind=8), dimension(number_of_bins_z12) :: bin_center_z1, bin_center_z2
   real(kind=8), dimension(number_of_bins_zpe) :: bin_center_zp, bin_center_ze

   real(kind=8), dimension(:,:),   allocatable :: time
   real(kind=8), dimension(:,:),   allocatable :: z1, z2
   real(kind=8), dimension(:,:),   allocatable :: zp, ze
   real(kind=8), dimension(:,:),   allocatable :: ekin, efe
   real(kind=8), dimension(:,:),   allocatable :: w1a, w1b, w2a, w2b
   real(kind=8), dimension(:,:),   allocatable :: histogram_z12
   real(kind=8), dimension(:,:),   allocatable :: histogram_zpe
   real(kind=8), dimension(:,:,:), allocatable :: state_histogram_z12
   real(kind=8), dimension(:,:,:), allocatable :: state_histogram_zpe
   real(kind=8), dimension(:),     allocatable :: z1_mean, z2_mean, zp_mean, ze_mean
   real(kind=8), dimension(:),     allocatable :: z1_var, z2_var, z12_var, zp_var, ze_var, zpe_var
   real(kind=8), dimension(:),     allocatable :: z11_tcf, z22_tcf, z12_tcf, zpp_tcf, zee_tcf, zpe_tcf
   real(kind=8), dimension(:),     allocatable :: efe_mean, ekin_mean
   real(kind=8), dimension(:),     allocatable :: w1a_mean, w1b_mean, w2a_mean, w2b_mean
   real(kind=8), dimension(:,:),   allocatable :: cw_cross
   real(kind=8), dimension(:,:),   allocatable :: pop_ad

   integer, dimension(:,:), allocatable :: istate
   integer, dimension(:), allocatable :: all_states
   integer, dimension(:), allocatable :: istate_occ

   character(len=20), dimension(100) :: carr
   integer,           dimension(100) :: iarr
   real(kind=8),      dimension(100) :: rarr

!--------------------------------------------------------------------------------------

   total_time_start = secondi()

   number_of_traj = iargc()
   if (number_of_traj == 0) then
      write(*,*) "No arguments... Come again..."
      stop
   endif

   !-- scan the first trajectory file
   !   to get the number of time steps

   call getarg(1,filename)
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

   if (number_of_timesteps.eq.0) then
      write(*,*)
      write(*,*) "The first trajectory file is empty, Something is wrong..."
      write(*,*)
      stop
   endif

   !-- allocate the data arrays

   allocate(time(number_of_traj,number_of_timesteps))
   allocate(z1(number_of_traj,number_of_timesteps))
   allocate(z2(number_of_traj,number_of_timesteps))
   allocate(zp(number_of_traj,number_of_timesteps))
   allocate(ze(number_of_traj,number_of_timesteps))
   allocate(ekin(number_of_traj,number_of_timesteps))
   allocate(efe(number_of_traj,number_of_timesteps))
   allocate(istate(number_of_traj,number_of_timesteps))
   allocate(w1a(number_of_traj,number_of_timesteps))
   allocate(w1b(number_of_traj,number_of_timesteps))
   allocate(w2a(number_of_traj,number_of_timesteps))
   allocate(w2b(number_of_traj,number_of_timesteps))

   write(*,*)
   write(*,'(1x,"All data arrays allocated, now reading data files... ",$)')
   time_start = secondi()

   !=============================================================
   !-- reading all trajectory files and fill the data arrays
   !=============================================================

   loop_over_trajectories: do itraj=1,number_of_traj

      call getarg(itraj,filename)
      filename = trim(filename)

      open(1,file=filename)

      nsteps = 0

      loop_over_timesteps: do

         read(1,'(a)',iostat=eof) record

         !-- end-of-file condition
         if (eof.ne.0) exit

         !-- skip comments and empty lines
         if (record(1:1).eq."#".or.record.eq."") cycle

         call readrec(record,"rrrrrrrrrrrirrrr",carr,iarr,rarr)

         nsteps = nsteps + 1

         if (nsteps.gt.number_of_timesteps) then
            call deallocate_all_arrays
            write(*,*)
            write(*,*) "Trajectory file ",trim(filename)," has inconsistent number of timesteps."
            write(*,*) "Something is wrong..."
            write(*,*)
            stop
         endif

         time(itraj,nsteps) = rarr(1)

         z1(itraj,nsteps) = rarr(2)
         z2(itraj,nsteps) = rarr(3)

         zp(itraj,nsteps) = rarr(6)
         ze(itraj,nsteps) = rarr(7)

         ekin(itraj,nsteps) = rarr(10)
         efe(itraj,nsteps) = rarr(11)

         istate(itraj,nsteps) = iarr(1)
 
         w1a(itraj,nsteps) = rarr(12)
         w1b(itraj,nsteps) = rarr(13)
         w2a(itraj,nsteps) = rarr(14)
         w2b(itraj,nsteps) = rarr(15)

      enddo loop_over_timesteps

   enddo loop_over_trajectories

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !=============================================================
   !-- Write the time values (ps) to the external file
   !=============================================================

   open(1,file="timesteps.dat",form="formatted")
   do i=1,number_of_timesteps
      write(1,'(i10,2x,f20.6)') i, time(1,i)
   enddo
   close(1)

   !=============================================================
   !-- Scan the occupied state array and determine the highest
   !   occupied state for the whole set of trajectories.
   !=============================================================

   write(*,*)
   write(*,'(1x,"Scanning MDQT trajectories to identify unique occupied states... ")')
   time_start = secondi()

   number_of_states = 0

   do itraj=1,number_of_traj
      do istep=1,number_of_timesteps
         iocc = istate(itraj,istep)
         if (iocc.gt.number_of_states) number_of_states = iocc
      enddo
   enddo

   !=============================================================
   !-- Scan the occupied state array and create a list of states
   !   occupied at least once for the whole set of trajectories.
   !=============================================================

   number_of_occ_states = 0

   !-- unwrap 2D array and store in the temporary array

   allocate (all_states(number_of_traj*number_of_timesteps))

   k = 0

   do itraj=1,number_of_traj
      do istep=1,number_of_timesteps
         k = k + 1
         all_states(k) = istate(itraj,istep)
      enddo
   enddo

   !-- find the number of unique occupied states

   do i=1,number_of_traj*number_of_timesteps
      found = .false.
      do k=i+1,number_of_traj*number_of_timesteps
         if (all_states(i).eq.all_states(k)) found = .true.
      enddo
      if (.not.found) then
         number_of_occ_states = number_of_occ_states + 1
      endif
   enddo

   !-- fill the array of unique occupied states

   allocate(istate_occ(number_of_occ_states))
   istate_occ = 0

   ii = 0
   do i=1,number_of_traj*number_of_timesteps
      found = .false.
      do k=i+1,number_of_traj*number_of_timesteps
         if (all_states(i).eq.all_states(k)) found = .true.
      enddo
      if (.not.found) then
         ii = ii + 1
         istate_occ(ii) = all_states(i)
      endif
   enddo

   !-- release temporary array
   deallocate(all_states)

   !-- sort the array istate_occ

   call shellsort(number_of_occ_states,istate_occ)

   write(*,'(1x,"Highest occupied adiabatic state: ",i3/)') number_of_states
   write(*,'(1x,"Unique occupied states:",/,(10(1x,i3)))') istate_occ

   time_end = secondi()
   write(*,'(/1x,"Done in ",f12.3," sec"/)') time_end-time_start

   !======================================
   !-- Start calculating the observables
   !======================================

   !---------------------------------------------------------------
   !--(1)-- Time-dependent distributions of solvent coordinates
   !---------------------------------------------------------------

   !-- find the global ranges for solvent coordinates

   z1_min =  999999.d0
   z1_max = -999999.d0
   z2_min =  999999.d0
   z2_max = -999999.d0
   zp_min =  999999.d0
   zp_max = -999999.d0
   ze_min =  999999.d0
   ze_max = -999999.d0

   do itraj=1,number_of_traj
      do istep=1,number_of_timesteps

         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)

         if (z1_curr.lt.z1_min) z1_min = z1_curr
         if (z1_curr.gt.z1_max) z1_max = z1_curr
         if (z2_curr.lt.z2_min) z2_min = z2_curr
         if (z2_curr.gt.z2_max) z2_max = z2_curr

         if (zp_curr.lt.zp_min) zp_min = zp_curr
         if (zp_curr.gt.zp_max) zp_max = zp_curr
         if (ze_curr.lt.ze_min) ze_min = ze_curr
         if (ze_curr.gt.ze_max) ze_max = ze_curr

      enddo
   enddo

   !-- shift min and max by very small amounts to ensure that
   !   binning will always put the values within bounds

   z1_min = z1_min - 1.d-10
   z2_min = z2_min - 1.d-10
   zp_min = zp_min - 1.d-10
   ze_min = ze_min - 1.d-10

   z1_max = z1_max + 1.d-10
   z2_max = z2_max + 1.d-10
   zp_max = zp_max + 1.d-10
   ze_max = ze_max + 1.d-10

   !-- define the histogram parameters

   bin_width_z1 = (z1_max - z1_min)/number_of_bins_z12
   bin_width_z2 = (z2_max - z2_min)/number_of_bins_z12

   bin_width_zp = (zp_max - zp_min)/number_of_bins_zpe
   bin_width_ze = (ze_max - ze_min)/number_of_bins_zpe

   do i=1,number_of_bins_z12
      bin_center_z1(i) = z1_min + (i-1)*bin_width_z1 + bin_width_z1/2.d0
      bin_center_z2(i) = z2_min + (i-1)*bin_width_z2 + bin_width_z2/2.d0
   enddo

   do i=1,number_of_bins_zpe
      bin_center_zp(i) = zp_min + (i-1)*bin_width_zp + bin_width_zp/2.d0
      bin_center_ze(i) = ze_min + (i-1)*bin_width_ze + bin_width_ze/2.d0
   enddo

   !---------------------------------------------------------------
   !-- build histograms for occupied adiabatic states
   !---------------------------------------------------------------

   write(*,*)
   write(*,'(1x,"Building state-resolved histograms for solvent coordinates... ",$)')
   time_start = secondi()
   
   allocate(state_histogram_z12(number_of_states,number_of_bins_z12,number_of_bins_z12))
   allocate(state_histogram_zpe(number_of_states,number_of_bins_zpe,number_of_bins_zpe))

   !open(21,file="z12_distribution_state.dat",form="formatted")
   !open(22,file="zpe_distribution_state.dat",form="formatted")

   !-- open a separate file for each occupied state

   do i=1,number_of_occ_states
      ichannel_z12 = 10 + 2*i - 1
      ichannel_zpe = 10 + 2*i
      write(state_suffix,'(i3.3)') istate_occ(i)
      open(ichannel_z12,file="z12_distribution_"//state_suffix//".dat",form="unformatted")
      open(ichannel_zpe,file="zpe_distribution_"//state_suffix//".dat",form="unformatted")
   enddo

   do istep=1,number_of_timesteps
   
      state_histogram_z12 = 0.d0
      state_histogram_zpe = 0.d0
   
      do itraj=1,number_of_traj
   
         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)
         
         iocc = istate(itraj,istep)
         
         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         ibin_2 = nint((z2_curr-z2_min)/bin_width_z2 + 0.5d0)

         !--DEBUG start--
         !if (ibin_1.gt.number_of_bins_z12) then
         !   write(*,*)
         !   write(*,*) "Problem with binning: ibin_1 =",ibin_1
         !   write(*,*) "z1_curr-z1_min =",z1_curr-z1_min
         !   write(*,*) "z1_max-z1_min  =",z1_max-z1_min
         !   call deallocate_all_arrays
         !   stop
         !endif
         !if (ibin_2.gt.number_of_bins_z12) then
         !   write(*,*)
         !   write(*,*) "Problem with binning: ibin_2 =",ibin_2
         !   write(*,*) "z2_curr-z2_min =",z2_curr-z2_min
         !   write(*,*) "z2_max-z2_min  =",z2_max-z2_min
         !   call deallocate_all_arrays
         !   stop
         !endif
         !--DEBUG end--

         state_histogram_z12(iocc,ibin_1,ibin_2) = state_histogram_z12(iocc,ibin_1,ibin_2) + 1.d0
    
         ibin_p = nint((zp_curr-zp_min)/bin_width_zp + 0.5d0)
         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)

         !--DEBUG start--
         !if (ibin_p.gt.number_of_bins_zpe) then
         !   write(*,*)
         !   write(*,*) "Problem with binning: ibin_p =",ibin_p
         !   write(*,*) "zp_curr-zp_min =",zp_curr-zp_min
         !   write(*,*) "zp_max-zp_min  =",zp_max-zp_min
         !   call deallocate_all_arrays
         !   stop
         !endif
         !if (ibin_e.gt.number_of_bins_zpe) then
         !   write(*,*)
         !   write(*,*) "Problem with binning: ibin_e =",ibin_e
         !   write(*,*) "ze_curr-ze_min =",ze_curr-ze_min
         !   write(*,*) "ze_max-ze_min  =",ze_max-ze_min
         !   call deallocate_all_arrays
         !   stop
         !endif
         !--DEBUG end--

         state_histogram_zpe(iocc,ibin_p,ibin_e) = state_histogram_zpe(iocc,ibin_p,ibin_e) + 1.d0
   
      enddo
   
      !-- normalize distributions
   
      do i1=1,number_of_bins_z12
         do i2=1,number_of_bins_z12
            do k=1,number_of_states
               state_histogram_z12(k,i1,i2) = state_histogram_z12(k,i1,i2)/number_of_traj
            enddo
         enddo
      enddo
   
      do i1=1,number_of_bins_zpe
         do i2=1,number_of_bins_zpe
            do k=1,number_of_states
               state_histogram_zpe(k,i1,i2) = state_histogram_zpe(k,i1,i2)/number_of_traj
            enddo
         enddo
      enddo
   
      !-- output to the external files for visualization
      !
      !do i1=1,number_of_bins_z12
      !   do i2=1,number_of_bins_z12
      !      write(21,'(f10.6,1x,f10.6,100(1x,i3,1x,g12.6))') bin_center_z1(i1), bin_center_z2(i2), (istate_occ(k),state_histogram_z12(istate_occ(k),i1,i2),k=1,number_of_occ_states)
      !      !write(21) bin_center_z1(i1), bin_center_z2(i2), (istate_occ(k),state_histogram_z12(istate_occ(k),i1,i2),k=1,number_of_occ_states)
      !   enddo
      !enddo
      !
      !do i1=1,number_of_bins_zpe
      !   do i2=1,number_of_bins_zpe
      !      write(22,'(f10.6,1x,f10.6,100(1x,i3,1x,g12.6))') bin_center_z1(i1), bin_center_z2(i2), (istate_occ(k),state_histogram_zpe(istate_occ(k),i1,i2),k=1,number_of_occ_states)
      !      !write(22) bin_center_z1(i1), bin_center_z2(i2), (istate_occ(k),state_histogram_zpe(istate_occ(k),i1,i2),k=1,number_of_occ_states)
      !   enddo
      !enddo


      !-- output to the external files for visualization

      do i=1,number_of_occ_states

         ichannel_z12 = 10 + 2*i - 1
         ichannel_zpe = 10 + 2*i

         do i1=1,number_of_bins_z12
            do i2=1,number_of_bins_z12
               !write(ichannel_z12,'(f10.6,1x,f10.6,1x,g15.8)') bin_center_z1(i1), bin_center_z2(i2), state_histogram_z12(istate_occ(i),i1,i2)
               write(ichannel_z12) bin_center_z1(i1), bin_center_z2(i2), state_histogram_z12(istate_occ(i),i1,i2)
            enddo
         enddo

         do i1=1,number_of_bins_zpe
            do i2=1,number_of_bins_zpe
               !write(ichannel_zpe,'(f10.6,1x,f10.6,1x,g15.8)') bin_center_zp(i1), bin_center_ze(i2), state_histogram_zpe(istate_occ(i),i1,i2)
               write(ichannel_zpe) bin_center_zp(i1), bin_center_ze(i2), state_histogram_zpe(istate_occ(i),i1,i2)
            enddo
         enddo

      enddo
  
   enddo
   
   !close(21)
   !close(22)

   do i=1,number_of_occ_states
      ichannel_z12 = 10 + 2*i - 1
      ichannel_zpe = 10 + 2*i
      close(ichannel_z12)
      close(ichannel_zpe)
   enddo
   
   deallocate (state_histogram_z12,state_histogram_zpe)
   deallocate (istate_occ)
   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !-- build global histograms
   !---------------------------------------------------------------

   write(*,*)
   write(*,'(1x,"Building global histograms for solvent coordinates... ",$)')
   time_start = secondi()
   
   allocate(histogram_z12(number_of_bins_z12,number_of_bins_z12))
   allocate(histogram_zpe(number_of_bins_zpe,number_of_bins_zpe))
   
   open(21,file="z12_distribution_global.dat",form="unformatted")
   open(22,file="zpe_distribution_global.dat",form="unformatted")
   
   do istep=1,number_of_timesteps
   
      histogram_z12 = 0.d0
      histogram_zpe = 0.d0
   
      do itraj=1,number_of_traj
   
         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)
         
         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         ibin_2 = nint((z2_curr-z2_min)/bin_width_z2 + 0.5d0)
         histogram_z12(ibin_1,ibin_2) = histogram_z12(ibin_1,ibin_2) + 1.d0
   
         ibin_p = nint((zp_curr-zp_min)/bin_width_zp + 0.5d0)
         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)
         histogram_zpe(ibin_p,ibin_e) = histogram_zpe(ibin_p,ibin_e) + 1.d0
   
      enddo
   
      !-- normalize distributions
   
      do i1=1,number_of_bins_z12
         do i2=1,number_of_bins_z12
            histogram_z12(i1,i2) = histogram_z12(i1,i2)/number_of_traj
         enddo
      enddo
   
      do i1=1,number_of_bins_zpe
         do i2=1,number_of_bins_zpe
            histogram_zpe(i1,i2) = histogram_zpe(i1,i2)/number_of_traj
         enddo
      enddo
   
      !-- output to the external files for visualization
   
      do i1=1,number_of_bins_z12
         do i2=1,number_of_bins_z12
            !write(21,'(3g15.6)') bin_center_z1(i1), bin_center_z2(i2), histogram_z12(i1,i2)
            write(21) bin_center_z1(i1), bin_center_z2(i2), histogram_z12(i1,i2)
         enddo
      enddo
   
      do i1=1,number_of_bins_zpe
         do i2=1,number_of_bins_zpe
            !write(22,'(3g15.6)') bin_center_z1(i1), bin_center_z2(i2), histogram_zpe(i1,i2)
            write(22) bin_center_z1(i1), bin_center_z2(i2), histogram_zpe(i1,i2)
         enddo
      enddo
   
   enddo
   
   close(21)
   close(22)
   
   deallocate (histogram_z12,histogram_zpe)
   
   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(2)-- Time-dependent averages of solvent coordinates
   !---------------------------------------------------------------

   write(*,'(/1x,"Building time-dependent averages... ",$)')
   time_start = secondi()

   allocate (z1_mean(number_of_timesteps))
   allocate (z2_mean(number_of_timesteps))
   allocate (zp_mean(number_of_timesteps))
   allocate (ze_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      z1_mean(istep) = 0.d0
      z2_mean(istep) = 0.d0
      zp_mean(istep) = 0.d0
      ze_mean(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)

         z1_mean(istep) = z1_mean(istep) + z1_curr
         z2_mean(istep) = z2_mean(istep) + z2_curr
         zp_mean(istep) = zp_mean(istep) + zp_curr
         ze_mean(istep) = ze_mean(istep) + ze_curr

      enddo

      z1_mean(istep) = z1_mean(istep)/number_of_traj
      z2_mean(istep) = z2_mean(istep)/number_of_traj
      zp_mean(istep) = zp_mean(istep)/number_of_traj
      ze_mean(istep) = ze_mean(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_mean.dat")
   do istep=1,number_of_timesteps
       write(2,'(5g20.10)') time(1,istep), z1_mean(istep), z2_mean(istep), zp_mean(istep), ze_mean(istep)
   enddo
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(3)-- Time-dependent second moments of solvent distribution
   !---------------------------------------------------------------

   write(*,'(/1x,"Building time-dependent variances... ",$)')
   time_start = secondi()

   allocate (z1_var(number_of_timesteps))
   allocate (z2_var(number_of_timesteps))
   allocate (z12_var(number_of_timesteps))
   allocate (zp_var(number_of_timesteps))
   allocate (ze_var(number_of_timesteps))
   allocate (zpe_var(number_of_timesteps))

   do istep=1,number_of_timesteps

      z1_var(istep) = 0.d0
      z2_var(istep) = 0.d0
      z12_var(istep) = 0.d0
      zp_var(istep) = 0.d0
      ze_var(istep) = 0.d0
      zpe_var(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_curr = z1(itraj,istep) - z1_mean(istep)
         z2_curr = z2(itraj,istep) - z2_mean(istep)
         zp_curr = zp(itraj,istep) - zp_mean(istep)
         ze_curr = ze(itraj,istep) - ze_mean(istep)

         z1_var(istep) = z1_var(istep) + z1_curr*z1_curr
         z2_var(istep) = z2_var(istep) + z2_curr*z2_curr
         z12_var(istep) = z12_var(istep) + z1_curr*z2_curr

         zp_var(istep) = zp_var(istep) + zp_curr*zp_curr
         ze_var(istep) = ze_var(istep) + ze_curr*ze_curr
         zpe_var(istep) = zpe_var(istep) + zp_curr*ze_curr

      enddo

      z1_var(istep) = z1_var(istep)/number_of_traj
      z2_var(istep) = z2_var(istep)/number_of_traj
      z12_var(istep) = z12_var(istep)/number_of_traj

      zp_var(istep) = zp_var(istep)/number_of_traj
      ze_var(istep) = ze_var(istep)/number_of_traj
      zpe_var(istep) = zpe_var(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_var.dat")
   do istep=1,number_of_timesteps
       write(2,'(7g20.10)') time(1,istep), z1_var(istep), z2_var(istep), z12_var(istep), zp_var(istep), ze_var(istep), zpe_var(istep)
   enddo
   close(2)

   deallocate (z1_mean, z2_mean, zp_mean, ze_mean, z1_var, z2_var, z12_var, zp_var, ze_var, zpe_var)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(4)-- Time-correlation functions of solvent coordinates
   !        (not clear what it means...)
   !---------------------------------------------------------------

   write(*,'(/1x,"Building time-correlation functions... ",$)')
   time_start = secondi()

   allocate (z11_tcf(number_of_timesteps))
   allocate (z22_tcf(number_of_timesteps))
   allocate (z12_tcf(number_of_timesteps))
   allocate (zpp_tcf(number_of_timesteps))
   allocate (zee_tcf(number_of_timesteps))
   allocate (zpe_tcf(number_of_timesteps))

   do istep=1,number_of_timesteps

      z11_tcf(istep) = 0.d0
      z22_tcf(istep) = 0.d0
      z12_tcf(istep) = 0.d0
      zpp_tcf(istep) = 0.d0
      zee_tcf(istep) = 0.d0
      zpe_tcf(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_0 = z1(itraj,1)
         z2_0 = z2(itraj,1)
         zp_0 = zp(itraj,1)
         ze_0 = ze(itraj,1)

         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)

         z11_tcf(istep) = z11_tcf(istep) + z1_0*z1_curr
         z22_tcf(istep) = z22_tcf(istep) + z2_0*z2_curr
         z12_tcf(istep) = z12_tcf(istep) + z1_0*z2_curr

         zpp_tcf(istep) = zpp_tcf(istep) + zp_0*zp_curr
         zee_tcf(istep) = zee_tcf(istep) + ze_0*ze_curr
         zpe_tcf(istep) = zpe_tcf(istep) + zp_0*ze_curr

      enddo

      z11_tcf(istep) = z11_tcf(istep)/number_of_traj
      z22_tcf(istep) = z22_tcf(istep)/number_of_traj
      z12_tcf(istep) = z12_tcf(istep)/number_of_traj

      zpp_tcf(istep) = zpp_tcf(istep)/number_of_traj
      zee_tcf(istep) = zee_tcf(istep)/number_of_traj
      zpe_tcf(istep) = zpe_tcf(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_tcf.dat")
   do istep=1,number_of_timesteps
       write(2,'(7g20.10)') time(1,istep), z11_tcf(istep), z22_tcf(istep), z12_tcf(istep), zpp_tcf(istep), zee_tcf(istep), zpe_tcf(istep)
   enddo
   close(2)

   deallocate (z11_tcf)
   deallocate (z22_tcf)
   deallocate (z12_tcf)
   deallocate (zpp_tcf)
   deallocate (zee_tcf)
   deallocate (zpe_tcf)

   deallocate (z1, z2, zp, ze)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(5)-- Time-dependent average free energy and kinetic energy
   !---------------------------------------------------------------

   write(*,'(/1x,"Building average energies... ",$)')
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
   do istep=1,number_of_timesteps
      write(2,'(3g20.10)') time(1,istep), efe_mean(istep), ekin_mean(istep)
   enddo
   close(2)

   deallocate (efe_mean)
   deallocate (ekin_mean)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(6)-- Time-dependent EVB weights
   !---------------------------------------------------------------

   write(*,'(/1x,"Building average EVB weights... ",$)')
   time_start = secondi()

   allocate (w1a_mean(number_of_timesteps))
   allocate (w1b_mean(number_of_timesteps))
   allocate (w2a_mean(number_of_timesteps))
   allocate (w2b_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      w1a_mean(istep) = 0.d0
      w1b_mean(istep) = 0.d0
      w2a_mean(istep) = 0.d0
      w2b_mean(istep) = 0.d0

      do itraj=1,number_of_traj
         w1a_mean(istep) = w1a_mean(istep) + w1a(itraj,istep)
         w1b_mean(istep) = w1b_mean(istep) + w1b(itraj,istep)
         w2a_mean(istep) = w2a_mean(istep) + w2a(itraj,istep)
         w2b_mean(istep) = w2b_mean(istep) + w2b(itraj,istep)
      enddo

      w1a_mean(istep) = w1a_mean(istep)/number_of_traj
      w1b_mean(istep) = w1b_mean(istep)/number_of_traj
      w2a_mean(istep) = w2a_mean(istep)/number_of_traj
      w2b_mean(istep) = w2b_mean(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="weights_mean.dat")
   do istep=1,number_of_timesteps
      write(2,'(5g20.10)') time(1,istep), w1a_mean(istep), w1b_mean(istep), w2a_mean(istep), w2b_mean(istep)
   enddo
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   write(*,'(/1x,"Building time-dependent cross-correlation matrix of EVB weights... ",$)')
   time_start = secondi()

   allocate (cw_cross(6,number_of_timesteps))  !-- 1a-1b, 1a-2a, 1a-2b, 1b-2a, 1b-2b, 2a-2b

   do istep=1,number_of_timesteps

      cw_cross(:,istep) = 0.d0
      cd1a1a = 0.d0
      cd1b1b = 0.d0
      cd2a2a = 0.d0
      cd2b2b = 0.d0
      cc1a1b = 0.d0
      cc1a2a = 0.d0
      cc1a2b = 0.d0
      cc1b2a = 0.d0
      cc1b2b = 0.d0
      cc2a2b = 0.d0

      do itraj=1,number_of_traj
         cd1a1a = cd1a1a + w1a(itraj,istep)*w1a(itraj,istep)
         cd1b1b = cd1b1b + w1b(itraj,istep)*w1b(itraj,istep)
         cd2a2a = cd2a2a + w2a(itraj,istep)*w2a(itraj,istep)
         cd2b2b = cd2b2b + w2b(itraj,istep)*w2b(itraj,istep)
         cc1a1b = cc1a1b + w1a(itraj,istep)*w1b(itraj,istep)
         cc1a2a = cc1a2a + w1a(itraj,istep)*w2a(itraj,istep)
         cc1a2b = cc1a2b + w1a(itraj,istep)*w2b(itraj,istep)
         cc1b2a = cc1b2a + w1b(itraj,istep)*w2a(itraj,istep)
         cc1b2b = cc1b2b + w1b(itraj,istep)*w2b(itraj,istep)
         cc2a2b = cc2a2b + w2a(itraj,istep)*w2b(itraj,istep)
      enddo

      cd1a1a = cd1a1a/number_of_traj
      cd1b1b = cd1b1b/number_of_traj
      cd2a2a = cd2a2a/number_of_traj
      cd2b2b = cd2b2b/number_of_traj
      cc1a1b = cc1a1b/number_of_traj
      cc1a2a = cc1a2a/number_of_traj
      cc1a2b = cc1a2b/number_of_traj
      cc1b2a = cc1b2a/number_of_traj
      cc1b2b = cc1b2b/number_of_traj
      cc2a2b = cc2a2b/number_of_traj

      sq1a = sqrt(cd1a1a - w1a_mean(istep)*w1a_mean(istep))
      sq1b = sqrt(cd1b1b - w1b_mean(istep)*w1b_mean(istep))
      sq2a = sqrt(cd2a2a - w2a_mean(istep)*w2a_mean(istep))
      sq2b = sqrt(cd2b2b - w2b_mean(istep)*w2b_mean(istep))

      cw_cross(1,istep) = (cc1a1b - w1a_mean(istep)*w1b_mean(istep))/sq1a/sq1b
      cw_cross(2,istep) = (cc1a2a - w1a_mean(istep)*w2a_mean(istep))/sq1a/sq2a
      cw_cross(3,istep) = (cc1a2b - w1a_mean(istep)*w2b_mean(istep))/sq1a/sq2b
      cw_cross(4,istep) = (cc1b2a - w1b_mean(istep)*w2a_mean(istep))/sq1b/sq2a
      cw_cross(5,istep) = (cc1b2b - w1b_mean(istep)*w2b_mean(istep))/sq1b/sq2b
      cw_cross(6,istep) = (cc2a2b - w2a_mean(istep)*w2b_mean(istep))/sq2a/sq2b

   enddo

   !-- output to the external file for visualization

   open(2,file="weights_cross_corr.dat")
   do istep=1,number_of_timesteps
      write(2,'("#",t10,"time(ps)",t30,"1a-1b",t50,"1a-2a",t70,"1a-2b",t90,"1b-2a",t110,"1b-2b",t130,"2a-2b")')
      write(2,'(7g20.10)') time(1,istep), (cw_cross(k,istep),k=1,6)
   enddo
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   deallocate (w1a,w1a_mean)
   deallocate (w1b,w1b_mean)
   deallocate (w2a,w2a_mean)
   deallocate (w2b,w2b_mean)
   deallocate (cw_cross)

   !---------------------------------------------------------------
   !--(7)-- Time-dependent adiabatic populations
   !---------------------------------------------------------------

   write(*,'(/1x,"Building time-dependent adiabatic populations... ",$)')
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
   do istep=1,number_of_timesteps
      write(2,'(101f10.6)') time(1,istep), (pop_ad(k,istep),k=1,number_of_states)
   enddo
   close(2)

   deallocate (istate, pop_ad)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   total_time_end = secondi()
   write(*,'(/"===> All Done in ",f12.3," sec"/)') total_time_end-total_time_start

contains

   subroutine deallocate_all_arrays
      if (allocated(time)) deallocate(time)
      if (allocated(z1)) deallocate(z1)
      if (allocated(z2)) deallocate(z2)
      if (allocated(zp)) deallocate(zp)
      if (allocated(ze)) deallocate(ze)
      if (allocated(ekin)) deallocate(ekin)
      if (allocated(efe)) deallocate(efe)
      if (allocated(istate)) deallocate(istate)
      if (allocated(istate_occ)) deallocate(istate_occ)
      if (allocated(all_states)) deallocate(all_states)
      if (allocated(w1a)) deallocate(w1a)
      if (allocated(w1b)) deallocate(w1b)
      if (allocated(w2a)) deallocate(w2a)
      if (allocated(w2b)) deallocate(w2b)
      if (allocated(histogram_z12)) deallocate(histogram_z12)
      if (allocated(histogram_zpe)) deallocate(histogram_zpe)
      if (allocated(state_histogram_z12)) deallocate(state_histogram_z12)
      if (allocated(state_histogram_zpe)) deallocate(state_histogram_zpe)
      if (allocated(z1_mean)) deallocate(z1_mean)
      if (allocated(z2_mean)) deallocate(z2_mean)
      if (allocated(zp_mean)) deallocate(zp_mean)
      if (allocated(ze_mean)) deallocate(ze_mean)
      if (allocated(z1_var)) deallocate(z1_var)
      if (allocated(z2_var)) deallocate(z2_var)
      if (allocated(z12_var)) deallocate(z12_var)
      if (allocated(zp_var)) deallocate(zp_var)
      if (allocated(ze_var)) deallocate(ze_var)
      if (allocated(zpe_var)) deallocate(zpe_var)
      if (allocated(z11_tcf)) deallocate(z11_tcf)
      if (allocated(z22_tcf)) deallocate(z22_tcf)
      if (allocated(z12_tcf)) deallocate(z12_tcf)
      if (allocated(zpp_tcf)) deallocate(zpp_tcf)
      if (allocated(zee_tcf)) deallocate(zee_tcf)
      if (allocated(zpe_tcf)) deallocate(zpe_tcf)
      if (allocated(ekin_mean)) deallocate(ekin_mean)
      if (allocated(efe_mean)) deallocate(efe_mean)
      if (allocated(w1a_mean)) deallocate(w1a_mean)
      if (allocated(w1b_mean)) deallocate(w1b_mean)
      if (allocated(w2a_mean)) deallocate(w2a_mean)
      if (allocated(w2b_mean)) deallocate(w2b_mean)
      if (allocated(cw_cross)) deallocate(cw_cross)
      if (allocated(pop_ad)) deallocate(pop_ad)
   end subroutine deallocate_all_arrays

end program analyze_trajectories
