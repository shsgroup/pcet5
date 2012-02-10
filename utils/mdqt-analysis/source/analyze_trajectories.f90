!======================================================================================
!
!  Analysis of the dynamical trajectories
!  Alexander V. Soudackov, Penn State University
!  March 14, 2011
!
!  $Author: souda $
!  $Date: 2012-02-10 22:03:44 $
!  $Revision: 1.15 $
!  $Log: not supported by cvs2svn $
!  Revision 1.14  2011/06/20 22:04:39  souda
!  added mean components of kinetic energy (global)
!
!  Revision 1.13  2011/05/07 04:38:03  souda
!  Several silly bugs fixed (sorry, analysis to be rerun...)
!
!  Revision 1.12  2011/05/07 03:06:28  souda
!  Added combined mean weights 1a/2a and 1b/2b (Sharons suggestion)
!
!  Revision 1.11  2011/05/07 00:15:44  souda
!  Various fixes; EVB weights are probab. units; added cross-corr. for time averaged weights.
!
!  Revision 1.10  2011/05/06 21:45:54  souda
!  added marginal distributions for ground state distributions;
!  added correlation coefficients for running time averages of solvent coordinates;
!  added mean EVB weights calculated by assignment of state with largest contribution;
!  (correlation coefficients for mean EVB weights are to be added in the next revision).
!
!  Revision 1.9  2011/05/05 23:51:57  souda
!  now outputs correlation coefficients for solvent coordinates instead of covariances
!
!  Revision 1.8  2011/05/05 03:56:48  souda
!  Added marginal distributions for excited states distributions;
!  fixed the bug causing a memory leak;
!  fixed the bug in calculation of state-resolved distributions.
!
!======================================================================================

program analyze_trajectories

   use string_utilities
   use sorting
   use timers
   implicit none

   integer, parameter :: number_of_bins_z12=50
   integer, parameter :: number_of_bins_zpe=50
   integer, parameter :: ndt=27
   integer, parameter :: ndt2 = (ndt-1)/2

   character(len=60)  :: filename
   character(len=200) :: record
   character(len=3)   :: state_suffix

   logical :: found

   integer :: number_of_traj, number_of_timesteps, number_of_states
   integer :: number_of_occ_states
   integer :: itraj, istep, eof=0
   integer :: i, ilargest, ii, iocc, k, i1, i2, iargc, nsteps
   integer :: ibin_1, ibin_2, ibin_p, ibin_e
   integer :: ichannel_z12, ichannel_zpe

   real(kind=8) :: z1_min, z1_max, z2_min, z2_max
   real(kind=8) :: zp_min, zp_max, ze_min, ze_max
   real(kind=8) :: z1_curr, z2_curr, zp_curr, ze_curr
   real(kind=8) :: z1_0, z2_0, zp_0, ze_0, efe_curr, ekin_curr, ekin_z1_curr, ekin_z2_curr
   real(kind=8) :: time_start, time_end, total_time_start, total_time_end

   real(kind=8) :: bin_width_z1, bin_width_z2, bin_width_zp, bin_width_ze

   real(kind=8) :: cd1a1a, cd1b1b, cd2a2a, cd2b2b
   real(kind=8) :: cc1a1b, cc1a2a, cc1a2b, cc1b2a, cc1b2b, cc2a2b
   real(kind=8) :: sq1a, sq1b, sq2a, sq2b

   real(kind=8) :: z1tav, z2tav, z11tav, z22tav, z12tav, c12tav
   real(kind=8) :: zptav, zetav, zpptav, zeetav, zpetav, cpetav

   real(kind=8) :: w1atav, w1btav, w2atav, w2btav
   real(kind=8) :: w1a1atav, w1b1btav, w2a2atav, w2b2btav
   real(kind=8) :: w1a1btav, w1a2atav, w1a2btav, w1b2atav, w1b2btav, w2a2btav
   real(kind=8) :: c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   real(kind=8), dimension(number_of_bins_z12) :: bin_center_z1, bin_center_z2
   real(kind=8), dimension(number_of_bins_zpe) :: bin_center_zp, bin_center_ze

   real(kind=8), dimension(:,:),   allocatable :: time
   real(kind=8), dimension(:,:),   allocatable :: z1, z2
   real(kind=8), dimension(:,:),   allocatable :: zp, ze
   real(kind=8), dimension(:,:),   allocatable :: ekin, ekin_z1, ekin_z2, efe
   real(kind=8), dimension(:,:),   allocatable :: w1a, w1b, w2a, w2b
   real(kind=8), dimension(:,:),   allocatable :: histogram_z12
   real(kind=8), dimension(:,:),   allocatable :: histogram_zpe
   real(kind=8), dimension(:),     allocatable :: histogram_marg_z1, histogram_marg_z2
   real(kind=8), dimension(:),     allocatable :: histogram_marg_zp, histogram_marg_ze
   real(kind=8), dimension(:,:,:), allocatable :: state_histogram_z12
   real(kind=8), dimension(:,:,:), allocatable :: state_histogram_zpe
   real(kind=8), dimension(:,:),   allocatable :: excited_states_histogram_z12
   real(kind=8), dimension(:,:),   allocatable :: excited_states_histogram_zpe
   real(kind=8), dimension(:),     allocatable :: z1_mean, z2_mean, zp_mean, ze_mean
   real(kind=8), dimension(:),     allocatable :: z1_var, z2_var, z12_var, zp_var, ze_var, zpe_var
   real(kind=8), dimension(:),     allocatable :: z11_tcf, z22_tcf, z12_tcf, z21_tcf
   real(kind=8), dimension(:),     allocatable :: zpp_tcf, zee_tcf, zpe_tcf, zep_tcf
   real(kind=8), dimension(:),     allocatable :: efe_mean, ekin_mean, ekin_z1_mean, ekin_z2_mean
   real(kind=8), dimension(:),     allocatable :: w1a_mean, w1b_mean, w2a_mean, w2b_mean
   real(kind=8), dimension(:),     allocatable :: wh1a_mean, wh1b_mean, wh2a_mean, wh2b_mean
   real(kind=8), dimension(:),     allocatable :: w1ab_mean, w2ab_mean, w12a_mean, w12b_mean
   real(kind=8), dimension(:,:),   allocatable :: cw_cross
   real(kind=8), dimension(:,:),   allocatable :: pop_ad

   integer, dimension(:,:), allocatable :: istate
   integer, dimension(:), allocatable :: all_states
   integer, dimension(:), allocatable :: istate_occ

   character(len=20), dimension(100) :: carr
   integer,           dimension(100) :: iarr
   real(kind=8),      dimension(100) :: rarr
   real(kind=8),      dimension(4)   :: tmparray

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
   allocate(ekin_z1(number_of_traj,number_of_timesteps))
   allocate(ekin_z2(number_of_traj,number_of_timesteps))
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

         call readrec(record,"rrrrrrrrrrrirrrrrr",carr,iarr,rarr)

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
 
         w1a(itraj,nsteps) = rarr(12)/100.d0
         w1b(itraj,nsteps) = rarr(13)/100.d0
         w2a(itraj,nsteps) = rarr(14)/100.d0
         w2b(itraj,nsteps) = rarr(15)/100.d0

         ekin_z1(itraj,nsteps) = rarr(16)
         ekin_z2(itraj,nsteps) = rarr(17)

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
   allocate(excited_states_histogram_z12(number_of_bins_z12,number_of_bins_z12))
   allocate(excited_states_histogram_zpe(number_of_bins_zpe,number_of_bins_zpe))

   !-- open a separate file for each occupied state

   do i=1,number_of_occ_states
      ichannel_z12 = 10 + 2*istate_occ(i) - 1
      ichannel_zpe = 10 + 2*istate_occ(i)
      write(state_suffix,'(i3.3)') istate_occ(i)
      open(ichannel_z12,file="z12_distribution_"//state_suffix//".dat",form="unformatted")
      open(ichannel_zpe,file="zpe_distribution_"//state_suffix//".dat",form="unformatted")
   enddo

   !-- open files for combined excited states

   open(571,file="z12_distribution_excited.dat",form="unformatted")
   open(572,file="zpe_distribution_excited.dat",form="unformatted")

   !-- open files for combined marginal distributions for excited states

   open(581,file="z1_marg_distribution_excited.dat",form="formatted")
   open(582,file="z2_marg_distribution_excited.dat",form="formatted")
   open(583,file="zp_marg_distribution_excited.dat",form="formatted")
   open(584,file="ze_marg_distribution_excited.dat",form="formatted")

   open(681,file="z1_marg_distribution_ground.dat",form="formatted")
   open(682,file="z2_marg_distribution_ground.dat",form="formatted")
   open(683,file="zp_marg_distribution_ground.dat",form="formatted")
   open(684,file="ze_marg_distribution_ground.dat",form="formatted")

   do istep=1,number_of_timesteps
   
      state_histogram_z12 = 0.d0
      state_histogram_zpe = 0.d0
      excited_states_histogram_z12 = 0.d0
      excited_states_histogram_zpe = 0.d0
   
      do itraj=1,number_of_traj
   
         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)
         
         iocc = istate(itraj,istep)
         
         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         ibin_2 = nint((z2_curr-z2_min)/bin_width_z2 + 0.5d0)

         state_histogram_z12(iocc,ibin_1,ibin_2) = state_histogram_z12(iocc,ibin_1,ibin_2) + 1.d0
    
         ibin_p = nint((zp_curr-zp_min)/bin_width_zp + 0.5d0)
         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)

         state_histogram_zpe(iocc,ibin_p,ibin_e) = state_histogram_zpe(iocc,ibin_p,ibin_e) + 1.d0
   
      enddo
   
      !-- normalize distributions
   
      do i1=1,number_of_bins_z12
         do i2=1,number_of_bins_z12
            do k=1,number_of_occ_states
               state_histogram_z12(istate_occ(k),i1,i2) = state_histogram_z12(istate_occ(k),i1,i2)/number_of_traj
            enddo
         enddo
      enddo
   
      do i1=1,number_of_bins_zpe
         do i2=1,number_of_bins_zpe
            do k=1,number_of_occ_states
               state_histogram_zpe(istate_occ(k),i1,i2) = state_histogram_zpe(istate_occ(k),i1,i2)/number_of_traj
            enddo
         enddo
      enddo

      !-- build combined histograms for all excited states

      do i1=1,number_of_bins_z12
         do i2=1,number_of_bins_z12
            do k=1,number_of_occ_states
               if (istate_occ(k).gt.1) &
               & excited_states_histogram_z12(i1,i2) = excited_states_histogram_z12(i1,i2) &
               & + state_histogram_z12(istate_occ(k),i1,i2)
            enddo
         enddo
      enddo

      do i1=1,number_of_bins_zpe
         do i2=1,number_of_bins_zpe
            do k=1,number_of_occ_states
               if (istate_occ(k).gt.1) &
               & excited_states_histogram_zpe(i1,i2) = excited_states_histogram_zpe(i1,i2) &
               & + state_histogram_zpe(istate_occ(k),i1,i2)
            enddo
         enddo
      enddo

      !-- output to the external files for visualization

      !-- state-resolved distributions

      do i=1,number_of_occ_states

         ichannel_z12 = 10 + 2*istate_occ(i) - 1
         ichannel_zpe = 10 + 2*istate_occ(i)

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

      !-- combined excited states distributions (including marginal distributions)

      do i1=1,number_of_bins_z12
         write(581,'(3g15.6)') time(1,istep), bin_center_z1(i1), sum(excited_states_histogram_z12(i1,:))
         write(582,'(3g15.6)') time(1,istep), bin_center_z2(i1), sum(excited_states_histogram_z12(:,i1))
         do i2=1,number_of_bins_z12
            !write(571,'(f10.6,1x,f10.6,1x,g15.8)') bin_center_z1(i1), bin_center_z2(i2), excited_states_histogram_z12(i1,i2)
            write(571) bin_center_z1(i1), bin_center_z2(i2), excited_states_histogram_z12(i1,i2)
         enddo
      enddo
      write(581,*)
      write(582,*)

      do i1=1,number_of_bins_zpe
         write(583,'(3g15.6)') time(1,istep), bin_center_zp(i1), sum(excited_states_histogram_zpe(i1,:))
         write(584,'(3g15.6)') time(1,istep), bin_center_ze(i1), sum(excited_states_histogram_zpe(:,i1))
         do i2=1,number_of_bins_zpe
            !write(572,'(f10.6,1x,f10.6,1x,g15.8)') bin_center_zp(i1), bin_center_ze(i2), excited_states_histogram_zpe(i1,i2)
            write(572) bin_center_zp(i1), bin_center_ze(i2), excited_states_histogram_zpe(i1,i2)
         enddo
      enddo
      write(583,*)
      write(584,*)

      !-- ground state marginal distributions

      do i=1,number_of_bins_z12
         write(681,'(3g15.6)') time(1,istep), bin_center_z1(i), sum(state_histogram_z12(1,i,:))
         write(682,'(3g15.6)') time(1,istep), bin_center_z2(i), sum(state_histogram_z12(1,:,i))
      enddo
      write(681,*)
      write(682,*)

      do i=1,number_of_bins_zpe
         write(683,'(3g15.6)') time(1,istep), bin_center_zp(i), sum(state_histogram_zpe(1,i,:))
         write(684,'(3g15.6)') time(1,istep), bin_center_ze(i), sum(state_histogram_zpe(1,:,i))
      enddo
      write(683,*)
      write(684,*)
  
   enddo
   
   do i=1,number_of_occ_states
      ichannel_z12 = 10 + 2*istate_occ(i) - 1
      ichannel_zpe = 10 + 2*istate_occ(i)
      close(ichannel_z12)
      close(ichannel_zpe)
   enddo
   close(571)
   close(572)
   close(581)
   close(582)
   close(583)
   close(584)
   close(681)
   close(682)
   close(683)
   close(684)
   
   deallocate (state_histogram_z12,state_histogram_zpe)
   deallocate (excited_states_histogram_z12,excited_states_histogram_zpe)
   deallocate (istate_occ)
   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !-- build global histograms
   !---------------------------------------------------------------

   write(*,*)
   write(*,'(1x,"Building global histograms for solvent coordinates... ",$)')
   time_start = secondi()
   
   !-- global 2D-distributions
   allocate(histogram_z12(number_of_bins_z12,number_of_bins_z12))
   allocate(histogram_zpe(number_of_bins_zpe,number_of_bins_zpe))

   open(21,file="z12_distribution_global.dat",form="unformatted")
   open(22,file="zpe_distribution_global.dat",form="unformatted")

   !-- marginal 1D-distributions
   allocate(histogram_marg_z1(number_of_bins_z12))
   allocate(histogram_marg_z2(number_of_bins_z12))
   allocate(histogram_marg_zp(number_of_bins_zpe))
   allocate(histogram_marg_ze(number_of_bins_zpe))

   open(31,file="z1_marg_distribution_global.dat",form="formatted")
   open(32,file="z2_marg_distribution_global.dat",form="formatted")
   open(33,file="zp_marg_distribution_global.dat",form="formatted")
   open(34,file="ze_marg_distribution_global.dat",form="formatted")
   
   do istep=1,number_of_timesteps
   
      histogram_z12 = 0.d0
      histogram_zpe = 0.d0
      histogram_marg_z1 = 0.d0
      histogram_marg_z2 = 0.d0
      histogram_marg_zp = 0.d0
      histogram_marg_ze = 0.d0
   
      do itraj=1,number_of_traj
   
         z1_curr = z1(itraj,istep)
         z2_curr = z2(itraj,istep)
         zp_curr = zp(itraj,istep)
         ze_curr = ze(itraj,istep)
         
         ibin_1 = nint((z1_curr-z1_min)/bin_width_z1 + 0.5d0)
         ibin_2 = nint((z2_curr-z2_min)/bin_width_z2 + 0.5d0)
         histogram_z12(ibin_1,ibin_2) = histogram_z12(ibin_1,ibin_2) + 1.d0
         histogram_marg_z1(ibin_1) = histogram_marg_z1(ibin_1) + 1.d0
         histogram_marg_z2(ibin_2) = histogram_marg_z2(ibin_2) + 1.d0
   
         ibin_p = nint((zp_curr-zp_min)/bin_width_zp + 0.5d0)
         ibin_e = nint((ze_curr-ze_min)/bin_width_ze + 0.5d0)
         histogram_zpe(ibin_p,ibin_e) = histogram_zpe(ibin_p,ibin_e) + 1.d0
         histogram_marg_zp(ibin_p) = histogram_marg_zp(ibin_p) + 1.d0
         histogram_marg_ze(ibin_e) = histogram_marg_ze(ibin_e) + 1.d0
   
      enddo
   
      !-- normalize distributions
   
      do i1=1,number_of_bins_z12
         histogram_marg_z1(i1) = histogram_marg_z1(i1)/number_of_traj
         histogram_marg_z2(i1) = histogram_marg_z2(i1)/number_of_traj
         do i2=1,number_of_bins_z12
            histogram_z12(i1,i2) = histogram_z12(i1,i2)/number_of_traj
         enddo
      enddo
   
      do i1=1,number_of_bins_zpe
         histogram_marg_zp(i1) = histogram_marg_zp(i1)/number_of_traj
         histogram_marg_ze(i1) = histogram_marg_ze(i1)/number_of_traj
         do i2=1,number_of_bins_zpe
            histogram_zpe(i1,i2) = histogram_zpe(i1,i2)/number_of_traj
         enddo
      enddo
   
      !-- output to the external files for visualization
   
      do i1=1,number_of_bins_z12
         write(31,'(3g15.6)') time(1,istep), bin_center_z1(i1), histogram_marg_z1(i1)
         write(32,'(3g15.6)') time(1,istep), bin_center_z2(i1), histogram_marg_z2(i1)
         do i2=1,number_of_bins_z12
            !write(21,'(3g15.6)') bin_center_z1(i1), bin_center_z2(i2), histogram_z12(i1,i2)
            write(21) bin_center_z1(i1), bin_center_z2(i2), histogram_z12(i1,i2)
         enddo
      enddo
      write(31,*)
      write(32,*)
   
      do i1=1,number_of_bins_zpe
         write(33,'(3g15.6)') time(1,istep), bin_center_zp(i1), histogram_marg_zp(i1)
         write(34,'(3g15.6)') time(1,istep), bin_center_ze(i1), histogram_marg_ze(i1)
         do i2=1,number_of_bins_zpe
            !write(22,'(3g15.6)') bin_center_z1(i1), bin_center_z2(i2), histogram_zpe(i1,i2)
            write(22) bin_center_zp(i1), bin_center_ze(i2), histogram_zpe(i1,i2)
         enddo
      enddo
      write(33,*)
      write(34,*)
   
   enddo
   
   close(21)
   close(22)
   close(31)
   close(32)
   close(33)
   close(34)
   
   deallocate (histogram_z12,histogram_zpe)
   deallocate (histogram_marg_z1,histogram_marg_z2)
   deallocate (histogram_marg_zp,histogram_marg_ze)
   
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
   write(2,'("#",t10,"time(ps)",t30,"<z1>",t50,"<z2>",t70,"<zp>",t90,"<ze>")')
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
      z12_var(istep) = z12_var(istep)/sqrt(z1_var(istep)*z2_var(istep))

      zp_var(istep) = zp_var(istep)/number_of_traj
      ze_var(istep) = ze_var(istep)/number_of_traj
      zpe_var(istep) = zpe_var(istep)/number_of_traj
      zpe_var(istep) = zpe_var(istep)/sqrt(zp_var(istep)*ze_var(istep))

   enddo

   !-- output to the external file for visualization

   open(2,file="z_var.dat")
   write(2,'("#",t10,"time(ps)",t30,"sigma(z1)",t50,"sigma(z2)",t68,"corr(z1,z2)",t90,"sigma(zp)",t110,"sigma(ze)",t128,"corr(zp,ze)")')
   do istep=1,number_of_timesteps
       write(2,'(7g20.10)') time(1,istep), &
       & sqrt(z1_var(istep)), sqrt(z2_var(istep)), &
       & z12_var(istep)/sqrt(z1_var(istep))/sqrt(z2_var(istep)), &
       & sqrt(zp_var(istep)), sqrt(ze_var(istep)), &
       & zpe_var(istep)/sqrt(zp_var(istep))/sqrt(ze_var(istep))
   enddo
   close(2)

   deallocate (z1_var, z2_var, z12_var, zp_var, ze_var, zpe_var)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !--------------------------------------------------------------------
   !--(4)-- Cross-correlations of time-averages of solvent coordinates
   !--------------------------------------------------------------------

   write(*,'(/1x,"Building cross-correlations of time-averages for solvent coordinates... ",$)')
   time_start = secondi()

   z1tav = z1_mean(1)
   z2tav = z2_mean(1)
   z11tav = z1_mean(1)*z1_mean(1)
   z22tav = z2_mean(1)*z2_mean(1)
   z12tav = z1_mean(1)*z2_mean(1)

   zptav = zp_mean(1)
   zetav = ze_mean(1)
   zpptav = zp_mean(1)*zp_mean(1)
   zeetav = ze_mean(1)*ze_mean(1)
   zpetav = zp_mean(1)*ze_mean(1)

   c12tav = 0.d0
   cpetav = 0.d0

   open(2,file="z_corr_tav.dat")
   write(2,'("#",t10,"time(ps)",t30,"<z1tav>",t50,"<z2tav>",t70,"<c12tav>",t90,"<zptav>",t110,"<zetav>",t130,"<cpetav>")')
   write(2,'(7g20.10)') time(1,1), z1tav, z2tav, c12tav, zptav, zetav, cpetav

   do istep=2,number_of_timesteps

      z1tav = z1tav + z1_mean(istep)
      z2tav = z2tav + z2_mean(istep)
      z11tav = z11tav + z1_mean(istep)*z1_mean(istep)
      z22tav = z22tav + z2_mean(istep)*z2_mean(istep)
      z12tav = z12tav + z1_mean(istep)*z2_mean(istep)

      zptav = zptav + zp_mean(istep)
      zetav = zetav + ze_mean(istep)
      zpptav = zpptav + zp_mean(istep)*zp_mean(istep)
      zeetav = zeetav + ze_mean(istep)*ze_mean(istep)
      zpetav = zpetav + zp_mean(istep)*ze_mean(istep)

      c12tav = (z12tav/istep - z1tav*z2tav/istep/istep) / &
      & sqrt((z11tav/istep - z1tav*z1tav/istep/istep)*(z22tav/istep - z2tav*z2tav/istep/istep))

      cpetav = (zpetav/istep - zptav*zetav/istep/istep) / &
      & sqrt((zpptav/istep - zptav*zptav/istep/istep)*(zeetav/istep - zetav*zetav/istep/istep))

      write(2,'(7g20.10)') time(1,istep), z1tav/istep, z2tav/istep, c12tav, zptav/istep, zetav/istep, cpetav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !--------------------------------------------------------------------
   !--(5)-- Cross-correlations of running (local) time-averages
   !        of solvent coordinates
   !--------------------------------------------------------------------

   write(*,'(/1x,"Building cross-correlations of running time-averages for solvent coordinates... ",$)')
   time_start = secondi()

   !-- time interval for averaging (number of timesteps),
   !   must be an odd number

   open(2,file="z_corr_rtav.dat")
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time(1,2)-time(1,1))
   write(2,'("#",t10,"time(ps)",t30,"<z1tav>",t50,"<z2tav>",t70,"<c12tav>",t90,"<zptav>",t110,"<zetav>",t130,"<cpetav>")')

   do istep=1+ndt2,number_of_timesteps-ndt2

      z1tav  = 0.d0
      z2tav  = 0.d0
      z11tav = 0.d0
      z22tav = 0.d0
      z12tav = 0.d0

      zptav  = 0.d0
      zetav  = 0.d0
      zpptav = 0.d0
      zeetav = 0.d0
      zpetav = 0.d0

      do i=istep-ndt2,istep+ndt2

         z1tav = z1tav + z1_mean(i)
         z2tav = z2tav + z2_mean(i)
         z11tav = z11tav + z1_mean(i)*z1_mean(i)
         z22tav = z22tav + z2_mean(i)*z2_mean(i)
         z12tav = z12tav + z1_mean(i)*z2_mean(i)

         zptav = zptav + zp_mean(i)
         zetav = zetav + ze_mean(i)
         zpptav = zpptav + zp_mean(i)*zp_mean(i)
         zeetav = zeetav + ze_mean(i)*ze_mean(i)
         zpetav = zpetav + zp_mean(i)*ze_mean(i)

      enddo

      z1tav  = z1tav/ndt
      z2tav  = z2tav/ndt
      z11tav = z11tav/ndt
      z22tav = z22tav/ndt
      z12tav = z12tav/ndt

      zptav  = zptav/ndt
      zetav  = zetav/ndt
      zpptav = zpptav/ndt
      zeetav = zeetav/ndt
      zpetav = zpetav/ndt

      c12tav = (z12tav - z1tav*z2tav) / &
      & sqrt((z11tav - z1tav*z1tav)*(z22tav - z2tav*z2tav))

      cpetav = (zpetav - zptav*zetav) / &
      & sqrt((zpptav - zptav*zptav)*(zeetav - zetav*zetav))

      write(2,'(7g20.10)') time(1,istep), z1tav, z2tav, c12tav, zptav, zetav, cpetav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(6)-- Non-equilibrium time-correlation functions
   !        of solvent coordinates
   !        (not clear what it means...)
   !---------------------------------------------------------------

   write(*,'(/1x,"Building time-correlation functions... ",$)')
   time_start = secondi()

   allocate (z11_tcf(number_of_timesteps))
   allocate (z22_tcf(number_of_timesteps))
   allocate (z12_tcf(number_of_timesteps))
   allocate (z21_tcf(number_of_timesteps))
   allocate (zpp_tcf(number_of_timesteps))
   allocate (zee_tcf(number_of_timesteps))
   allocate (zpe_tcf(number_of_timesteps))
   allocate (zep_tcf(number_of_timesteps))

   do istep=1,number_of_timesteps

      z11_tcf(istep) = 0.d0
      z22_tcf(istep) = 0.d0
      z12_tcf(istep) = 0.d0
      z21_tcf(istep) = 0.d0
      zpp_tcf(istep) = 0.d0
      zee_tcf(istep) = 0.d0
      zpe_tcf(istep) = 0.d0
      zep_tcf(istep) = 0.d0

      do itraj=1,number_of_traj

         z1_0 = z1(itraj,1) - z1_mean(1)
         z2_0 = z2(itraj,1) - z2_mean(1)
         zp_0 = zp(itraj,1) - zp_mean(1)
         ze_0 = ze(itraj,1) - ze_mean(1)

         z1_curr = z1(itraj,istep) - z1_mean(istep)
         z2_curr = z2(itraj,istep) - z2_mean(istep)
         zp_curr = zp(itraj,istep) - zp_mean(istep)
         ze_curr = ze(itraj,istep) - ze_mean(istep)

         z11_tcf(istep) = z11_tcf(istep) + z1_0*z1_curr
         z22_tcf(istep) = z22_tcf(istep) + z2_0*z2_curr
         z12_tcf(istep) = z12_tcf(istep) + z1_0*z2_curr
         z21_tcf(istep) = z21_tcf(istep) + z2_0*z1_curr

         zpp_tcf(istep) = zpp_tcf(istep) + zp_0*zp_curr
         zee_tcf(istep) = zee_tcf(istep) + ze_0*ze_curr
         zpe_tcf(istep) = zpe_tcf(istep) + zp_0*ze_curr
         zep_tcf(istep) = zep_tcf(istep) + ze_0*zp_curr

      enddo

      z11_tcf(istep) = z11_tcf(istep)/number_of_traj
      z22_tcf(istep) = z22_tcf(istep)/number_of_traj
      z12_tcf(istep) = z12_tcf(istep)/number_of_traj
      z21_tcf(istep) = z21_tcf(istep)/number_of_traj

      zpp_tcf(istep) = zpp_tcf(istep)/number_of_traj
      zee_tcf(istep) = zee_tcf(istep)/number_of_traj
      zpe_tcf(istep) = zpe_tcf(istep)/number_of_traj
      zep_tcf(istep) = zep_tcf(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="z_tcf.dat")
   write(2,'("#",t10,"time(ps)",t25,"<z1(0)z1(t)>",t45,"<z2(0)z2(t)>",t65,"<z1(0)z2(t)>",t85,"<z2(0)z1(t)>",t105,"<zp(0)zp(t)>",t125,"<ze(0)ze(t)>",t145,"<zp(0)ze(t)>",t165,"<ze(0)zp(t)>")')
   do istep=1,number_of_timesteps
       write(2,'(9g20.10)') time(1,istep), z11_tcf(istep), z22_tcf(istep), z12_tcf(istep), z21_tcf(istep), zpp_tcf(istep), zee_tcf(istep), zpe_tcf(istep), zep_tcf(istep)
   enddo
   close(2)

   deallocate (z11_tcf)
   deallocate (z22_tcf)
   deallocate (z12_tcf)
   deallocate (z21_tcf)
   deallocate (zpp_tcf)
   deallocate (zee_tcf)
   deallocate (zpe_tcf)
   deallocate (zep_tcf)

   deallocate (z1, z2, zp, ze)
   deallocate (z1_mean, z2_mean, zp_mean, ze_mean)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(7)-- Time-dependent average free energy and kinetic energy
   !---------------------------------------------------------------

   write(*,'(/1x,"Building average energies... ",$)')
   time_start = secondi()

   allocate (efe_mean(number_of_timesteps))
   allocate (ekin_mean(number_of_timesteps))

   allocate (ekin_z1_mean(number_of_timesteps))
   allocate (ekin_z2_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      efe_mean(istep) = 0.d0
      ekin_mean(istep) = 0.d0
      ekin_z1_mean(istep) = 0.d0
      ekin_z2_mean(istep) = 0.d0

      do itraj=1,number_of_traj
         efe_curr = efe(itraj,istep)
         ekin_curr = ekin(itraj,istep)
         ekin_z1_curr = ekin_z1(itraj,istep)
         ekin_z2_curr = ekin_z2(itraj,istep)
         efe_mean(istep) = efe_mean(istep) + efe_curr
         ekin_mean(istep) = ekin_mean(istep) + ekin_curr
         ekin_z1_mean(istep) = ekin_z1_mean(istep) + ekin_z1_curr
         ekin_z2_mean(istep) = ekin_z2_mean(istep) + ekin_z2_curr
      enddo

      efe_mean(istep) = efe_mean(istep)/number_of_traj
      ekin_mean(istep) = ekin_mean(istep)/number_of_traj
      ekin_z1_mean(istep) = ekin_z1_mean(istep)/number_of_traj
      ekin_z2_mean(istep) = ekin_z2_mean(istep)/number_of_traj

   enddo

   !-- output to the external file for visualization

   open(2,file="ene_mean.dat")
   write(2,'("#",79("-"))')
   write(2,'("#",t10,"time",t25,"Free energy",t45,"Kinetic energy",t65,"Kinetic (z1)",t85,"Kinetic (z2)")')
   write(2,'("#",t10,"(ps)",t25,"(kcal/mol) ",t45,"  (kcal/mol)  ",t65,"  (kcal/mol)  ",t85,"  (kcal/mol)  ")')
   write(2,'("#",79("-"))')
   do istep=1,number_of_timesteps
      write(2,'(5g20.10)') time(1,istep), efe_mean(istep), ekin_mean(istep), ekin_z1_mean(istep), ekin_z2_mean(istep)
   enddo
   close(2)

   deallocate (efe_mean)
   deallocate (ekin_mean, ekin_z1_mean, ekin_z2_mean)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(8)-- Time-dependent EVB weights
   !---------------------------------------------------------------

   write(*,'(/1x,"Building average EVB weights... ",$)')
   time_start = secondi()

   allocate (w1a_mean(number_of_timesteps))
   allocate (w1b_mean(number_of_timesteps))
   allocate (w2a_mean(number_of_timesteps))
   allocate (w2b_mean(number_of_timesteps))

   allocate (wh1a_mean(number_of_timesteps))
   allocate (wh1b_mean(number_of_timesteps))
   allocate (wh2a_mean(number_of_timesteps))
   allocate (wh2b_mean(number_of_timesteps))

   allocate (w1ab_mean(number_of_timesteps))
   allocate (w2ab_mean(number_of_timesteps))

   allocate (w12a_mean(number_of_timesteps))
   allocate (w12b_mean(number_of_timesteps))

   do istep=1,number_of_timesteps

      w1a_mean(istep) = 0.d0
      w1b_mean(istep) = 0.d0
      w2a_mean(istep) = 0.d0
      w2b_mean(istep) = 0.d0

      wh1a_mean(istep) = 0.d0
      wh1b_mean(istep) = 0.d0
      wh2a_mean(istep) = 0.d0
      wh2b_mean(istep) = 0.d0

      w1ab_mean(istep) = 0.d0
      w2ab_mean(istep) = 0.d0

      w12a_mean(istep) = 0.d0
      w12b_mean(istep) = 0.d0

      do itraj=1,number_of_traj

         w1a_mean(istep) = w1a_mean(istep) + w1a(itraj,istep)
         w1b_mean(istep) = w1b_mean(istep) + w1b(itraj,istep)
         w2a_mean(istep) = w2a_mean(istep) + w2a(itraj,istep)
         w2b_mean(istep) = w2b_mean(istep) + w2b(itraj,istep)

         !-- assign 1 for the largest weight (whXX arrays)
         tmparray(1) = w1a(itraj,istep)
         tmparray(2) = w1b(itraj,istep)
         tmparray(3) = w2a(itraj,istep)
         tmparray(4) = w2b(itraj,istep)

         ilargest = 1
         do i=2,4
            if (tmparray(i).gt.tmparray(ilargest)) ilargest = i
         enddo

         select case(ilargest)
            case(1)
               wh1a_mean(istep) = wh1a_mean(istep) + 1.d0
            case(2)
               wh1b_mean(istep) = wh1b_mean(istep) + 1.d0
            case(3)
               wh2a_mean(istep) = wh2a_mean(istep) + 1.d0
            case(4)
               wh2b_mean(istep) = wh2b_mean(istep) + 1.d0
         end select

         if (w1a(itraj,istep)+w1b(itraj,istep).gt.w2a(itraj,istep)+w2b(itraj,istep)) then
            w1ab_mean(istep) = w1ab_mean(istep) + 1.d0
         else
            w2ab_mean(istep) = w2ab_mean(istep) + 1.d0
         endif

         if (w1a(itraj,istep)+w2a(itraj,istep).gt.w1b(itraj,istep)+w2b(itraj,istep)) then
            w12a_mean(istep) = w12a_mean(istep) + 1.d0
         else
            w12b_mean(istep) = w12b_mean(istep) + 1.d0
         endif

      enddo

      w1a_mean(istep) = w1a_mean(istep)/number_of_traj
      w1b_mean(istep) = w1b_mean(istep)/number_of_traj
      w2a_mean(istep) = w2a_mean(istep)/number_of_traj
      w2b_mean(istep) = w2b_mean(istep)/number_of_traj

      wh1a_mean(istep) = wh1a_mean(istep)/number_of_traj
      wh1b_mean(istep) = wh1b_mean(istep)/number_of_traj
      wh2a_mean(istep) = wh2a_mean(istep)/number_of_traj
      wh2b_mean(istep) = wh2b_mean(istep)/number_of_traj

      w1ab_mean(istep) = w1ab_mean(istep)/number_of_traj
      w2ab_mean(istep) = w2ab_mean(istep)/number_of_traj

      w12a_mean(istep) = w12a_mean(istep)/number_of_traj
      w12b_mean(istep) = w12b_mean(istep)/number_of_traj

   enddo

   !-- output to the external files for visualization

   open(2,file="weights_mean.dat")
   write(2,'("#",179("-"))')
   write(2,'("#",t10,"time",t30,"<1a>",t50,"<1b>",t70,"<2a>",t90,"<2b>",t110,"<1a+1b>",t130,"<2a+2b>",t150,"<1a+2a>",t170,"<1b+2b>")')
   write(2,'("#",179("-"))')
   do istep=1,number_of_timesteps
      write(2,'(9g20.10)') time(1,istep), w1a_mean(istep), w1b_mean(istep),  w2a_mean(istep), w2b_mean(istep), &
                                        & w1a_mean(istep) + w1b_mean(istep), w2a_mean(istep) + w2b_mean(istep), &
                                        & w1a_mean(istep) + w2a_mean(istep), w1b_mean(istep) + w2b_mean(istep)
   enddo
   close(2)

   open(2,file="weights_assigned_mean.dat")
   write(2,'("#",179("-"))')
   write(2,'("#",t10,"time",t30,"<1a>",t50,"<1b>",t70,"<2a>",t90,"<2b>",t110,"<1a/1b>",t130,"<2a/2b>",t150,"<1a/2a>",t170,"<1b/2b>")')
   write(2,'("#",179("-"))')
   do istep=1,number_of_timesteps
      write(2,'(9g20.10)') time(1,istep), wh1a_mean(istep), wh1b_mean(istep), wh2a_mean(istep), wh2b_mean(istep), &
                                        & w1ab_mean(istep), w2ab_mean(istep), w12a_mean(istep), w12b_mean(istep)
   enddo
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !---------------------------------------------------------------
   !--(9)-- Time-dependent cross-correlation matrix of EVB weights
   !---------------------------------------------------------------

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
   write(2,'("#",139("-"))')
   write(2,'("#",t10,"time(ps)",t30,"1a-1b",t50,"1a-2a",t70,"1a-2b",t90,"1b-2a",t110,"1b-2b",t130,"2a-2b")')
   write(2,'("#",139("-"))')
   do istep=1,number_of_timesteps
      write(2,'(8g20.10)') time(1,istep), (cw_cross(k,istep),k=1,6)
   enddo
   write(2,'("#",139("-"))')
   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   deallocate (w1a)
   deallocate (w1b)
   deallocate (w2a)
   deallocate (w2b)

   !------------------------------------------------------------------------------
   !--(10)-- Time-dependent correlation coefficients of time-averaged EVB weights
   !------------------------------------------------------------------------------

   write(*,'(/1x,"Building time-dependent cross-correlation matrix of EVB weights... ",$)')
   time_start = secondi()

   w1atav = w1a_mean(1)
   w1btav = w1b_mean(1)
   w2atav = w2a_mean(1)
   w2btav = w2b_mean(1)

   w1a1atav = w1a_mean(1)*w1a_mean(1)
   w1b1btav = w1b_mean(1)*w1b_mean(1)
   w2a2atav = w2a_mean(1)*w2a_mean(1)
   w2b2btav = w2b_mean(1)*w2b_mean(1)

   w1a1btav = w1a_mean(1)*w1b_mean(1)
   w1a2atav = w1a_mean(1)*w2a_mean(1)
   w1a2btav = w1a_mean(1)*w2b_mean(1)
   w1b2atav = w1b_mean(1)*w2a_mean(1)
   w1b2btav = w1b_mean(1)*w2b_mean(1)
   w2a2btav = w2a_mean(1)*w2b_mean(1)

   c1a1btav = 0.d0
   c1a2atav = 0.d0
   c1a2btav = 0.d0
   c1b2atav = 0.d0
   c1b2btav = 0.d0
   c2a2btav = 0.d0

   open(2,file="weights_corr_tav.dat")
   write(2,'("#",t10,"time(ps)",t30,"<w1atav>",t50,"<w1btav>",t70,"<w2atav>",t90,"<w2btav>",t110,"<c1a1btav>",t130,"<c1a2atav>",t150,"<c1a2btav>",t170,"<c1b2atav>",t190,"<c1b2btav>",t210,"<c2a2btav>")')
   write(2,'(11g20.10)') time(1,1), w1atav, w1btav, w2atav, w2btav, c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   do istep=2,number_of_timesteps

      w1atav = w1atav + w1a_mean(istep)
      w1btav = w1btav + w1b_mean(istep)
      w2atav = w2atav + w2a_mean(istep)
      w2btav = w2btav + w2b_mean(istep)

      w1a1atav = w1a1atav + w1a_mean(istep)*w1a_mean(istep)
      w1b1btav = w1b1btav + w1b_mean(istep)*w1b_mean(istep)
      w2a2atav = w2a2atav + w2a_mean(istep)*w2a_mean(istep)
      w2b2btav = w2b2btav + w2b_mean(istep)*w2b_mean(istep)

      w1a1btav = w1a1btav + w1a_mean(istep)*w1b_mean(istep)
      w1a2atav = w1a2atav + w1a_mean(istep)*w2a_mean(istep)
      w1a2btav = w1a2btav + w1a_mean(istep)*w2b_mean(istep)
      w1b2atav = w1b2atav + w1b_mean(istep)*w2a_mean(istep)
      w1b2btav = w1b2btav + w1b_mean(istep)*w2b_mean(istep)
      w2a2btav = w2a2btav + w2a_mean(istep)*w2b_mean(istep)

      c1a1btav = (w1a1btav/istep - w1atav*w1btav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w1b1btav/istep - w1btav*w1btav/istep/istep))

      c1a2atav = (w1a2atav/istep - w1atav*w2atav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w2a2atav/istep - w2atav*w2atav/istep/istep))

      c1a2btav = (w1a2btav/istep - w1atav*w2btav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      c1b2atav = (w1b2atav/istep - w1btav*w2atav/istep/istep) / &
      & sqrt((w1b1btav/istep - w1btav*w1btav/istep/istep)* &
      &      (w2a2atav/istep - w2atav*w2atav/istep/istep))

      c1b2btav = (w1b2btav/istep - w1btav*w2btav/istep/istep) / &
      & sqrt((w1b1btav/istep - w1btav*w1btav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      c2a2btav = (w2a2btav/istep - w2atav*w2btav/istep/istep) / &
      & sqrt((w2a2atav/istep - w2atav*w2atav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      write(2,'(11g20.10)') time(1,istep), w1atav/istep, w1btav/istep, w2atav/istep, w2btav/istep, &
                                     & c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   enddo

   close(2)

   !-- for mean weights calculated by assignment

   w1atav = wh1a_mean(1)
   w1btav = wh1b_mean(1)
   w2atav = wh2a_mean(1)
   w2btav = wh2b_mean(1)

   w1a1atav = wh1a_mean(1)*wh1a_mean(1)
   w1b1btav = wh1b_mean(1)*wh1b_mean(1)
   w2a2atav = wh2a_mean(1)*wh2a_mean(1)
   w2b2btav = wh2b_mean(1)*wh2b_mean(1)

   w1a1btav = wh1a_mean(1)*wh1b_mean(1)
   w1a2atav = wh1a_mean(1)*wh2a_mean(1)
   w1a2btav = wh1a_mean(1)*wh2b_mean(1)
   w1b2atav = wh1b_mean(1)*wh2a_mean(1)
   w1b2btav = wh1b_mean(1)*wh2b_mean(1)
   w2a2btav = wh2a_mean(1)*wh2b_mean(1)

   c1a1btav = 0.d0
   c1a2atav = 0.d0
   c1a2btav = 0.d0
   c1b2atav = 0.d0
   c1b2btav = 0.d0
   c2a2btav = 0.d0

   open(2,file="weights_assigned_corr_tav.dat")
   write(2,'("#",t10,"time(ps)",t30,"<w1atav>",t50,"<w1btav>",t70,"<w2atav>",t90,"<w2btav>",t110,"<c1a1btav>",t130,"<c1a2atav>",t150,"<c1a2btav>",t170,"<c1b2atav>",t190,"<c1b2btav>",t210,"<c2a2btav>")')
   write(2,'(11g20.10)') time(1,1), w1atav, w1btav, w2atav, w2btav, c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   do istep=2,number_of_timesteps

      w1atav = w1atav + wh1a_mean(istep)
      w1btav = w1btav + wh1b_mean(istep)
      w2atav = w2atav + wh2a_mean(istep)
      w2btav = w2btav + wh2b_mean(istep)

      w1a1atav = w1a1atav + wh1a_mean(istep)*wh1a_mean(istep)
      w1b1btav = w1b1btav + wh1b_mean(istep)*wh1b_mean(istep)
      w2a2atav = w2a2atav + wh2a_mean(istep)*wh2a_mean(istep)
      w2b2btav = w2b2btav + wh2b_mean(istep)*wh2b_mean(istep)

      w1a1btav = w1a1btav + wh1a_mean(istep)*wh1b_mean(istep)
      w1a2atav = w1a2atav + wh1a_mean(istep)*wh2a_mean(istep)
      w1a2btav = w1a2btav + wh1a_mean(istep)*wh2b_mean(istep)
      w1b2atav = w1b2atav + wh1b_mean(istep)*wh2a_mean(istep)
      w1b2btav = w1b2btav + wh1b_mean(istep)*wh2b_mean(istep)
      w2a2btav = w2a2btav + wh2a_mean(istep)*wh2b_mean(istep)

      c1a1btav = (w1a1btav/istep - w1atav*w1btav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w1b1btav/istep - w1btav*w1btav/istep/istep))

      c1a2atav = (w1a2atav/istep - w1atav*w2atav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w2a2atav/istep - w2atav*w2atav/istep/istep))

      c1a2btav = (w1a2btav/istep - w1atav*w2btav/istep/istep) / &
      & sqrt((w1a1atav/istep - w1atav*w1atav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      c1b2atav = (w1b2atav/istep - w1btav*w2atav/istep/istep) / &
      & sqrt((w1b1btav/istep - w1btav*w1btav/istep/istep)* &
      &      (w2a2atav/istep - w2atav*w2atav/istep/istep))

      c1b2btav = (w1b2btav/istep - w1btav*w2btav/istep/istep) / &
      & sqrt((w1b1btav/istep - w1btav*w1btav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      c2a2btav = (w2a2btav/istep - w2atav*w2btav/istep/istep) / &
      & sqrt((w2a2atav/istep - w2atav*w2atav/istep/istep)* &
      &      (w2b2btav/istep - w2btav*w2btav/istep/istep))

      write(2,'(11g20.10)') time(1,istep), w1atav/istep, w1btav/istep, w2atav/istep, w2btav/istep, &
                                     & c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   !--------------------------------------------------------------------
   !--(11)-- Cross-correlations of running (local) time-averages
   !         of EVB weights
   !--------------------------------------------------------------------

   write(*,'(/1x,"Building cross-correlations of running time-averages for EVB weights... ",$)')
   time_start = secondi()

   !-- time interval for averaging (number of timesteps),
   !   must be an odd number

   open(2,file="weights_corr_rtav.dat")
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time(1,2)-time(1,1))
   write(2,'("#",t10,"time(ps)",t30,"<w1atav>",t50,"<w1btav>",t70,"<w2atav>",t90,"<w2btav>",t110,"<c1a1btav>",t130,"<c1a2atav>",t150,"<c1a2btav>",t170,"<c1b2atav>",t190,"<c1b2btav>",t210,"<c2a2btav>")')

   do istep=1+ndt2,number_of_timesteps-ndt2

      w1atav  = 0.d0
      w1btav  = 0.d0
      w2atav  = 0.d0
      w2btav  = 0.d0

      w1a1atav = 0.d0
      w1b1btav = 0.d0
      w2a2atav = 0.d0
      w2b2btav = 0.d0

      w1a1btav = 0.d0
      w1a2atav = 0.d0
      w1a2btav = 0.d0
      w1b2atav = 0.d0
      w1b2btav = 0.d0
      w2a2btav = 0.d0

      c1a1btav = 0.d0
      c1a2atav = 0.d0
      c1a2btav = 0.d0
      c1b2atav = 0.d0
      c1b2btav = 0.d0
      c2a2btav = 0.d0

      do i=istep-ndt2,istep+ndt2

         w1atav = w1atav + w1a_mean(i)
         w1btav = w1btav + w1b_mean(i)
         w2atav = w2atav + w2a_mean(i)
         w2btav = w2btav + w2b_mean(i)

         w1a1atav = w1a1atav + w1a_mean(i)*w1a_mean(i)
         w1b1btav = w1b1btav + w1b_mean(i)*w1b_mean(i)
         w2a2atav = w2a2atav + w2a_mean(i)*w2a_mean(i)
         w2b2btav = w2b2btav + w2b_mean(i)*w2b_mean(i)

         w1a1btav = w1a1btav + w1a_mean(i)*w1b_mean(i)
         w1a2atav = w1a2atav + w1a_mean(i)*w2a_mean(i)
         w1a2btav = w1a2btav + w1a_mean(i)*w2b_mean(i)
         w1b2atav = w1b2atav + w1b_mean(i)*w2a_mean(i)
         w1b2btav = w1b2btav + w1b_mean(i)*w2b_mean(i)
         w2a2btav = w2a2btav + w2a_mean(i)*w2b_mean(i)

      enddo

      w1atav = w1atav/ndt
      w1btav = w1btav/ndt
      w2atav = w2atav/ndt
      w2btav = w2btav/ndt

      w1a1atav = w1a1atav/ndt
      w1b1btav = w1b1btav/ndt
      w2a2atav = w2a2atav/ndt
      w2b2btav = w2b2btav/ndt

      w1a1btav = w1a1btav/ndt
      w1a2atav = w1a2atav/ndt
      w1a2btav = w1a2btav/ndt
      w1b2atav = w1b2atav/ndt
      w1b2btav = w1b2btav/ndt
      w2a2btav = w2a2btav/ndt

      c1a1btav = (w1a1btav - w1atav*w1btav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w1b1btav - w1btav*w1btav))

      c1a2atav = (w1a2atav - w1atav*w2atav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w2a2atav - w2atav*w2atav))

      c1a2btav = (w1a2btav - w1atav*w2btav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w2b2btav - w2btav*w2btav))

      c1b2atav = (w1b2atav - w1btav*w2atav) / &
      & sqrt((w1b1btav - w1btav*w1btav)*(w2a2atav - w2atav*w2atav))

      c1b2btav = (w1b2btav - w1btav*w2btav) / &
      & sqrt((w1b1btav - w1btav*w1btav)*(w2b2btav - w2btav*w2btav))

      c2a2btav = (w2a2btav - w2atav*w2btav) / &
      & sqrt((w2a2atav - w2atav*w2atav)*(w2b2btav - w2btav*w2btav))

      write(2,'(11g20.10)') time(1,istep), w1atav, w1btav, w2atav, w2btav, c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   enddo

   close(2)


   !-- for mean weights calculated by assignment

   open(2,file="weights_assigned_corr_rtav.dat")
   write(2,'("#--- Time interval for averaging (ps): ",f12.6)') (ndt-1)*(time(1,2)-time(1,1))
   write(2,'("#",t10,"time(ps)",t30,"<w1atav>",t50,"<w1btav>",t70,"<w2atav>",t90,"<w2btav>",t110,"<c1a1btav>",t130,"<c1a2atav>",t150,"<c1a2btav>",t170,"<c1b2atav>",t190,"<c1b2btav>",t210,"<c2a2btav>")')

   do istep=1+ndt2,number_of_timesteps-ndt2

      w1atav  = 0.d0
      w1btav  = 0.d0
      w2atav  = 0.d0
      w2btav  = 0.d0

      w1a1atav = 0.d0
      w1b1btav = 0.d0
      w2a2atav = 0.d0
      w2b2btav = 0.d0

      w1a1btav = 0.d0
      w1a2atav = 0.d0
      w1a2btav = 0.d0
      w1b2atav = 0.d0
      w1b2btav = 0.d0
      w2a2btav = 0.d0

      c1a1btav = 0.d0
      c1a2atav = 0.d0
      c1a2btav = 0.d0
      c1b2atav = 0.d0
      c1b2btav = 0.d0
      c2a2btav = 0.d0

      do i=istep-ndt2,istep+ndt2

         w1atav = w1atav + wh1a_mean(i)
         w1btav = w1btav + wh1b_mean(i)
         w2atav = w2atav + wh2a_mean(i)
         w2btav = w2btav + wh2b_mean(i)

         w1a1atav = w1a1atav + wh1a_mean(i)*wh1a_mean(i)
         w1b1btav = w1b1btav + wh1b_mean(i)*wh1b_mean(i)
         w2a2atav = w2a2atav + wh2a_mean(i)*wh2a_mean(i)
         w2b2btav = w2b2btav + wh2b_mean(i)*wh2b_mean(i)

         w1a1btav = w1a1btav + wh1a_mean(i)*wh1b_mean(i)
         w1a2atav = w1a2atav + wh1a_mean(i)*wh2a_mean(i)
         w1a2btav = w1a2btav + wh1a_mean(i)*wh2b_mean(i)
         w1b2atav = w1b2atav + wh1b_mean(i)*wh2a_mean(i)
         w1b2btav = w1b2btav + wh1b_mean(i)*wh2b_mean(i)
         w2a2btav = w2a2btav + wh2a_mean(i)*wh2b_mean(i)

      enddo

      w1atav = w1atav/ndt
      w1btav = w1btav/ndt
      w2atav = w2atav/ndt
      w2btav = w2btav/ndt

      w1a1atav = w1a1atav/ndt
      w1b1btav = w1b1btav/ndt
      w2a2atav = w2a2atav/ndt
      w2b2btav = w2b2btav/ndt

      w1a1btav = w1a1btav/ndt
      w1a2atav = w1a2atav/ndt
      w1a2btav = w1a2btav/ndt
      w1b2atav = w1b2atav/ndt
      w1b2btav = w1b2btav/ndt
      w2a2btav = w2a2btav/ndt

      c1a1btav = (w1a1btav - w1atav*w1btav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w1b1btav - w1btav*w1btav))

      c1a2atav = (w1a2atav - w1atav*w2atav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w2a2atav - w2atav*w2atav))

      c1a2btav = (w1a2btav - w1atav*w2btav) / &
      & sqrt((w1a1atav - w1atav*w1atav)*(w2b2btav - w2btav*w2btav))

      c1b2atav = (w1b2atav - w1btav*w2atav) / &
      & sqrt((w1b1btav - w1btav*w1btav)*(w2a2atav - w2atav*w2atav))

      c1b2btav = (w1b2btav - w1btav*w2btav) / &
      & sqrt((w1b1btav - w1btav*w1btav)*(w2b2btav - w2btav*w2btav))

      c2a2btav = (w2a2btav - w2atav*w2btav) / &
      & sqrt((w2a2atav - w2atav*w2atav)*(w2b2btav - w2btav*w2btav))

      write(2,'(11g20.10)') time(1,istep), w1atav, w1btav, w2atav, w2btav, c1a1btav, c1a2atav, c1a2btav, c1b2atav, c1b2btav, c2a2btav

   enddo

   close(2)

   time_end = secondi()
   write(*,'("Done in ",f12.3," sec"/)') time_end-time_start

   deallocate (w1a_mean,wh1a_mean)
   deallocate (w1b_mean,wh1b_mean)
   deallocate (w2a_mean,wh2a_mean)
   deallocate (w2b_mean,wh2b_mean)
   deallocate (w1ab_mean,w2ab_mean)
   deallocate (w12a_mean,w12b_mean)

   !---------------------------------------------------------------
   !--(10)-- Time-dependent adiabatic populations
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
      if (allocated(ekin_z1)) deallocate(ekin_z1)
      if (allocated(ekin_z2)) deallocate(ekin_z2)
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
      if (allocated(histogram_marg_z1)) deallocate(histogram_marg_z1)
      if (allocated(histogram_marg_z2)) deallocate(histogram_marg_z2)
      if (allocated(histogram_marg_zp)) deallocate(histogram_marg_zp)
      if (allocated(histogram_marg_ze)) deallocate(histogram_marg_ze)
      if (allocated(state_histogram_z12)) deallocate(state_histogram_z12)
      if (allocated(state_histogram_zpe)) deallocate(state_histogram_zpe)
      if (allocated(excited_states_histogram_z12)) deallocate(excited_states_histogram_z12)
      if (allocated(excited_states_histogram_zpe)) deallocate(excited_states_histogram_zpe)
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
      if (allocated(ekin_z1_mean)) deallocate(ekin_z1_mean)
      if (allocated(ekin_z2_mean)) deallocate(ekin_z2_mean)
      if (allocated(efe_mean)) deallocate(efe_mean)
      if (allocated(w1a_mean)) deallocate(w1a_mean)
      if (allocated(w1b_mean)) deallocate(w1b_mean)
      if (allocated(w2a_mean)) deallocate(w2a_mean)
      if (allocated(w2b_mean)) deallocate(w2b_mean)
      if (allocated(w1a_mean)) deallocate(wh1a_mean)
      if (allocated(w1b_mean)) deallocate(wh1b_mean)
      if (allocated(w2a_mean)) deallocate(wh2a_mean)
      if (allocated(w2b_mean)) deallocate(wh2b_mean)
      if (allocated(w1ab_mean)) deallocate(w1ab_mean)
      if (allocated(w2ab_mean)) deallocate(w2ab_mean)
      if (allocated(w12a_mean)) deallocate(w12a_mean)
      if (allocated(w12b_mean)) deallocate(w12b_mean)
      if (allocated(cw_cross)) deallocate(cw_cross)
      if (allocated(pop_ad)) deallocate(pop_ad)
   end subroutine deallocate_all_arrays

end program analyze_trajectories
