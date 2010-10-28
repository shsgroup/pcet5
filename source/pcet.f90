program pcet

!======================================================================!
!                                                                      !
!     " 2b, or not 2b: that is the question "                          !
!                                                                      !
!                 William Shakespear, Hamlet (Act 3, Scene 1)          !
!                                                                      !
!======================================================================!
!
!     I. Introduction
!     ---------------
!
!     The program calculates free energies and nonadiabatic
!     rates for general proton-coupled electron transfer (PCET)
!     reaction in polar solvent. The solvent is described
!     in dielectric continuum approximation, the solute is
!     described in terms of the four-state EVB model, the proton
!     and electrons are treated quantum mechanically on equal
!     footing.
!
!     II. References
!     --------------
!
!     1. A. V. Soudackov and S. Hammes-Schiffer.
!        Multistate Continuum Theory for Multiple Charge
!        Transfer Reactions in Solution.
!        J. Chem. Phys. 111 (1999) 4672-4687.
!
!     2. A. V. Soudackov and S. Hammes-Schiffer.
!        Theoretical Study of Photoinduced Proton-Coupled
!        Electron Transfer through Asymmetric Salt Bridges.
!        J. Am. Chem. Soc. 121 (1999) 10598-10607.
!
!
!     III. Versions history
!     --------------------
!
!     1998 - Program for calculation of free energy surfaces
!            as functions of two scalar solvent variables.
!            Solvent - electrostatic ellipsoidal model;
!            Solute  - 4-state EVB with parameters fit
!            to the quantum-chemical calculations of linear
!            D-NH...O-A complexes.
!            (by Alexander Soudackov, University of Notre Dame)
!
!     1999 - Program for calculation of nonadiabatic rates
!            of nonradiative transitions between ET
!            diabatic free energy surfaces for PCET.
!            Solvent - electrostatic ellipsoidal model;
!            Solute  - 4-state EVB with parameters fit
!            to the quantum-chemical calculations of linear
!            D-NH...O-A complexes.
!            (by Alexander Soudackov, University of Notre Dame)
!
!     1999 - Incorporation of the Voth's EVB model for water chain.
!            (by Helene Decornez, University of Notre Dame)
!
!     1999 - Incorporation of the FRCM model for calculation
!            of reorganization energy matrices.
!            (by Ivan Rostov, University of Notre Dame)
!
!     2000 - The common interface combining all previous
!     v2.0   programs for one quantum proton (deuterium)
!     v2.1   and four-state EVB description of the solute.
!            (by Alexander Soudackov, University of Notre Dame)
!
!     2003 - Quantization for gating coordinate (proton donor-acceptor
!     v3.0   distance) added along with the LEPS gas-phase potential.
!            (by Alexander Soudackov, Penn State)
!
!     2003 - Fortran-90 adaptation + modular structure
!     v4.0   (by Alexander Soudackov, Penn State)
!
!     2010 - added functionality: MDQT dynamics in the space
!     v5.0   of two solvent coordinates for four state model
!            (by Alexander Soudackov, Penn State)
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!======================================================================!

   use pardim
   use keys
   use strings
   use cst
   use geometry

   implicit none

   character( 8) :: cdum
   character(10) :: curdat, curtim
   character(160) :: input, line
   
   real*8  :: tstart, time_s, time0, time1, time2
   real*8  :: telapsed, tend, ttotal
   integer :: n1, lfile, ijob, ijobn, ispa

   interface
      function second() result(tic)
      implicit none
      real(8) :: tic
      end function second
   end interface

   !---------------------------------------------------
   ! Uncomment two lines below for debugging
   ! (otherwise for some debuggers you won't have time
   !  to set a breakpoint)
   !---------------------------------------------------
   !write(*,*) ' Debug dummy input...'
   !read(*,*) idum

   tstart = second()

   !---------------------------------------------
   ! Initialize maximum dimensions
   !---------------------------------------------
   call init_pardim

   !---------------------------------------------------
   ! Initialization of constants and conversion factors
   !---------------------------------------------------
   call init_cst

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read the main INPUT file
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !-------------------------------------------------------------------
   ! extract the name of the input file from the command line argument
   !-------------------------------------------------------------------
   n1 = 1
   call getarg(n1,input)
   lfile = nblen(input)

   if (lfile.eq.0) then
      write(*,'(/''*** Please specify the name of the input'',&
                   &'' file as the command line argument ***''/)')
      stop 'Sorry...'
   endif

   !---------------------------------------------
   ! Open the file and parse keywords
   !---------------------------------------------
   open(5,file=input(1:lfile),status='old')

   ijob = 0

   !-----------------------------------------------------
   ! MAIN LOOP OVER THE JOBS SPECIFIED IN THE INPUT FILE
   !-----------------------------------------------------
   do

      read(5,'(a)') line

      if (index(line,'END').ne.0.or.index(line,'end').ne.0) exit

      !-------------------------------------------------
      ! Extract the name of the job from the first line
      ! of the input. If a JOBNAME is not specified use
      ! the consequtive number of the job instead.
      ! JOB  - name of the job (no spaces or special symbols)
      ! LJOB - length of the JOB string
      !---------------------------------------------
      ijobn = index(line,'JOBNAME=')
      if (ijobn.ne.0) then
         ispa = index(line(ijobn+8:),' ')
         job = line(ijobn+8:ijobn+ispa+6)
         ljob = ispa - 1
      else
         ijob = ijob + 1
         write(job,'(i2.2)') ijob
         ljob = 2
      endif

      !-------------------------------------------------
      ! Create a subdirectory for the current job.
      ! The name is given in the string JOB.
      !
      ! ATTENTION!!! The following line is machine
      !              dependent and is a call to a
      !              system routine. It might be
      !              different on different platforms.
      !---------------------------------------------
      call system('mkdir '//job(1:ljob))
      call system('cp '//input(1:lfile)//' '//job(1:ljob))

      !-------------------------------------------------
      ! Read titles and keywords for the current job
      !---------------------------------------------
      backspace 5
      call read0(5)

      !-------------------------------------------------
      ! Initiate timing and set date
      !---------------------------------------------
      time_s = second()
      time0  = time_s
      call date_and_time(curdat,curtim)
      strdat=curdat(5:6)//"/"//curdat(7:8)//"/"//curdat(1:4)//&
       " at "//curtim(1:2)//":"//curtim(3:4)//":"//curtim(5:6)

      !-------------------------------------------------
      ! Print a banner for the calculation and
      ! a list of specified keywords with definitions
      !---------------------------------------------
      call printb

      !-------------------------------------------------
      ! Set options and parameters for the current job
      !---------------------------------------------
      call setjob

      !-------------------------------------------------
      !     ET2 - one-dimensional free energy profiles and
      !           non-adiabatic rates for single ET reaction
      !           (standard Marcus two-state model)
      !---------------------------------------------
      if (index(keywrd,' ET2').ne.0) then
         call et2
         if (index(keywrd,' QUANTUM').eq.0) cycle
      endif

      !-------------------------------------------------
      ! Initialize matrices on the grid
      !---------------------------------------------
      call initmat

      !-------------------------------------------------
      ! Calculate matrices on the grid
      ! along the coordinate of the quantum particle
      ! To be done: parallelization over grid points
      !---------------------------------------------
      time1 = second()
      call setmat
      time2 = second()

      telapsed = time2 - time1
      WRITE(6,'(/1X,70(''-''))')
      WRITE(6,'( 1X,''Time for calculation of matrices on the grid: '',F15.3,'' sec.'')') telapsed
      WRITE(6,'( 1X,70(''-'')/)')

      !-------------------------------------------------
      ! Print out gas-phase and electronically solvated
      ! energy profiles along the proton coordinate
      ! (if PTGAS keyword is specified)
      !-------------------------------------------------
      if (index(keywrd,' PTGAS').ne.0) call prngas

      !-------------------------------------------------
      ! Print out solvated free energy profiles along the
      ! proton coordinate (if PTSOL keyword is specified)
      !--------------------------------------------------
      if (index(keywrd,' PTSOL').ne.0) call prnsol

      !==================================================
      !  NOW START THE PRODUCTION:
      !==================================================

      !-------------------------------------------------
      ! PT2 - one-dimensional free energy profiles and
      !       non-adiabatic rates for single PT reaction
      !       (standard two-state model is utilized)
      !-------------------------------------------------
      if (index(keywrd,' PT2').ne.0) then
         write(*,*) 'PT2 is not implemented in this version...'
         !call pt2
      endif

      !-------------------------------------------------
      ! SURF3 - three-dimensional free energy surfaces
      !         If only one grid point along the gating
      !         coordinate is defined output two-dimensional
      !         surfaces at fixed gating distance
      !-----------------------------------------------------
      if (index(keywrd,' SURF3').ne.0) call surface3

      !-------------------------------------------------
      ! SURF2 - two-dimensional free energy surfaces
      !         Gating coordinate should be treated
      !         quantum-mechanically (GQUANT=.true.)
      !-------------------------------------------------
      if (index(keywrd,' SURF2').ne.0) call surface2

      !-------------------------------------------------
      ! PATH - free energy profiles along a given path
      !        (quantum gating)
      !-------------------------------------------------
      if (index(keywrd,' PATH2').ne.0) call path2

      !-------------------------------------------------
      ! PATH - free energy profiles along a given path
      !        (classical gating: 2D surfaces along the path)
      !------------------------------------------------------
      if (index(keywrd,' PATH3').ne.0) call path3

      !-------------------------------------------------
      ! MIN2  - find minimums on the 2D free energy surfaces
      !------------------------------------------------------
      if (index(keywrd,' MIN2').ne.0) call minima2

      !-------------------------------------------------
      ! MIN3  - find minimums on the 3D free energy surfaces
      !-----------------------------------------------------
      if (index(keywrd,' MIN3').ne.0) call minima3

      !-------------------------------------------------
      ! WAVEFUN3 - calculate and write out the proton
      !            and vibrational wavefunctions
      !            depending on the gating coordinate
      !-------------------------------------------------
      if (index(keywrd,' WAVEFUN3').ne.0) call wavef3

      !-------------------------------------------------
      ! WAVEFUN2 - calculate and write out the proton
      !            and gating vibrational wavefunctions
      !-------------------------------------------------
      if (index(keywrd,' WAVEFUN2').ne.0) call wavef2

      !-------------------------------------------------
      ! RATE3 - estimate non-adiabatic rates by averaging
      !         over the gating coordinate
      !-------------------------------------------------
      if (index(keywrd,' RATE3').ne.0) call rate3

      !-------------------------------------------------
      ! RATEB - estimate non-adiabatic rates using
      !         exact flux expressions
      !-------------------------------------------------
      if (index(keywrd,' RATEB').ne.0) call rateb

      !-------------------------------------------------
      ! RATE2 - estimate non-adiabatic rates with
      !         quantization over the gating coordinate
      !-------------------------------------------------
      if (index(keywrd,' RATE2').ne.0) call rate2

      !--------------------------------------------------------------
      !                  ######################
      !__________________# New in version 5.x #______________________
      !                  ######################
      ! DYNAMICS - solvent dynamics in the 2D space of solvent
      !            coordinates (optional Surface Hopping dynamics)
      !            and fixed gating distance (gating coordinate
      !            will be dynamica in the later version)
      !--------------------------------------------------------------
      if (index(keywrd,' DYNAMICS3').ne.0) call dynamics3

      !-------------------------------------------------
      ! deallocate arrays and prepare for the next job
      !-------------------------------------------------
      call deinitmat

   enddo
   !---------------------------------------------------------------
   ! END OF THE MAIN LOOP OVER THE JOBS SPECIFIED IN THE INPUT FILE
   !---------------------------------------------------------------

   close(5)

   tend = second()
   ttotal = tend - tstart
   write(*,'(/1x,72(''=''))')
   write(*,'( 1x,''all jobs done ...'')')
   write(*,'( 1x,''total time elapsed: '',f15.3,'' sec.'')') ttotal
   write(*,'( 1x,72(''='')/)')

end program pcet
