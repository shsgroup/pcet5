program reorg_frcm

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
!     The program calculates solvent reorganization energy matrices
!     for four fixed charge distributions represented as collections
!     of point charges. The solvent is described in dielectric
!     continuum approximation within the Frequency Resolved Cavity
!     Model (FRCM) or ellipsoidal cavity model (ELCM).
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
!======================================================================!

   use pardim
   use keys
   use strings
   use cst
   use timers
   use geometry

   implicit none

   character(len=  8) :: cdum
   character(len= 10) :: curdat, curtim
   character(len=160) :: input, line, lsline, timestamp

   logical :: dir_exists=.false.

   real(8) :: tstart, time_s, time0, time1, time2
   real(8) :: telapsed, tend, ttotal
   integer :: n1, lfile, ijob, ijobn, ispa, io_status

   !interface
   !   function second() result(tic)
   !   implicit none
   !   real(8) :: tic
   !   end function second
   !end interface


   !---(AVS)---> add interface statements for all called routines


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
   !-now called in setjob
   !call init_pardim

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
      write(*,'(/''*** Please specify the name of the main input'',&
                   &'' file as the command line argument ***''/)')
      stop
   endif

   !---------------------------------------------
   ! Open the file and parse keywords
   !---------------------------------------------
   open(5,file=input(1:lfile),status='old')

   read(5,'(a)') line

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
         ijob = 1
         write(job,'(i2.2)') ijob
         ljob = 2
      endif

      !-------------------------------------------------
      ! Create a subdirectory for the current job.
      ! The name is given in the string JOB.
      !
      ! ATTENTION!!! The following lines are machine
      !              dependent and are calls to a
      !              system routine. It might be
      !              different on different platforms.
      !---------------------------------------------

      !-- check if directory already exists

      call system('ls -p > ls.scr')
      open(1,file='ls.scr')

      io_status = 0
      dir_exists = .false.

      loop_over_ls: do while (io_status.eq.0)
         read(1,'(a)',iostat=io_status) lsline
         if (index(trim(lsline),job(1:ljob)//slash).ne.0) then
            dir_exists = .true.
            exit loop_over_ls
         endif
      enddo loop_over_ls

      close(1)
      call system('rm -f ls.scr')

      if (dir_exists) then

         write(*,'(/1x,"****************************** W A R N I N G !!! *******************************")')
         write(*,'( 1x,"The output directory for the current job <",a,"> already exists.")') job(1:ljob)
         write(*,'( 1x,"It might be too important to be overwritten and thus the directory will be")')
         write(*,'( 1x,"backed up (timestamp will be added to its original name).")')
         write(*,'( 1x,"********************************************************************************")')

         !-------------------------------------------------
         ! Create timestamp
         !---------------------------------------------
         call date_and_time(curdat,curtim)
         timestamp=curdat(5:6)//"-"//curdat(7:8)//"-"//curdat(1:4)//&
         &"-"//curtim(1:2)//":"//curtim(3:4)//":"//curtim(5:6)

         call system("mv "//job(1:ljob)//" "//job(1:ljob)//"."//trim(timestamp))

      endif

      !-- create the output directory and copy the input files

      call system('mkdir '//job(1:ljob))
      call system('cp '//input(1:lfile)//' '//job(1:ljob))

      !-------------------------------------------------
      ! Initiate timing and set date
      !---------------------------------------------
      time_s = second()
      time0  = time_s
      call date_and_time(curdat,curtim)
      strdat=curdat(5:6)//"/"//curdat(7:8)//"/"//curdat(1:4)//&
       " at "//curtim(1:2)//":"//curtim(3:4)//":"//curtim(5:6)

      !-------------------------------------------------
      ! Read titles and keywords for the current job
      !---------------------------------------------
      backspace 5
      call read0(5)

      !-------------------------------------------------
      ! Print a banner for the calculation and
      ! a list of specified keywords with definitions
      !---------------------------------------------
      call printb

      !-------------------------------------------------
      ! Set options and parameters for the current job
      !---------------------------------------------
      call setjob

      !-----------------------------------------------------
      ! EREORG - reorganization energy matrices and
      !          one-dimensional free energy profiles
      !          for a standard Marcus two-state ET model
      !----------------------------------------------------
      if (index(keywrd,' EREORG').ne.0) then
         call ereorg
      endif

      !-------------------------------------------------
      ! deallocate arrays and prepare for the next job
      !-------------------------------------------------
      call deinitmat

   close(5)

   tend = second()
   ttotal = tend - tstart
   write(*,'(/1x,72(''=''))')
   write(*,'( 1x,''job done ...'')')
   write(*,'( 1x,''total time elapsed: '',f15.3,'' sec.'')') ttotal
   write(*,'( 1x,72(''='')/)')

end program reorg_frcm
