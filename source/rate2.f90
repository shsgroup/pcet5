subroutine rate2
!=======================================================================
!     Calculates the nonadiabatic rate of the PCET reaction in solution
!
!     Reference:
!     [*] Alexander Soudackov and Sharon Hammes-Schiffer,
!         Derivation of rate expressions for nonadiabatic
!         proton-coupled electron transfer reactions in solution,
!         J. Chem. Phys. 113, 2385 (2000).
!
!     OPTIONS:
!
!     T=<value>[/<value>/<nsteps>] - absolute temperature
!                                    [range (init/fin/npoints)]
!
!     T=<value1>-<value2>-...-<valueN> - absolute temperature list
!                                        (up to 20 values)
!
!     PREC=I - precursor state I is the only initial state.
!              No Boltzmann averaging over the initial states
!              is performed.
!
!     NPREC=N - N precursor states are taken into account.
!               The total rate is averaged over initial
!               precursor states with Boltzmann weights.
!
!     NSUCC=N - N successor states are taken into account
!
!     ZPREC=<ZP0>/<ZE0> - guess for the reactant minimum on the
!                         precursor free energy surface (kcal/mol)
!
!     ZSUCC=<ZP0>/<ZE0> - guess for the product minimum on the
!                         successor free energy surface (kcal/mol)
!
!     ERANALY - the reorganization free energy is evaluated analytically
!               using the expression (55) from Ref.[*].
!               Otherwise (if not specified) it is evaluated
!               numerically according to Marcus definition
!
!     EANUMER - the activation free energy is evaluated numerically
!               as the difference between the free energies at the
!               intersection point and at the reactant minimum.
!               Otherwise (if not specified) the activation free
!               energy is evaluated analytically using the Marcus
!               relation E_act=(DeltaG+Er)**2/(4*Er)
!
!     WEIGHTS - total EVB weights are printed out for all minima
!               and intersection points
!
!     LOG=<filename> - all RATE output is written to the external file
!                      <filename>. If the LOG option is not specified
!                      the output is directed to the standard output
!                      (unit 6).
!
!     ACC=<VALUE> - the desired accuracy of the minimization
!                   during location of reactants and products
!                   minima (kcal/mol)
!
!     SLIM=<VALUE> - the maximum length of the Newton-Raphson step
!                    in the minimization (kcal/mol)
!
!     XACC=<VALUE> - the desired accuracy in the location
!                    of the crossing point between one-dimensional
!                    slices of the reactant and product free
!                    energy surfaces (kcal/mol)
!
!     IPRINT=N - printing level in LBFGS minimization routine (default=-1)
!
!     MCORR=m  - number of corrections used in the limited memory matrix.
!             Values of m < 3  are not recommended, and large values
!             of m can result in excessive computing time. The range
!             3 <= m <= 20 is recommended (default=5)

!     FACTR=<value> - accuracy of LBFGS minimization (factor of machine precision)
!                  1.d+12 for low accuracy;
!                  1.d+7  for moderate accuracy; 
!                  1.d+1  for extremely high accuracy
!                  Default value is 1.d+6
!
!     PGTOL=<value> - projected gradient tolerance in LBFGS minimization (positive)
!                  The iteration will stop when max{|proj g_i|, i=1,...,n} <= pgtol
!                  where pg_i is the ith component of the projected gradient.
!
!     MAXIT=N  - maximum number of Newton-Raphson steps allowed
!
!     MAXITX=N  - maximum number of steps in locating a crossing point
!
!     ZPS=<SCALE> - scaling factor for ZP coordinate in the minimization
!
!     ZES=<SCALE> - scaling factor for ZE coordinate in the minimization
!
!     ALIN=<VALUE> - inner sphere reorganization energy in kcal/mol
!
!-----------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:37
!  4.1
!  Exp
!  rate2.f90,v 4.1 2010/06/25 20:02:37 souda Exp
!  rate2.f90,v
!  Revision 4.1  2010/06/25 20:02:37  souda
!  Release 4.1
!
!  Revision 1.11  2008/04/11 00:07:20  souda
!  length of string OPTIONS increased to 1024
!  to accomodate more options (not critical)
!
!  Revision 1.10  2007/03/12 23:08:04  souda
!  Modifications related to using LBFGS minimization method.
!
!  Revision 1.9  2004/06/14 15:58:24  souda
!  decrease amount of out put in verbose mode in rate*
!
!  Revision 1.8  2004/06/07 17:33:05  souda
!  changed output formats
!
!  Revision 1.7  2004/06/07 16:58:06  souda
!  corrected filenames
!
!  Revision 1.6  2004/02/10 21:55:59  souda
!  Added the temperature list specification
!
!  Revision 1.5  2004/02/04 16:09:59  souda
!  Extended output to Arrhenius.dat
!
!  Revision 1.4  2004/01/21 22:45:34  souda
!  innersphere contrib. added to output of reorg. energies
!
!  Revision 1.3  2004/01/16 17:27:07  souda
!  Apparent activation energy and prefactor added
!
!  Revision 1.2  2004/01/15 23:30:24  souda
!  minor bug fixes/formats
!
!  Revision 1.1.1.1  2004/01/13 19:09:42  souda
!  Initial PCET-4.0 Release
!
!
!=======================================================================
   use pardim
   use keys
   use strings
   use cst
   use control
   use quantum
   use gasmat
   use solmat
   use fesmin_2d
   use feszz_2d, only: feszz2

   implicit none
   character(1024) :: options
   character(  40) :: fname
   character(   5) :: mode

   logical :: weights, outlog,eanumer, eranaly, vrpout
   logical :: filemin, rdmin, trange, tlist

   integer :: ikey, ilog, lenf, ispa, log
   integer :: iiprec, inprec, iprec, nprec, insucc, nsucc
   integer :: irdm, itemp, itemp_start, itemp_end, islash, islash1, islash2
   integer :: idash, idashn, n_t
   integer :: izprec, izsucc, izprec1, izsucc1
   integer :: islim, iacc, ixacc, imaxit, maxitx, imaxitx
   integer :: izps, izes, ilin
   integer :: ielst, nz, nzdim, ndabf, nfedim, npsiga, i, j, np, ng
   integer :: nvibprec, nvibsucc, nstates
   integer :: it, kzf, ierr
   integer :: iprin, imcorr, ifactr, ipgtol

   ! hbarps is the Planck constant in (kcal/mole)*picosecond (defined in cst module)
   
   real(8) :: t_start, t_end, ti(20)
   real(8) :: zpprec0, zeprec0, zpsucc0, zesucc0
   real(8) :: zpprec, zeprec, zpsucc, zesucc, feprec, fesucc
   real(8) :: dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin
   real(8) :: xacc, zps, zes, alin
   real(8) :: denom, elxx, elyy, elxy
   real(8) :: zpij, zeij, feij, elpt, elet, erija, erijn
   real(8) :: eaijn, eaija, fcad
   real(8) :: vij
   real(8) :: range, dtemp, temp, prefac, zf, fezero
   real(8) :: wtotal, defe, roi, wi, eaexp
   real(8) :: beta, betafirst, ratefirst, betalast, ratelast
   real(8) :: eaeff, aeff, deltarate, deltabeta

   logical, allocatable, dimension(:) :: min_succ
   integer, allocatable, dimension(:) :: istel, istpr, istga
   integer, allocatable, dimension(:) :: istvb_px, istvb_sx
   integer, allocatable, dimension(:) :: istga_px, istga_sx

   real(8), allocatable, dimension(:)       :: fe, z_px, z_sx
   real(8), allocatable, dimension(:)       :: psiga_px, psiga_sx
   real(8), allocatable, dimension(:,:,:,:) :: psiel_px, psiel_sx
   real(8), allocatable, dimension(:,:,:,:) :: psipr_px, psipr_sx
   real(8), allocatable, dimension(:,:,:)   :: enel, envib, envibg
   real(8), allocatable, dimension(:)       :: fei,fej
   real(8), allocatable, dimension(:)       :: zpi, zei, zpj, zej
   real(8), allocatable, dimension(:,:)     :: fex, zpx, zex
   real(8), allocatable, dimension(:,:)     :: de, er, ea
   real(8), allocatable, dimension(:,:)     :: vc2, wij

   !real(8), allocatable, dimension(:,:)     :: vrpkk

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Check whether gating quantization is specified
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (.not.gquant) then
      write(*,'(/1x,"*** input error (in rate2): ",&
      &"you must specify quantum option for gating keyword ***"/)')
      stop 'try again (rate2)...'
   endif

   if (method.eq.2) then
      write(*,'(/1x,"*** method=2 is not implemented for ",&
      &"rate calculations (coupling is too complicated) ***"/)')
      stop 'try again (rate)...'
   endif

   vrpout=.true.

   ! default
   rdmin   = .false.
   filemin = .false.
   trange  = .false.
   tlist   = .false.

   !====================================================
   ! Extract options
   !====================================================
   ikey = index(keywrd,' RATE2(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** input error (in rate): ",&
      &"you must specify options for the rate keyword ***"/)')
      stop 'try again (rate)...'
   else
      call getopt(keywrd,ikey+7,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Log file
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   outlog = index(options,' LOG').ne.0

   if (outlog) then

      ilog = index(options,' LOG=')

      if (ilog.eq.0) then

         fname = job(1:ljob)//'/rates.log'
         lenf = ljob + 10

      else

         ispa = index(options(ilog+5:),' ')
         fname = options(ilog+5:ilog+ispa+3)
         lenf = ispa - 1
         call locase(fname,lenf)
         fname = job(1:ljob)//'/'//fname(1:lenf)
         lenf = lenf + ljob + 1

      endif

      write(6,'(/1x,"rates output is stored in the external file <",a,">")') fname(1:lenf)

      log = 1
      open(log,file=fname(1:lenf),status='new')

   else

      log = 6

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! BANNER
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   write(log,'(/1x,70("="))')
   write(log,'( 1x,"nonadiabatic rate calculation")')
   write(log,'( 1x,"with quantum gating mode")')
   write(log,'( 1x,70("=")/)')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Precursor states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   iiprec = index(options,' IPREC=')
   inprec = index(options,' NPREC=')

   if (iiprec.ne.0.and.inprec.ne.0) then
      write(*,'(/1x,"*** input error (in rate): ",/,&
      &"you must specify either iprec or nprec, not both!",/,&
      &"it is confusing, so read manuals and make up your mind"/)')
      stop 'try again (rate)...'
   endif

   if (iiprec.ne.0) then

      iprec = reada(options,iiprec+7)
      nprec = iprec
      write(log,'(/1x,"initial precursor state:",i3)') iprec
      write(log,'( 1x,"no Boltzmann averaging over initial states")')

   elseif (inprec.ne.0) then

      iprec = 1
      nprec = reada(options,inprec+7)
      write(log,'(/1x,"initial precursor state:",i3)') iprec
      write(log,'( 1x,"Boltzmann averaging over lowest ",i2," initial states")') nprec

   else

      iprec = 1
      nprec = 1
      write(log,'(/1x,"initial precursor state:",i3)') iprec
      write(log,'( 1x,"no Boltzmann averaging over initial states")')

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of successor states to include
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   insucc = index(options,' NSUCC=')
   if (insucc.ne.0) then
      nsucc = reada(options,insucc+7)
   else
      nsucc = 10
   endif
   write(log,'(/1x,"number of successor states:",i3)') nsucc

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Allocate arrays
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   allocate (fei(iprec:nprec))
   allocate (zpi(iprec:nprec))
   allocate (zei(iprec:nprec))
   allocate (fej(1:nsucc))
   allocate (zpj(1:nsucc))
   allocate (zej(1:nsucc))
   allocate (fex(iprec:nprec,1:nsucc))
   allocate (zpx(iprec:nprec,1:nsucc))
   allocate (zex(iprec:nprec,1:nsucc))
   allocate ( de(iprec:nprec,1:nsucc))
   allocate ( er(iprec:nprec,1:nsucc))
   allocate ( ea(iprec:nprec,1:nsucc))
   allocate (vc2(iprec:nprec,1:nsucc))
   allocate (wij(iprec:nprec,1:nsucc))

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! File containing minima for calculation of the rates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   filemin = index(options,' READMIN').ne.0
   irdm = index(options,' READMIN=')

   if (.not.filemin) then

      write(log,'(/1x,"minima will be determined by Newton-Raphson Minimization"/)')
      rdmin = .false.

   else

      rdmin=.true.
      ispa = index(options(irdm+9:),' ')
      fname = options(irdm+9:irdm+ispa+7)
      lenf = ispa -1
      call locase (fname,lenf)
      open(2,file=fname(1:lenf),status='old')
      call getmin(2,nprec-iprec+1,nsucc,zpi,zei,zpj,zej)
      close(2)
      write(log,'(/1x,"minima will be read in from  the external file <",a,">")') fname(1:lenf)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Temperature (or temperature range)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   itemp = index(options,' T=')

   if (itemp.ne.0) then

      itemp_start = itemp + 3
      itemp_end = itemp_start + index(options(itemp_start:),space)-1
      islash = index(options(itemp_start:itemp_end),'/')
      idash  = index(options(itemp_start:itemp_end),'-')

      if (islash.eq.0.and.idash.eq.0) then

         t_start = reada(options,itemp+3)
         t_end = t_start
         n_t = 1
         trange = .false.
	 tlist  = .false.

      elseif (islash.ne.0) then

         islash1 = itemp_start + islash - 1
         islash2 = islash1 + index(options(islash1+1:itemp_end),'/')

         if (islash2.eq.0) then
            write(*,*) " error in temperature range specification..."
            write(*,*) " (from rate2)"
            stop
         endif

         t_start = reada(options(itemp_start:islash1-1),1)
         t_end   = reada(options(islash1+1:islash2-1),1)
         n_t = reada(options,islash2+1)
         if (n_t.le.1) then
            write(*,*) " error in temperature range specification..."
            write(*,*) " you should specify more than one point"
            stop
         endif
         trange = .true.

      elseif (idash.ne.0) then

         idashn = idash
	 n_t = 0
         do while (idashn.ne.0)
	    n_t = n_t + 1
	    ti(n_t) = reada(options(itemp_start:itemp_start+idashn-1),1)
	    itemp_start = itemp_start + idashn
	    idashn = index(options(itemp_start:itemp_end),'-')
	 enddo
	 n_t = n_t + 1
         ti(n_t) = reada(options(itemp_start:itemp_end),1)
         tlist = .true.

      endif

   else

      t_start = 298.15d0
      t_end = t_start
      n_t = 1
      trange = .false.
      tlist  = .false.

   endif

   if (trange) then
      write(log,'(/1x,"Temperature range (K):",i3," steps from",f8.3," to",f8.3)')&
      n_t, t_start, t_end
   elseif (tlist) then
      write(log,'(/1x,"Temperature list (K):",/,(5f8.3))') (ti(i),i=1,n_t)
   else
      write(log,'(/1x,"Temperature (K):",f8.3)') t_start
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Initial guesses for reactant and product minima
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izprec = index(options,' ZPREC=')
   izsucc = index(options,' ZSUCC=')

   if (izprec.ne.0.and.izsucc.ne.0) then

      ! reactant minimum guess (zpprec0,zeprec0)

      izprec1 = izprec+7
      islash = index(options(izprec1:),'/')
      zpprec0 = reada(options(izprec1:izprec1+islash-2),1)
      zeprec0 = reada(options,izprec1+islash)

      ! product minimum (zpsucc0,zesucc0)

      izsucc1 = izsucc+7
      islash = index(options(izsucc1:),'/')
      zpsucc0 = reada(options(izsucc1:izsucc1+islash-2),1)
      zesucc0 = reada(options,izsucc1+islash)

      write(log,'(/1x,"initial guesses for minima (kcal/mol):")')
      write(log,'( 1x,"reactant=(",f8.3,",",f8.3,") and product=(",f8.3,",",f8.3,")")')&
      zpprec0, zeprec0, zpsucc0, zesucc0

   else

      write(*,'(/1x,"*** input error (in rate): ",/,&
      &"you must specify both zprec and zsucc options",&
      &" for rate keyword ***"/)')
      stop 'try again (rate2)...'

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for the activation energies (EANUMER)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eanumer = index(options,' EANUMER').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for the reorganization energies (ERANALY)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eranaly = index(options,' ERANALY').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Flag for output of the total EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   weights = index(options,' WEIGHTS').ne.0

   if (iminim.eq.1) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Maximum length of the Newton-Raphson step
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      islim = index(options,' SLIM=')
      if (islim.ne.0) then
         slim = reada(options,islim+6)
      else
         slim = 1.d-1
      endif
      write(log,'(/1X,"Maximum length of the Newton-Raphson step: ",g12.6)') slim

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Accuracy of the minimization
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      iacc = index(options,' ACC=')
      if (iacc.ne.0) then
         acc = reada(options,iacc+5)
      else
         acc = 1.d-6
      endif
      write(log,'(1x,"Accuracy of the minimization: ",g12.6)') acc

   elseif (iminim.eq.2) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Printing level in LBFGS routine
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      iprin = index(options,' IPRINT=')
      if (iprin.ne.0) then
         iprint = reada(options,iprin+8)
      else
         iprint = -1
      endif
      write(6,'(/1x,"Printing level in LBFGS minimization: ",i10)') iprint

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Printing level in LBFGS routine
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      imcorr = index(options,' MCORR=')
      if (imcorr.ne.0) then
         mcorr = reada(options,imcorr+7)
      else
         mcorr = 5
      endif
      write(6,'(/1x,"Number of corrections used in the limited memory matrix (LBFGS): ",i10)') mcorr

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! LBFGS tolerance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ifactr = index(options,' FACTR=')
      if (ifactr.ne.0) then
         factr = reada(options,ifactr+7)
      else
         factr = 1.d+6
      endif
      write(6,'(/1X,"LBFGS tolerance: ",G12.6)') factr

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! LBFGS projected gradient tolerance
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ipgtol = index(options,' PGTOL=')
      if (ipgtol.ne.0) then
         pgtol = reada(options,ipgtol+7)
      else
         pgtol = 1.d-6
      endif
      write(6,'(/1X,"LBFGS gradient tolerance: ",G12.6)') pgtol

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Accuracy of the crossing point search
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ixacc = index(options,' XACC=')
   if (ixacc.ne.0) then
      xacc = reada(options,ixacc+6)
   else
      xacc = 1.d-3
   endif
   write(log,'(1x,"Accuracy of the crossing point location: ",g12.6)') xacc

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations in searching minima
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxit = index(options,' MAXIT=')
   if (imaxit.ne.0) then
      maxit = reada(options,imaxit+7)
   else
      maxit = 100
   endif
   write(log,'(1x,"Maximum number of Newton-Raphson steps: ",i5)') maxit

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations in searching a crossing point
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxitx = index(options,' MAXITX=')
   if (imaxitx.ne.0) then
      maxitx = reada(options,imaxitx+8)
   else
      maxitx = 100
   endif
   write(log,'(1x,"maximum number of steps in locating a crossing point: ",i5)') maxitx

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Scaling factors for the solvent coordinates (ZPS,ZES)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   izps = index(options,' ZPS=')
   if (izps.ne.0.and.iminim.eq.1) then
      zps = reada(options,izps+5)
      write(log,'(1x,"Scaling factor for ZP coordinate: ",g12.6)') zps
   else
      zps = 1.d0
   endif

   izes = index(options,' ZES=')
   if (izes.ne.0.and.iminim.eq.1) then
      zes = reada(options,izes+5)
      write(log,'(1x,"Scaling factor for ZE coordinate: ",g12.6)') zes
   else
      zes = 1.d0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Innersphere reorganization energy
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ilin = index(options,' ALIN=')
   if (ilin.ne.0) then
      alin = reada(options,ilin+6)
      write(log,'(1x,"Inner sphere reorg E added (kcal/mol): ",G12.6)') alin
   else
      alin = 0.d0
   endif


   !====================================================
   ! START THE CALCULATION
   !====================================================

   ! The following is fixed in this version
   mode = 'DIAB2'
   ielst = 2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Some quantities needed for analytical
   ! evaluation of reorganization energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   denom = 2.d0*(erpt*eret-erx*erx)
   elxx  = eret/denom
   elyy  = erpt/denom
   elxy  = -erx/denom

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Allocate arrays for wavefunctions and energies;
   ! calculate indexing arrays for electronic and
   ! vibrational quantum numbers
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call alloc_psi(ielst,iprec,nprec,nsucc,nzdim,nz,nvibprec,nvibsucc,nfedim,npsiga,ierr)
   if (ierr.ne.0) then
      write(*,*) '*** allocation error in rate2:alloc_psi ***'
      stop
   endif

   allocate (min_succ(1:nsucc))

   ! call flush(log)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! In case of MGQUANT=3 reorder the free energy
   ! quantum states according to their energies;
   ! fill out the pointer arrays for precursor
   ! and successor states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !---------------- here ----------------------------
   ! not sure it is necessary...
   !---------------- here ----------------------------
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over the initial states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   min_succ = .false.

   do i=iprec,nprec

      ! Locate a minimum on the reactant surface I
      ! or calculate free energy if minima read in

      if (.not.rdmin) then

         write(log,'(/"  search of minumum for reactant state ",i3/)') i
         if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"h(1,1)",t91,"h(2,2)",t106,"h(1,2)")')
         call fes2min(mode,1,i,zpprec0,zps,zeprec0,zes,&
                      ierr,zpprec,zeprec,feprec,&
                      dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)

         if (ierr.eq.0) then
            write(6,'(1x,a4,2x,8f15.6)') "done",&
            &zpprec,zeprec,feprec,sqrt(dzpmin*dzpmin+dzemin*dzemin)   !,d2zpmin,d2zemin,d2zpzemin
         else
            write(6,'(1x,a4,2x,8f15.6)') "FAIL",&
            &zpprec,zeprec,feprec,sqrt(dzpmin*dzpmin+dzemin*dzemin)   !,d2zpmin,d2zemin,d2zpzemin
         endif

         if (ierr.eq.0) then
            zpprec0 = zpprec
            zeprec0 = zeprec
         endif
         zpi(i) = zpprec
         zei(i) = zeprec
         fei(i) = feprec

      else

         fe       = 0.d0
         z_px     = 0.d0
         enel     = 0.d0
         envib    = 0.d0
         envibg   = 0.d0
         psiel_px = 0.d0
         psipr_px = 0.d0
         psiga_px = 0.d0

         if (mgquant.eq.1.or.mgquant.eq.2) then
            nstates = i
         elseif (mgquant.eq.3) then
            nstates = nvibprec
         endif

         call feszz2(mode,1,zpi(i),zei(i),nstates,nfedim,fe,&
                     nz,z_px,ndabf,ielst,enel,envib,envibg,&
                     psiel_px,psipr_px,npsiga,psiga_px)

         zpprec = zpi(i)
         zeprec = zei(i)
         fei(i) = fe(i)

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate and print out the EVB weights
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (weights) then
         write(log,'(/1x,70("="))')
         write(log,'(1x,"EVB weights at the minimum of the ",i2,a3," reactant state")') i,th(i)
         write(log,'(1x,"zp(precursor)=",f9.3,5x,"ze(precursor)=",f9.3/)') zpprec,zeprec
         call weight2(log,zpprec,zeprec,mode,1,i)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Loop over the successor (product) states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do j=1,nsucc

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
         ! Locate a minimum on the product surface J  !
         ! or calculate free energy if minima read in !
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

         if (.not.min_succ(j)) then

            if (.not.rdmin) then

               write(log,'(/"  search of minumum for product state ",i3/)') j
               if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"h(1,1)",t91,"h(2,2)",t106,"h(1,2)")')

               call fes2min(mode,2,j,zpsucc0,zps,zesucc0,zes,&
                            ierr,zpsucc,zesucc,fesucc,&
                            dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)

               if (ierr.eq.0) then
                  write(6,'(1x,a4,2x,8f15.6)') "done",&
                  &zpsucc,zesucc,fesucc,sqrt(dzpmin*dzpmin+dzemin*dzemin)   !,d2zpmin,d2zemin,d2zpzemin
               else
                  write(6,'(1x,a4,2x,8f15.6)') "FAIL",&
                  &zpsucc,zesucc,fesucc,sqrt(dzpmin*dzpmin+dzemin*dzemin)   !,d2zpmin,d2zemin,d2zpzemin
               endif

               if (ierr.eq.0) then
                  zpsucc0 = zpsucc
                  zesucc0 = zesucc
               endif
               zpj(j) = zpsucc
               zej(j) = zesucc
               fej(j) = fesucc

            else

               fe       = 0.d0
               z_sx     = 0.d0
               enel     = 0.d0
               envib    = 0.d0
               envibg   = 0.d0
               psiel_px = 0.d0
               psipr_px = 0.d0
               psiga_sx = 0.d0

               if (mgquant.eq.1.or.mgquant.eq.2) then
                  nstates = j
               elseif (mgquant.eq.3) then
                  nstates = nvibsucc
               endif

               call feszz2(mode,2,zpj(j),zej(j),nstates,nfedim,fe,&
                           nz,z_sx,ndabf,ielst,enel,envib,envibg,&
                           psiel_px,psipr_px,npsiga,psiga_sx)

               zpsucc = zpj(j)
               zesucc = zej(j)
               fej(j) = fe(j)

            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! calculate and print out the evb weights
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (weights) then
               write(log,'(/1x,70("="))')
               write(log,'(1x,"EVB weights at the minimum of the ",i2,a3," product state")') j,th(j)
               write(log,'(1x,"zp(successor)=",f9.3,5x,"ze(successor)=",f9.3/)') zpsucc,zesucc
               call weight2(log,zpsucc,zesucc,mode,2,j)
            endif

            min_succ(j) = .true.

         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate reaction free energy for the pair (I,J)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         de(i,j) = fej(j) - fei(i)

         write(log,'(/1x,"the Reaction Free Energy",&
         &" for the pair of states ",i2,"-",i2,": ",f12.6/)') i,j,de(i,j)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate two-dimensional reorganization energy
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         zpij = zpj(j) - zpi(i)
         zeij = zej(j) - zei(i)
         elpt = elxx*zpij + elxy*zeij
         elet = elxy*zpij + elyy*zeij
         erija = elpt*elpt*erpt + elet*elet*eret + 2.d0*elpt*elet*erx + alin

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate the free energy of the reactant state I
         ! at the minimum of the product surface J
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         fe = 0.d0
         z_px = 0.d0
         enel = 0.d0
         envib = 0.d0
         envibg = 0.d0
         psiel_px = 0.d0
         psipr_px = 0.d0
         psiga_px = 0.d0

         if (mgquant.eq.1.or.mgquant.eq.2) then
            nstates = i
         elseif (mgquant.eq.3) then
            nstates = nvibprec
         endif

         call feszz2(mode,1,zpsucc,zesucc,nstates,nfedim,fe,&
                     nz,z_px,ndabf,ielst,enel,envib,envibg,&
                     psiel_px,psipr_px,npsiga,psiga_px)
         erijn = fe(i) - fei(i) + alin

         write(log,'(/1x,"the reorganization energy",&
         &" for the pair of states ",i2,"-",i2,":",/,&
         &" numerical  estimate (from the surfaces)  : ",f12.6,/,&
         &" analytical estimate (from our expression): ",f12.6/)')&
         i,j,erijn,erija

         if (eranaly) then
            er(i,j) = erija
         else
            er(i,j) = erijn
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Find the crossing point between electronically
         ! diabatic one-dimensional (slices along the straight
         ! line connecting two minima) free-energy curves
         ! for the states I and J
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call crossp2(i,j,zpprec,zeprec,zpsucc,zesucc,xacc,maxitx,zpij,zeij,feij)
         zpx(i,j) = zpij
         zex(i,j) = zeij
         fex(i,j) = feij

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !Calculate and print out the EVB weights
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (weights) then
            write(log,'(/1x,70("="))')
            write(log,'(1x,"EVB weights at the crossing point ",&
                         &"between ",i2,a3," reactant and ",&
                         &i2,a3," product states")') i,th(i),j,th(j)
            write(log,'(1x,"zp(xing point)=",f9.3,5x,"ze(xing point)=",f9.3/)') zpij,zeij
            call weight2(log,zpij,zeij,mode,1,i)
            call weight2(log,zpij,zeij,mode,2,j)
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate activation free energy for the pair (I,J)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         eaijn = fex(i,j) - fei(i)
         eaija = (er(i,j)+de(i,j))*(er(i,j)+de(i,j))/4.d0/er(i,j)

         WRITE(LOG,'(/1X,"The activation energy",&
      &" for the pair of states ",I2,"-",I2,":",/,&
      &" numerical estimate (from the crossing point): ",F12.6,/,&
      &" analytical estimate (from Marcus expression): ",F12.6/)') I,J,EAIJN,EAIJA

         if (eanumer) then
            ! numerical estimate
            ea(i,j) = eaijn
         else
            ! Marcus relation
            ea(i,j) = eaija
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate coupling matrix element at the crossing point
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !fe       = 0.d0
         !z_px     = 0.d0
         !enel     = 0.d0
         !envib    = 0.d0
         !envibg   = 0.d0
         !psiel_px = 0.d0
         !psipr_px = 0.d0
         !psiga_px = 0.d0

         if (mgquant.eq.1.or.mgquant.eq.2) then
            nstates = i
         elseif (mgquant.eq.3) then
            nstates = nvibprec
         endif

         call feszz2(mode,1,zpij,zeij,nstates,nfedim,fe,&
                     nz,z_px,ndabf,ielst,enel,envib,envibg,&
                     psiel_px,psipr_px,npsiga,psiga_px)

         !fe       = 0.d0
         !z_sx     = 0.d0
         !enel     = 0.d0
         !envib    = 0.d0
         !envibg   = 0.d0
         !psiel_sx = 0.d0
         !psipr_sx = 0.d0
         !psiga_sx = 0.d0

         if (mgquant.eq.1.or.mgquant.eq.2) then
            nstates = j
         elseif (mgquant.eq.3) then
            nstates = nvibsucc
         endif

         call feszz2(mode,2,zpij,zeij,nstates,nfedim,fe,&
                     nz,z_sx,ndabf,ielst,enel,envib,envibg,&
                     psiel_sx,psipr_sx,npsiga,psiga_sx)

         ! Actual coupling evalulation
         ! Note that METHOD=2 is not implemented for this (sorry).

         vij = coupling(i,j)
         vc2(i,j) = vij*vij

         !call flush(log)

      enddo         !loop over i - precursor states

   enddo         !loop over j - successor states

   !----------------------------------------------------------
   ! Deallocate arrays for wavefunctions
   !----------------------------------------------------------
   call dealloc_psi(ierr)
   if (ierr.ne.0) then
      write(*,*) '*** deallocation error in feszz2:dealloc_psi ***'
      stop
   endif

   !====================================================
   ! Initial states characteristics
   !====================================================
   write(log,'(/1x,"INITIAL STATES CHARACTERISTICS (kcal/mol)")')
   write(log,'(1x,70("="))')
   write(log,'(1x,t4,"state",t17,"zp0",t32,"ze0",t43,"fe(zp0,ze0)")')
   write(log,'(1x,70("-"))')
   do i=iprec,nprec
      write(log,'(1x,i5,5x,3(f10.3,5x))') i,zpi(i),zei(i),fei(i)
   enddo
   write(log,'(1x,70("-"))')

   !====================================================
   ! Crossing points
   !====================================================
   write(log,'(/1x,"CROSSING POINTS CHARACTERISTICS (kcal/mol)")')
   write(log,'(1x,70("="))')
   write(log,'(1x,t3,"pair",t14,"zpx",t29,"zex",t39,"fex(zp0,ze0)",t59,"v^2")')
   write(log,'(1x,70("-"))')
   do i=iprec,nprec
      do j=1,nsucc
         write(log,'(1x,i2,"-",i2,3(f10.3,5x),g12.6)')&
         i, j, zpx(i,j), zex(i,j), fex(i,j), vc2(i,j)
      enddo
   enddo
   write(log,'(1x,70("-"))')

   !====================================================
   ! Final states characteristics
   !====================================================
   write(log,'(/1x,"FINAL STATES CHARACTERISTICS (kcal/mol)")')
   write(log,'(1x,70("="))')
   write(log,'(1x,t4,"state",t17,"zp0",t32,"ze0",t43,"fe(zp0,ze0)")')
   write(log,'(1x,70("-"))')
   do j=1,nsucc
      write(log,'(1x,i5,5x,3(f10.3,5x))') j,zpj(j),zej(j),fej(j)
   enddo
   write(log,'(1x,70("-"))')

   !==================================================================
   ! Start the loop over temperatures
   !==================================================================

   if (trange.or.tlist) then
      open(99,file=job(1:ljob)//"/arrhenius_q.dat")
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("# Arrhenius plots: quantum rates (1/sec)")')
      write(99,'("#")')
      write(99,'("# 6 Columns:")')
      write(99,'("#")')
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("#    1           2         3          4                    5                   6      ")')
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("#   T(K)       100/T     1/kT        k_tot             log10(k_tot)         ln(k_tot) ")')
      write(99,'("#-------------------------------------------------------------------------------------")')
   endif

   if (trange) then
      range = t_end - t_start
      dtemp = range/(n_t - 1)
   else
      dtemp = 0.d0
   endif

   do it=1,n_t

      if (tlist) then
         temp = ti(it)
      else      
         temp = t_start + (it-1)*dtemp
      endif
      beta = 1.d0/(kb*temp)

      write(log,'(/1x,70("*"))')
      write(log,'( 1x,"Rate analysis at ",f8.3," K")') temp
      write(log,'(/1x,70("*"))')

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Preexponential factor (constant) in the rate expression
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      prefac = dsqrt(pi*beta)/hbarps

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! zf - partition function
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zf = 0.d0
      kzf = 0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate partial rates I-->J (1/picoseconds)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Loop over the initial states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=iprec,nprec

         kzf = kzf + 1

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Sum up partition function
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (kzf.eq.1) fezero = fei(i)
         zf = zf + dexp(-beta*(fei(i) - fezero))

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Loop over the successor (product) states
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do j=1,nsucc
            wij(i,j) = prefac*vc2(i,j)*dexp(-beta*ea(i,j))/dsqrt(er(i,j))
         enddo

      enddo

      !====================================================
      ! Total rate (WTOTAL) and final printing
      !====================================================

      wtotal = 0.d0
      write(log,'(/1x,"RATE CONTRIBUTIONS FOR T=",f8.3,"K")') temp
      write(log,'(1x,79("="))')
      write(log,'(1x,"column headers:",/,&
      &1x,"J   - successor state",/,&
      &1x,"DeG - reaction free energy (kcal/mol)",/,&
      &1x,"Er  - reorganization energy (kcal/mol)",/,&
      &1x,"Ea  - activation energy (kcal/mol)",/,&
      &1x,"Ea* - exp(-Ea/kT) factor",/,&
      &1x,"V^2 - squared nonadiabatic coupling (kcal/mol)**2",/,&
      &1x,"Wij - partial rate (1/sec)")')
      write(log,'(1x,79("="))')

      do i=iprec,nprec

         !~~~~~~~~~~~~~~~~~
         ! Boltzmann weight
         !~~~~~~~~~~~~~~~~~
         defe = fei(i) - fezero
         roi = dexp(-beta*defe)/zf

         write(log,'(1x,"Precursor state:",I2,5X,"Boltzmann weight:",E20.9)') i,roi
         write(log,'(1x,"Minimum: ",3F12.3)') ZPI(I),ZEI(I),FEI(I)
         write(log,'(1x,79("-"))')
         write(log,'(1x,T2,"J",T9,"DeG",T20,"Er",T31,"Ea",T45,"Ea*",T60,"V^2",T75,"Wij")')

         wi = 0.d0
         do j=1,nsucc
            wi = wi + roi*wij(i,j)
            eaexp = dexp(-beta*ea(i,j))
            wtotal = wtotal + roi*wij(i,j)
            write(log,'(i2,3f11.6,3e15.6)') j,de(i,j),er(i,j),ea(i,j),eaexp,vc2(i,j),wij(i,j)*1.d12
         enddo

         write(log,'(1x,79("-"))')
         write(log,'(1x,"TOTAL RATE FROM STATE ",i2," (1/sec): ",e20.9)') i, wi*1.d12
         write(log,'(1x,79("="))')

      enddo

      write(log,'(1x,"TOTAL NONADIABATIC RATE  (1/sec): ",e20.9)') wtotal*1.d12
      write(log,'(1x,79("=")/)')

      !======================================================
      ! Store the rates for the first and last temperature to
      ! estimate an effective (apparent) activation energy
      !======================================================

      if (trange.or.tlist)&
         &write(99,'(3f10.3,3g20.10)')&
	 &temp,100.d0/temp,beta,wtotal*1.d12,dlog10(wtotal*1.d12),dlog(wtotal*1.d12)

      if (trange) then

         if (it.eq.1) then
            betalast = beta
            ratelast = dlog(wtotal*1.d12)
         endif

         if (it.eq.n_t) then
            betafirst = beta
            ratefirst = dlog(wtotal*1.d12)
         endif

      endif

   enddo ! loop over temperature range

   if (trange.or.tlist) close(99)

   if (trange) then

      !=====================================================
      ! Estimate an effective (apparent) activation energy
      ! and a prefactor
      !=====================================================
      eaeff = -(ratelast - ratefirst)/(betalast - betafirst)
      deltarate = dabs(ratelast-ratefirst)
      deltabeta = dabs(betalast-betafirst)
      aeff = dexp(ratefirst + betafirst*deltarate/deltabeta)
      write(log,'(1x,79("="))')
      write(log,'(1x,"Apparent activation energy (kcal/mol): ",f20.6)') eaeff
      write(log,'(1x,"Apparent prefactor (1/sec): ",g20.6)') aeff
      write(log,'(1x,79("=")/)')

   endif

   deallocate (fei,zpi,zei,fej,zpj,zej,fex,zpx,zex)
   deallocate (de,er,ea)
   deallocate (vc2,wij)
   deallocate (min_succ)

   !close (73)

   return

   contains

   !============================================================================
   subroutine alloc_psi&
             (ielst,iprec,nprec,nsucc,nzdim,nz,nvibprec,nvibsucc,nfedim,npsiga,ierr)
   !============================================================================
   ! Allocate arrays for wavefunctions and energies;
   ! calculate indexing arrays for electronic and vibrational quantum numbers
   !============================================================================
      integer, intent(in)  :: ielst, iprec, nprec, nsucc
      integer, intent(out) :: nzdim, nz, nvibprec, nvibsucc, nfedim, npsiga
      integer, intent(out) :: ierr

      integer, parameter :: nalloc = 19
      integer, dimension(nalloc) :: ierror
      integer :: i

      if (mgquant.eq.1.or.mgquant.eq.2) then
         nzdim = ielst*nprst*ngast
         nz = nzdim*nzdim
         nvibprec = nprec
         nvibsucc = nsucc
      elseif (mgquant.eq.3) then
         nzdim = ielst*nprst
         nz = nzdim*nzdim*npntsg
         nvibprec = (nprec-1)/ngast + 1
         nvibsucc = (nsucc-1)/ngast + 1
      endif

      nfedim = max(nprec,nsucc)
      npsiga = ielst*nprst*ngast*npntsg
      
      ierror = 0

      allocate (z_px(nz),                           stat=ierror(1))
      allocate (z_sx(nz),                           stat=ierror(2))
      allocate (fe(nfedim),                         stat=ierror(3))
      allocate (psiel_px(ielst,npnts,npntsg,ielst), stat=ierror(4))
      allocate (psipr_px(ielst,nprst,npnts,npntsg), stat=ierror(5))
      allocate (psiel_sx(ielst,npnts,npntsg,ielst), stat=ierror(6))
      allocate (psipr_sx(ielst,nprst,npnts,npntsg), stat=ierror(7))
      allocate (psiga_px(npsiga),                   stat=ierror(8))
      allocate (psiga_sx(npsiga),                   stat=ierror(9))
      allocate (enel  (ielst,npnts,npntsg),         stat=ierror(10))
      allocate (envib (ielst,nprst,npntsg),         stat=ierror(11))
      allocate (envibg(ielst,nprst,ngast),          stat=ierror(12))
      allocate (istel(nzdim),                       stat=ierror(13))
      allocate (istpr(nzdim),                       stat=ierror(14))

      if (mgquant.eq.1.or.mgquant.eq.2) then

         allocate (istga(nzdim),stat=ierror(15))
         ndabf = 0
         do i=1,ielst
            do np=1,nprst
               do ng=1,ngast
                  ndabf = ndabf + 1
                  istel(ndabf) = i
                  istpr(ndabf) = np
                  istga(ndabf) = ng
               enddo
            enddo
         enddo

      elseif (mgquant.eq.3) then

         allocate (istvb_px(iprec:nprec),stat=ierror(16))
         allocate (istvb_sx(1:nsucc),stat=ierror(17))
         allocate (istga_px(iprec:nprec),stat=ierror(18))
         allocate (istga_sx(1:nsucc),stat=ierror(19))
         ndabf = 0
         do i=1,ielst
            do np=1,nprst
               ndabf = ndabf + 1
               istel(ndabf) = i
               istpr(ndabf) = np
            enddo
         enddo

         do i=iprec,nprec
            istvb_px(i) = (i-1)/ngast + 1
            istga_px(i) = i - (istvb_px(i)-1)*ngast
         enddo

         do i=1,nsucc
            istvb_sx(i) = (i-1)/ngast + 1
            istga_sx(i) = i - (istvb_sx(i)-1)*ngast
         enddo

      endif

      ierr = 0
      do i=1,nalloc
         if (ierror(i).ne.0) ierr = ierr + 1
      enddo

      return
      
   end subroutine alloc_psi

   !============================================================================
   subroutine dealloc_psi(ierr)
   !============================================================================
   ! Allocate arrays for wavefunctions and energies;
   ! calculate indexing arrays for electronic and vibrational quantum numbers
   !============================================================================
      integer, intent(out) :: ierr

      integer, parameter :: nalloc = 19
      integer, dimension(nalloc) :: ierror
      integer :: i

      ierror = 0

      deallocate (z_px,     stat=ierror(1))
      deallocate (z_sx,     stat=ierror(2))
      deallocate (fe,       stat=ierror(3))
      deallocate (psiel_px, stat=ierror(4))
      deallocate (psipr_px, stat=ierror(5))
      deallocate (psiel_sx, stat=ierror(6))
      deallocate (psipr_sx, stat=ierror(7))
      deallocate (psiga_px, stat=ierror(8))
      deallocate (psiga_sx, stat=ierror(9))
      deallocate (enel,     stat=ierror(10))
      deallocate (envib,    stat=ierror(11))
      deallocate (envibg,   stat=ierror(12))
      deallocate (istel,    stat=ierror(13))
      deallocate (istpr,    stat=ierror(14))

      if (mgquant.eq.1.or.mgquant.eq.2) then
         deallocate (istga, stat=ierror(15))
      elseif (mgquant.eq.3) then
         deallocate (istvb_px, stat=ierror(16))
         deallocate (istvb_sx, stat=ierror(17))
         deallocate (istga_px, stat=ierror(18))
         deallocate (istga_sx, stat=ierror(19))
      endif

      ierr = 0
      do i=1,nalloc
         if (ierror(i).ne.0) ierr = ierr + 1
      enddo

      return
      
   end subroutine dealloc_psi

   !============================================================================
   subroutine crossp2(m_,n_,pzp0_,pze0_,szp0_,sze0_,xacc_,maxit_,zpx_,zex_,fex_)
   !============================================================================
   ! Locates a crossing point between two free energy surfaces (M,N)
   ! along the straight line connecting two points (PZP0,PZE0) and
   ! (SZP0,SZE0).
   ! The method of secants is used. Algorithm is based on the
   ! original program RTSEC from "Fortran recipes",
   ! (C) Copr. 1986-92 Numerical Recipes Software $$$.
   !============================================================================
      integer, intent(in)  :: m_, n_, maxit_
      real(8), intent(in)  :: pzp0_, pze0_, szp0_, sze0_, xacc_
      real(8), intent(out) :: zpx_, zex_, fex_

      integer :: nstates_px, nstates_sx, j_
      real(8)  :: deltazp, deltaze, x1, zp, ze, fe1, fe2
      real(8)  :: fl, x2, f, rtsec, xl, swap, deltax

      if (mgquant.eq.1.or.mgquant.eq.2) then
         nstates_px = m_
         nstates_sx = n_
      else
         nstates_px = (m_ - 1)/ngast + 1
         nstates_sx = (n_ - 1)/ngast + 1
      endif

      deltazp = szp0_ - pzp0_
      deltaze = sze0_ - pze0_

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! evaluate the function (energy splitting)
      ! at the endpoints (pzp0,pze0) and (szp0,sze0)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      x1 = 0.d0
      zp = pzp0_
      ze = pze0_
      call feszz2('DIAB2',1,zp,ze,nstates_px,nfedim,fe,nz,z_px,ndabf,&
                ielst,enel,envib,envibg,psiel_px,psipr_px,npsiga,psiga_px)
      fe1 = fe(m_)
      call feszz2('DIAB2',2,zp,ze,nstates_sx,nfedim,fe,nz,z_sx,ndabf,&
                ielst,enel,envib,envibg,psiel_sx,psipr_sx,npsiga,psiga_sx)
      fe2 = fe(n_)
      fl = fe2 - fe1

      x2 = 1.d0
      zp = szp0_
      ze = sze0_
      call feszz2('DIAB2',1,zp,ze,nstates_px,nfedim,fe,nz,z_px,ndabf,&
                ielst,enel,envib,envibg,psiel_px,psipr_px,npsiga,psiga_px)
      fe1 = fe(m_)
      call feszz2('DIAB2',2,zp,ze,nstates_sx,nfedim,fe,nz,z_sx,ndabf,&
                ielst,enel,envib,envibg,psiel_sx,psipr_sx,npsiga,psiga_sx)
      fe2 = fe(n_)
      f = fe2 - fe1

      if (dabs(fl).lt.dabs(f)) then
         rtsec = x1
         xl = x2
         swap = fl
         fl = f
         f = swap
      else
         xl = x1
         rtsec = x2
      endif

      !~~~~~~~~~~~~~~~~~
      ! Begin iterations
      !~~~~~~~~~~~~~~~~~

      do j_=1,maxit_

         deltax = (xl-rtsec)*f/(f-fl)

         xl = rtsec
         fl = f
         rtsec = rtsec + deltax
         zp = pzp0_ + rtsec*deltazp
         ze = pze0_ + rtsec*deltaze
         call feszz2('DIAB2',1,zp,ze,nstates_px,nfedim,fe,nz,z_px,ndabf,&
                   ielst,enel,envib,envibg,psiel_px,psipr_px,npsiga,psiga_px)
         fe1 = fe(m_)
         call feszz2('DIAB2',2,zp,ze,nstates_sx,nfedim,fe,nz,z_sx,ndabf,&
                   ielst,enel,envib,envibg,psiel_sx,psipr_sx,npsiga,psiga_sx)
         fe2 = fe(n_)
         f = fe2 - fe1

         if (dabs(deltax).lt.xacc_.or.f.eq.0.d0) then
            zpx_ = pzp0_ + rtsec*deltazp
            zex_ = pze0_ + rtsec*deltaze
            fex_ = (fe1 + fe2)/2.d0
            return
         endif

      enddo

      write(*,'(/1x,''rtsec exceeded maximum iterations ('',i3.3,'')''/)') maxit_
      zpx_ = pzp0_ + rtsec*deltazp
      zex_ = pze0_ + rtsec*deltaze
      fex_ = (fe1 + fe2)/2.d0

      return

   end subroutine crossp2

   !============================================================================
   function coupling(i_,j_) result(vij_)
   !============================================================================
   ! calculates coupling between two wavefunctions
   !============================================================================
      integer, intent(in) :: i_, j_
      real(8) :: vij_

      integer :: option
      integer :: k, l, mu, nu, m, n, ll, levb, mevb, km, ln
      integer :: kg, kp, ivb_px, jvb_sx, iga_px, jga_sx, kmu, lnu
      integer :: iel, imu, img, jel, jnu, jng
      
      real(8)  :: vklmunumn, vklmunumn_kg, psiga_ovl, psipr_ovl, el_ovl
      real(8)  :: vij_kg, vklmunu, h0ij, h0kl

      real(8), allocatable, dimension(:,:)     :: z2_px, z2_sx
      real(8), allocatable, dimension(:,:,:)   :: z3_px, z3_sx
      real(8), allocatable, dimension(:,:,:,:) :: psiga4_px, psiga4_sx
      real(8), allocatable, dimension(:,:,:)   :: psiga3_px, psiga3_sx

      select case(mgquant)
         case(1,2)
            allocate (z2_px(nzdim,nzdim))
            allocate (z2_sx(nzdim,nzdim))
            allocate (psiga4_px(ielst,nprst,ngast,npntsg))
            allocate (psiga4_sx(ielst,nprst,ngast,npntsg))
            z2_px = reshape(z_px,(/ndabf,ndabf/))
            z2_sx = reshape(z_sx,(/ndabf,ndabf/))
            psiga4_px = reshape(psiga_px,(/ielst,nprst,ngast,npntsg/))
            psiga4_sx = reshape(psiga_sx,(/ielst,nprst,ngast,npntsg/))
         case(3)
            allocate (z3_px(nzdim,nzdim,npntsg))
            allocate (z3_sx(nzdim,nzdim,npntsg))
            allocate (psiga3_px(nvibprec,ngast,npntsg))
            allocate (psiga3_sx(nvibsucc,ngast,npntsg))
            z3_px = reshape(z_px,(/ndabf,ndabf,npntsg/))
            z3_sx = reshape(z_sx,(/ndabf,ndabf,npntsg/))
            psiga3_px = reshape(psiga_px,(/nstates,ngast,npntsg/))
            psiga3_sx = reshape(psiga_sx,(/nstates,ngast,npntsg/))
      end select

      option = method*10 + mgquant

      vij_ = 0.d0

      select case(option)

         case(11,12)   ! method=1, mgquant=1/2

            do k=1,2    ! over electronic basis states 1a,1b
               do mu=1,nprst  ! over proton vibrational states within k
                  do m=1,ngast    ! over gating vibrational states within mu
                     do l=1,2          !over electronic basis states 2a,2b
                        do nu=1,nprst       ! over proton vibrational states within l
                           do n=1,ngast         ! over gating vibrational states within nu

                              vklmunumn = 0.d0

                              do kg=1,npntsg   ! integrate over gating coordinate
                                 vklmunumn_kg = 0.d0
                                 do kp=1,npnts   ! integrate over proton coordinate
                                    psipr_ovl = psipr_px(k,mu,kp,kg)*psipr_sx(l,nu,kp,kg)
                                    vklmunumn_kg = vklmunumn_kg + psipr_ovl*h0(k,l+2,kp,kg)
                                 enddo
                                 psiga_ovl = psiga4_px(k,mu,m,kg)*psiga4_sx(l,nu,n,kg)
                                 vklmunumn = vklmunumn + psiga_ovl*vklmunumn_kg
                              enddo

                              km = (k-1)*nprst*ngast + (mu-1)*ngast + m
                              ln = (l-1)*nprst*ngast + (nu-1)*ngast + n
                              vij_ = vij_ + z2_px(km,i_)*z2_sx(ln,j_)*vklmunumn

                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo

            deallocate (z2_px)
            deallocate (z2_sx)
            deallocate (psiga4_px)
            deallocate (psiga4_sx)

         case(13)   ! method=1, mgquant=3

            ivb_px = istvb_px(i_)
            jvb_sx = istvb_sx(j_)
            iga_px = istga_px(i_)
            jga_sx = istga_sx(j_)

            do kg=1,npntsg   ! integrate over gating coordinate

               vij_kg = 0.d0

               do k=1,2    ! over electronic basis states 1a,1b
                  do mu=1,nprst  ! over proton vibrational states within "k"
                     do l=1,2    ! over electronic basis states 2a,2b
                        do nu=1,nprst  ! over proton vibrational states within "l"

                           vklmunu = 0.d0
                           do kp=1,npnts   ! integrate over proton coordinate
                              psipr_ovl = psipr_px(k,mu,kp,kg)*psipr_sx(l,nu,kp,kg)
                              vklmunu = vklmunu + psipr_ovl*h0(k,l+2,kp,kg)
                           enddo

                           kmu = (k-1)*nprst + mu
                           lnu = (l-1)*nprst + nu
                           vij_kg = vij_kg + z3_px(kmu,ivb_px,kg)*z3_sx(lnu,jvb_sx,kg)*vklmunu

                        enddo
                     enddo
                  enddo
               enddo

               psiga_ovl = psiga3_px(ivb_px,iga_px,kg)*psiga3_sx(jvb_sx,jga_sx,kg)
               vij_ = vij_ + psiga_ovl*vij_kg

            enddo

            deallocate (z3_px)
            deallocate (z3_sx)
            deallocate (psiga3_px)
            deallocate (psiga3_sx)

         case(21:23)   ! method=2

            write(*,'(/1x,"*** METHOD=2 is not implemented (coupling is too complicated) ***"/)')
            stop

          case(31)   ! method=3, mgquant=1

            iel = istel(i_)
            jel = istel(j_)

            do k=1,2
               do mu=1,nprst
                  do m=1,ngast

                     do l=1,2
                        do nu=1,nprst
                           do n=1,ngast

                              vklmunumn = 0.d0
                              do kg=1,npntsg   ! integrate over gating coordinate

                                 vklmunumn_kg = 0.d0
                                 do kp=1,npnts   ! integrate over proton coordinate

                                    h0kl = 0.d0
                                    do levb=1,2     ! integrate over electronic coordinate
                                       do mevb=1,2
                                          el_ovl = psiel_px(iel,kp,kg,levb)*psiel_sx(jel,kp,kg,mevb)
                                          h0kl = h0kl + el_ovl*h0(levb,mevb+2,kp,kg)
                                       enddo
                                    enddo  ! loop over EVB states

                                    psipr_ovl = psipr_px(k,mu,kp,kg)*psipr_sx(l,nu,kp,kg)
                                    vklmunumn_kg = vklmunumn_kg + psipr_ovl*h0kl

                                 enddo  ! loop over proton grid
                                 
                                 psiga_ovl = psiga4_px(k,mu,m,kg)*psiga4_sx(l,nu,n,kg)
                                 vklmunumn = psiga_ovl*vklmunumn_kg

                              enddo  ! loop over gating grid

                              km = (k-1)*nprst*ngast + (mu-1)*ngast + m
                              ln = (l-1)*nprst*ngast + (nu-1)*ngast + n
                              vij_ = vij_ + z2_px(km,i_)*z2_sx(ln,j_)*vklmunumn

                           enddo
                        enddo
                     enddo

                  enddo
               enddo
            enddo

            deallocate (z2_px)
            deallocate (z2_sx)
            deallocate (psiga4_px)
            deallocate (psiga4_sx)

         case(32)   ! method=3, mgquant=2

            iel = istel(i_)
            imu = istpr(i_)
            img = istga(i_)
            jel = istel(j_)
            jnu = istpr(j_)
            jng = istga(j_)

            do kg=1,npntsg  ! integrate over gating grid

               vij_kg = 0.d0
               do kp=1,npnts  ! integrate over proton grid

                  h0ij = 0.d0
                  do levb=1,2      ! integrate over electronic coordinate
                     do mevb=1,2
                        el_ovl = psiel_px(iel,kp,kg,levb)*psiel_sx(jel,kp,kg,mevb)
                        h0ij = h0ij + el_ovl*h0(levb,mevb+2,kp,kg)
                     enddo
                  enddo     ! loop over EVB states

                  psipr_ovl = psipr_px(iel,imu,kp,kg)*psipr_sx(jel,jnu,kp,kg)
                  vij_kg = vij_kg + psipr_ovl*h0ij

               enddo  ! loop over proton grid
               
               psiga_ovl = psiga4_px(iel,imu,img,kg)*psiga4_sx(jel,jnu,jng,kg)
               vij_ = vij_ + psiga_ovl*vij_kg

            enddo   ! loop over gating grid

            deallocate (z2_px)
            deallocate (z2_sx)
            deallocate (psiga4_px)
            deallocate (psiga4_sx)

         case(33)   ! method=3, mgquant=3

            ivb_px = istvb_px(i_)
            jvb_sx = istvb_sx(j_)
            iel = istel(ivb_px)
            imu = istpr(ivb_px)
            img = istga_px(i_)
            jel = istel(jvb_sx)
            jnu = istpr(jvb_sx)
            jng = istga_sx(j_)

            do kg=1,npntsg  ! integrate over gating grid

               vij_kg = 0.d0
               do kp=1,npnts  ! integrate over proton grid

                  h0ij = 0.d0
                  do levb=1,2      ! integrate over electronic coordinate
                     do mevb=1,2
                        el_ovl = psiel_px(iel,kp,kg,levb)*psiel_sx(jel,kp,kg,mevb)
                        h0ij = h0ij + el_ovl*h0(levb,mevb+2,kp,kg)
                     enddo
                  enddo            ! loop over EVB states

                  psipr_ovl = psipr_px(iel,imu,kp,kg)*psipr_sx(jel,jnu,kp,kg)
                  vij_kg = vij_kg + psipr_ovl*h0ij
                  
               enddo  ! loop over proton grid
               
               psiga_ovl = psiga3_px(ivb_px,img,kg)*psiga3_sx(jvb_sx,jng,kg)
               vij_ = vij_ + psiga_ovl*vij_kg

            enddo   ! loop over gating grid

            deallocate (z3_px)
            deallocate (z3_sx)
            deallocate (psiga3_px)
            deallocate (psiga3_sx)

      end select

      return

   end function coupling

end subroutine rate2
