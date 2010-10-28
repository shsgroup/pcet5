subroutine rate3

!=======================================================================
!
!     Calculates the nonadiabatic rate of the PCET reaction in solution
!
!     Reference:
!     [*] Alexander Soudackov and Sharon Hammes-Schiffer,
!         Derivation of rate expressions for nonadiabatic
!         proton-coupled electron transfer reactions in solution,
!         J. Chem. Phys. 113, 2385 (2000).
!
!     The gating is taken into account according to the following
!     prescription:
!
!     k(T) = Int[k(R,T)g(R,T)dR]
!
!     where k(R,T) is a non-adiabatic rate constant at the fixed
!     gating distance R, and g(R) is a distribution function
!     for gating coordinate (classical or exact quantum depending
!     on the options specified).
!
!     If there is only one gating grid point then the standard
!     rate calculation will be performed for the gating distance
!     corresponding to the input geometry.
!
!     OPTIONS:
!
!     T=<value>[/<value>/<nsteps>] - absolute temperature
!                                    [range (init/fin/npoints)]
!
!     T=<value1>-<value2>-...-<valueN> - absolute temperature list
!                                        (up to 20 values)
!
!     GRQUANT - exact quantum distribution function for
!               harmonic gating motion will be used
!
!     GRCLASS - classical (Boltzmann) distribution function for
!               harmonic gating motion will be used
!
!     GRAVER  - classical (Boltzmann) explicit averaging over all
!               gating distances will be performed
!
!     GRHARM  - dynamical expressions with harmonic R-mode
!               (both low and high-T limits) will be used
!               for rate calculations
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
!     CROSS00 - use the crossing point for the pair of ground states
!               for all pairs when calculating the vibronic couplings
!
!     IPRINT=N - printing level in LBFGS minimization routine (default=-1)
!
!     MCORR=m  - number of corrections used in the limited memory matrix.
!             Values of m < 3  are not recommended, and large values
!             of m can result in excessive computing time. The range
!             3 <= m <= 20 is recommended (default=5)
!
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
!     VRPOUT  - output of Franck-Condon overlaps for all pairs of states
!               (Helene Decornez) into separate files
!
!     VERBOSE - print detailed rate info for each gating grid point.
!               If not specified only the total rate is printed out.
!               Not also that not specifying VERBOSE will cancel the
!               effects of WEIGHTS option.
!
!---- keywords for analytical rate calculations including R-mode effects
!
!     GMASS=<VALUE> - reduced mass for R-mode (Daltons)
!
!     GFREQ=<VALUE> - R-mode harmonic frequency (in 1/cm)
!
!     DELTAR=<VALUE> - Shift of equilibrium R-distance relative to the reactant states
!                      ( DeltaR = R0(product) - R0(reactant) )
!
!     ALPHA=<VALUE> - alpha parameter (common for all pairs) (in 1/A)
!
!     RATELIMIT=HIGH - calculate
!
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
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
   use geogas
   use fesmin_3d
   use feszz_3d, only: feszz3

   implicit real*8 (a-h,o-z)            !----- FIX THIS! Change to implicit none

   character(1024) :: options
   character( 40)  :: fname
   character(  5)  :: mode
   character(  4)  :: ctemp
   logical         :: weights, outlog, eanumer, eranaly, vrpout
   logical         :: filemin, rdmin, verbose, trange, tlist, cross00

   real(8), dimension(20) :: ti(20)

   logical, allocatable, dimension(:)    :: min_succ
   integer, allocatable, dimension(:)    :: istel, istpr
   integer :: iprin, imcorr, ifactr, ipgtol

   real*8, allocatable, dimension(:)     :: fe
   real*8, allocatable, dimension(:,:)   :: z_px, z_sx
   real*8, allocatable, dimension(:,:)   :: enel, envib
   real*8, allocatable, dimension(:,:,:) :: psiel_px, psipr_px
   real*8, allocatable, dimension(:,:,:) :: psiel_sx, psipr_sx
   real*8, allocatable, dimension(:,:)   :: zpi, zei, zpj, zej
   real*8, allocatable, dimension(:,:,:) :: zpx, zex
   real*8, allocatable, dimension(:,:)   :: fei, fej
   real*8, allocatable, dimension(:,:,:) :: fex, de, ea, ea_low, er, v2
   real*8, allocatable, dimension(:)     :: vrpkk

   real(8) :: zpij00, zeij00, zpxij, zexij

   real(8) :: alphamunu=0.d0, deltar=0.d0


   if (method.eq.2) then
      write(*,'(/1x,"*** method=2 is not implemented for rate calculations ***"/)')
      write(*,'( 1x,"    try again (rate3)..."/)')
      stop
   endif

   ! default
   rdmin   = .false.
   filemin = .false.
   trange  = .false.

   !====================================================
   ! Extract options
   !====================================================

   ikey = index(keywrd,' RATE3(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** in rate3: you must specify options for the RATE3 keyword ***"/)')
      stop
   else
      call getopt(keywrd,ikey+7,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VRPOUT output
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   vrpout = index(options,' VRPOUT').ne.0
   if (vrpout) then
      write(6,'(/1x,"two-dimensional overlap factors output")')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! VERBOSE output
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   verbose = index(options,' VERBOSE').ne.0

   if (verbose) then
      write(6,'(/1x,"rates output is verbose (for each gating grid point)")')
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

      write(6,'(/1x,"Rates output is in the external file <",a,">")') fname(1:lenf)

      log = 1
      open(log,file=fname(1:lenf),status='new')

   else

      log = 6

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! BANNER
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   write(log,'(/1x,70("="))')
   write(log,'( 1x,"Nonadiabatic Rate Calculation")')
   write(log,'( 1x,70("=")/)')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate distribution function
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   igrq = index(options,' GRQUANT')
   igrc = index(options,' GRCLASS')
   igra = index(options,' GRAVER')
   igrh = index(options,' GRHARM')

   if (igrq*igrc.ne.0) then
      write(*,'(/1x,"*** INPUT WARNING (in RATE3): ",/,&
      &"You MUST specify either GRQUANT or GRCLASS, NOT BOTH!",/,&
      &"It is confusing, so read manuals and make up your mind"/)')
      stop 'try again (rate3)...'
   endif

   if (igrc.ne.0) then
      igr = 1
      write(log,'(/1x,"Classical Boltzmann distribution for the harmonic gating")')
   elseif (igrq.ne.0) then
      igr = 2
      write(log,'(/1x,"Exact Quantum distribution for the harmonic gating")')
   elseif (igra.ne.0) then
      igr = 3
      write(log,'(/1x,"Classical Boltzmann distribution (explicit treatment)")')
   elseif (igrh.ne.0) then
      igr = 4
      write(log,'(/1x,"Dynamical expressions for harmonic R-mode (High- and Low-T limits)")')
      if (npntsg.gt.1) then
         write(*,'(/1x,"*** INCONSISTENT INPUT ERROR (in RATE3): ",/,&
         &"If you specify GRHARM, you MUST have ONLY A SINGLE GATING DISTANCE",/,&
         &"corresponding to the equilibrium distance in reactant.",/,&
         &"(In GATING group NGRIDS must be 1)."/)')
         stop 'read non-existing chapter in the manual and try again (rate3)...'
      endif
   else
      igr = 3     !default is the classical Boltzmann distribution
      write(log,'(/1x,"Classical Boltzmann distribution (explicit treatment)")')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate effective mass and frequency
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   igmass = index(options,' GMASS=')

   if (igmass.ne.0) then
      gmass = reada(options,igmass+7)
      gmass = gmass*dalton
   else
      gmass = dm*am/(dm+am)
   endif

   igfreq = index(options,' GFREQ=')

   if (igfreq.ne.0) then

      ! inverse centimeters
      gfreq = reada(options,igfreq+7)

   elseif (igr.ne.3) then

      write(*,'(/1x,"*** INPUT WARNING (in RATE3): ",/,&
      &"You MUST specify GFREQ for use with GRQUANT or GRCLASS.",/,&
      &"It is confusing, so read manuals and make up your mind"/)')
      stop

   else

      gfreq = 0.d0

   endif
   

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gating coordinate average value
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   igaver = index(options,' GAVER=')

   if (npntsg.le.1) then
      midpoint = 1
   else
      midpoint = npntsg/2
   endif

   if (igaver.ne.0) then
      gaver = reada(options,igaver+7)
   else
      gaver = glist(midpoint)*bohr2a
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Shift of the equilibrium R-distance (Product minus Reactant)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ikey = index(options,' DELTAR=')
   if (ikey.ne.0) then
      deltar = reada(options,ikey+8)
   else
      deltar = 0.d0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Coupling parameter alpha
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ikey = index(options,' ALPHA=')
   if (ikey.ne.0) then
      alphamunu = reada(options,ikey+7)
   else
      alphamunu = 0.d0
   endif

   !-- calculate lambda_R and lambda_alpha (in case of GRHARM)

   if (igr.eq.4) then

      !-- lambda_R
      gmass_au = gmass
      gfreq_au = gfreq*cm2ev*ev2au
      deltar_au = deltar*a2bohr
      er_r = (0.5*gmass_au*gfreq_au*gfreq_au*deltar_au*deltar_au)*au2cal
      write(log,*)
      write(log,*) "R-mode reorganization energy (kcal/mol): ",er_r

      !-- lambda_alpha
      alphamunu_au = alphamunu/a2bohr
      er_alpha = (0.5d0*alphamunu_au*alphamunu_au/gmass_au)*au2cal
      write(log,*)
      write(log,*) "lambda_alpha (coupling) reorganization energy (kcal/mol): ",er_alpha

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Precursor states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   iiprec = index(options,' IPREC=')
   inprec = index(options,' NPREC=')

   if (iiprec.ne.0.and.inprec.ne.0) then
      write(*,'(/1x,"*** INPUT ERROR (in RATE3): ",/,&
      &"You MUST specify either IPREC or NPREC, NOT BOTH!",/,&
      &"It is confusing, so read manuals and make up your mind"/)')
      stop
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
   write(log,'(/1x,"Number of successor states:",i3)') nsucc

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Allocate some arrays
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   allocate (zpi(npntsg,iprec:nprec), zei(npntsg,iprec:nprec))                   ! precursors minima
   allocate (zpj(npntsg,1:nsucc), zej(npntsg,1:nsucc))                           ! successors minima
   allocate (zpx(npntsg,iprec:nprec,1:nsucc), zex(npntsg,iprec:nprec,1:nsucc))   ! crossing points

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! File containing minima for calculation of the rates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   filemin = index(options,' READMIN').ne.0
   irdm = index(options,' READMIN=')

   if (.not.filemin) then

      write(*,'(/1x,"Minima will be determined in calculation"/)')
      rdmin=.false.

   else

      rdmin=.true.
      ispa = index(options(irdm+9:),' ')
      fname = options(irdm+9:irdm+ispa+7)
      lenf = ispa -1
      call locase (fname,lenf)
      open(2,file=fname(1:lenf),status='old')
      call getmin(2,nprec-iprec+1,nsucc,zpi(1,iprec:nprec),zei(1,iprec:nprec),zpj(1,1:nsucc),zej(1,1:nsucc))
      close(2)
      do kg=2,npntsg
         zpi(kg,iprec:nprec) = zpi(1,iprec:nprec)
         zei(kg,iprec:nprec) = zei(1,iprec:nprec)
         zpj(kg,1:nsucc) = zpj(1,1:nsucc)
         zej(kg,1:nsucc) = zej(1,1:nsucc)
      enddo
      write(6,'(/1x,"Minima will be read from the external file <",a,">")') fname(1:lenf)

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
            write(*,*) " Error in temperature range specification..."
            write(*,*) " (from rate3)"
            stop
         endif

         t_start = reada(options(itemp_start:islash1-1),1)
         t_end   = reada(options(islash1+1:islash2-1),1)
         n_t = reada(options,islash2+1)
         if (n_t.le.1) then
            write(*,*) " Error in temperature range specification..."
            write(*,*) " You should specify more than one point"
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

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Reactant minimum guess (ZPPREC0,ZEPREC0)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      izprec1 = izprec+7
      islash = index(options(izprec1:),'/')
      zpprec0 = reada(options(izprec1:izprec1+islash-2),1)
      zeprec0 = reada(options,izprec1+islash)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Product minimum (ZPSUCC0,ZESUCC0)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      izsucc1 = izsucc+7
      islash = index(options(izsucc1:),'/')
      zpsucc0 = reada(options(izsucc1:izsucc1+islash-2),1)
      zesucc0 = reada(options,izsucc1+islash)

      write(log,'(/1x,"Initial guesses for minima (kcal/mol):")')
      write(log,'( 1x,"Reactant=(",f8.3,",",f8.3,") and Product=(",f8.3,",",f8.3,")")')&
      zpprec0,zeprec0,zpsucc0,zesucc0

   else

      write(*,'(/1x,"*** INPUT ERROR (in RATE): ",/,&
      &"You MUST specify both ZPREC and ZSUCC options",&
      &" for RATE3 keyword ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for the activation energies (EANUMER)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eanumer = index(options,' EANUMER').ne.0
   if (igr.eq.4) eanumer=.false.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for the reorganization energies (ERANALY)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eranaly = index(options,' ERANALY').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Common crossing point for all pairs (CROSS00)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   cross00 = index(options,' CROSS00').ne.0

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

      write(log,'(/1x,"Maximum length of the Newton-Raphson step: ",g12.6)') slim

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

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations in searching a crossing point
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxitx = index(options,' MAXITX=')
   if (imaxitx.ne.0) then
      maxitx = reada(options,imaxitx+8)
   else
      maxitx = 100
   endif
   write(log,'(1x,"Maximum number of steps in locating a crossing point: ",i5)') maxitx

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Scaling factors for the solvent coordinates (ZPS,ZES)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

   ilin = index(options,' ALIN=')
   if (ilin.ne.0) then
      alin = reada(options,ilin+6)
      write(log,'(1x,"Inner sphere reorganization energy added (kcal/mol): ",g12.6)') alin
   else
      alin = 0.d0
   endif

   !====================================================
   ! START THE CALCULATION
   !====================================================

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! The following is fixed in this calculation.
   ! It corresponds to the diabatic representation
   ! with respect to the ET reaction
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   mode = 'DIAB2'
   ielst = 2
   nzdim = ielst*nprst

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Some quantities needed for analytical
   ! evaluation of reorganization energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   denom = 2.d0*(erpt*eret-erx*erx)
   elxx = eret/denom
   elyy = erpt/denom
   elxy = -erx/denom

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Indices of electronic and vibrational basis states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   allocate (istel(nzdim),istpr(nzdim))
   istel = 0
   istpr = 0

   ndabf = 0
   do i=1,ielst
      do np=1,nprst
         ndabf = ndabf + 1
         istel(ndabf) = i
         istpr(ndabf) = np
      enddo
   enddo

   !--------------------------------------
   ! Allocate arrays
   !--------------------------------------
   allocate (min_succ(nsucc))

   if (method.eq.3.and.vrpout) allocate (vrpkk(npnts))

   allocate (fe(nzdim),z_px(nzdim,nzdim),z_sx(nzdim,nzdim))
   allocate (psiel_px(ielst,npnts,ielst),psipr_px(ielst,nprst,npnts))
   allocate (psiel_sx(ielst,npnts,ielst),psipr_sx(ielst,nprst,npnts))
   allocate (enel(ielst,npnts),envib(ielst,nprst))

   allocate (fei(npntsg,iprec:nprec))           ! precursors equilibrium free energies
   allocate (fej(npntsg,1:nsucc))               ! successors equilibrium free energies
   allocate (fex(npntsg,iprec:nprec,1:nsucc))   ! free energies at the intersection points
   allocate ( de(npntsg,iprec:nprec,1:nsucc))   ! reaction free energies (delta G's)
   allocate ( er(npntsg,iprec:nprec,1:nsucc))   ! reorganization energies
   allocate ( ea(npntsg,iprec:nprec,1:nsucc))   ! activation energies
   if (igr.eq.4) allocate ( ea_low(npntsg,iprec:nprec,1:nsucc))   ! activation energies for low-T rate expression
   allocate ( v2(npntsg,iprec:nprec,1:nsucc))   ! squared couplings

   !---------------------------------------
   ! Start the loop over gating grid points
   !---------------------------------------

   totalrate = 0.d0
   if (npntsg.le.1) then
      gint = 1.d0
   else
      gint = dabs(glist(2) - glist(1))
   endif

   open(75,file=job(1:ljob)//'/fminprec.dat')
   write(75,'("#------------------------------------------------------------------------")')
   write(75,'("# Minimum energy of the precursor states as a function of gating distance")')
   write(75,'("#")')
   write(75,'("# Columns:")')
   write(75,'("#")')
   write(75,'("#   R    Energy(i)   zp_min(i)    ze_min(i),   i=1,nprec")')
   write(75,'("#------------------------------------------------------------------------")')
   write(75,'("#")')

   open(76,file=job(1:ljob)//'/fminsucc.dat')
   write(76,'("#------------------------------------------------------------------------")')
   write(76,'("# Minimum energy of the successor states as a function of gating distance")')
   write(76,'("#")')
   write(76,'("# Columns:")')
   write(76,'("#")')
   write(76,'("#   R    Energy(i)   Zp(i)    Ze(i), i=1,NSUCC")')
   write(76,'("#------------------------------------------------------------------------")')
   write(76,'("#")')

   open(77,file=job(1:ljob)//'/fmincros.dat')
   write(77,'("#------------------------------------------------------------------------")')
   write(77,'("# Energies and couplings at the crossing points")')
   write(77,'("#")')
   write(77,'("# Columns:")')
   write(77,'("#")')
   write(77,'("#   R  Energy(i,j) V^2(i,j) zp(i,j) ze(i,j),"," j=1,nsucc,i=1,nprec")')
   write(77,'("#------------------------------------------------------------------------")')
   write(77,'("#")')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Start loop over the gating coordinate grid
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do kg=1,npntsg

      rr = glist(kg)*bohr2a

      write(log,'(/1x,70("*"))')
      write(log,'( 1x,"Gating distance (A): ",f10.6)') rr
      write(log,'( 1x,70("*"))')

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Loop over the initial states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      min_succ = .false.

      do i=iprec,nprec

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Locate a minimum on the reactant surface I
         ! or calculate free energy if minima read in
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (.not.rdmin) then

            write(6,'(/"  search of the minumum for reactant state ",i3/)') i
            if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"h(1,1)",t91,"h(2,2)",t106,"h(1,2)")')
            call fes3gmin(mode,1,i,kg,zpprec0,zps,zeprec0,zes,&
                         &ierr,zpprec,zeprec,&
                         &feprec,dzpmin,dzemin,d2zpmin,d2zemin,&
                         &d2zpzemin)

            if (ierr.eq.0) then
               write(6,'(1x,a4,2x,8f15.6)') "done",&
	       &zpprec,zeprec,feprec,sqrt(dzpmin*dzpmin+dzemin*dzemin)  !,d2zpmin,d2zemin,d2zpzemin
            else
               write(6,'(1x,a4,2x,8f15.6)') "FAIL",&
	       &zpprec,zeprec,feprec,sqrt(dzpmin*dzpmin+dzemin*dzemin)  !,d2zpmin,d2zemin,d2zpzemin
	    endif

            if (ierr.eq.0) then
               zpprec0 = zpprec
               zeprec0 = zeprec
            endif
            zpi(kg,i) = zpprec
            zei(kg,i) = zeprec
            fei(kg,i) = feprec

         else

            fe = 0.d0
            z_px = 0.d0
            enel = 0.d0
            envib = 0.d0
            psiel_px = 0.d0
            psipr_px = 0.d0
            call feszz3(mode,1,kg,zpi(kg,i),zei(kg,i),nzdim,fe,&
                       &nzdim,z_px,ndabf,ielst,enel,envib,&
                       &psiel_px,psipr_px)
            zpprec = zpi(kg,i)
            zeprec = zei(kg,i)
            fei(kg,i) = fe(i)

         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate and print out the EVB weights
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (weights.and.verbose) then
            write(log,'(/1x,70("="))')
            write(log,'(1x,"EVB weights at the minimum of the ",i2,a3," reactant state")') i,th(i)
            write(log,'(1x,"zp(precursor)=",f9.3,5x,"ze(precursor)=",f9.3/)') zpprec,zeprec
            call weight3(log,kg,zpprec,zeprec,mode,1,i)
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Loop over the successor (product) states
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do j=1,nsucc

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Locate a minimum on the product surface J
            ! or calculate free energy if minima read in
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (.not.min_succ(j)) then

               if (.not.rdmin) then

                  write(6,'(/" Search of minumum for product state ",i3/)') j
                  if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"h(1,1)",t91,"h(2,2)",t106,"h(1,2)")')
                  call fes3gmin(mode,2,j,kg,zpsucc0,zps,zesucc0,zes,&
                               &ierr,zpsucc,zesucc,&
                               &fesucc,dzpmin,dzemin,d2zpmin,d2zemin,&
                               &d2zpzemin)

                  if (ierr.eq.0) then
                     write(6,'(1x,a4,2x,8f15.6)') "done",&
	             &zpsucc,zesucc,fesucc,sqrt(dzpmin*dzpmin+dzemin*dzemin)  !,d2zpmin,d2zemin,d2zpzemin
                  else
                     write(6,'(1x,a4,2x,8f15.6)') "FAIL",&
	             &zpsucc,zesucc,fesucc,sqrt(dzpmin*dzpmin+dzemin*dzemin)  !,d2zpmin,d2zemin,d2zpzemin
	          endif

                  if (ierr.eq.0) then
                     zpsucc0 = zpsucc
                     zesucc0 = zesucc
                  endif
                  zpj(kg,j) = zpsucc
                  zej(kg,j) = zesucc
                  fej(kg,j) = fesucc

               else

                  fe = 0.d0
                  z_px = 0.d0
                  enel = 0.d0
                  envib = 0.d0
                  psiel_px = 0.d0
                  psipr_px = 0.d0
                  call feszz3(mode,2,kg,zpj(kg,j),zej(kg,j),&
                             &nzdim,fe,nzdim,z_px,ndabf,ielst,enel,envib,&
                             &psiel_px,psipr_px)
                  zpsucc = zpj(kg,j)
                  zesucc = zej(kg,j)
                  fej(kg,j) = fe(j)

               endif

               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! Calculate and print out the EVB weights
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if (weights.and.verbose) then
                  write(log,'(/1x,70("="))')
                  write(log,'(1x,"EVB weights at the minimum of the ",i2,a3," product state")') j,th(j)
                  write(log,'(1x,"zp(successor)=",f9.3,5x,"ze(successor)=",f9.3/)') zpsucc,zesucc
                  call weight3(log,kg,zpsucc,zesucc,mode,2,j)
               endif

               min_succ(j) = .true.

            endif


            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate reaction free energy for the pair (I,J)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            de(kg,i,j) = fej(kg,j) - fei(kg,i)

            if (verbose)&
            write(log,'(/1x,"The reaction free energy (kcal/mol) for the pair ",2i2,":",f20.6)')&
            i, j, de(kg,i,j)

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate two-dimensional reorganization energy
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            ! Analytical reorganization energy

            zpij = zpj(kg,j) - zpi(kg,i)
            zeij = zej(kg,j) - zei(kg,i)
            elpt = elxx*zpij + elxy*zeij
            elet = elxy*zpij + elyy*zeij
            erij_a = elpt*elpt*erpt + elet*elet*eret + 2.d0*elpt*elet*erx
            erij_a = erij_a + alin

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Numerical reorganization energy
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate the free energy of the reactant state I
            ! at the minimum of the product surface J
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            fe = 0.d0
            z_px = 0.d0
            enel = 0.d0
            envib = 0.d0
            psiel_px = 0.d0
            psipr_px = 0.d0
            call feszz3(mode,1,kg,zpsucc,zesucc,nzdim,fe,nzdim,z_px,&
                       &ndabf,ielst,enel,envib,psiel_px,psipr_px)
            erij_n = fe(i) - fei(kg,i) + alin

            if (verbose)&
            write(log,'(/1x,"The reorganization energy (kcal/mol) for the pair ",2i2,":",/,&
            &1x,"analytical:",f20.6,/,&
            &1x,"numerical: ",f20.6)') i, j, erij_a, erij_n

            if (eranaly) then
               er(kg,i,j) = erij_a
               if (verbose) write(log,'(/1x,"analytical value will be used")')
            else
               er(kg,i,j) = erij_n
               if (verbose) write(log,'(/1x,"numerical value will be used")')
            endif

            if (er(kg,i,j).lt.0.d0) then
               write(*,'(/1x,"the reorg. energy for the pair ",2i2," is negative:",f20.6)')&
               &i,j,er(kg,i,j)
               write(*,'("this is very suspicious, so the program terminates...")')
               stop
            endif


            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Find the crossing point between electronically
            ! diabatic one-dimensional (slices along the straiht
            ! line connecting two minima) free-energy curves
            ! for the states I and J
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            call crossp3(i,j,kg,zpprec,zeprec,zpsucc,zesucc,xacc,maxitx,zpij,zeij,feij)

            !-- store the crossing point for the 0-0 pair
            if (i.eq.1.and.j.eq.1) then
               zpij00 = zpij
               zeij00 = zeij
            endif

            zpx(kg,i,j) = zpij
            zex(kg,i,j) = zeij
            fex(kg,i,j) = feij

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate and print out the EVB weights
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (weights.and.verbose) then
               write(log,'(/1x,70("="))')
               write(log,'(1x,"EVB weights at the crossing point between ",i2,a3,&
               &" reactant and ",i2,a3," product states")') i,th(i),j,th(j)
               write(log,'(1x,"zp(x-ing point)=",f9.3,5x,"ze(xing point)=",f9.3/)') zpij, zeij
               call weight3(log,kg,zpij,zeij,mode,1,i)
               call weight3(log,kg,zpij,zeij,mode,2,j)
            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate activation free energy for the pair (I,J)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            eaijn = fex(kg,i,j) - fei(kg,i)

            if (igr.ne.4) then
               eaija = (er(kg,i,j)+de(kg,i,j))*(er(kg,i,j)+de(kg,i,j))/4.d0/er(kg,i,j)
	    else
	       er_tot = er(kg,i,j)+er_r+er_alpha
	       eaija_high = (er_tot + de(kg,i,j))**2/4.d0/er_tot           !- withour temperature dependent terms
	       eaija_low  = (er(kg,i,j) + de(kg,i,j))**2/4.d0/er(kg,i,j)
	    endif

            if (verbose) then
	       if (igr.ne.4) then
                  write(log,'(/1x,"The activation energy for the pair of states ",i2,"-",i2,":",/,&
                  &" numerical (from the crossing point): ",f12.6,/,&
                  &" analytical (from Marcus expression): ",f12.6/)') i, j, eaijn, eaija
               else
                  write(log,'(/1x,"The activation energy for the pair of states ",i2,"-",i2,":",/,&
                  &" numerical (from the crossing point): ",f12.6,/,&
                  &" analytical (Low temperature limit) : ",f12.6,/,&
                  &" analytical (High temperature limit): ",f12.6/)') i, j, eaijn, eaija_high, eaija_low
	       endif
	    endif

            if (eanumer) then
               ea(kg,i,j) = eaijn   ! numerical estimate
            else
	       if (igr.ne.4) then
                  ea(kg,i,j) = eaija   ! Marcus relation
	       else
                  ea(kg,i,j) = eaija_high      ! High-T for R-mode
                  ea_low(kg,i,j) = eaija_low   ! Low-T for R-mode
	       endif
            endif

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate coupling matrix element at the crossing point
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            !-- CROS00 option: reset the crossing point to 0-0 pair

            if (cross00) then
               zpxij = zpij00
               zexij = zeij00
            endif

            call feszz3(mode,1,kg,zpxij,zexij,nzdim,&
                       &fe,nzdim,z_px,ndabf,ielst,enel,envib,&
                       &psiel_px,psipr_px)

            call feszz3(mode,2,kg,zpxij,zexij,nzdim,&
                       &fe,nzdim,z_sx,ndabf,ielst,enel,envib,&
                       &psiel_sx,psipr_sx)

            !=============HYD OUTPUT V(R_P)=================

            if (method.eq.3.and.vrpout) then

               write(fname,'("vrpout-",i2.2,"-",i2.2)') i,j
               nij = 100+10*i+j
               open(nij,file=job(1:ljob)//'/'//trim(fname),position='append')

               iel = istel(i)
               imu = istpr(i)
               jel = istel(j)
               jnu = istpr(j)

               fcad = 0.d0

               do kk=1,npnts

                  vrpkk(kk) = 0.d0

                  do levb=1,2
                     do mevb=1,2
                        psiel_ovl = psiel_px(iel,kk,levb)*psiel_sx(jel,kk,mevb)
                        vrpkk(kk) = vrpkk(kk) + psiel_ovl*h0(levb,mevb+2,kk,kg)
                     enddo
                  enddo

                  if (kk.eq.npnts) then
                     write(*,*) 'coefficients 1a,1b,2a,2b',&
                     &psiel_px(iel,kk,1),psiel_px(iel,kk,2),&
                     &psiel_sx(jel,kk,1),psiel_sx(jel,kk,2)
                  endif
                  
                  psipr_ovl = psipr_px(iel,imu,kk)*psipr_sx(jel,jnu,kk)

                  write(nij,'(6g20.6)') rlist(kk)*bohr2a,&
                                      & glist(kg)*bohr2a,&
                                      & vrpkk(kk),&
                                      & psipr_px(iel,imu,kk),&
                                      & psipr_sx(jel,jnu,kk),&
                                      & vrpkk(kk)*psipr_ovl

                  fcad = fcad + psipr_ovl

               enddo

               write(nij,*)
               close(nij)

               write(log,'(/1x,"FC overlap for the pair ",i2,"-",i2,": ",g20.10/)') i, j, fcad

            endif

            !==============END HYD OUTPUT V(R_P)==========Dec 6 2000======


            ! Actual coupling calculation

            vij = coupling(i,j)
            v2(kg,i,j) = vij*vij

         enddo         ! loop over i (precursors)

      enddo          ! loop over j (successors)

      !====================================================
      ! Initial states characteristics
      !====================================================
      if (verbose) then
         write(log,'(/1x,"INITIAL STATES CHARACTERISTICS (kcal/mol)")')
         write(log,'(1x,70("="))')
         write(log,'(1x,t4,"state",t17,"zp0",t32,"ze0",t43,"fe(zp0,ze0)")')
         write(log,'(1x,70("-"))')
         do i=iprec,nprec
           write(log,'(1x,i5,5x,3(f10.3,5x))') i,zpi(kg,i),zei(kg,i),fei(kg,i)
         enddo
         write(log,'(1x,70("-"))')
      endif
      write(75,'(30g20.6:)') rr,(fei(kg,i),zpi(kg,i),zei(kg,i),i=iprec,nprec)

      !====================================================
      ! Crossing points
      !====================================================
      if (verbose) then
         write(log,'(/1x,"CROSSING POINTS CHARACTERISTICS (kcal/mol)")')
         write(log,'(1x,70("="))')
         write(log,'(1x,t3,"pair",t14,"zpx",t29,"zex",t39,"fex(zp0,ze0)",t59,"v^2")')
         write(log,'(1x,70("-"))')
         do i=iprec,nprec
            do j=1,nsucc
               write(log,'(1x,i2,"-",i2,3(f10.3,5x),g12.6)')&
               &i,j,zpx(kg,i,j),zex(kg,i,j),fex(kg,i,j),v2(kg,i,j)
            enddo
         enddo
         write(log,'(1x,70("-"))')
      endif

      write(77,'(100g15.6)')&
      &rr,((fex(kg,i,j),v2(kg,i,j),zpx(kg,i,j),zex(kg,i,j),j=1,nsucc),i=iprec,nprec)

      !====================================================
      ! Final states characteristics
      !====================================================
      if (verbose) then
         write(log,'(/1x,"FINAL STATES CHARACTERISTICS (kcal/mol)")')
         write(log,'(1x,70("="))')
         write(log,'(1x,t4,"state",t17,"zp0",t32,"ze0",t43,"fe(zp0,ze0)")')
         write(log,'(1x,70("-"))')
         do j=1,nsucc
            write(log,'(1x,i5,5x,3(f10.3,5x))') j,zpj(kg,j),zej(kg,j),fej(kg,j)
         enddo
         write(log,'(1x,70("-"))')
      endif
      write(76,'(30g20.6:)') rr,(fej(kg,j),zpj(kg,j),zej(kg,j),j=1,nsucc)

   enddo !loop over gating grid

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! End of loop over gating coordinate grid
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !==================================================================
   ! Start loop over temperatures
   !==================================================================

   if (trange.or.tlist) then
      open(99,file=job(1:ljob)//"/arrhenius_c.dat")
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("# Arrhenius plots: averaged rates in 1/sec")')
      write(99,'("#")')
      write(99,'("# 6 Columns:")')
      write(99,'("#")')
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("#    1           2        3            4                   5                    6     ")')
      write(99,'("#-------------------------------------------------------------------------------------")')
      write(99,'("#   T(K)       100/T     1/kT      k_tot(1/sec)        log10(k_tot)          ln(k_tot)")')
      write(99,'("#-------------------------------------------------------------------------------------")')
   endif


   if (igr.eq.4) then

      if (trange.or.tlist) then

         open(98,file=job(1:ljob)//"/arrhenius_rharm_high.dat")
         write(98,'("#-------------------------------------------------------------------------------------")')
         write(98,'("# Arrhenius plots: Dynamical rates (R-harm, High-T) in 1/sec")')
         write(98,'("#")')
         write(98,'("# 6 Columns:")')
         write(98,'("#")')
         write(98,'("#-------------------------------------------------------------------------------------")')
         write(98,'("#    1           2        3            4                   5                    6     ")')
         write(98,'("#-------------------------------------------------------------------------------------")')
         write(98,'("#   T(K)       100/T     1/kT      k_tot(1/sec)        log10(k_tot)          ln(k_tot)")')
         write(98,'("#-------------------------------------------------------------------------------------")')

         open(97,file=job(1:ljob)//"/arrhenius_rharm_low.dat")
         write(97,'("#-------------------------------------------------------------------------------------")')
         write(97,'("# Arrhenius plots: Dynamical rates (R-harm, Low-T) in 1/sec")')
         write(97,'("#")')
         write(97,'("# 6 Columns:")')
         write(97,'("#")')
         write(97,'("#-------------------------------------------------------------------------------------")')
         write(97,'("#    1           2        3            4                   5                    6     ")')
         write(97,'("#-------------------------------------------------------------------------------------")')
         write(97,'("#   T(K)       100/T     1/kT      k_tot(1/sec)        log10(k_tot)          ln(k_tot)")')
         write(97,'("#-------------------------------------------------------------------------------------")')

      endif

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

      write(ctemp,'(i4.4)') int(temp)
      open(74,file=job(1:ljob)//"/rateofr"//ctemp//".dat")

      if (igr.eq.1.or.igr.eq.2) then

         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("# non-adiabatic rate as a function of gating distance")')
         write(74,'("#")')
         write(74,'("# harmonic approximation for the gating potential")')
         write(74,'("# k_tot = int[dr*g(r)*k(r)]")')
         write(74,'("# g(r) is a quantum distribution function for harmonic oscillator:")')
         write(74,'("# oscillator mass (daltons)  :",f8.3)') gmass/dalton
         write(74,'("# oscillator frequency (1/cm):",f8.3)') gfreq
         write(74,'("# equilibrium r distance (a) :",f8.3)') gaver
         write(74,'("#")')
         write(74,'("# 5 columns:")')
         write(74,'("#")')
         write(74,'("#   r    g(r)   k(r)  g(r)k(r)dr   k_tot")')
         write(74,'("#-----------------------------------------------------------------")')

      elseif (igr.eq.3) then

         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("# Non-adiabatic rate as a function of gating distance")')
         write(74,'("#")')
         write(74,'("# explicit averaging over gating coordinate")')
         write(74,'("# k_tot = sum_i [ int [ dr*p(i,r)*sum_f [ k_if(r) ] ] ]")')
         write(74,'("# p(i,r) is a distribution function for reactant state i:")')
         write(74,'("#")')
         write(74,'("# 3 columns:")')
         write(74,'("#")')
         write(74,'("# r sum_i[p(i,r)*sum_f[k_if(r)]] k(tot)")')
         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("#")')

      elseif (igr.eq.4) then

         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("# Non-adiabatic rate as a function of gating distance (harmonic R-mode)")')
         write(74,'("#")')
         write(74,'("# 3 columns:")')
         write(74,'("#")')
         write(74,'("# r  k(tot,High-T)  k(tot,low-T)")')
         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("#")')

      else

         write(74,'("#")')

      endif

      write(log,'(/1x,70("*"))')
      write(log,'( 1x,"RATE ANALYSIS AT ",f8.3," K")') temp
      write(log,'(/1x,70("*"))')

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Preexponential factor (constant) in the rate expression
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      prefac = dsqrt(pi*beta)/hbarps

      if (igr.eq.4) then
         gfreq_kcal = gfreq_au*au2cal
         prefac_high = prefac*dexp(4.d0*er_alpha/beta/gfreq_kcal/gfreq_kcal)
	 prefac_low  = prefac*dexp(-alphamunu*deltar + (er_alpha-er_r)/gfreq_kcal)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ZFT - total partition function summed over gating
      !       coordinate grid and all the initial states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zft = 0.d0

      !--------------------------------------
      ! find lowest free energy and
      ! calculate total partition function
      !--------------------------------------

      !fezerot = fei(1,1)
      fezerot = minval(fei)

      do kg=1,npntsg
         do i=iprec,nprec
            !if (fei(kg,i).lt.fezerot) fezerot = fei(kg,i)
            zft = zft + dexp(-beta*(fei(kg,i)-fezerot))
         enddo
      enddo

      !--------------------------------------
      ! adjust total partition function
      !--------------------------------------
      !zft = zft*dexp(beta*fezerot)

      !----------------------------------------
      ! Total rate over the entire gating range
      !----------------------------------------
      totalrate = 0.d0

      !====================================================
      ! Rate analysis for fixed gating distances
      !====================================================

      do kg=1,npntsg

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! ZFR - partition function for proton vibrational states
         !      at the fixed gating distance
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         zfr = 0.d0
         kzfr = 0
	 fezero = minval(fei(kg,iprec:nprec))

         do i=iprec,nprec
            kzfr = kzfr + 1
            !if (kzfr.eq.1) fezero = fei(kg,i)
            zfr = zfr + dexp(-beta*(fei(kg,i)-fezero))
         enddo
         !zfr = zfr*dexp(beta*fezero)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Analytical distribution function g(R)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (npntsg.le.1) then
            rr = dabs(xyzgas(1,iptgas(1))-xyzgas(1,iptgas(3)))
            gr = 1.d0
         else
            rr = glist(kg)*bohr2a
            gr = gofr(igr,gmass,gfreq,temp,gaver,rr)
         endif

         wtotal = 0.d0
	 wtotal_low = 0.d0
         wr = 0.d0
	 wr_low = 0.d0

         if (verbose) then
            write(log,'(/1x,79("="))')
            write(log,'(1x,"RATE CONTRIBUTIONS, R=",f8.3," A")') rr
            write(log,'(1x,79("="))')
            write(log,'(1x,"column headers:",/,&
            &1X,"J   - successor state",/,&
            &1X,"DeG - reaction free energy (kcal/mol)",/,&
            &1X,"Er  - reorganization energy (kcal/mol)",/,&
            &1X,"Ea  - activation energy (kcal/mol)",/,&
            &1X,"Ea* - exp(-Ea/kT) factor",/,&
            &1X,"V^2 - squared nonadiabatic coupling (kcal/mol)**2",/,&
            &1X,"Wij - partial rate (1/sec)")')
            write(log,'(1x,79("="))')
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Averaging over initial states
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do i=iprec,nprec

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Boltzmann weight at fixed gating coordinate
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            defe = fei(kg,i) - fezero
            roi = dexp(-beta*defe)/zfr

            !~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Boltzmann weight (total)
            !~~~~~~~~~~~~~~~~~~~~~~~~~
            defet = fei(kg,i) - fezerot
            roit = dexp(-beta*defet)/zft

            if (verbose) then
               write(log,'(1x,"Precursor state:",i2,5x,"Boltzmann weight (at fixed R):",e20.9)') i,roi
               write(log,'(1x,t24,"Boltzmann weight (total):     ",e20.9)') roit
               write(log,'(1x,"minimum: ",3f12.3)') zpi(kg,i),zei(kg,i),fei(kg,i)
               write(log,'(1x,79("-"))')
               if (igr.ne.4) then
                  write(log,'(1x,t2,"j",t9,"deg",t20,"Er",t31,"Ea",t45,"Ea*",t60,"V^2",t75,"Wij")')
	       else
                  write(log,'(1x,t2,"j",t9,"deg",t20,"Er(tot)",t31,"Ea",t45,"Ea*",t60,"V^2",t75,"Wij-highT",t90,"Wij-lowT")')
	       endif
            endif

            wi = 0.d0
	    wi_low = 0.d0
            do j=1,nsucc

               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! Calculate partial rate I-->J
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               if (igr.ne.4) then
                  wij = prefac*v2(kg,i,j)*dexp(-beta*ea(kg,i,j))/dsqrt(er(kg,i,j))
	       else
	          er_total = er(kg,i,j) + er_r + er_alpha
	          ea_high_total = ea(kg,i,j) + alphamunu*deltar*(er_total+de(kg,i,j))/er_total/beta &
		                           & + alphamunu*alphamunu*deltar*deltar/er_total/beta/beta
                  wij = prefac_high*v2(kg,i,j)*dexp(-beta*ea_high_total)/dsqrt(er_total)
                  wij_low = prefac_low*v2(kg,i,j)*dexp(-beta*ea_low(kg,i,j))/dsqrt(er(kg,i,j))
	       endif

               wi = wi + roi*wij
	       if (igr.eq.4) wi_low = wi_low + roi*wij_low
               eaexp = dexp(-beta*ea(kg,i,j))

               wtotal = wtotal + roi*wij
               wr = wr + roit*wij
	       
	       if (igr.eq.4) then
                  wtotal_low = wtotal_low + roi*wij
                  wr_low = wr_low + roit*wij_low
	       endif

               if (verbose) write(log,'(i2,3f11.6,4e15.6)')&
               &j,de(kg,i,j),er(kg,i,j),ea(kg,i,j),eaexp,v2(kg,i,j),wij*1.d12,wij_low*1.d12

            enddo

            if (verbose) then
               write(log,'(1x,79("-"))')
	       if (igr.eq.4) then
                  write(log,'(1x,"Total High-T rate from state ",i2," (1/sec): ",e20.9)') i,wi*1.d12
                  write(log,'(1x,"Total Low-T  rate from state ",i2," (1/sec): ",e20.9)') i,wi_low*1.d12
	       else
                  write(log,'(1x,"Total rate from state ",i2," (1/sec): ",e20.9)') i,wi*1.d12
	       endif
               write(log,'(1x,79("="))')
            endif

         enddo

         if (igr.eq.1.or.igr.eq.2) then

            totalrate = totalrate + wtotal*gr*gint
            write(74,'(5g20.6)') rr,gr/bohr2a,wtotal*1.d12,wtotal*gr*gint*1.d12,totalrate*1.d12
            write(log,'(1x,"(igr=1/2) Rate (1/sec) at R=",f6.3,": ",e20.9)') rr,wtotal*1.d12

         elseif (igr.eq.3) then

            totalrate = totalrate + wr
            totalout = totalrate*1.d12
            wrout = wr*1.d12
            if (totalout.lt.1.d-20) totalout = 0.d0
            if (wrout.lt.1.d-20) wrout = 0.d0
            write(74,'(5g20.6)') rr,wrout,totalout
            write(log,'(1x,"(igr=3)   Rate (1/sec) at R=",f6.3,": ",e20.9)') rr,wr*1.d12

	 elseif (igr.eq.4) then

            totalrate = totalrate + wr
	    totalrate_low = totalrate_low + wr_low
            totalout = totalrate*1.d12
            totalout_low = totalrate_low*1.d12
            wrout = wr*1.d12
	    wrout_low = wr_low*1d12
            if (totalout.lt.1.d-20) totalout = 0.d0
            if (totalout_low.lt.1.d-20) totalout_low = 0.d0
            if (wrout.lt.1.d-20) wrout = 0.d0
            if (wrout_low.lt.1.d-20) wrout_low = 0.d0
            write(74,'(5g20.6)') rr,totalout,totalout_low
            write(log,'(1x,"(igr=4) High-T Rate (1/sec) at R=",f6.3,": ",e20.9)') rr,wr*1.d12
            write(log,'(1x,"(igr=4) Low-T  Rate (1/sec) at R=",f6.3,": ",e20.9)') rr,wr_low*1.d12

         endif

         write(log,'(1x,79("=")/)')

      enddo ! loop over gating coordinate grid

      close (74)
      close (75)
      close (76)
      close (77)

      write(log,'(1x,79("=")/)')
      if (igr.ne.4) then
         write(log,'(1x,"TOTAL NONADIABATIC RATE (1/sec): ",e20.9)') totalrate*1.d12
      else
         write(log,'(1x,"TOTAL NONADIABATIC RATE (High-T) (1/sec): ",e20.9)') totalrate*1.d12
         write(log,'(1x,"TOTAL NONADIABATIC RATE (Low-T)  (1/sec): ",e20.9)') totalrate_low*1.d12
      endif
      write(log,'(1x,79("=")/)')

      !======================================================
      ! Store the rates for the first and last temperature to
      ! estimate an effective (apparent) activation energy
      !======================================================

      if (trange.or.tlist) then

         if (igr.ne.4) then
            write(99,'(3f10.3,3g20.10)')&
            &temp,100.d0/temp,beta,totalrate*1.d12,dlog10(totalrate*1.d12),dlog(totalrate*1.d12)
	 else
            write(98,'(3f10.3,3g20.10)')&
            &temp,100.d0/temp,beta,totalrate*1.d12,dlog10(totalrate*1.d12),dlog(totalrate*1.d12)
            write(97,'(3f10.3,3g20.10)')&
            &temp,100.d0/temp,beta,totalrate_low*1.d12,dlog10(totalrate_low*1.d12),dlog(totalrate_low*1.d12)
	 endif

      endif

      if (trange) then

         if (it.eq.1) then
            betalast = beta
            ratelast = dlog(totalrate*1.d12)
         endif

         if (it.eq.n_t) then
            betafirst = beta
            ratefirst = dlog(totalrate*1.d12)
         endif

      endif

      !=====================================================
      ! Total rate analysis for each precursor state (IGR=3)
      !=====================================================

      if (igr.eq.3) then

         open(74,file=job(1:ljob)//"/rate_i_ofr"//ctemp//".dat")
         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("# Non-adiabatic rate for each precursor state")')
         write(74,'("#")')
         write(74,'("# Explicit averaging over gating coordinate")')
         write(74,'("# k_tot = Sum_i [ Int [ dR*P(i,R)*Sum_f [ k_if(R) ] ] ]")')
         write(74,'("# P(i,R) is a distribution function for reactant state i:")')
         write(74,'("#")')
         write(74,'("# 5 Columns:")')
         write(74,'("#")')
         write(74,'("# R   P(i,R)   Sum_f[k_if(R)]] P(i,R)*Sum_f[k_if(R)]] k_i(tot)")')
         write(74,'("#-----------------------------------------------------------------")')
         write(74,'("#")')

         !--------------------------------------
         ! calculate total rate
         !--------------------------------------
         totalrate = 0.d0

         do i=iprec,nprec

            write(74,'("#------ Initial reactant state: ",i2)') i

            wi = 0.d0
            do kg=1,npntsg

               rr = glist(kg)*bohr2a

               defet = fei(kg,i) - fezerot
               roit = dexp(-beta*defet)/zft

               wijsum = 0.d0

               do j=1,nsucc

                  wij = prefac*v2(kg,i,j)*dexp(-beta*ea(kg,i,j))/dsqrt(er(kg,i,j))
                  wijsum = wijsum + wij
                  wi = wi + roit*wij
                  totalrate = totalrate + roit*wij

               enddo

               write(74,'(4g20.6)') rr,roi,wijsum*1.d12,wi*1.d12

            enddo

            write(74,*)
            write(74,*)
            write(74,*)

         enddo

         close (74)

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

   !--------------------------
   ! deallocate working arrays
   !--------------------------
   deallocate (istel,istpr)
   if (method.eq.3.and.vrpout) deallocate (vrpkk)
   deallocate (fe,z_px,z_sx,psiel_px,psipr_px,psiel_sx,psipr_sx)
   deallocate (enel,envib)
   deallocate (zpi,zei,zpj,zej,zpx,zex)
   deallocate (min_succ)
   deallocate (fei,fej,fex,er,ea,de,v2)
   if (igr.eq.4) deallocate (ea_low)

   return


   !============================================================================
   contains

   !============================================================================
   real*8 function gofr(igr,gmass,gfreq,temp,gaver,rr)
   !============================================================================
   !  Calculates the distribution function for gating coordinate
   !
   !  IGR = 1 - classical Boltzmann distribution for Harmonic potential
   !        2 - exact quantum distribution for Harmonic potential
   !        3 - classical Boltzmann distribution for real potential
   !            (explicit treatment)
   !============================================================================
      integer, intent(in)  :: igr
      real*8,  intent(in)  :: gmass, gfreq, temp, gaver, rr

      real*8 :: gm, gf, tangent, pref, drr, expon

      gm = gmass
      gf = gfreq*cm2ev*ev2au

      if (igr.eq.2.or.igr.eq.1) then

         ! quantum distribution function for a harmonic oscillator

         tangent = dtanh(hbar*gfreq/(2.d0*kb*temp*cal2au))
         pref = pi*hbar/(gmass*gfreq*tangent)
         pref = 1.d0/dsqrt(pref)

         drr = (rr - gaver)*a2bohr
         expon = dexp(-gmass*gfreq*drr*drr*tangent/hbar)

         gofr = pref*expon

      elseif (igr.eq.3) then

         gofr = 1.d0

      else

         write(6,*) ' error in gofr: invalid igr key:',igr
         stop

      endif

      return

   end function gofr

   !============================================================================
   subroutine crossp3(m_,n_,kg_,pzp0_,pze0_,szp0_,sze0_,xacc_,maxit_,zpx_,zex_,fex_)
   !============================================================================
   !  Locates a crossing point between two free energy surfaces (M,N)
   !  along the straight line connecting two points (PZP0,PZE0) and
   !  (SZP0,SZE0) (gating distance is fixed).
   !
   !  The method of secants is used. Algorithm is based on the
   !  original program RTSEC from "Fortran recipes",
   !  (C) Copr. 1986-92 Numerical Recipes Software $$$.
   !============================================================================
      integer, intent(in)  :: m_, n_, kg_, maxit_
      real*8,  intent(in)  :: pzp0_, pze0_, szp0_, sze0_, xacc_
      real*8,  intent(out) :: zpx_, zex_, fex_

      integer :: j_
      real(8) :: deltazp, deltaze, x1, zp, ze, fe1, fe2, fl
      real(8) :: x2, f, rtsec, xl, swap, deltax

      deltazp = szp0_ - pzp0_
      deltaze = sze0_ - pze0_

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Evaluate the function (energy splitting)
      ! at the endpoints (PZP0,PZE0) and (SZP0,SZE0)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      x1 = 0.d0
      zp = pzp0_
      ze = pze0_
      call feszz3('DIAB2',1,kg_,zp,ze,nzdim,fe,nzdim,z_px,ndabf,ielst,enel,envib,psiel_px,psipr_px)
      fe1 = fe(m_)
      call feszz3('DIAB2',2,kg_,zp,ze,nzdim,fe,nzdim,z_sx,ndabf,ielst,enel,envib,psiel_sx,psipr_sx)
      fe2 = fe(n_)
      fl = fe2 - fe1

      x2 = 1.d0
      zp = szp0_
      ze = sze0_
      call feszz3('DIAB2',1,kg_,zp,ze,nzdim,fe,nzdim,z_px,ndabf,ielst,enel,envib,psiel_px,psipr_px)
      fe1 = fe(m_)
      call feszz3('DIAB2',2,kg_,zp,ze,nzdim,fe,nzdim,z_sx,ndabf,ielst,enel,envib,psiel_sx,psipr_sx)
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

      ! begin iterations

      do j_=1,maxit_

         deltax = (xl-rtsec)*f/(f-fl)

         xl = rtsec
         fl = f
         rtsec = rtsec + deltax
         zp = pzp0_ + rtsec*deltazp
         ze = pze0_ + rtsec*deltaze
         call feszz3('DIAB2',1,kg_,zp,ze,nzdim,fe,nzdim,z_px,ndabf,ielst,enel,envib,psiel_px,psipr_px)
         fe1 = fe(m_)
         call feszz3('DIAB2',2,kg_,zp,ze,nzdim,fe,nzdim,z_sx,ndabf,ielst,enel,envib,psiel_sx,psipr_sx)
         fe2 = fe(n_)
         f = fe2 - fe1

         if (dabs(deltax).lt.xacc_.or.f.eq.0.d0) then
            zpx_ = pzp0_ + rtsec*deltazp
            zex_ = pze0_ + rtsec*deltaze
            fex_ = (fe1 + fe2)/2.d0
            return
         endif

      enddo

      write(*,'(/1x,''RTSEC exceeded maximum iterations ('',i3.3,'')''/)') maxit_
      zpx_ = pzp0_ + rtsec*deltazp
      zex_ = pze0_ + rtsec*deltaze
      fex_ = (fe1 + fe2)/2.d0

      return

   end subroutine crossp3

   !============================================================================
   function coupling(i_,j_) result(vij_)
   !============================================================================
   ! calculates coupling between two wavefunctions
   !============================================================================
      integer, intent(in) :: i_, j_
      real*8 :: vij_

      vij_ = 0.d0

      select case(method)

         case(1)

            do k=1,2
               do mu=1,nprst
                  do l=1,2
                     do nu=1,nprst
                        vklmunu = 0.d0
                        do kk=1,npnts
                           psipr_ovl = psipr_px(k,mu,kk)*psipr_sx(l,nu,kk)
                           vklmunu = vklmunu + psipr_ovl*h0(k,l+2,kk,kg)
                        enddo
                        kmu = (k-1)*nprst+mu
                        lnu = (l-1)*nprst+nu
                        vij_ = vij_ + z_px(kmu,i_)*z_sx(lnu,j_)*vklmunu
                     enddo
                  enddo
               enddo
            enddo

         case(2)

            write(*,'(/1x,"*** method=2 is not implemented for rate3 calculations"/)')
            stop

         case(3)

            iel = istel(i_)
            imu = istpr(i_)
            jel = istel(j_)
            jnu = istpr(j_)

            do levb=1,2
               do mevb=1,2
                  do kk=1,npnts
                     psiel_ovl = psiel_px(iel,kk,levb)*psiel_sx(jel,kk,mevb)
                     psipr_ovl = psipr_px(iel,imu,kk)*psipr_sx(jel,jnu,kk)
                     vij_ = vij_ + psiel_ovl*psipr_ovl*h0(levb,mevb+2,kk,kg)
                  enddo
               enddo
            enddo

      end select

      return
      
   end function coupling

end subroutine rate3
