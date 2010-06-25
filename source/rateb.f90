subroutine rateb
!=======================================================================
!     Calculates the nonadiabatic rate of the PCET reaction in solution
!
!     Reference:
!     [*] Alexander Soudackov and Sharon Hammes-Schiffer,
!         Gating effects on the nonadiabatic rate of
!         proton-coupled electron transfer reactions in solution,
!         [in preparation]
!
!     The dynamical gating effects are taken into account according
!     to the treatment similar to the Borgis-Lee-Hynes treatment
!     of proton transfer reactions, see
!
!     [**] ...
!
!     OPTIONS:
!
!     T=<value>[/<value>/<step>] - absolute temperature
!                                  [range (init/fin/npoints)]
!
!     T=<value1>-<value2>-...-<valueN> - absolute temperature list
!                                        (up to 20 values)
!
!     PREC=<I> - precursor state I is the only initial state.
!                No Boltzmann averaging over the initial states
!                is performed.
!
!     NPREC=<N> - N precursor states are taken into account.
!                 The total rate is averaged over initial
!                 precursor states with Boltzmann weights.
!
!     NSUCC=<N> - N successor states are taken into account
!
!     ZPREC=<ZP0>/<ZE0>/<R0> - guess for the reactant minimum on the
!                              precursor free energy surface (kcal/mol)
!
!     ZSUCC=<ZP0>/<ZE0>/<R0> - guess for the product minimum on the
!                              successor free energy surface (kcal/mol)
!
!     ALPHA=<value> - coefficient (in A^-1) in the exponential dependence
!                     of the coupling on the gating distance.
!                     If omitted then it will be approximately determined
!                     from the difference of couplings calculated
!                     at the equilibrium gating distance and at the
!                     gating distance corresponding to the crossing point
!                     along a straight line between reactant and product
!                     minima. The following relations
!                     can be used:
!
!                     (a) alpha = -(1/dR) Log[1+dV/V0], where
!
!                         dR - displacement from equilibrium gating distance (1 A)
!                         R = R0 + dR
!                         dV = V(R) - V(R0)
!                         V0 = V(R0)
!
!     ERANALY - the reorganization free energy is evaluated analytically
!               using the expression (55) from Ref.[*].
!               Otherwise (if not specified) it is evaluated
!               numerically from two intersecting three-dimensional
!               free energy surfaces (parabaloids).
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
!     FLUXOUT=   - the time dependent flux in the general
!                  rate expression for each pair of states
!                  will be written to the external files.
!
!     NINTEGRATE=<N> - method for numerical integration of the flux
!               =1 - DQSF method from SSPLIB
!               =2 - DQH64 method from SSPLIB
!               =3 - DQAGS method from QADPACK
!               =4 - simple Monte-Carlo method
!
!     VERBOSE - print detailed rate info for each gating grid point.
!               If not specified only the total rate is printed out.
!               Not also that not specifying VERBOSE will cancel the
!               effects of WEIGHTS option.
!
!-----------------------------------------------------------------------
!
!     souda
!     2010/06/25 20:02:37
!     4.1
!     Exp
!     rateb.f90,v 4.1 2010/06/25 20:02:37 souda Exp
!     rateb.f90,v
!     Revision 4.1  2010/06/25 20:02:37  souda
!     Release 4.1
!
!     Revision 1.13  2008/04/11 00:07:20  souda
!     length of string OPTIONS increased to 1024
!     to accomodate more options (not critical)
!
!     Revision 1.12  2007/03/12 23:08:04  souda
!     Modifications related to using LBFGS minimization method.
!
!     Revision 1.11  2004/06/21 16:39:48  souda
!     some beautification...
!
!     Revision 1.10  2004/06/14 17:36:16  souda
!     testing loginfo
!
!     Revision 1.9  2004/06/14 16:07:19  souda
!     testing loginfo
!
!     Revision 1.8  2004/06/11 23:20:58  souda
!     fixed phase issue in alpha calculation in rateb
!
!     Revision 1.7  2004/06/09 21:44:40  souda
!     output to Arrhenius changed
!
!     Revision 1.6  2004/06/07 17:33:05  souda
!     changed output formats
!
!     Revision 1.5  2004/05/17 15:03:24  souda
!     minor change in rateb...
!
!     Revision 1.4  2004/05/15 03:32:45  souda
!     Added Borgis-Hynes rate routine
!
!     Revision 1.3  2004/05/13 18:32:01  souda
!     finalized rateb. routine
!
!     Revision 1.2  2004/01/14 23:54:04  souda
!     Start of modifying rate3
!
!     Revision 1.1.1.1  2004/01/05 00:31:12  souda
!     Initial PCET-4.0 Release
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
   use geogas
   use fesmin_3d
   use feszz_3d, only: feszz3
   use rate_flux

   implicit none

   character(1024) :: options
   character(  40) :: fname
   character(   5) :: mode
   character(   3) :: ctemp

   logical :: weights, outlog, eanumer, eranaly, alpha_fixed
   logical :: filemin, rdmin, verbose, trange, tlist

   logical, allocatable, dimension(:)    :: min_succ

   integer :: ikey, ilog, lenf, ispa, log_channel, nintegrate
   integer :: iiprec, inprec, insucc, kgprec, kgsucc
   integer :: itemp, itemp_start, itemp_end
   integer :: islash, idash, n_t, islash1, islash2, idashn
   integer :: izprec, izsucc, izprec1, izprec1_end, izsucc1, izsucc1_end
   integer :: islim, iacc, ixacc, imaxit, imaxitx, maxitx, izps, izes, ilin
   integer :: i, j, it, ielst, np, ndabf, nzdim, iprec, nprec, nsucc, ierr
   integer :: iprin, imcorr, ifactr, ipgtol
  
   integer, allocatable, dimension(:) :: istel, istpr
   integer, allocatable, dimension(:) :: kgi, kgj

   real(kind=8) :: mr, zpprec0, zeprec0, rrprec0, zpsucc0, zesucc0, rrsucc0
   real(kind=8) :: xacc, zps, zes, alin, denom, elxx, elyy, elxy
   real(kind=8) :: zpprec, zeprec, rrprec, feprec, zpsucc, zesucc, rrsucc, fesucc, fej_tmp
   real(kind=8) :: dzpij, dzpijp, dzeij, dzeijp, drrij, drrijp
   real(kind=8) :: erij_za, fij_ra, erij_ra, erij_a, erij_n, erij_zn, erij_rn, omegar_av
   real(kind=8) :: zpij, zeij, feij, eaijn, eaija, vij0, vij1, dvij, rrstep
   real(kind=8) :: t_start, t_end, range, dtemp, temp, beta, zft, fezerot
   real(kind=8) :: totalrate_exact, totalrate_high, totalrate_low
   real(kind=8) :: wi_exact, wi_high, wi_low
   real(kind=8) :: wij_exact, wij_high, wij_low
   real(kind=8) :: betalast, betafirst, ratelast, ratefirst
   real(kind=8) :: defet, roit, eaeff, aeff, deltarate, deltabeta, alpha_inp

   real(kind=8), dimension(20)  :: ti
   real(kind=8), dimension(3)   :: gradprec, gradsucc
   real(kind=8), dimension(3,3) :: hessprec, hesssucc

   real(kind=8), allocatable, dimension(:)     :: fe
   real(kind=8), allocatable, dimension(:,:)   :: z_px, z_sx
   real(kind=8), allocatable, dimension(:,:)   :: enel, envib
   real(kind=8), allocatable, dimension(:,:,:) :: psiel_px, psipr_px
   real(kind=8), allocatable, dimension(:,:,:) :: psiel_sx, psipr_sx

   real(kind=8), allocatable, dimension(:) :: zpi, zei, rri, zpj, zej, rrj
   real(kind=8), allocatable, dimension(:) :: elxr, elyr, omegar, frr

   real(kind=8), allocatable, dimension(:,:) :: el_x, el_y, el_r
   real(kind=8), allocatable, dimension(:,:) :: el0_xy, el0_r, el0_zr
   real(kind=8), allocatable, dimension(:,:) :: er, er_z, er_r, erq

   real(kind=8), allocatable, dimension(:,:) :: zpx, zex, rrx, alpha
   real(kind=8), allocatable, dimension(:)   :: fei, fej
   real(kind=8), allocatable, dimension(:,:) :: fex, de, ea, v2

   ! numerical integration variables
   integer :: ndim_dqsf
   real(kind=8) :: error_dqsf
   integer :: limit_qags, last_qags, neval_qags, ier_qags
   real(kind=8) :: epsabs_qags, epsrel_qags, abserr_qags
   real(kind=8), allocatable, dimension(:) :: alist_qags, blist_qags, rlist_qags, elist_qags
   integer, allocatable, dimension(:) :: iord_qags
   integer :: nsteps_mc, nsteps_accepted
   real(kind=8) :: maxstep_mc, tau_min, tau_max

   !============================================================================ 

   if (method.eq.2) then
      write(*,'(/1x,"*** method=2 is not implemented for rate calculations ***"/)')
      write(*,'( 1x,"    try again (RATEB)..."/)')
      stop
   endif

   ! default
   rdmin   = .false.
   filemin = .false.
   trange  = .false.
   alpha_fixed = .false.

   !====================================================
   ! Extract options
   !====================================================

   ikey = index(keywrd,' RATEB(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** in rateb: you must specify options for the RATEB keyword ***"/)')
      stop
   else
      call getopt(keywrd,ikey+7,options)
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

      log_channel = 1
      open(log_channel,file=fname(1:lenf),status='new')

   else

      log_channel = 6

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! BANNER
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   write(log_channel,'(/1x,70("="))')
   write(log_channel,'( 1x,"Nonadiabatic Rate Calculation")')
   write(log_channel,'( 1x,"(dynamical treatment of the gating motion)")')
   write(log_channel,'( 1x,70("=")/)')

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
      write(log_channel,'(/1x,"initial precursor state:",i3)') iprec
      write(log_channel,'( 1x,"no Boltzmann averaging over initial states")')

   elseif (inprec.ne.0) then

      iprec = 1
      nprec = reada(options,inprec+7)
      write(log_channel,'(/1x,"initial precursor state:",i3)') iprec
      write(log_channel,'( 1x,"Boltzmann averaging over lowest ",i2," initial states")') nprec

   else

      iprec = 1
      nprec = 1
      write(log_channel,'(/1x,"initial precursor state:",i3)') iprec
      write(log_channel,'( 1x,"no Boltzmann averaging over initial states")')

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
   write(log_channel,'(/1x,"Number of successor states:",i3)') nsucc

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Allocate some arrays
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   allocate (zpi(iprec:nprec),zei(iprec:nprec),rri(iprec:nprec))  ! precursors minima
   allocate (zpj(1:nsucc),zej(1:nsucc),rrj(1:nsucc))              ! successors minima
   allocate (kgi(iprec:nprec), kgj(iprec:nprec))
   allocate (zpx(iprec:nprec,1:nsucc),&
           & zex(iprec:nprec,1:nsucc),&
           & rrx(iprec:nprec,1:nsucc))                            ! crossing points
   allocate (alpha(iprec:nprec,1:nsucc))                          ! coupling exp. parameters

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! File containing minima for calculation of the rates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(*,'(/1x,"Minima will be determined in calculation"/)')
   rdmin = .false.

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
      write(log_channel,'(/1x,"Temperature range (K):",i3," steps from",f8.3," to",f8.3)')&
      n_t, t_start, t_end
   elseif (tlist) then
      write(log_channel,'(/1x,"Temperature list (K):",/,(5f8.3))') (ti(i),i=1,n_t)
   else
      write(log_channel,'(/1x,"Temperature (K):",f8.3)') t_start
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
      izprec1 = izprec + 7
      izprec1_end = izprec1 + index(options(izprec1:),space)-1
      islash = index(options(izprec1:izprec1_end),'/')
      islash1 = izprec1 + islash - 1
      islash2 = islash1 + index(options(islash1+1:izprec1_end),'/')
      zpprec0 = reada(options(izprec1:islash1-1),1)
      zeprec0 = reada(options(islash1+1:islash2-1),1)
      rrprec0 = reada(options,islash2+1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Product minimum (ZPSUCC0,ZESUCC0)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      izsucc1 = izsucc + 7
      izsucc1_end = izsucc1 + index(options(izsucc1:),space)-1
      islash = index(options(izsucc1:izsucc1_end),'/')
      islash1 = izsucc1 + islash - 1
      islash2 = islash1 + index(options(islash1+1:izsucc1_end),'/')
      zpsucc0 = reada(options(izsucc1:islash1-1),1)
      zesucc0 = reada(options(islash1+1:islash2-1),1)
      rrsucc0 = reada(options,islash2+1)

      write(log_channel,'(/1x,"Initial guesses for minima (kcal/mol):")')
      write(log_channel,'( 1x,"Reactant=(",f8.3,",",f8.3,",",f8.3,")",&
                 &" and Product=(",f8.3,",",f8.3,",",f8.3,")")')&
                 & zpprec0,zeprec0,rrprec0, &
                 & zpsucc0,zesucc0,rrsucc0

   else

      write(*,'(/1x,"*** INPUT ERROR (in RATEB): ",/,&
      &"You MUST specify both ZPREC and ZSUCC options",&
      &" for RATEB keyword ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for the total reorganization energy (ERANALY)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eranaly = index(options,' ERANALY').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation mode for activation energy (EANUMER)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   eanumer = index(options,' EANUMER').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Flag for output of the total EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   weights = index(options,' WEIGHTS').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Flag for output of the total EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   weights = index(options,' WEIGHTS').ne.0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Coupling parameter alpha
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ikey = index(options,' ALPHA=')
   if (ikey.ne.0) then
      alpha_inp = reada(options,ikey+7)
      alpha_fixed = .true.
   else
      alpha_fixed = .false.
   endif

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

      write(log_channel,'(/1x,"Maximum length of the Newton-Raphson step: ",g12.6)') slim

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Accuracy of the minimization
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      iacc = index(options,' ACC=')
      if (iacc.ne.0) then
         acc = reada(options,iacc+5)
      else
         acc = 1.d-6
      endif
      write(log_channel,'(1x,"Accuracy of the minimization: ",g12.6)') acc

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
   write(log_channel,'(1x,"Accuracy of the crossing point location: ",g12.6)') xacc

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations in searching minima
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxit = index(options,' MAXIT=')
   if (imaxit.ne.0) then
      maxit = reada(options,imaxit+7)
   else
      maxit = 100
   endif
   write(log_channel,'(1x,"Maximum number of Newton-Raphson steps: ",i5)') maxit

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Maximum number of iterations in searching a crossing point
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   imaxitx = index(options,' MAXITX=')
   if (imaxitx.ne.0) then
      maxitx = reada(options,imaxitx+8)
   else
      maxitx = 100
   endif
   write(log_channel,'(1x,"Maximum number of steps in locating a crossing point: ",i5)') maxitx

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Scaling factors for the solvent coordinates (ZPS,ZES)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   izps = index(options,' ZPS=')
   if (izps.ne.0.and.iminim.eq.1) then
      zps = reada(options,izps+5)
      write(log_channel,'(1x,"Scaling factor for ZP coordinate: ",g12.6)') zps
   else
      zps = 1.d0
   endif

   izes = index(options,' ZES=')
   if (izes.ne.0.and.iminim.eq.1) then
      zes = reada(options,izes+5)
      write(log_channel,'(1x,"Scaling factor for ZE coordinate: ",g12.6)') zes
   else
      zes = 1.d0
   endif

   ilin = index(options,' ALIN=')
   if (ilin.ne.0) then
      alin = reada(options,ilin+6)
      write(log_channel,'(1x,"Inner sphere reorganization energy added (kcal/mol): ",g12.6)') alin
   else
      alin = 0.d0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Method for numerical integration of the flux
   ! =1 - DQSF method from SSPLIB
   ! =2 - DQH64 method from SSPLIB
   ! =3 - DQAGS method from QUADPACK
   ! =4 - simple Monte-Carlo integration
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ikey = index(options,' NINTEGRATE=')
   if (ikey.ne.0) then
      nintegrate = reada(options,ikey+12)
   else
      nintegrate = 3
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Parameters for numerical integration
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! DQSF

   ikey = index(options,' NDIM_DQSF=')
   if (ikey.ne.0) then
      ndim_dqsf = reada(options,ikey+11)
   else
      ndim_dqsf = 5000
   endif
   
   ! DQH64
   ! none
   
   ! DQAGS

   ikey = index(options,' EPSABS_QAGS=')
   if (ikey.ne.0) then
      epsabs_qags = reada(options,ikey+13)
   else
      epsabs_qags = 1.d-10
   endif

   ikey = index(options,' EPSREL_QAGS=')
   if (ikey.ne.0) then
      epsrel_qags = reada(options,ikey+13)
   else
      epsrel_qags = 1.d-10
   endif

   ikey = index(options,' LIMIT_QAGS=')
   if (ikey.ne.0) then
      limit_qags = reada(options,ikey+12)
   else
      limit_qags  = 500
   endif

   ! Monte-Carlo

   ikey = index(options,' NSTEPS_MC=')
   if (ikey.ne.0) then
      nsteps_mc = reada(options,ikey+11)
   else
      nsteps_mc = 1000000
   endif

   ikey = index(options,' MAXSTEP_MC=')
   if (ikey.ne.0) then
      maxstep_mc = reada(options,ikey+12)
   else
      maxstep_mc = 0.1d0
   endif

   select case(nintegrate)
      case(1)
         write(log_channel,'(1x,"DQSF method will be used for numerical integration of the flux")')
         write(log_channel,'(1x,"- number of integration points: ",i10)') ndim_dqsf
      case(2)
         write(log_channel,'(1x,"DQH64 method will be used for numerical integration of the flux")')
      case(3)
         write(log_channel,'(1x,"DQAGS method will be used for numerical integration of the flux")')
         write(log_channel,'(1x,"- absolute accuracy requested: ",e15.6)') epsabs_qags
         write(log_channel,'(1x,"- relative accuracy requested: ",e15.6)') epsrel_qags
         write(log_channel,'(1x,"- max. number of subdivisions: ",i6)')    limit_qags
      case(4)
         write(log_channel,'(1x,"Monte-Carlo method will be used for numerical integration of the flux")')
         write(log_channel,'(1x,"- max. number of steps: ",i6)')    nsteps_mc
         write(log_channel,'(1x,"- max. step           : ",e15.6)') maxstep_mc
   end select


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
   denom = two*(erpt*eret-erx*erx)
   elxx = eret/denom
   elyy = erpt/denom
   elxy = -erx/denom

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! reduced mass for proton donor-acceptor mode
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   mr = dm*am/(dm+am)

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

   allocate (fe(nzdim),z_px(nzdim,nzdim),z_sx(nzdim,nzdim))
   allocate (psiel_px(ielst,npnts,ielst),psipr_px(ielst,nprst,npnts))
   allocate (psiel_sx(ielst,npnts,ielst),psipr_sx(ielst,nprst,npnts))
   allocate (enel(ielst,npnts),envib(ielst,nprst))

   allocate (fei(iprec:nprec))           ! precursors equilibrium free energies
   allocate (fej(1:nsucc))               ! successors equilibrium free energies
   allocate (fex(iprec:nprec,1:nsucc))   ! free energies at the crossing points (at Rmin)

   allocate (elxr(iprec:nprec),elyr(iprec:nprec),omegar(iprec:nprec),frr(iprec:nprec))

   allocate ( de  (iprec:nprec,1:nsucc))    ! reaction free energies (delta G's)
   allocate ( er  (iprec:nprec,1:nsucc))    ! total reorganization energies
   allocate ( er_z(iprec:nprec,1:nsucc))    ! solvent reorganization energies
   allocate ( er_r(iprec:nprec,1:nsucc))    ! R-mode reorganization energies
   allocate ( erq (iprec:nprec,1:nsucc))    ! coupling reorganization energies
   allocate ( ea  (iprec:nprec,1:nsucc))    ! activation energies
   allocate ( v2  (iprec:nprec,1:nsucc))    ! squared couplings

   allocate ( el_x(iprec:nprec,1:nsucc),&
           &  el_y(iprec:nprec,1:nsucc),&
           &  el_r(iprec:nprec,1:nsucc))

   allocate ( el0_xy(iprec:nprec,1:nsucc),&
           &  el0_r (iprec:nprec,1:nsucc),&
           &  el0_zr(iprec:nprec,1:nsucc))

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over the initial states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   min_succ = .false.

   do i=iprec,nprec

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Locate a minimum on the reactant surface I
      ! or calculate free energy if minima read in
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/"  Search of the minumum for reactant state ",i3/)') i
      if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"f(1)",t91,"f(2)")')
      call fes3min(mode,1,i,zpprec0,zps,zeprec0,zes,rrprec0,&
                   &ierr,zpprec,zeprec,kgprec,rrprec,&
                   &feprec,gradprec,hessprec)

      if (ierr.eq.0) then
         zpprec0 = zpprec
         zeprec0 = zeprec
         rrprec0 = rrprec
      endif

      zpi(i) = zpprec
      zei(i) = zeprec
      rri(i) = rrprec
      kgi(i) = kgprec
      fei(i) = feprec

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Quantities related to the calculation
      ! of reorganization energies
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      elxr(i) = hessprec(1,3)
      elyr(i) = hessprec(2,3)
      omegar(i) = sqrt(hessprec(3,3)*cal2au/mr)*au2ev*ev2cm/a2bohr  ! R-mode frequency (1/cm)
      frr(i) = hessprec(3,3) - two*erpt*elxr(i)*elxr(i) &
                             - two*eret*elyr(i)*elyr(i) &
                             - four*erx*elxr(i)*elyr(i)    ! R-mode force constant (kcal/mol/A^2)
      if (verbose) then
         write(log_channel,'(/1x,"R-mode frequency (1/cm)             : ",f12.6)') omegar(i)
         write(log_channel,'( 1x,"R-mode force constant (kcal/mol/A^2): ",f12.6)') frr(i)
         write(log_channel,'( 1x,"Lambda_xR (1/A^2): ",f12.6)') elxr(i)
         write(log_channel,'( 1x,"Lambda_yR (1/A^2): ",f12.6)') elyr(i)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Calculate and print out the EVB weights
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (weights.and.verbose) then
         write(log_channel,'(/1x,70("="))')
         write(log_channel,'(1x,"EVB weights at the minimum of the ",i2,a3," reactant state")') i,th(i)
         write(log_channel,'(1x,"zp(precursor)=",f9.3,5x,"ze(precursor)=",f9.3/)') zpprec,zeprec
         call weight3(log_channel,kgprec,zpprec,zeprec,mode,1,i)
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

            write(6,'(/" Search of minumum for product state ",i3/)') j
            if(iminim.eq.1) write(6,'(t2,"iter",t18,"zp",t32,"ze",t46,"e",t61,"grad",t76,"f(1)",t91,"f(2)")')
            call fes3min(mode,2,j,zpsucc0,zps,zesucc0,zes,rrsucc0,&
                         &ierr,zpsucc,zesucc,kgsucc,rrsucc,&
                         &fesucc,gradsucc,hesssucc)
            if (ierr.eq.0) then
               zpsucc0 = zpsucc
               zesucc0 = zesucc
               rrsucc0 = rrsucc
            endif
            zpj(j) = zpsucc
            zej(j) = zesucc
            rrj(j) = rrsucc
            kgj(j) = kgsucc
            fej(j) = fesucc

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate and print out the EVB weights
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (weights.and.verbose) then
               write(log_channel,'(/1x,70("="))')
               write(log_channel,'(1x,"EVB weights at the minimum of the ",i2,a3," product state")') j,th(j)
               write(log_channel,'(1x,"zp(successor)=",f9.3,5x,"ze(successor)=",f9.3/)') zpsucc,zesucc
               call weight3(log_channel,kgsucc,zpsucc,zesucc,mode,2,j)
            endif

            min_succ(j) = .true.

         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Quantities related to the calculation
         ! of reorganization energies
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         dzpij  = zpj(j) - zpi(i)
         dzpijp = zpj(j) + zpi(i)
         dzeij  = zej(j) - zei(i)
         dzeijp = zej(j) + zei(i)
         drrij  = rrj(j) - rri(i)
         drrijp = rrj(j) + rri(i)

         el_x(i,j) = -elxx*dzpij - elxy*dzeij - elxr(i)*drrij
         el_y(i,j) = -elxy*dzpij - elyy*dzeij - elyr(i)*drrij
         el_r(i,j) = -hessprec(3,3)*drrij - elxr(i)*dzpij - elyr(i)*dzeij

         el0_xy(i,j) = half*elxx*dzpij*dzpijp + &
                     & half*elyy*dzeij*dzeijp + &
                     &      elxy*(zpj(j)*zej(j)-zpi(i)*zei(i))

         el0_r(i,j)  = half*hessprec(3,3)*drrij*drrijp

         el0_zr(i,j) = elxr(i)*(rrj(j)*zpj(j)-rri(i)*zpi(i)) + &
                     & elyr(i)*(rrj(j)*zej(j)-rri(i)*zei(i))


         if (verbose) then

            write(log_channel,'(/1x,"Quantities related to reorganization energies")')
            write(log_channel,'( 1x,"---------------------------------------------")')

            write(log_channel,'(/1x,"Delta R (A): ",f12.6)') drrij

            write(log_channel,'(/1x,"Lambda_x (          ): ",f12.6)') el_x(i,j)
            write(log_channel,'( 1x,"Lambda_y (          ): ",f12.6)') el_y(i,j)
            write(log_channel,'( 1x,"Lambda_R (kcal/mol/A): ",f12.6)') el_y(i,j)

            write(log_channel,'(/1x,"Lambda0_xy (kcal/mol): ",f12.6)') el0_xy(i,j)
            write(log_channel,'( 1x,"Lambda0_R  (kcal/mol): ",f12.6)') el0_r(i,j)
            write(log_channel,'( 1x,"Lambda0_zR (kcal/mol): ",f12.6)') el0_zr(i,j)

	 endif


         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate reaction free energy for the pair (I,J)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         de(i,j) = fej(j) - fei(i)

         if (verbose)&
         write(log_channel,'(/1x,"The reaction free energy (kcal/mol) for the pair ",2i2,":",f20.6)')&
         i, j, de(i,j)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate two-dimensional reorganization energy
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         ! Analytical solvent reorganization energy

         erij_za = el_x(i,j)*el_x(i,j)*erpt + el_y(i,j)*el_y(i,j)*eret + &
                &  two*el_x(i,j)*el_y(i,j)*erx

         ! Analytical R-mode reorganization energy

         fij_ra = frr(i) -  two*erpt*elxr(i)*elxr(i) &
                       & -  two*eret*elyr(i)*elyr(i) &
                       & - four*erx* elxr(i)*elyr(i)
         erij_ra = half*fij_ra*drrij*drrij


         ! Analytical total reorganization energy

         erij_a = el0_xy(i,j) + el0_r(i,j) + el0_zr(i,j) + &
                & el_x(i,j)*zpi(i) + el_y(i,j)*zei(i) + el_r(i,j)*rri(i)

         !erij_a = erij_za + erij_ra
         erij_a = erij_a + alin

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Numerical total reorganization energy
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 ! calculate UI(zpj,zej,rrj), and
	 ! lambda_t_num = UI(zpj,zej,rrj) - UI(zpi,zei,rri)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         fe = 0.d0
         z_px = 0.d0
         enel = 0.d0
         envib = 0.d0
         psiel_px = 0.d0
         psipr_px = 0.d0
         call feszz3(mode,1,kgsucc,zpsucc,zesucc,nzdim,fe,nzdim,z_px,&
                   & ndabf,ielst,enel,envib,psiel_px,psipr_px)
         erij_n = fe(i) - fei(i) + alin

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Numerical solvent reorganization energy
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	 ! calculate UI(zpj,zej,rri), and
	 ! lambda_z_num = UI(zpj,zej,rri) - UI(zpi,zei,rri)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         fe = 0.d0
         z_px = 0.d0
         enel = 0.d0
         envib = 0.d0
         psiel_px = 0.d0
         psipr_px = 0.d0
         call feszz3(mode,1,kgprec,zpsucc,zesucc,nzdim,fe,nzdim,z_px,&
                   & ndabf,ielst,enel,envib,psiel_px,psipr_px)
         erij_zn = fe(i) - fei(i)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Numerical R-mode reorganization energy
	 ! calculate UI(zpi,zei,rrj), and
	 ! lambda_r_num = UI(zpi,zei,rrj) - UI(zpi,zei,rri)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         fe = 0.d0
         z_px = 0.d0
         enel = 0.d0
         envib = 0.d0
         psiel_px = 0.d0
         psipr_px = 0.d0
         call feszz3(mode,1,kgsucc,zpprec,zeprec,nzdim,fe,nzdim,z_px,&
                   & ndabf,ielst,enel,envib,psiel_px,psipr_px)
         erij_rn = fe(i) - fei(i)

         if (verbose)&
            write(log_channel,'(/1x,"The reorganization energies (kcal/mol) for the pair ",2i2,":",/,&
            &1x,"analytical solvent:",f20.6,/,&
            &1x,"analytical R-mode :",f20.6,/,&
            &1x,"analytical total  :",f20.6,//,&
            &1x,"numerical solvent :",f20.6,/,&
            &1x,"numerical R-mode  :",f20.6,/,&
            &1x,"numerical total   :",f20.6)') i, j, erij_za, erij_ra, erij_a, &
                                                   & erij_zn, erij_rn, erij_n

         if (eranaly) then
            er_z(i,j) = erij_za
            er_r(i,j) = erij_ra
            er  (i,j) = erij_a
            if (verbose) write(log_channel,'(/1x,"analytical values will be used")')
         else
            er_z(i,j) = erij_zn
            er_r(i,j) = erij_rn
            er  (i,j) = erij_n
            if (verbose) write(log_channel,'(/1x,"numerical values will be used")')
         endif

         if (er(i,j).lt.0.d0) then
            write(*,'(/1x,"The reorg. energy for the pair ",2i2," is negative:",f20.6)')&
            &i,j,er(i,j)
            write(*,'("this is very suspicious, so the program terminates...")')
            stop
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Find the crossing point between electronically
         ! diabatic one-dimensional (slices along the straiht
         ! line connecting two minima) free-energy curves
         ! for the states I and J at R = Rmin(i)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call crossp3(i,j,kgprec,zpprec,zeprec,zpsucc,zesucc,xacc,maxitx,zpij,zeij,feij)

         zpx(i,j) = zpij
         zex(i,j) = zeij
         fex(i,j) = feij

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate and print out the EVB weights
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (weights.and.verbose) then
            write(log_channel,'(/1x,70("="))')
            write(log_channel,'(1x,"EVB weights at the crossing point between ",i2,a3,&
            &" reactant and ",i2,a3," product states")') i,th(i),j,th(j)
            write(log_channel,'(1x,"zp(x-ing point)=",f9.3,5x,"ze(xing point)=",f9.3/)') zpij, zeij
            call weight3(log_channel,kgprec,zpij,zeij,mode,1,i)
            call weight3(log_channel,kgprec,zpij,zeij,mode,2,j)
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate activation free energy for the pair (I,J)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         eaijn = fex(i,j) - fei(i)
         eaija = (er(i,j)+de(i,j))*(er(i,j)+de(i,j))/4.d0/er(i,j)

         if (verbose)&
            &write(log_channel,'(/1x,"The activation energy for the pair of states ",i2,"-",i2,":",/,&
            &" numerical (from the crossing point): ",f12.6,/,&
            &" analytical (from Marcus expression): ",f12.6/)') i, j, eaijn, eaija

         if (eanumer) then
            ea(i,j) = eaijn   ! numerical estimate
         else
            ea(i,j) = eaija   ! Marcus relation
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate coupling matrix element at the crossing point
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         ! First, at R = Rmin(i)

         call feszz3(mode,1,kgprec,zpij,zeij,nzdim,&
                    &fe,nzdim,z_px,ndabf,ielst,enel,envib,&
                    &psiel_px,psipr_px)

         call feszz3(mode,2,kgprec,zpij,zeij,nzdim,&
                    &fe,nzdim,z_sx,ndabf,ielst,enel,envib,&
                    &psiel_sx,psipr_sx)

         vij0 = coupling(i,j,kgprec)
         v2(i,j) = vij0*vij0


         ! Second, at R = Rmin(i) + dR

         if (alpha_fixed) then

            alpha(i,j) = alpha_inp

         else

            call feszz3(mode,1,kgprec+1,zpij,zeij,nzdim,&
                       &fe,nzdim,z_px,ndabf,ielst,enel,envib,&
                       &psiel_px,psipr_px)

            call feszz3(mode,2,kgprec+1,zpij,zeij,nzdim,&
                       &fe,nzdim,z_sx,ndabf,ielst,enel,envib,&
                       &psiel_sx,psipr_sx)

            vij1 = coupling(i,j,kgprec+1)

            ! Estimate "alpha" parameter

            rrstep = (glist(2) - glist(1))*bohr2a
            dvij = (vij1*vij1 - vij0*vij0)/rrstep
            alpha(i,j) = -dvij/(two*vij0*vij0)

         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Calculate quantum coupling reorganization energy
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

          erq(i,j) = au2cal*alpha(i,j)*alpha(i,j)/(a2bohr*a2bohr*two*mr)

      enddo         ! loop over i (precursors)

   enddo          ! loop over j (successors)

   !====================================================
   ! Initial states characteristics
   !====================================================
   if (verbose) then
      write(log_channel,'(/1x,"INITIAL STATES CHARACTERISTICS (kcal/mol)")')
      write(log_channel,'(1x,70("="))')
      write(log_channel,'(1x,t4,"state",t17,"zp0",t32,"ze0",t47,"rr0",t58,"fe(zp0,ze0,rr0)")')
      write(log_channel,'(1x,70("-"))')
      do i=iprec,nprec
         write(log_channel,'(1x,i5,5x,4(f10.3,5x))') i,zpi(i),zei(i),rri(i),fei(i)
      enddo
      write(log_channel,'(1x,70("-"))')
   endif

   !====================================================
   ! Crossing points
   !====================================================
   if (verbose) then
      write(log_channel,'(/1x,"CROSSING POINTS CHARACTERISTICS AT R=Rmin (kcal/mol)")')
      write(log_channel,'(1x,70("="))')
      write(log_channel,'(1x,t3,"pair",t14,"zpx",t29,"zex",t39,"fex(zp0,ze0)",t59,"v^2")')
      write(log_channel,'(1x,70("-"))')
      do i=iprec,nprec
         do j=1,nsucc
            write(log_channel,'(1x,i2,"-",i2,3(f10.3,5x),g12.6)')&
            &i,j,zpx(i,j),zex(i,j),fex(i,j),v2(i,j)
         enddo
      enddo
      write(log_channel,'(1x,70("-"))')
   endif

   !====================================================
   ! Alpha's for pairs of states
   !====================================================
   if (verbose) then
      write(log_channel,'(/1x,"Quantum coupling reorganization energies")')
      write(log_channel,'(1x,70("="))')
      write(log_channel,'(1x,t3,"pair",t14,"alpha",t29,"lambda(alpha)")')
      write(log_channel,'(1x,70("-"))')
      do i=iprec,nprec
         do j=1,nsucc
            write(log_channel,'(1x,i2,"-",i2,3(f10.3,5x),g12.6)')&
            &i,j,alpha(i,j),erq(i,j)
         enddo
      enddo
      write(log_channel,'(1x,70("-"))')
   endif

   !====================================================
   ! Final states characteristics
   !====================================================
   if (verbose) then
      write(log_channel,'(/1x,"FINAL STATES CHARACTERISTICS (kcal/mol)")')
      write(log_channel,'(1x,70("="))')
      write(log_channel,'(1x,t4,"state",t17,"zp0",t32,"ze0",t47,"rr0",t58,"fe(zp0,ze0,rr0)")')
      write(log_channel,'(1x,70("-"))')
      do j=1,nsucc
         write(log_channel,'(1x,i5,5x,4(f10.3,5x))') j,zpj(j),zej(j),rrj(j),fej(j)
      enddo
      write(log_channel,'(1x,70("-"))')
   endif

   !==================================================================
   ! Start loop over temperatures
   !==================================================================

   if (trange.or.tlist) then
      open(99,file=job(1:ljob)//"/arrhenius_b.dat")
      write(99,'("#",205("-"))')
      write(99,'("# Arrhenius plots: dynamical rates in 1/sec")')
      write(99,'("#")')
      write(99,'("# 12 Columns:")')
      write(99,'("#",205("-"))')
      write(99,'("#",t5,"1",t16,"2",t26,"3",t37,"4",t57,"5",t77,"6",t95,"7",t115,"8",t135,"9",&
                   & t156,"10",t176,"11",t195,"12")')
      write(99,'("#",205("-"))')
      write(99,'("#",t5,"T(K)",t16,"100/T",t26,"hO/kT",t37,"k_exact",t57,"k_high",t77,"k_low",&
                   & t95,"log10(k_exact)",t115,"log10(k_high)",t135,"log10(k_low)",&
		   & t156,"ln(k_exact)",t176,"ln(k_high)",t195,"ln(k_low)")')
      write(99,'("#",205("-"))')
   endif

   !==================================================================
   ! Average R-mode frequency (1/cm)
   !==================================================================

   omegar_av = zero
   do i=iprec,nprec
      omegar_av = omegar_av + omegar(i)
   enddo
   omegar_av = omegar_av/(nprec-iprec+1)


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

      write(log_channel,'(/1x,70("*"))')
      write(log_channel,'( 1x,"RATE ANALYSIS AT ",f8.3," K")') temp
      write(log_channel,'(/1x,70("*"))')

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! ZFT - total partition function
      !       summed over all the initial states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      zft = 0.d0

      !--------------------------------------
      ! find lowest free energy and
      ! calculate total partition function
      !--------------------------------------

      fezerot = minval(fei)

      do i=iprec,nprec
         zft = zft + dexp(-beta*(fei(i)-fezerot))
      enddo

      !---------------
      ! Total rate
      !---------------
      totalrate_exact = zero
      totalrate_high  = zero
      totalrate_low   = zero

      if (verbose) then
         write(log_channel,'(/1x,79("="))')
         write(log_channel,'(1x,"RATE CONTRIBUTIONS")')
         write(log_channel,'(1x,79("="))')
         write(log_channel,'(1x,"column headers:",/,&
         &1X,"J     - successor state",/,&
         &1X,"DeG   - reaction free energy (kcal/mol)",/,&
         &1X,"Er    - total reorganization energy (kcal/mol)",/,&
         &1X,"Er_Z  - solvent reorganization energy (kcal/mol)",/,&
         &1X,"Er_R  - R-mode reorganization energy (kcal/mol)",/,&
         &1X,"alpha - coupling parameter (1/A)",/,&
         &1X,"Er_q  - coupling reorganization energy (kcal/mol)",/,&
         &1X,"V^2   - squared nonadiabatic coupling (kcal/mol)**2",/,&
         &1X,"Wij   - partial rate (1/sec)")')
         write(log_channel,'(1x,79("="))')
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Averaging over initial states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=iprec,nprec

         !~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Boltzmann weight (total)
         !~~~~~~~~~~~~~~~~~~~~~~~~~
         defet = fei(i) - fezerot
         roit = dexp(-beta*defet)/zft

         if (verbose) then
            write(log_channel,'(1x,"Precursor state: ",i2,5x,"Boltzmann weight (total): ",e20.9)') i,roit
            write(log_channel,'(1x,"minimum: ",3f12.3)') zpi(i),zei(i),rri(i),fei(i)
            write(log_channel,'(1x,79("-"))')
            write(log_channel,'(1x,t2,"j",t9,"DeG",t20,"Er",t31,"Er_Z",t41,"Er_R",t51,"alpha",&
	                       &t63,"Erq",t74,"V^2",t90,"Wij(exact)",t105,"Wij(high)",t119,"Wij(low)")')
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Sum over final states
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         wi_exact = zero
         wi_high  = zero
         wi_low   = zero

         do j=1,nsucc

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Calculate partial rate I-->J
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

             call set_flux_pars&
                &(omegar(i),mr/dalton,er_r(i,j),alpha(i,j),er_z(i,j),de(i,j),v2(i,j),temp)

             ! calculate exact rate (integrating flux)
             
             select case(nintegrate)

                case(1)
                    call rate_exact(ndim_dqsf, error_dqsf, wij_exact)

                case(2)
                    call rate_exact(wij_exact)

                case(3)
                    allocate (alist_qags(limit_qags),&
                            & blist_qags(limit_qags),&
                            & rlist_qags(limit_qags),&
                            & elist_qags(limit_qags))
                    allocate (iord_qags(limit_qags))
                    call rate_exact(epsabs_qags,epsrel_qags,limit_qags,&
                                  & abserr_qags,neval_qags,ier_qags,alist_qags,&
                                  & blist_qags,rlist_qags,elist_qags,iord_qags,&
                                  & last_qags,wij_exact)

                    deallocate (alist_qags,blist_qags,rlist_qags,elist_qags)
                    deallocate (iord_qags)

                case(4)
                    call rate_exact(nsteps_mc,maxstep_mc,nsteps_accepted,tau_min,tau_max,wij_exact)

             end select

             ! calculate asymptotic rates
             call rate_high(wij_high)
             call rate_low (wij_low)

             ! rates from state i
             wi_exact = wi_exact + roit*wij_exact
             wi_high  = wi_high  + roit*wij_high
             wi_low   = wi_low   + roit*wij_low

            if (verbose) write(log_channel,'(i2,6f11.6,4e15.6)')&
              & j,de(i,j),er(i,j),er_z(i,j),er_r(i,j),alpha(i,j),erq(i,j),v2(i,j),&
              & wij_exact,wij_high,wij_low
              
            ! total rates
            totalrate_exact = totalrate_exact + roit*wij_exact
            totalrate_high  = totalrate_high  + roit*wij_high
            totalrate_low   = totalrate_low   + roit*wij_low

         enddo

         if (verbose) then
            write(log_channel,'(1x,79("-"))')
            write(log_channel,'(1x,"Total rates from state ",i2," (1/sec)")') i
            write(log_channel,'(1x,"exact       : ",e20.9)') wi_exact
            write(log_channel,'(1x,"high-T limit: ",e20.9)') wi_high
            write(log_channel,'(1x,"low-T  limit: ",e20.9)') wi_low
            write(log_channel,'(1x,79("="))')
         endif

      enddo

      write(log_channel,'(1x,79("=")/)')

      close (74)
      close (75)
      close (76)
      close (77)

      write(log_channel,'(1x,79("=")/)')
      write(log_channel,'(1x,"TOTAL NONADIABATIC RATES (1/sec)")')
      write(log_channel,'(1x,"exact       : ",e20.9)') totalrate_exact
      write(log_channel,'(1x,"high-T limit: ",e20.9)') totalrate_high
      write(log_channel,'(1x,"low-T  limit: ",e20.9)') totalrate_low
      write(log_channel,'(1x,79("=")/)')

      !======================================================
      ! Store the rates for the first and last temperature to
      ! estimate an effective (apparent) activation energy
      !======================================================

      if (trange) then

         write(99,'(3f10.3,10g20.10)') temp,100.d0/temp,omegar_av*cm2ev*ev2cal*beta,&
             & totalrate_exact,totalrate_high,totalrate_low,&
             & log10(abs(totalrate_exact)),log10(abs(totalrate_high)),log10(abs(totalrate_low)),&
             & log(abs(totalrate_exact)),log(abs(totalrate_high)),log(abs(totalrate_low))

         if (it.eq.1) then
            betalast = beta
            ratelast = log(abs(totalrate_exact))
         endif

         if (it.eq.n_t) then
            betafirst = beta
            ratefirst = log(abs(totalrate_exact))
         endif

      endif

   enddo ! loop over temperature range

   if (trange) then

      close(99)

      !=====================================================
      ! Estimate an effective (apparent) activation energy
      ! and a prefactor
      !=====================================================
      eaeff = -(ratelast - ratefirst)/(betalast - betafirst)
      deltarate = dabs(ratelast-ratefirst)
      deltabeta = dabs(betalast-betafirst)
      aeff = dexp(ratefirst + betafirst*deltarate/deltabeta)
      write(log_channel,'(1x,79("="))')
      write(log_channel,'(1x,"Apparent activation energy (kcal/mol): ",f20.6)') eaeff
      write(log_channel,'(1x,"Apparent prefactor (1/sec): ",g20.6)') aeff
      write(log_channel,'(1x,79("=")/)')

   endif

   !--------------------------
   ! deallocate working arrays
   !--------------------------
   deallocate (istel,istpr)
   deallocate (fe,z_px,z_sx,psiel_px,psipr_px,psiel_sx,psipr_sx)
   deallocate (enel,envib)
   deallocate (zpi,zei,rri,zpj,zej,rrj,zpx,zex,rrx)
   deallocate (alpha)
   deallocate (kgi, kgj)
   deallocate (min_succ)
   deallocate (fei,fej,fex,elxr,elyr,omegar,frr)
   deallocate (de,er,er_z,er_r,erq,ea,v2)
   deallocate (el_x,el_y,el_r,el0_xy,el0_r,el0_zr)

   return

   !============================================================================
   contains

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
   function coupling(i_,j_,kg_) result(vij_)
   !============================================================================
   ! calculates coupling between two wavefunctions
   !============================================================================
      integer, intent(in) :: i_, j_, kg_
      real*8 :: vij_

      ! local variables      
      integer :: k, mu, l, nu, kk, kmu, lnu, iel, imu, jel, jnu, levb, mevb
      real(8) :: vklmunu, psipr_ovl, psiel_ovl

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
                           vklmunu = vklmunu + psipr_ovl*h0(k,l+2,kk,kg_)
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
                     vij_ = vij_ + psiel_ovl*psipr_ovl*h0(levb,mevb+2,kk,kg_)
                  enddo
               enddo
            enddo

      end select

      return
      
   end function coupling


end subroutine rateb

