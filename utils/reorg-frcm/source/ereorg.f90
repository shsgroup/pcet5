subroutine ereorg
!======================================================================C
!
!     Calculates reorganization energy matrices (electronic, inertial,
!     and reduced) for four basis charge distributions given as
!     collections of partial charges on atoms.
!
!     Also calculates the free energy along the ET solvent coordinate
!     for single ET reaction in the framework of the standard
!     Marcus-like two-state model. Any two states (out of four
!     given in the input) can be specified as
!     reactant and product states in the ET reaction.
!
!     Output data are written to the specified external
!     file in the format
!
!     ZE   FE_R   FE_P   FE1   FE2
!
!     where ZE          is the solvent (Zusman) coordinate;
!           FE_R/FE_P   are the diabatic reactant/product
!                       free energies;
!           FE1/FE2     are the adiabatic ground/excited
!                       state free energies.
!
!     The nonadiabatic rate is also calculated using the
!     standard Marcus-Levich expression.
!
!     OPTIONS:
!     --------
!
!     REAC=<N> - the EVB state N (1/2/3/4) will be considered
!                as a reactant state in the ET reaction. If this
!                option is not specified then the state 1 (1a)
!                is assumed to be the reactant state.
!
!     PROD=<M> - the EVB state M (1/2/3/4) will be considered
!                as a product state in the ET reaction. If this
!                option is not specified then the state 3 (2a)
!                is assumed to be the product state.
!
!     ZE1=<ZE1> - the left limit on the grid along the solvent
!                 coordinate ZE (kcal/mol). The default is
!                 -3*Er, where Er is the reorganization energy.
!
!     ZE2=<ZE2> - the right limit on the grid along the solvent
!                 coordinate ZE (kcal/mol). The default is
!                 +3*Er, where Er is the reorganization energy.
!
!     NZE=<N> - number of points along the grid. The default is 100.
!
!     ALIN=<ERIN> - innersphere reorganization energy in kcal/mol
!                   (default ERIN=0)
!
!     PLOT=<filename> - name of the output file with the free energy
!                       profiles (default filename is "et.dat").
!
!     EBIAS=<BIAS> - gas phase bias in kcal/mol (prod-reac)
!
!     VEL=<VEL> - electronic coupling in kcal/mol
!
!     RATE=<T> - the nonadiabatic rate at T=<T>K will be calculated.
!                The default temperature is 298.15 K.
!
!======================================================================C

   use pardim
   use keys
   use strings
   use cst
   use control
   use geosol
   use frcm
   use elcm

   implicit none
   character(1024) :: options
   character( 40)  :: fname

   integer :: ikey, ireac, iprod, nze, lenf
   integer :: ixp, ir, ip, ize1, ize2, inze, ispa
   integer :: ilin, iplot, i, ize, irate

   real(8) :: ze1, ze2, alin, hreac, hprod, ebias, beta
   real(8) :: zestep, ze, se, ureac, uprod, uplus, uminu
   real(8) :: uground, uexcite, temp, prefac, dg0, er, ea, wps, gsreac, gsprod

   real(8), dimension(nelst,nelst) :: tk, tinfk, trk, trinfk
   real(8), dimension(2,2) :: hg2, erm

   write(6,'(/1x,''======================================================'',/,&
             &1x,''Reorganization energy matrices and Marcus ET reaction:'',/,&
             &1x,''======================================================'')')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' EREORG(')

   if (ikey.eq.0) then

      !-- All options are set to the default values

      ireac = 1
      iprod = 2
      ze1 = 0.d0
      ze2 = 0.d0
      nze = 100
      alin = 0.d0
      fname = job(1:ljob)//'/et.dat'
      lenf = ljob + 7

   else

      call getopt(keywrd,ikey+5,options)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Reactant and product EVB states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ir = index(options,' REAC=')
      if (ir.ne.0) then
         ireac = reada(options,ir+6)
      else
         ireac = 1
      endif

      ip = index(options,' PROD=')
      if (ip.ne.0) then
         iprod = reada(options,ip+6)
      else
         iprod = 3
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Gas-phase bias and coupling for ET pair of states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ir = index(options,' EBIAS=')
      if (ir.ne.0) then
         ebias = reada(options,ir+7)
      else
         ebias = 0.d0
      endif

      ip = index(options,' VEL=')
      if (ip.ne.0) then
         beta = reada(options,ip+5)
      else
         beta = 0.d0
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Left and right limits for the solvent coordinate
      ! and a number of points
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ize1 = index(options,' ZE1=')
      if (ize1.ne.0) then
         ze1 = reada(options,ize1+5)
      else
         ze1 = 0.d0
      endif

      ize2 = index(options,' ZE2=')
      if (ize2.ne.0) then
         ze2 = reada(options,ize2+5)
      else
         ze2 = 0.d0
      endif

      inze = index(options,' NZE=')
      if (inze.ne.0) then
         nze = reada(options,inze+5)
      else
         nze = 100
      endif

      ilin = index(options,' ALIN=')
      if (ilin.ne.0) then
         alin = reada(options,ilin+6)
      else
         alin = 0.d0
      endif

      write(6,'(1x,''Inner sphere reorganization energy (kcal/mol): '',g13.6)') alin

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Output file for free energy plots
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      iplot = index(options,' PLOT=')

      if (iplot.eq.0) then
         fname = job(1:ljob)//'/et.dat'
         lenf = ljob + 7
      else
         ispa = index(options(iplot+6:),' ')
         fname = options(iplot+6:iplot+ispa+4)
         lenf = ispa - 1
         call locase(fname,lenf)
         fname = job(1:ljob)//'/'//fname(1:lenf)
         lenf = lenf + ljob + 1
      endif

   endif

   write(6,'( 1x,''Reactant EVB state: '',I3)') ireac
   write(6,'( 1x,''Product  EVB state: '',I3)') iprod
   write(6,'(/1x,''Gas-phase bias (kcal/mol):      '',f8.3)') ebias
   write(6,'( 1x,''Electronic coupling (kcal/mol): '',f8.3/)') beta

   if (ze1.eq.0.d0.and.ze2.eq.0.d0) then
      write(6,'( 1x,''Left  limit for ZE: -3*Er'')')
      write(6,'( 1x,''Right limit for ZE:  3*Er'')')
   else
      write(6,'( 1x,''Left  limit for ZE: '',F10.3)') ze1
      write(6,'( 1x,''Right limit for ZE: '',F10.3)') ze2
   endif

   write(6,'(/1x,''Number of points:   '',I5)') nze
   write(6,'( 1x,''Free energy data are are in the external file <'',a,''>'')') fname(1:lenf)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Form the gas phase 2x2 Hamiltonian matrix for the ET
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hreac = 0.d0
   hprod = ebias
   hg2(1,1) = hreac
   hg2(2,2) = hprod
   hg2(1,2) = beta
   hg2(2,1) = beta
   call primat(hg2,2,2,2,5,6,0,'Gas phase ET Hamiltonian')


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate [4x4] reorganization energy matrices
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (isolv.eq.1) then

      !-- Ellipsoidal model
      call tmat(tk,tinfk,trk,trinfk)

   elseif (isolv.eq.2) then

      !-- FRCM model
      call solint(tk,tinfk,trk,trinfk)

   else

      write(*,'(/1x,'' ERROR in SETMAT: ISOLV='',I2)')
      stop

   endif

   erm(1,1) = tk(ireac,ireac)
   erm(2,2) = tk(iprod,iprod)
   erm(1,2) = tk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Full inertial reorganization energy matrix')

   erm(1,1) = tinfk(ireac,ireac)
   erm(2,2) = tinfk(iprod,iprod)
   erm(1,2) = tinfk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Full electronic reorganization energy matrix')

   erm(1,1) = trk(ireac,ireac)
   erm(2,2) = trk(iprod,iprod)
   erm(1,2) = trk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Full inertial reduced reorganization energy matrix')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Electronic solvation
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hreac = hreac - tinfk(ireac,ireac)/2.d0
   hprod = hprod - tinfk(iprod,iprod)/2.d0
   hg2(1,1) = hreac
   hg2(2,2) = hprod
   call primat(hg2,2,2,2,5,6,0,'Gas phase ET Hamiltonian with electronic solvation')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate ET reorganization energy
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   er = tk(ireac,ireac) + tk(iprod,iprod) - 2.d0*tk(ireac,iprod)
   er = er/2.d0 + alin
   write(6,'(/1x,"ET reorganization energy      : ",T35,F13.6," kcal/mol")') er
   write(6,'( 1x,"Stokes shift (linear response): ",T35,F13.6," cm^-1")') 2.d0*er*cal2ev*ev2cm

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate ET solvation energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   gsreac = -tinfk(ireac,ireac)/2.d0 - tk(ireac,ireac)/2.d0
   gsprod = -tinfk(iprod,iprod)/2.d0 - tk(iprod,iprod)/2.d0
   write(6,'(/1x,"Total solvation free energy of the reactant state: ",F13.6," kcal/mol")') gsreac
   write(6,'( 1x,"Total solvation free energy of the product  state: ",F13.6," kcal/mol")') gsprod

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Free energy curves
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Open output file for free energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   open(1,file=fname(1:lenf),status='new')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over the solvent coordinate points
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (ze1.eq.0.d0.and.ze2.eq.0.d0) then
      ze1 = -3.d0*er
      ze2 =  3.d0*er
   endif

   if (nze.gt.1) then
      zestep = (ze2-ze1)/(nze-1)
   else
      zestep = 0.d0
   endif

   do ize=1,nze

      ze = ze1 + (ize-1)*zestep

      !-- self energy

      se = ze + tk(iprod,ireac) - tk(ireac,ireac)
      ! The line below was edited by Helene Decornez
      ! se = se*se/4.d0/er
      se = se*se/4.d0/(er-alin)

      !~~~~~~~~~~~~~~~~~~~~~~~
      ! diabatic free energies
      !~~~~~~~~~~~~~~~~~~~~~~~
      ureac = se + hreac - tk(ireac,ireac)/2.d0
      uprod = se + hprod - tk(ireac,ireac)/2.d0 + ze

      !~~~~~~~~~~~~~~~~~~~~~~~~
      ! adiabatic free energies
      !~~~~~~~~~~~~~~~~~~~~~~~~
      uplus = (ureac + uprod)/2.d0
      uminu = (ureac - uprod)/2.d0
      uground = uplus - dsqrt(uminu*uminu + beta*beta)
      uexcite = uplus + dsqrt(uminu*uminu + beta*beta)

      !~~~~~~~~~~~~~~~~~~~~~~~~
      ! Print out to FNAME file
      !~~~~~~~~~~~~~~~~~~~~~~~~
      write(1,'(5g13.6)') ze,ureac,uprod,uground,uexcite

   enddo

   close(1)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Nonadiabatic rate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (index(options,' RATE').ne.0) then

      irate = index(options,' RATE=')
      if (irate.ne.0) then
         temp = reada(options,irate+6)
      else
         temp = 298.15d0
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Preexponential factor (constant)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      prefac = dsqrt(pi/(kb*temp))/hbarps

      !~~~~~~~~~~~~~~~~~~~~~
      ! Reaction free energy
      !~~~~~~~~~~~~~~~~~~~~~
      dg0 = hprod - hreac - 0.5d0*(tk(iprod,iprod)-tk(ireac,ireac))
      write(6,'(/1x,''ET reaction free energy: '',T35,F13.6,'' kcal/mol'')') dg0

      !~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Marcus activation energy
      !~~~~~~~~~~~~~~~~~~~~~~~~~
      ea = (er + dg0)*(er + dg0)/er/4.d0
      write(6,'(1x,''ET activation free energy: '',T35,F13.6,'' kcal/mol'')') ea

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Nonadiabatic rate (1/picoseconds)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wps = prefac*beta*beta*dexp(-ea/(kb*temp))/dsqrt(er)
      write(6,'(/1x,''ET nonadiabatic rate  (1/sec) at '',F10.3,'' K: '',E20.9)') temp,wps*1.d12

   endif

   return

end subroutine ereorg
