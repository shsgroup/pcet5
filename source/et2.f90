subroutine et2
!======================================================================C
!     Calculates the free energy along the ET solvent coordinate
!     for single ET reaction in the framework of the standard
!     Marcus-like two-state model. Any two states (out of four
!     given in input for the PCET reaction) can be specified as
!     reactant and product states in the ET reaction. Proton
!     is kept fixed at the specified location.
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
!     XP=<XP> - X coordinate of the proton (it is assumed that the
!               PT interface is linear and aligned along the X axis).
!               Actually this is just a bond length between the
!               proton donor and hydrogen. If XP option is not
!               specified then the location of the proton is taken
!               from the geometry input file.
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
!     PLOT=<filename> - name of the output file with the free energy
!                       profiles (default filename is "et.dat").
!
!     RATE=<T> - the nonadiabatic rate at T=<T>K will be calculated.
!                The default temperature is 298.15 K.
!
!----------------------------------------------------------------------C
!
!  $Author: souda $
!  $Date: 2010-12-15 21:24:55 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:35  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!======================================================================C
   use pardim
   use keys
   use strings
   use cst
   use control
   use geogas
   use geosol
   use frcm
   use elcm
   use potential
   use msevb_water

   implicit none
   character(1024) :: options
   character( 40)  :: fname

   integer :: nhgas, nhsol, ikey, ireac, iprod, nze, lenf
   integer :: ixp, ir, ip, ize1, ize2, inze, ispa
   integer :: ilin, iplot, i, ize, irate
   
   real(8) :: xxp, yyp, zzp, ze1, ze2, alin, hreac, hprod, beta
   real(8) :: zestep, ze, se, ureac, uprod, uplus, uminu
   real(8) :: uground, uexcite, temp, prefac, dg0, er, ea, wps

   real(8), dimension(4,4) :: h0k, dh0k, d2h0k, dgh0k, dg2h0k
   real(8), dimension(4,4) :: tk, tinfk, trk, trinfk
   real(8), dimension(2,2) :: hg2, erm

   write(6,'(/1x,''==================='',/,&
             &1x,''Marcus ET reaction:'',/,&
             &1x,''==================='')')

   nhgas = iptgas(2)
   nhsol = iptsol(2)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' ET2(')

   if (ikey.eq.0) then

      ! All options are set to the default values

      xxp = xyzgas(1,nhgas)
      yyp = xyzgas(2,nhgas)
      zzp = xyzgas(3,nhgas)
      ireac = 1
      iprod = 3
      ze1 = 0.d0
      ze2 = 0.d0
      nze = 100
      alin = 0.d0
      fname = job(1:ljob)//'/et.dat'
      lenf = ljob + 7

   else

      call getopt(keywrd,ikey+5,options)

      ! Proton X-coordinate

      ixp = index(options,' XP=')
      if (ixp.ne.0) then
         xxp = reada(options,ixp+4)
         yyp = 0.d0
         zzp = 0.d0
      else
         xxp = xyzgas(1,iptgas(2))
         yyp = xyzgas(2,iptgas(2))
         zzp = xyzgas(3,iptgas(2))
      endif

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

      write(6,'(1x,''Inner sphere reorganization energy (kcal/mol): '',g12.6)') alin

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

   write(6,'(/1x,''Proton coordinates: '',3F8.3)') xxp,yyp,zzp
   write(6,'( 1x,''Reactant EVB state: '',I3)') ireac
   write(6,'( 1x,''Product  EVB state: '',I3)') iprod

   if (ze1.eq.0.d0.and.ze2.eq.0.d0) then
      write(6,'( 1x,''Left  limit for ZE: -3*Er'')')
      write(6,'( 1x,''Right limit for ZE:  3*Er'')')
   else
      write(6,'( 1x,''Left  limit for ZE: '',F8.3)') ze1
      write(6,'( 1x,''Right limit for ZE: '',F8.3)') ze2
   endif

   write(6,'( 1x,''Number of points:   '',I5)') nze
   write(6,'( 1x,''Free energy data are are in the external file <'',a,''>'')') fname(1:lenf)


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the full 4x4 gas-phase Hamiltonian matrix
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   xyzgas(1,nhgas) = xxp
   xyzgas(2,nhgas) = yyp
   xyzgas(3,nhgas) = zzp

   if (igas.eq.1) then

      ! MM5 potential for Nocera-like systems
      call h0mat_mm5(h0k)

   elseif (igas.eq.2) then

      ! MS-EVB (Voth, Schmitt) potential for water clusters
      x0evb(nhgas) = xyzgas(1,nhgas)
      call msevb(h0k)

   elseif (igas.eq.3) then

      ! original LEPS 2D potential
      call h0mat_leps0(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

   elseif (igas.eq.4) then

      ! modified LEPS 2D potential
      call h0mat_leps2(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

   elseif (igas.eq.5) then

      ! Hybrid LEPS/Harmonic 2D potential
      call h0mat_hyb(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

   elseif (igas.eq.6) then

      ! MM5 potential for general O-H...O systems
      call h0mat_mm5gen(h0k)

   elseif (igas.eq.7) then

      ! Harmonic potential for general PT systems
      call h0mat_harm(h0k)

   else

      write(*,'(/1x,'' error in et2: igas='',i2)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Form the gas phase 2x2 Hamiltonian matrix for the ET
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hreac = h0k(ireac,ireac)
   hprod = h0k(iprod,iprod)
   beta  = h0k(ireac,iprod)
   hg2(1,1) = hreac
   hg2(2,2) = hprod
   hg2(1,2) = beta
   hg2(2,1) = beta
   call primat(hg2,2,2,2,5,6,0,'Gas phase ET Hamiltonian')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the full reorganization energy matrices
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   xyzsol(1,nhsol) = xxp
   xyzsol(2,nhsol) = yyp
   xyzsol(3,nhsol) = zzp

   if (isolv.eq.1) then

      ! Ellipsoidal model
      call tmat(tk,tinfk,trk,trinfk)

   elseif (isolv.eq.2) then

      ! FRCM model
      do i=1,3
         coord(i,nhsol) = xyzsol(i,nhsol)
      enddo
      xx(nhsol) = coord(1,nhsol)
      yy(nhsol) = coord(2,nhsol)
      zz(nhsol) = coord(3,nhsol)
      call solint(tk,tinfk,trk,trinfk)

   else

      write(*,'(/1x,'' ERROR in SETMAT: ISOLV='',I2)')
      stop

   endif

   erm(1,1) = tk(ireac,ireac)
   erm(2,2) = tk(iprod,iprod)
   erm(1,2) = tk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Inertial reorganization energy matrix')

   erm(1,1) = tinfk(ireac,ireac)
   erm(2,2) = tinfk(iprod,iprod)
   erm(1,2) = tinfk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Electronic reorganization energy matrix')

   erm(1,1) = trk(ireac,ireac)
   erm(2,2) = trk(iprod,iprod)
   erm(1,2) = trk(ireac,iprod)
   erm(2,1) = erm(1,2)
   call primat(erm,2,2,2,5,6,0,'Inertial reduced reorganization energy matrix')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Electronic solvation
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   hreac = hreac - tinfk(ireac,ireac)/2.d0
   hprod = hprod - tinfk(iprod,iprod)/2.d0
   hg2(1,1) = hreac
   hg2(2,2) = hprod
   call primat(hg2,2,2,2,5,6,0,'Gas phase ET Hamiltonian with electronic solvation')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the ET reorganization energy
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   er = tk(ireac,ireac) + tk(iprod,iprod) - 2.d0*tk(ireac,iprod)
   er = er/2.d0 + alin
   write(6,'(/1x,"ET reorganization energy      : ",T35,F12.6," kcal/mol")') er
   write(6,'( 1x,"Stokes shift (linear response): ",T35,F12.6," cm^-1")') 2.d0*er*cal2ev*ev2cm

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

      ! self energy

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
      write(1,'(5g12.6)') ze,ureac,uprod,uground,uexcite

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
      write(6,'(1x,''ET reaction free energy: '',T35,F12.6,'' kcal/mol'')') dg0

      !~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Marcus activation energy
      !~~~~~~~~~~~~~~~~~~~~~~~~~
      ea = (er + dg0)**2/er/4.d0
      write(6,'(1x,''ET activation free energy: '',T35,F12.6,'' kcal/mol'')') ea

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Nonadiabatic rate (1/picoseconds)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      wps = prefac*beta*beta*dexp(-ea/(kb*temp))/dsqrt(er)
      write(6,'(/1x,''ET nonadiabatic rate  (1/sec) at '',F8.3,'' K: '',E20.9)') temp,wps*1.d12

   endif

   return

end subroutine et2
