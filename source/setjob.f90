subroutine setjob

!====================================================================C
!   Analyses specified keywords and sets options
!   and parameters for the current job
!---------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-03-13 22:06:35 $
!  $Revision: 5.8 $
!  $Log: not supported by cvs2svn $
!  Revision 5.7  2011/03/01 23:52:15  souda
!  insignificant change in output
!
!  Revision 5.6  2011/02/25 19:11:25  souda
!  Now using a separate set of dielectric constant for solvent dynamics.
!
!  Revision 5.5  2011/01/04 19:59:14  souda
!  added: now the custom delta and kappa FRCM parameters can be
!         specified even if the solvent name is explicitly given.
!
!  Revision 5.4  2010/12/15 21:24:56  souda
!  various fixes (non-critical)
!
!  Revision 5.3  2010/11/10 21:14:21  souda
!  Last addition/changes. IMPORTANT: Solvation keyword is back to SOLV( for compatibility reasons
!
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!====================================================================C

   use pardim
   use cst
   use keys
   use strings
   use control
   use quantum
   use parsol
   use geogas
   use geosol
   use geometry
   use potential
   use frcm
   use msevb_water

   implicit none

   character(1024) :: options
   character(  40) :: fname
   logical :: ok
   logical :: quantum_proton

   integer :: ikey, ingrid, ingrids, npntsi, inprst, ilims, ialim, islash, itset
   integer :: irhmin, npntspow2, iarlim, igquant, ingast, icoulomb
   integer :: idist, ideltag, ivept, ipars, ispa, lenf
   integer :: igeo, ixyz, imajor, iminor, ikappa, idelta, igeom, itread
   integer :: natsol_, natgas_, nhpt, k, ieps0, ieps8, npntssol2
   
   real(8) :: x, ba, ba2

   !==================================================
   ! Set control parameters
   !==================================================

   !==================================================
   ! Method of calculation of vibronic states:
   !==================================================
   !
   ! 1 - using the DIABATIC approach: diagonalization
   !     of the total Hamiltonian in the basis of the
   !     products of DIABATIC electronic states (EVB)
   !     and vibrational states (in the DIABATIC
   !     EVB potentials).
   !
   ! 2 - using the ADIABATIC approach: diagonalization
   !     of the total Hamiltonian in the basis of the
   !     products of ADIABATIC electronic states
   !     and vibrational states (in the ADIABATIC
   !     potentials). Note that this approach involves
   !     the calculation of non-adiabatic "d" and "g"
   !     derivative terms.
   !
   ! 3 - using the DOUBLE ADIABATIC approach: the
   !     non-adiabatic "d" and "g" derivative terms
   !     are assumed to be zero.
   !==================================================

   ikey = index(keywrd,' METHOD=')

   if (ikey.ne.0) then

      method = reada(keywrd,ikey+8)
      if (method.lt.1.or.method.gt.3) then
         write(*,'(/1x,"*** (in setjob): method value is invalid ***"/)')
         stop
      endif
      write(6,'(/1x,"Method ",i1," will be used for calculation of the vibronic states")') method

   else

      !-- default value: 1
      method = 1
      write(6,'(/1x,"Default method=1 will be used for calculation of the vibronic states (if needed)")')
      !write(*,'(/1x,"*** (in setjob): you must specify the METHOD keyword ***"/)')
      !stop

   endif

   !====================================================================================
   ! Method of minimization on the 2D and 3D
   ! free energy surfaces
   !====================================================================================
   !
   ! 1 - standard Newton-Rafson method
   !
   ! 2 - Limited memory BFGS method (LBFGS)
   !
   !     References
   !
   !   R. H. Byrd, P. Lu and J. Nocedal.
   !      A Limited Memory Algorithm for Bound Constrained Optimization, (1995),
   !      SIAM Journal on Scientific and Statistical Computing, 16, 5, pp. 1190-1208.
   !   C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B,
   !      FORTRAN routines for large scale bound constrained optimization (1997),
   !      ACM Transactions on Mathematical Software, Vol 23, Num. 4, pp. 550 - 560.
   !
   !====================================================================================

   ikey = index(keywrd,' MINIM=')

   if (ikey.ne.0) then

      iminim = reada(keywrd,ikey+7)
      if (method.lt.1.or.iminim.gt.2) then
         write(*,'(/1x,"*** (in setjob): minimization method is invalid ***"/)')
         stop
      endif

   else

      !-- default value: Newton-Rafson
      iminim = 1

   endif

   if (iminim.eq.1) then
      write(6,'(/1x,"Newton-Rafson optimization routines will be used for locating the minimna on 2D and 3D free energy surfaces")')
   elseif (iminim.eq.2) then
      write(6,'(/1x,"LBFGS optimization routines will be used for locating the minimna on 2D and 3D free energy surfaces")')
   endif

   !==================================================
   ! Total charge of the solute
   !==================================================

   ikey = index(keywrd,' CHARGE=')
   if (ikey.ne.0) then
      charge = reada(keywrd,ikey+8)
   else
      charge = 0.d0
   endif
   write(6,'(/1x,"Total charge of the solute: ",f13.6)') charge

   !====================================================
   ! Parameters of the grid for 1D-Schroedinger equation
   ! for proton/deuterium
   !====================================================
   ! PM     - mass of the quantum particle
   ! DM     - mass of the donor
   ! AM     - mass of the acceptor
   ! ALIM   - left  integration limit for the proton (A)
   ! BLIM   - right integration limit for the proton (A)
   ! RHMIN  - minimum distance between the proton and
   !          donor/acceptor
   ! STEP   - integration step for the proton (a.u.)
   ! RLIST  - gridpoints along the quantum coordinate (a.u.)
   !====================================================

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! default parameters
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   alim     = -1.d0
   blim     =  1.d0
   rhmin    = 0.5d0
   npnts    = 128
   npntssol = npnts
   nprst    = 40
   pm       = hmass

   ikey = index(keywrd,' QUANTUM(')

   quantum_proton = .false.

   if (ikey.ne.0) then

      quantum_proton = .true.

      call getopt(keywrd,ikey+9,options)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Mass of the quantum particle (PM) and
      ! masses of donor (DM) and acceptor (AM) atoms
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (index(options,' PROTON').ne.0) then
         pm = hmass
      else if (index(options,' DEUTERIUM').ne.0) then
         pm = dmass
      else if (index(options,' MASS=').ne.0) then
         pm = reada(options,index(options,' MASS=')+6)
      else
         pm = hmass
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of grid points for the proton coordinate
      ! (must be a powers of two)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ingrid = index(options,' NGRID=')

      if (ingrid.ne.0) then
         npnts = reada(options,ingrid+7)
      else
         npnts = 128
      endif

      npntsi = npnts

      if (.not.power2(npntsi)) then

         write(*,'(/1x,"*** The number of grid points (",I3,") is not a power of two.")') npntsi

         npntspow2 = 1
         do while (npntsi.gt.npntspow2.and.npntspow2.lt.maxpnt)
            npntspow2 = npntspow2 * 2
         enddo
         npntsi = 2 * npntspow2
         npnts = npntsi
         write(*,'(" NGRID is set to ",i3/)') npntsi

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of grid points for solvation calculations.
      ! (must be a power of two and smaller that NPNTS)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ingrids = index(options,' NGRIDS=')

      if (ingrids.ne.0) then
         npntssol = reada(options,ingrids+8)
      else
         npntssol = npnts
      endif

      write(*,'(/1x,"Solvation calculations on the proton coordinate grid",/)')

      if (npntssol.lt.0) then

         write(*,'(/1x,"Solvation calculations will be performed only at the proton position",/,&
         &" specified in the geometry input file",/)')

      else if (npntssol.eq.0) then

         write(*,'(/1x,"Solvation calculations will be performed",/,&
         &" only at one configuration of the solute",/,&
         &" with a hydrogen atom being placed at the center",/,&
         &" of mass of the PT interface",/)')

      else if (.not.power2(npntssol).or.npntssol.gt.npnts) then

         write(*,'(/1x,"*** The number of grid points for solvation calculations (",I3,")",/,&
         &" is either not a power of two or greater then the total number of grid points",/,&
         &" NPNTSSOL is set equal to NPNTS (",i3,")",/)') npntssol,npnts
         npntssol = npnts

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of proton vibrational states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      inprst = index(options,' NPRST=')

      if (inprst.ne.0) then
         nprst = reada(options,inprst+7)
      endif

      if (nprst.gt.nprstmax) then

         WRITE(6,'(/" Number of protonic vibronic functions specified NPRST =",I5,/,&
         &" It is larger than maximum number allowed NPRSTMAX =",I5,/,&
         &" **  Program stopped  ** ",//,&
         &" Recommendation: Either reduce NPRST or recompile the",&
         &" program with larger NRPRSTMAX")') nprst,nprstmax
         stop

      elseif (nprst.lt.1) then

         write(6,'(/" Wrong number of protonic vibronic functions specified NPRST =",I5,/,&
         &" NPRST must be larger than 0 and smaller than NPRSTMAX =",I5,/,&
         &" **  Program stopped  ** ")') nprst,nprstmax
         stop

      ENDIF

      if (nprst.gt.npnts) then
         WRITE(6,'(/" Number of protonic vibronic functions specified NPRST =",I5,/,&
         &" It is larger than the number of grid points =",I5,/,&
         &" **  Program stopped  ** ")') nprst,npnts
         stop
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Grid (integration) limits for the proton
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ilims = index(options,' LIMITS=')
      if (ilims.ne.0) then
         ialim = ilims+8
         islash = index(options(ialim:),'/')
         alim = reada(options(ialim:ialim+islash-2),1)
         blim = reada(options,ialim+islash)
      endif

      write(6,'(/1x,"Quantum particle parameters:")')
      write(6,'(10x,"Mass of the particle (Daltons): ",f13.6)') pm/dalton
      write(6,'(10x,"Integration limits (Angstroms): ",2f13.6)') alim,blim
      write(6,'(10x,"Number of integration points:   ",i12)') npnts
      write(6,'(10x,"Number of proton basis states:  ",i12)') nprst

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Integration steps (STEP,STEPR) and
      ! arrays of grid points (RLIST)
      ! (now it's done in module_quantum:alloc_pquant)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !step = a2bohr*(blim-alim)/(npnts-1)
      !x = alim*a2bohr
      !do k=1,npnts+1
      !   rlist(k) = x
      !   x = x + step
      !enddo

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Minimum distance (Angstroms) between the proton
      ! and donor/acceptor
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      irhmin = index(options,' RHMIN=')

      if (irhmin.ne.0) then
         rhmin = reada(options,irhmin+7)
      else
         rhmin = 0.5d0
      endif

   endif

   !====================================================
   ! Parameters of the grid and options for gating
   ! coordinate
   !====================================================
   ! DM       - mass of the proton donor
   ! AM       - mass of the proton acceptor
   ! AGLIM    - left grid limit
   ! BGLIM    - right grid limit
   ! GLIST    - grid points along the gating coordinate
   ! NPNTSG   - number of grid points along the gating
   !            coordinate
   ! NPNTSSOLG- number of grid points along the gating
   !            coordinate for solvation calculations
   ! NGAST    - number of gating vibrational states
   !            per vibronic state
   ! GQUANT   - logical variable for quantization along
   !            the gating coordinate
   ! MGQUANT  - method of incorporation of quantum gating
   !====================================================

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! default parameters: one point along gating coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   aglim     = 0.d0
   bglim     = 0.d0
   npntsg    = 1
   npntssolg = 0
   ngast     = 0
   dm        = hmass
   am        = hmass
   gquant    = .false.
   mgquant   = 1
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' GATING(')

   if (ikey.ne.0) then

      call getopt(keywrd,ikey+8,options)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Masses of the proton donor (DM) and acceptor (AM)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (index(options,' DMASS=').ne.0) then
         dm = reada(options,index(options,' DMASS=')+7)
      else
         dm = 1.d0
      endif
      dm = dm*dalton

      if (index(options,' AMASS=').ne.0) then
         am = reada(options,index(options,' AMASS=')+7)
      else
         am = 1.d0
      endif
      am = am*dalton

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of grid points for the gating coordinate
      ! (must be a power of two)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ingrid = index(options,' NGRID=')

      if (ingrid.ne.0) then
         npntsg = reada(options,ingrid+7)
      else
         npntsg = 1
      endif

      if (npntsg.le.1) then

         write(*,'(/1x,"Calculations will be performed only for the gating distance",/,&
         &" specified in the geometry input file",/)')
         npntsg = 1

      elseif (.not.power2(npntsg)) then

         write(*,'(/1x,"*** Number of grid points (",I3,") is not a power of two.")') npntsg

         npntspow2 = 1
         do while (npntsi.gt.npntspow2.and.npntspow2.lt.maxpnt)
            npntspow2 = npntspow2 * 2
         enddo
         npntsg = 2 * npntspow2
         write(*,'(" NGRID is set to ",i3/)') npntsg

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of grid points for solvation calculations.
      ! (must be a power of two and smaller that NPNTSG)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ingrids = index(options,' NGRIDS=')

      if (ingrids.ne.0) then
         npntssolg = reada(options,ingrids+8)
      else
         npntssolg = npntsg
      endif

      npntssol2 = npntssolg

      write(*,'(/1x,"Solvation calculations on the gating coordinate grid",/)')

      if (npntssol2.le.0) then

         write(*,'(/1x,"Solvation calculations will be performed only for the gating distance",/,&
         &" specified in the geometry input file",/)')

      else if (.not.power2(npntssol2).or.npntssol2.gt.npntsg) then

         write(*,'(/1x,"*** Number of grid points for solvation calculations (",I3,")",/,&
         &" is either not a power of two or greater then the total number of grid points",/,&
         &" NPNTSSOLG is set equal to NPNTSG (",i3,")",/)') npntssol2,npntsg
         npntssolg = npntsg

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Grid (integration) limits for the gating coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ilims = index(options,' LIMITS=')
      if (ilims.ne.0) then
         iarlim = ilims+8
         islash = index(options(iarlim:),'/')
         aglim = reada(options(iarlim:iarlim+islash-2),1)
         bglim = reada(options,iarlim+islash)
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Quantization for the gating coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      igquant = index(options,' QUANTUM')
      gquant = igquant.ne.0

      if (gquant) then

         if (method.eq.2) then
            write(*,'(/1x,"*** Quantization of the gating is not compatible with METHOD=2 ***"/)')
            stop
         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Method of incorporation of gating quantum states
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ikey = index(options,' QUANTUM=')

         if (ikey.ne.0) then

            mgquant = reada(options,ikey+9)
            if (mgquant.lt.1.or.mgquant.gt.3) then
               write(*,'(/1x,"*** (in SETJOB): QUANTUM value is invalid ***"/)')
               stop
            endif

         else

            mgquant = 1

         endif

         if (mgquant.eq.1.and.(method.eq.2.or.method.eq.3)) then
            write(*,'(/1x,"*** MGQUANT=1 is not compatible with METHOD=2,3 ***"/)')
            stop
         endif

         write(6,'(/1x,"MGQUANT ",i1," will be used for quantization of the gating")') mgquant

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Number of quantum states per vibronic state
         ! for the gating coordinate
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ingast = index(options,' NGAST=')

         if (ingast.ne.0) then
            ngast = reada(options,ingast+7)
         endif

         if (ngast.gt.ngastmax) then

            write(6,'(/" Number of gating vibrational states NGAST =",I5,/,&
            &" It is larger than maximum number allowed NGASTMAX =",I5,/,&
            &" **  Program stopped  ** "//&
            &" Recommendation: Either reduce NGAST or recompile the",&
            &" program with larger NRGASTMAX")') ngast,ngastmax
            stop

         elseif (ngast.lt.1) then

            write(6,'(/" Wrong number of gating quantum states NGAST =",I5,/,&
            &" NGAST must be larger than 0 and smaller than NGASTMAX =",I5,/,&
            &" **  Program stopped  ** ")') ngast,ngastmax
            stop

         endif

      endif

      write(6,'(/1x,"Gating coordinate parameters:")')
      write(6,'(10x,"Mass of the donor (Daltons)   : ",f13.6)') dm/dalton
      write(6,'(10x,"Mass of the acceptor (Daltons): ",f13.6)') am/dalton
      write(6,'(10x,"Integration limits (Angstroms): ",2f13.6)') aglim,bglim
      write(6,'(10x,"Number of grid points:          ",i12)') npntsg
      if (gquant) then
         write(6,'(10x,"Quantization is turned on")')
         write(6,'(10x,"Number of vibrational states to include:",4x,i12)') ngast
      else
         write(6,'(10x,"Quantization is turned off")')
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Integration step (STEPG) and array of grid points (GLIST)
      ! for the gating coordinate
      ! (now it's done in module_quantum:alloc_gquant)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! IF (NPNTSG.GT.1) THEN
      !
      !    STEPG = A2BOHR*(BGLIM-AGLIM)/(NPNTSG-1)
      !    X = AGLIM*A2BOHR
      !    DO K=1,NPNTSG+1
      !       GLIST(K) = X
      !       X = X + STEPG
      !    ENDDO
      !
      ! ELSE
      !
      !    GLIST(1) = 0.D0
      !
      ! ENDIF

   endif

   !================================================================================
   ! Type of the gas-phase potential
   !================================================================================
   ! CONST  - constant potential for a fixed solute (no internal degrees of freedom)
   ! MM5    - MM EVB five site potential (for Nocera-like systems)
   ! WATER  - MM EVB potential for water cluster (Voth-Schmidt)
   ! LEPS0  - original five site LEPS 2D potential (with ET Coulomb)
   ! LEPS5  - modified five site LEPS 2D potential (without ET Coulomb)
   ! HYBRID - hybrid LEPS/harmonic 2D potential
   ! MMGEN  - MM EVB five site potential (for general O-H...O systems)
   ! GEOM=  - read external file with the geometry
   ! PARS=  - read external file with parameters
   !================================================================================

   ikey = index(keywrd,' HGAS(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in SETJOB): You MUST specify the HGAS keyword with options ***"/)')
      stop
   else
      call getopt(keywrd,ikey+6,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Type of the potential: MM5
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' CONSTANT').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Constant potential for a fixed solute
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Constant gas phase potential will be used (fixed solute)")')
      igas = 0

   elseif (index(options,' MM5').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! MM5 for Nocera-like systems
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"MM5 gas phase potential will be used")')
      igas = 1

      icoulomb = index(options,' NOCOULOMB')
      if (icoulomb.ne.0) then
         write(*,'(/1x,"Coulomb interactions in the gas phase are turned off"/)')
         coulomb = .false.
      endif

   elseif (index(options,' MMGEN').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! MM5 for Nocera-like systems
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"MMGEN general gas phase potential will be used")')
      igas = 6

      icoulomb = index(options,' NOCOULOMB')
      if (icoulomb.ne.0) then
         write(*,'(/1x,"Coulomb interactions in the gas phase are turned off"/)')
         coulomb = .false.
      endif

   elseif (index(options,' WATER').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Voth-Schmidt EVB potential for water clusters
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      write(6,'(/1x,"Voth-Schmidt EVB gas phase potential for water clusters")')
      igas = 2

      if (method.eq.2) then
         write(*,'(1x,"** You can not use METHOD=2 with WATER potential at this time ***"/)')
         stop
      endif

      idist = index(options,' DIST_ET=')
      if (idist.ne.0) then
         dist_et = reada(options,idist+9)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify DIST_ET option in GAS ***"/)')
         stop
      endif

      idist = index(options,' DIST_PT=')
      if (idist.ne.0) then
         dist_pt = reada(options,idist+9)
      endif

      ideltag = index(options,' DELTAG=')
      if (ideltag.ne.0) then
         deltag = reada(options,ideltag+8)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify DELTAG option in GAS ***"/)')
         stop
      endif

      ivept = index(options,' VEPTMETH=')
      if (ivept.ne.0) then
         veptmeth = reada(options,ivept+10)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify VEPTMETH option in GAS ***"/)')
         stop
      endif

      write(6,'(/1x,"Some water potential parameters:")')
      write(6,'(10x,"Electron transfer distance (A): ",f13.6)') dist_et
      write(6,'(10x,"Bias for ET states (kcal/mol):  ",f13.6)') deltag
      write(6,'(10x,"Method for VEPT estimate:       ",i12)')   veptmeth


   else if (index(options,' LEPS0').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Original LEPS potential (2-dimensional)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Original LEPS gas phase potential")')
      igas = 3


   else if (index(options,' LEPS5').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Modified LEPS potential (2-dimensional)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"LEPS gas phase potential will be used")')
      igas = 4

   else if (index(options,' HYBRID').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Hybrid LEPS/Harmonic potential (2-dimensional)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Hybrid LEPS/harmonic gas phase potential")')
      igas = 5

   else if (index(options,' HARM').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Harmonic potential (1-dimensional)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Harmonic gas phase potential")')
      igas = 7

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   else

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Put here your potential...
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(*,'( 1x,"*** (in SETJOB): CONSTANT, HARM, MM5, MMGEN, WATER, LEPS0, LEPS5, and HYBRID potentials are available ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read external file with the parameters
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ipars = index(options,' PARS=')

   if (ipars.eq.0) then

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the PARS= option for HGAS ***"/)')
      stop

   else

      ispa = index(options(ipars+6:),' ')
      fname = options(ipars+6:ipars+ispa+4)
      lenf = ispa - 1
      call locase(fname,lenf)
      open(1,file=fname(1:lenf),status='old')
      call set_potential(1,igas)
      close(1)
      write(6,'(/1x,"Parameters for the gas phase potential from external file <",a,">")') fname(1:lenf)
      call system("cp "//fname(1:lenf)//" "//job(1:ljob))

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read external file with the geometry of the complex
   ! and atomic charges for EVB states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   igeo = index(options,' GEOM=')

   if (igeo.eq.0) then

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the GEOM= option for HGAS ***"/)')
      stop

   else

      ! extract the name of input file
      ispa = index(options(igeo+6:),' ')
      fname = options(igeo+6:igeo+ispa+4)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for gas phase calculations from the file <",a,">")') fname(1:lenf)

      ! copy the input file to the output directory
      call system("cp "//fname(1:lenf)//" "//job(1:ljob))

      ! open the file
      open(1,file=fname(1:lenf),status='old')

      ! scan input file for number of atoms
      call scan_geo(1,natgas_)

      ! allocate the geometry/charges arrays
      call alloc_geogas(natgas_)

      ! read the gas phase geometry
      rewind 1
      call geoin(1,dm,am,natgas_,natgas,iptgas,ietgas,labgas,xyzgas,chrgas)
      close(1)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Write out the cartesian coordinates to the external
   ! file if specified in XYZOUT option
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ixyz = index(options,' XYZOUT=')

   if (ixyz.ne.0) then

      ispa = index(options(ixyz+8:),' ')
      fname = options(ixyz+8:ixyz+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for gas phase calculations is written to the file <",A,">",/,&
      &" in xyz format (use MOLDEN to view)")') fname(1:lenf)
      open(1,file=job(1:ljob)//'/'//fname(1:lenf),status='new')
      call xyzout(1,natgas,iptgas,labgas,xyzgas)
      close(1)

   endif

   !===============================================================
   ! Solvent model
   !===============================================================
   ! FRCM    - Frequency Resolved Cavity Model (Basilevsky, Rostov)
   ! ELLIPSE - Simple electrostatic model with ellipsoidal cavity
   ! TSET    - Reorganization energy matrix is reconstructed from
   !           the values given in input (new in version 5.2)
   !===============================================================

   ikey = index(keywrd,' SOLV(')

   if (ikey.eq.0) then

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the SOLV keyword with options ***"/)')
      stop

   else

      call getopt(keywrd,ikey+6,options)

   endif


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read values of dielectric constants (EPS0,EPS8) at 298.15 K
   ! [Y. Marcus, Ion Solvation, Wiley, NY, 1985, p.137-138]
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' H2O').ne.0.or.index(options,' WATER').ne.0) then

      eps0 = 78.39d0
      eps8 = 1.78d0

   elseif (index(options,' CH2CL2').ne.0.or.index(options,' DICHLOROETHANE').ne.0) then

      eps0 = 8.93d0
      eps8 = 2.02d0

   elseif (index(options,' MEOH').ne.0.or.index(options,' METHANOL').ne.0) then

      eps0 = 33.7d0
      eps8 = 5.6d0

   elseif (index(options,' ETOH').ne.0.or.index(options,' ETHANOL').ne.0) then

      eps0 = 24.6d0
      eps8 = 1.85d0

   elseif (index(options,' CH3CN').ne.0) then

      eps0 = 37.5d0
      eps8 = 1.80d0

   elseif (index(options,' DCE').ne.0) then

      eps0 = 10.4d0
      eps8 = 2.08d0

   elseif (index(options,' THF').ne.0) then

      eps0 = 7.58d0
      eps8 = 1.97d0

   elseif (index(options,' NBZ').ne.0) then

      eps0 = 34.82d0
      eps8 = 2.40d0

   elseif (index(options,' DMF').ne.0) then

      eps0 = 36.7d0
      eps8 = 2.04d0

   else

      ieps0 = index(options,' EPS0=')
      ieps8 = index(options,' EPS8=')

      if (ieps0.ne.0) then
         eps0 = reada(options,ieps0+6)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify EPS0 constant in SOLV ***"/)')
         stop
      endif

      if (ieps8.ne.0) then
         eps8 = reada(options,ieps8+6)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify EPS8 constant in SOLV ***"/)')
         stop
      endif

   endif

   !-- initialize inverse Pekar factor f_0
   f0 = four*pi*eps0*eps8/(eps0 - eps8)

   !-- symmetrization options for reorganization energy matrices

   if (index(options,' SYMT').ne.0) then
      symt=.true.
   else
      symt=.false.
   endif

   if (index(options,' SYMPT').ne.0) then
      sympt=.true.
   else
      sympt=.false.
   endif

   if (index(options,' SYMET').ne.0) then
      symet=.true.
   else
      symet=.false.
   endif

   if (index(options,' REDDENS').ne.0) then
      reddens=.true.
   else
      reddens=.false.
   endif

   if (index(options,' NOSYMD').ne.0) then
      nosymd=.true.
   else
      nosymd=.false.
   endif

   write(6,'(/1x,"Static solvent parameters:")')
   write(6,'(10x,"Static dielectric constant:  ",f13.6)') eps0
   write(6,'(10x,"Optical dielectric constant: ",f13.6)') eps8
   !write(6,'(10x,"Inverse Pekar factor:       ",f13.6)') f0
   write(6,'(10x,"Total charge of the solute:  ",f13.6)') charge

   if (symt) then
      write(6,'("Solvation T matrices symmetrized for PCET in a symmetric system")')
   endif

   if (sympt) then
      write(6,'("Solvation T matrices symmetrized for PT in a symmetric system")')
   endif

   if (symet) then
      write(6,'("Solvation T matrices symmetrized for ET in a symmetric system")')
   endif

   if (reddens) then
      write(6,'("Reduced basis used for solvation matrices")')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Solvation model
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' TSET').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Reconstruction of reorganization energy matrix from
      ! input values of partial reorganization energies
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Reorganization energy matrices will be reconstructed from partial reorganization energies specified in input")')
      isolv = 0


   elseif (index(options,' ELLIPSE').NE.0) THEN

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Simple electrostatic model with ellipsoidal cavity
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"Ellipsoidal model will be used for solvation calculations")')
      isolv = 1

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read values for ellipsoid parameters
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! A - major semiaxis
      ! B - minor semiaxis
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      imajor = index(options,' A=')
      if (imajor.ne.0) then
         a = reada(options,imajor+3)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify the major axis of ellipsoid ***"/)')
         stop
      endif

      iminor = index(options,' B=')
      if (iminor.ne.0) then
         b = reada(options,iminor+3)
      else
         write(*,'(/1x,"*** (in SETJOB): You MUST specify the minor axis of ellipsoid ***"/)')
         stop
      endif

      ba  = b/a
      ba2 = ba*ba
      l0  = dsqrt(1.d0/(1.d0-ba2))
      r   = 2.d0*a/l0

      write(6,'(/1x,"Ellipsoidal cavity parameters:")')
      write(6,'(10x,"Major semiaxis:      ",f13.6)') a
      write(6,'(10x,"Minor semiaxis:      ",f13.6)') b
      write(6,'(10x,"Interfocal distance: ",f13.6)') r

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   elseif (index(options,' FRCM').ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Frequency Resolved Cavity Model (Basilevsky, Rostov)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      write(6,'(/1x,"FRCM model will be used for solvation calculations")')
      isolv = 2

      if (method.eq.2) then
         write(*,'(1x,"** You can not use METHOD=2 with FRCM solvation model ***"/)')
         stop
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read values for cavity parameters in FRCM
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! KAPPA - factor for VdW radii
      ! DELTA - the width of the layer between two cavities
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      !-- default values for different solvents

      if (index(options,' H2O').ne.0.or.index(options,' WATER').ne.0) then
         kappa = 0.9d0
         delta = 1.0d0
      elseif (index(options,' CH2CL2').ne.0) then
         kappa = 0.9d0
         delta = 1.2d0
      elseif (index(options,' MEOH').ne.0) then
         kappa = 0.9d0
         delta = 1.1d0
      elseif (index(options,' ETOH').ne.0) then
         kappa = 0.9d0
         delta = 1.4d0
      elseif (index(options,' CH3CN').ne.0) then
         kappa = 0.9d0
         delta = 1.8d0
      elseif (index(options,' DCE').ne.0) then
         kappa = 0.9d0
         delta = 2.3d0
      elseif (index(options,' THF').ne.0) then
         kappa = 0.9d0
         delta = 2.3d0
      elseif (index(options,' NBZ').ne.0) then
         kappa = 0.9d0
         delta = 3.7d0
      elseif (index(options,' DMF').ne.0) then
         kappa = 0.9d0
         delta = 2.7d0
      else
         kappa = 0.9d0
         delta = 1.0d0
      endif

      ikappa = index(options,' KAPPA=')
      idelta = index(options,' DELTA=')

      if (ikappa.ne.0) then
         kappa = reada(options,ikappa+7)
         write(*,'(/1x,"The custom value of KAPPA parameter was specified in the input"/)')
      endif

      if (idelta.ne.0) then
         delta = reada(options,idelta+7)
         write(*,'(/1x,"The custom value of DELTA parameter was specified in the input"/)')
      endif

      write(6,'(/1x,"FRCM cavity parameters:")')
      write(6,'(10x,"Factor for VdW radii (KAPPA):     ",f13.6)') kappa
      write(6,'(10x,"Width of the inner layer (DELTA): ",f13.6," Angstroms")') delta

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   else

      write(*,'(/1x,"*** (in SETJOB): You MUST specify the solvation model (TSET/ELLIPSE/FRCM) with SOLV ***"/)')
      stop

   endif

   itread = index(options,' TREAD=')
   if (itread.ne.0) then
      write(6,'(/1x,"Reorganization energy matrices will be read from disk file in SETMAT")')
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read external file with the geometry of the complex
   ! for solvation calculations and charges for EVB states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   igeom = index(options,' GEOM=')

   if (igeom.ne.0) then

      ! extract the name of the geometry input file
      ispa = index(options(igeom+6:),' ')
      fname = options(igeom+6:igeom+ispa+4)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for solvation calculations from the file <",a,">")') fname(1:lenf)

      ! copy the input file to the output directory
      call system("cp "//fname(1:lenf)//" "//job(1:ljob))

      ! open the file
      open(1,file=fname(1:lenf),status='old')

      ! scan input file for number of atoms
      call scan_geo(1,natsol_)

      ! allocate the geometry/charges arrays
      call alloc_geosol(natsol_)

      ! read the gas phase geometry
      rewind 1
      call geoin(1,dm,am,natsol_,natsol,iptsol,ietsol,labsol,xyzsol,chrsol)

      if (isolv.eq.2) then

         if (npntssol.ge.0.and.quantum_proton) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! put the hydrogen atom at the center of mass
            ! of the pt interface
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            nhpt = iptsol(2)
            xyzsol(1,nhpt) = 0.d0
            xyzsol(2,nhpt) = 0.d0
            xyzsol(3,nhpt) = 0.d0
         endif

         rewind 1
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! In case of FRCM solvation model:
         ! - read specific FRCM keywords
         ! - initialize internal COMMON-blocks used in FRCM
         !   calculations of the reorganization energy matrices
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         call frcminit(1,eps0,eps8,kappa,delta,charge,natsol,iptsol,labsol,xyzsol,chrsol)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! - use specific FRCM keywords
         ! - build spheres around the atoms
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! call getsfe(1,xyzsol,natsol)
         ! now it is called from within FRCMINIT

      endif
      close(1)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! For ellipsoidal model check whether all the point
      ! charges are on the straight line (x-axis)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (isolv.eq.1) then

         ok = straight(natsol,xyzsol)
         if (.not.ok) then
            write(*,'(/1x,"*** (in SETJOB): ",/,&
            &"The geometry for solvation calculations using an ellipsoidal model",/,&
            &" should be LINEAR !"/)')
            stop
         endif

      endif

      if (isolv.eq.0) then
         if (npntssol.ge.0.and.quantum_proton) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! put the hydrogen atom at the center of mass
            ! of the pt interface
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            nhpt = iptsol(2)
            xyzsol(1,nhpt) = 0.d0
            xyzsol(2,nhpt) = 0.d0
            xyzsol(3,nhpt) = 0.d0
         endif
      endif

   else

      write(*,'(/1x,"*** (in SETJOB): You MUST specify file with the geometry (GEOM=) in SOLV keyword ***"/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Write out the cartesian coordinates to the external
   ! file if option XYZOUT is specified
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ixyz = index(options,' XYZOUT=')

   if (ixyz.ne.0.and.isolv.gt.0) then

      ispa = index(options(ixyz+8:),' ')
      fname = options(ixyz+8:ixyz+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,"Geometry for solvation calculations is written to the file <",A,">",/,&
      &"in xyz format (use MOLDEN to view)")') fname(1:lenf)
      open(1,file=job(1:ljob)//'/'//fname(1:lenf),status='new')
      call xyzout(1,natsol,iptsol,labsol,xyzsol)
      close(1)

   endif

   return
   
contains
!=======================================================================
!  Internal procedures
!=======================================================================
   
   subroutine scan_geo(iunit,n)

   ! scans geometry input file (iunit) and extracts
   ! number of atoms (n)

      implicit none
      integer, intent(in)  :: iunit
      integer, intent(out) :: n

      character(80) :: line
      integer :: io_status

      ! title records
      read(iunit,'(a)') line
      read(iunit,'(a)') line
      read(iunit,'(a)') line

      n = 0

      do
         read(iunit,'(a)',iostat=io_status) line
         if (io_status.ne.0) then
            write(*,*) "Error/eof while scanning geometry input file"
            close(iunit)
            stop
         endif
         if (line.eq. ' ') exit
         n = n + 1
      enddo

      if (n.le.0) then
         write(*,*) "No atoms in the geometry input file"
         close(iunit)
         stop
      endif

      return

   end subroutine scan_geo


   !=======================================================================
   logical function straight(nat,xyz)
   ! Checks whether all the atoms are on the straight line (x-axis)

      implicit none

      integer, intent(in) :: nat
      real*8,  intent(in), dimension(:,:) :: xyz

      integer :: i, j

      straight = .true.

      do i=1,nat
         do j=2,3
            if (dabs(xyz(j,i)).gt.1.d-7) then
               straight = .false.
               return
            endif
         enddo
      enddo

      return
   end function straight

   !=======================================================================
   logical function power2(n)
   ! Checks whether the given integer (N) is a power of two

      implicit none
      integer, intent(in) :: n

      integer, parameter :: itwo = 2
      integer :: nw

      nw = n
      if (nw.lt.2) then
         power2 = .false.
         return
      else
         power2 = .true.
      endif

      do while (nw.gt.1)
         if (mod(nw,itwo).ne.0) then
            power2 = .false.
            return
         else
            nw = nw/itwo
         endif
      enddo

      return
   end function power2

end subroutine setjob
