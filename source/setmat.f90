subroutine setmat

!===================================================================C
!  Calculates the matrices on the grid along the proton coordinate
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-11-04 22:49:28 $
!  $Revision: 5.4 $
!  $Log: not supported by cvs2svn $
!  Revision 5.3  2010/11/04 22:43:09  souda
!  Next iteration... and two additional Makefiles for building the code with debug options.
!
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!===================================================================C

   use pardim
   use keys
   use strings
   use cst
   use control
   use parsol
   use geogas
   use geosol
   use gasmat
   use solmat
   use quantum
   use frcm
   use elcm
   use msevb_water
   use potential

   implicit none

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! One COMMON-block from Voth-Schmitt related routines
   ! with transformed cartesian coordinates for the gas phase:
   ! to adjust proton position (X0,Y0,Z0)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! PARAMETER(NATMAX=31)
   ! REAL*8 X0,Y0,Z0
   ! COMMON /EVBATO000/ X0(NATMAX),Y0(NATMAX),Z0(NATMAX)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! One COMMON-block from FRCM related routines
   ! with cartesian coordinates for the solute:
   ! to adjust proton position (COORD)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! COMMON /GEOMXYZ/ COORD(3,NUMATM)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Local variables
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   character(1024) :: options
   character(  40) :: fname
   logical :: tread, deriv, derivg
   logical :: file_exists

   integer :: i, j, ndgas, nhgas, nagas, nkr, kr, kleft, krigh, k
   integer :: ikey, itread, ispa, lenf, kk, kg, npnts1, npnts2
   integer :: ndsol, nhsol, nasol, kratio, nsc, kgratio, ngsc
   integer :: kgsol, ksol, l, kp, it, np2, ng2, kpoint, kgpoint
   integer :: ii, jj, itwrite

   real(8) :: xint, xgint, rminv, totm, rr, q, rdh, rah
   real(8) :: h0fill, dh0fill, d2h0fill, dgh0fill, dg2h0fill
   real(8) :: alim1, blim1, alim2, blim2
   real(8) :: ratio, gratio, rl, rcoef
   real(8) :: talp1, talp2, tbeta, tdet, tr1a1b, tr1a2a
   real(8) :: ersum, erdif, ersq, tanth, sqtan

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Local arrays
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   integer, allocatable, dimension(:) :: ksta, kend
   real(8), dimension(4,4) :: h0k, dh0k, d2h0k, dgh0k, dg2h0k
   real(8), dimension(4,4) :: tk, tinfk, trk, trinfk
   real(8), dimension(2,2) :: t1k
   real(8), dimension(4)   :: dtk, dtinfk, d2tinfk

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   tread = .false.
   deriv  = method.eq.2.and.(igas.eq.1.or.igas.eq.3.or.igas.eq.4.or.igas.eq.5.or.igas.eq.6).and.isolv.eq.1
   derivg = gquant.and.mgquant.eq.1.and.(igas.eq.3.or.igas.eq.4.or.igas.eq.5)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation of matrices (kinetic energy and first derivative)
   ! in the basis of delta functions on the grid along the proton
   ! coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   xint = rlist(npnts) - rlist(1)
   write(*,*) ' '
   write(*,*) ' Kinetic energy matrix for proton...'
   call gridke(npnts,xint,hke)
   write(*,*) ' done...'

   do i=1,npnts
      do j=1,npnts
         hke(i,j) = hke(i,j)/pm
      enddo
   enddo
   write(6,*) ' done'

   if (deriv) then
      write(6,*) ' '
      write(6,*) ' First derivative matrix for proton...'
      call griddx(npnts,xint,dx)
      write(6,*) ' done'
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculation of matrices (kinetic energy and first derivative)
   ! in the basis of delta functions on the grid along the gating
   ! coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (gquant) then

      xgint = glist(npntsg) - glist(1)
      write(6,*) ' '
      write(6,*) ' kinetic energy matrix for gating mode...'
      call gridke(npntsg,xgint,hgke)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Inversed reduced mass for gating coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      rminv = 1.d0/dm + 1.d0/am

      do i=1,npntsg
         do j=1,npntsg
            hgke(i,j) = hgke(i,j)*rminv
         enddo
      enddo
      write(6,*) ' done'

      if (derivg) then
         write(6,*) ' '
         write(6,*) ' First derivative matrix...'
         call griddx(npntsg,xgint,dgx)
         write(6,*) ' done'
      endif

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate gas phase Hamiltonian matrix elements,
   ! reorganization energies matrix elements and its
   ! derivatives on the grid along the proton coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(6,*)
   write(6,*)' gas-phase hamiltonian matrix and its derivatives...'
   write(6,*)' reorganization energy matrices, their derivatives...'

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Start calculations on the grid
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Gas-phase matrix elements and its derivatives
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! In case of water potential initialize MS EVB
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (igas.eq.2) call initevb

   ndgas = iptgas(1)
   nhgas = iptgas(2)
   nagas = iptgas(3)
   totm = dm + am

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! External loop over grid points along the
   ! gating coordinate.
   ! Center of mass is fixed at the origin.
   ! Note that we assume that the electronic donor
   ! and acceptor are fixed in space.
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (npntsg.gt.0) then
      nkr = npntsg
   else
      nkr = 1
   endif

   allocate (ksta(nkr),kend(nkr))

   do kr=1,nkr

      rr = glist(kr)*bohr2a

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! In case RR=0 there is only one point along
      ! gating coordinate corresponding to the
      ! input geometry
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (rr.ne.0.d0) then
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Move donor and acceptor so that the center of mass
         ! stays at the coordinate system origin
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         xyzgas(1,ndgas) = -am*rr/totm
         xyzgas(1,nagas) =  dm*rr/totm
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Internal loop over grid points along the
      ! quantum proton coordinate.
      ! Center of mass is fixed at the system origin.
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      kleft = 1
      krigh = npnts

      do k=1,npnts

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Position of the proton (Angstroms) on the x-axis
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         q = rlist(k)*bohr2a

         xyzgas(1,nhgas) = q
         xyzgas(2,nhgas) = 0.d0
         xyzgas(3,nhgas) = 0.d0

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Distance between the proton
         ! and ET donor/acceptor (Angstroms)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         rdh = q - xyzgas(1,ndgas)
         rah = xyzgas(1,nagas) - q

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! The matrices are calculated only when the proton
         ! is located inside the safe (reasonable) interval
         ! between donor and acceptor
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (rdh.lt.rhmin) then
            kleft = kleft + 1
            cycle
         endif

         if (rah.lt.rhmin) then
            krigh = k - 1
            exit
         endif

         if (igas.eq.1) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! MM5 potential for Nocera-like systems
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_mm5(h0k)
            if (deriv) then
               call  dh0mat_mm5(dh0k)
               call d2h0mat_mm5(d2h0k)
            endif

         elseif (igas.eq.2) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! MS-EVB (Voth, Schmitt) potential for water clusters
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            x0evb(nhgas) = xyzgas(1,nhgas)
            call msevb(h0k)
            if (deriv) then
               write(*,'(/1x,'' No derivatives for MSEVB'')')
               stop
            endif

         elseif (igas.eq.3) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! LEPS 2D potential (original version)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_leps0(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

         elseif (igas.eq.4) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! LEPS 2D potential (simplified version)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_leps2(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

         elseif (igas.eq.5) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Hybrid LEPS/harmonic 2D potential
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_hyb(h0k,dh0k,d2h0k,dgh0k,dg2h0k)

         elseif (igas.eq.6) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! MM5 potential for general O-H...O systems
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_mm5gen(h0k)
            if (deriv) then
               call  dh0mat_mm5gen(dh0k)
               call d2h0mat_mm5gen(d2h0k)
            endif

         elseif (igas.eq.7) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Harmonic potential for general D-H...A systems
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call h0mat_harm(h0k)
            if (deriv) then
               call  dh0mat_harm(dh0k)
               call d2h0mat_harm(d2h0k)
            endif

         else

            write(*,'(/1x,'' ERROR in SETMAT: unknown IGAS='',i2)') igas
            stop

         endif

         do i=1,4
            do j=1,4
               h0(i,j,k,kr) =   h0k(i,j)
               if (deriv) then
                  dh0(i,j,k,kr) =  dh0k(i,j)
                 d2h0(i,j,k,kr) = d2h0k(i,j)
               endif
               if (derivg) then
                  dgh0(i,j,k,kr) =  dgh0k(i,j)
                 dg2h0(i,j,k,kr) = dg2h0k(i,j)
               endif
            enddo
         enddo

      enddo  !end of loop over proton grid

      ksta(kr) = kleft
      kend(kr) = krigh

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Fill the matrices for forbidden (excluded) k-points
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do k=1,kleft-1
         do i=1,4
            do j=1,4
               h0fill = h0(i,j,kleft,kr)
               h0(i,j,k,kr) = h0fill
               if (deriv) then
                  dh0fill = dh0(i,j,kleft,kr)
                  d2h0fill = d2h0(i,j,kleft,kr)
                  dh0(i,j,k,kr) =  dh0fill
                 d2h0(i,j,k,kr) = d2h0fill
               endif
               if (derivg) then
                  dgh0fill = dgh0(i,j,kleft,kr)
                  dg2h0fill = dg2h0(i,j,kleft,kr)
                  dgh0(i,j,k,kr) =  dgh0fill
                 dg2h0(i,j,k,kr) = dg2h0fill
               endif
            enddo
         enddo
      enddo

      do k=krigh+1,npnts
         do i=1,4
            do j=1,4
               h0fill = h0(i,j,krigh,kr)
               h0(i,j,k,kr) = h0fill
               if (deriv) then
                  dh0fill = dh0(i,j,krigh,kr)
                  d2h0fill = d2h0(i,j,krigh,kr)
                  dh0(i,j,k,kr) =  dh0fill
                 d2h0(i,j,k,kr) = d2h0fill
               endif
               if (derivg) then
                  dgh0fill = dgh0(i,j,krigh,kr)
                  dg2h0fill = dg2h0(i,j,krigh,kr)
                  dgh0(i,j,k,kr) =  dgh0fill
                 dg2h0(i,j,k,kr) = dg2h0fill
               endif
            enddo
         enddo
      enddo

   enddo  !end of loop over gating grid


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reorganization energy matrices and their derivatives
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Read external file with the reorganization energy matrices
   ! if TREAD option is specified for the keyword SOLVENT()
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' SOLVENT(')
   call getopt(keywrd,ikey,options)

   itread = index(options,' TREAD=')

   IF (ITREAD.NE.0) THEN

      ispa = index(options(itread+7:),' ')
      fname = options(itread+7:itread+ispa+5)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,''Reorganization energy matrices from the binary file <'',a,''>'')') fname(1:lenf)

      inquire(file=fname(1:lenf),exist=file_exists)

      if (.not.file_exists) then
         write(*,'(/1x,"File ",a," does not exist. Check the file name in TREAD option.")') fname(1:lenf)
         stop
      endif

      open(1,file=fname(1:lenf),status='old',form='unformatted')

      read(1) npnts1, alim1, blim1, npnts2, alim2, blim2

      if (npnts1.ne.npnts.or.alim1.ne.alim.or.blim1.ne.blim.or.&
         &npnts2.ne.npntsg.or.alim2.ne.aglim.or.blim2.ne.bglim) then
         write(6,'(/1x,''Reorganization energy matrices stored in the binary file <'',a,''>'',/,&
         &'' are not consistent with the current grid parameters.'',/,&
         &'' The full calculation will be performed...'')') fname(1:lenf)
         close(1)
         goto 10
      endif

      read(1) ((((    t(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      read(1) ((((   tr(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      read(1) (((( tinf(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      read(1) ((((trinf(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)

      if (deriv) then
         read(1) (((    dt(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         read(1) (((   dtr(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         read(1) ((( dtinf(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         read(1) (((d2tinf(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
      endif

      close(1)

      tread = .true.
      goto 20

   endif

   10 continue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Start calculations on the grid
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ndsol = iptsol(1)
   nhsol = iptsol(2)
   nasol = iptsol(3)

   if (npntssol.gt.0) then
      kratio = npnts/npntssol
      nsc = npntssol + 1
   elseif(npntssol.eq.0) then
      kratio = 0
      nsc = 1
   else
      kratio = 0
      nsc = 1
   endif

   if (npntssolg.gt.0) then
      kgratio = nkr/npntssolg
      ngsc = npntssolg + 1
   elseif(npntssolg.eq.0) then
      kgratio = 0
      ngsc = 1
   else
      kgratio = 0
      ngsc = 1
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over gating grid points
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do kgsol=1,ngsc

      kg = 1 + (kgsol-1)*kgratio
      rr = glist(kg)*bohr2a

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Gating coordinate value:
      ! If several points are to be calculated then
      ! move proton donor and acceptor so that the
      ! center of mass stays at the coordinate system
      ! origin.
      ! Otherwise use the input geometry
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (npntssolg.gt.0) then
         xyzsol(1,ndsol) = -am*rr/totm
         xyzsol(1,nasol) =  dm*rr/totm
      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Loop over proton grid points
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do ksol=1,nsc

         k = 1 + (ksol-1)*kratio

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Position of the proton (Angstroms)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if (npntssol.gt.0) then
            xyzsol(1,nhsol) = rlist(k)*bohr2a
            xyzsol(2,nhsol) = 0.d0
            xyzsol(3,nhsol) = 0.d0
         elseif (npntssol.eq.0) then
            xyzsol(1,nhsol) = 0.d0
            xyzsol(2,nhsol) = 0.d0
            xyzsol(3,nhsol) = 0.d0
         endif

         q = xyzsol(1,nhsol)
         rdh = q - xyzsol(1,ndsol)
         rah = xyzsol(1,nasol) - q

         if (rdh.lt.rhmin) cycle
         if (rah.lt.rhmin) exit

         if (isolv.eq.1) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Ellipsoidal model
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            call tmat(tk,tinfk,trk,trinfk)

         elseif (isolv.eq.2) then
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! FRCM model
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (npntssol.ge.0) then
               write(6,'(/'' proton at:'',f20.6/)')xyzsol(1,nhsol)
            else
               write(6,'(/'' proton position is fixed '')')
               write(6,'('' proton at: '',f20.6/)')xyzsol(1,nhsol)
            endif

            do i=1,3
               coord(i,ndsol) = xyzsol(i,ndsol)
               coord(i,nhsol) = xyzsol(i,nhsol)
               coord(i,nasol) = xyzsol(i,nasol)
            enddo
            call solint(tk,tinfk,trk,trinfk)

         else
            write(*,'(/1x,'' ERROR in SETMAT: ISOLV='',I2)')
            stop
         endif

         do i=1,4
            do j=1,4
               t    (i,j,k,kg) =  tk    (i,j)
               tr   (i,j,k,kg) =  trk   (i,j)
               tinf (i,j,k,kg) =  tinfk (i,j)
               trinf(i,j,k,kg) =  trinfk(i,j)
            enddo
         enddo

      enddo   !ksol: over proton grid

   enddo   !kgsol: over gating grid

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Fill out the intermediate grid points
   ! along the proton coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (kratio.gt.1) then

      ratio = dble(kratio)

      do kgsol=1,ngsc

         kg = 1+ (kgsol-1)*kgratio

         do k=1,npnts,kratio
            do l=1,kratio-1
               rl = dble(l)
               rcoef = rl/ratio
               do i=1,4
                  do j=1,4
                     t(i,j,k+l,kg)     = t(i,j,k,kg)     + rcoef*(t(i,j,k+kratio,kg)-t(i,j,k,kg))
                     tr(i,j,k+l,kg)    = tr(i,j,k,kg)    + rcoef*(tr(i,j,k+kratio,kg)-tr(i,j,k,kg))
                     tinf(i,j,k+l,kg)  = tinf(i,j,k,kg)  + rcoef*(tinf(i,j,k+kratio,kg)-tinf(i,j,k,kg))
                     trinf(i,j,k+l,kg) = trinf(i,j,k,kg) + rcoef*(trinf(i,j,k+kratio,kg)-trinf(i,j,k,kg))
                  enddo !j
               enddo !i
            enddo !l
         enddo !k

      enddo !kgsol

   elseif (kratio.eq.0) then

      do kgsol=1,ngsc

         kg = 1+ (kgsol-1)*kgratio

         do k=2,npnts
            do i=1,4
               do j=1,4
                  t    (i,j,k,kg) = t    (i,j,1,kg)
                  tr   (i,j,k,kg) = tr   (i,j,1,kg)
                  tinf (i,j,k,kg) = tinf (i,j,1,kg)
                  trinf(i,j,k,kg) = trinf(i,j,1,kg)
               enddo
            enddo
         enddo

      enddo

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Fill out the intermediate grid points
   ! along the gating coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (kgratio.gt.1) then

      gratio = dble(kgratio)

      do kp=1,npnts

         do k=1,nkr,kgratio
            do l=1,kgratio-1
               rl = dble(l)
               rcoef = rl/gratio
               do i=1,4
                  do j=1,4
                     t(i,j,kp,k+l)     = t(i,j,kp,k)     + rcoef*(t(i,j,kp,k+kgratio)-t(i,j,kp,k))
                     tr(i,j,kp,k+l)    = tr(i,j,kp,k)    + rcoef*(tr(i,j,kp,k+kgratio)-tr(i,j,kp,k))
                     tinf(i,j,kp,k+l)  = tinf(i,j,kp,k)  + rcoef*(tinf(i,j,kp,k+kgratio)-tinf(i,j,kp,k))
                     trinf(i,j,kp,k+l) = trinf(i,j,kp,k) + rcoef*(trinf(i,j,kp,k+kgratio)-trinf(i,j,kp,k))
                  enddo !j
               enddo !i
            enddo !l
         enddo !k

      enddo !kp

   elseif (kgratio.eq.0) then

      do kp=1,npnts

         do k=2,nkr
            do i=1,4
               do j=1,4
                  t    (i,j,kp,k) = t    (i,j,kp,1)
                  tr   (i,j,kp,k) = tr   (i,j,kp,1)
                  tinf (i,j,kp,k) = tinf (i,j,kp,1)
                  trinf(i,j,kp,k) = trinf(i,j,kp,1)
               enddo
            enddo
         enddo

      enddo

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Fill the matrices for forbidden (excluded) k-points
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   do kg=1,nkr

      kleft = ksta(kg)
      krigh = kend(kg)

      do k=1,kleft-1
         do i=1,4
            do j=1,4
               t    (i,j,k,kg) = t    (i,j,kleft,kg)
               tr   (i,j,k,kg) = tr   (i,j,kleft,kg)
               tinf (i,j,k,kg) = tinf (i,j,kleft,kg)
               trinf(i,j,k,kg) = trinf(i,j,kleft,kg)
            enddo
         enddo
      enddo

      do k=krigh+1,npnts
         do i=1,4
            do j=1,4
               t    (i,j,k,kg) = t    (i,j,krigh,kg)
               tr   (i,j,k,kg) = tr   (i,j,krigh,kg)
               tinf (i,j,k,kg) = tinf (i,j,krigh,kg)
               trinf(i,j,k,kg) = trinf(i,j,krigh,kg)
            enddo
         enddo
      enddo

   enddo

   20 continue

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Inverse truncated matrix [t'(trunc)]^{-1}
   ! and first and second derivatives
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   do kg=1,nkr

      kleft = ksta(kg)
      krigh = kend(kg)

      do k=1,npnts

         if (k.lt.kleft) cycle
         if (k.gt.krigh) exit

         do i=2,3
            do j=2,3
               ii = i-1
               jj = j-1
               t1k(ii,jj) = tr(i,j,k,kg)
            enddo
         enddo

         talp1 = t1k(1,1)
         talp2 = t1k(2,2)
         tbeta = t1k(1,2)
         tdet  = talp1*talp2-tbeta*tbeta
         if (dabs(tdet).lt.1.d-4) then
            write(*,'(/1x,''The truncated reorganization energy matrix'',/,&
            &'' at the gridpoint '',2I3,'' is singular (Det[TR]='',G12.6,'')'',/,&
            &'' Further calculations are meaningless...'',/,&
            &'' (Check charges for EVB states)''/)') k,kg,tdet
            stop
         endif
         t1(1,1,k,kg) = talp2/tdet
         t1(2,2,k,kg) = talp1/tdet
         t1(1,2,k,kg) = -tbeta/tdet
         t1(2,1,k,kg) = t1(1,2,k,kg)

         if (deriv.and..not.tread) then

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Derivatives of the reorganization energies matrices
            ! Note that the derivatives are available only for
            ! ellipsoidal model (ISOLV=1)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (npntssol+npntssolg.eq.0.or.isolv.ne.1) then

               do i=1,4
                   dtinf(i,k,kg) = 0.d0
                  d2tinf(i,k,kg) = 0.d0
               enddo
               do i=1,4
                  dt (i,k,kg) =  0.d0
                  dtr(i,k,kg) =  0.d0
               enddo

            else

               if (npntssol.gt.0) then
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Position of the proton (Angstroms)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  q = rlist(k)*bohr2a
                  xyzsol(1,nhsol) = q
                  xyzsol(2,nhsol) = 0.d0
                  xyzsol(3,nhsol) = 0.d0
               endif

               if (npntssolg.gt.0) then
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Positions of the proton donor and
                  ! acceptor (Angstroms)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  rr = glist(kg)*bohr2a
                  xyzsol(1,ndsol) = -am*rr/totm
                  xyzsol(1,nasol) =  dm*rr/totm
               endif

               call dtinfm(dtinfk)
               call d2tinfm(d2tinfk)
               call dtm(dtk)
               do i=1,4
                   dtinf(i,k,kg) =  dtinfk(i)
                  d2tinf(i,k,kg) = d2tinfk(i)
               enddo
               do i=1,4
                  dt(i,k,kg) =  dtk(i)
               enddo

               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               ! Derivatives of the reduced reorganization matrices
               !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

               do it=1,4
                  if (it.eq.1) then
                     dtr(it,k,kg) = dt(1,k,kg)
                  else
                     dtr(it,k,kg) = dt(it,k,kg) - dt(1,k,kg)
                  endif
               enddo

            endif

         endif

      enddo

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Fill the matrices for forbidden (excluded) k-points
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do k=1,kleft-1
         do i=1,2
            do j=1,2
               t1(i,j,k,kg) = t1(i,j,kleft,kg)
            enddo
         enddo
      enddo

      do k=krigh+1,npnts
         do i=1,2
            do j=1,2
               t1(i,j,k,kg) = t1(i,j,krigh,kg)
            enddo
         enddo
      enddo

      if (deriv.and..not.tread) then

         do k=1,kleft-1
            do i=1,4
               dtinf(i,k,kg) = dtinf(i,kleft,kg)
               d2tinf(i,k,kg) = d2tinf(i,kleft,kg)
               dt(i,k,kg) = dt(i,kleft,kg)
               dtr(i,k,kg) = dtr(i,kleft,kg)
            enddo
         enddo

         do k=krigh+1,npnts
            do i=1,4
               dtinf(i,k,kg) = dtinf(i,krigh,kg)
               d2tinf(i,k,kg) = d2tinf(i,krigh,kg)
               dt(i,k,kg) = dt(i,krigh,kg)
               dtr(i,k,kg) = dtr(i,krigh,kg)
            enddo
         enddo

      endif

   enddo

   deallocate (ksta,kend)
   write(*,*) ' Done'

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Write the reorganization energy matrices into an external file
   ! if TWRITE= option is specified for the keyword SOLVENT()
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   itwrite = index(options,' TWRITE=')

   if (itwrite.ne.0) then

      ispa = index(options(itwrite+8:),' ')
      fname = options(itwrite+8:itwrite+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1x,''Reorganization energy matrices are dumped to the binary file <'',a,''>'')') fname(1:lenf)

      inquire(file=fname(1:lenf),exist=file_exists)

      if (file_exists) then
         write(*,'(/1x,"File ",a," already exists and might be important. Consider changing the file name in TWRITE.")') fname(1:lenf)
         stop
      endif

      open(1,file=fname(1:lenf),status='new',form='unformatted')
      write(1) npnts, alim, blim, npntsg, aglim, bglim
      write(1) ((((    t(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      write(1) ((((   tr(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      write(1) (((( tinf(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)
      write(1) ((((trinf(i,j,kk,kg),i=1,4),j=1,4),kk=1,npnts),kg=1,nkr)

      if (deriv) then
         write(1) (((    dt(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         write(1) (((   dtr(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         write(1) ((( dtinf(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
         write(1) (((d2tinf(i,kk,kg),i=1,4),kk=1,npnts),kg=1,nkr)
      endif

      close(1)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !(HYD) CHANGED ABOVE SO THAT NICE TMATRICES WRITTEN AT NPNTS/2 and
   !      NPNTS/2 +1 that is at either point around r_p=0.
   !      NOTE that if NGRIDS=0 then all tmat elements have the
   !      same value at all points along rp.  That is because in the
   !      beginning, tmats are calculated once for r_p=0.
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !(AVS) Gating coordinate point added...
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   np2 = npnts/2
   if (nkr.gt.1) then
      ng2 = nkr/2
   else
      ng2 = 1
   endif

   open(72,file=job(1:ljob)//'/tmat.dat',status='new',form='formatted')

   write(72,'("Proton grid point: ",i4,";  Limits: from ",f8.3," to ",f8.3)') np2, alim, blim
   write(72,'("Gating grid point: ",i4,";  Limits: from ",f8.3," to ",f8.3)') ng2, aglim, bglim
   write(72,*)
   write(72,*) '---------- TMAT INNERTIAL (kcal/mol)'
   write(72,'(4(3x,f15.6))') (    t(1,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (    t(2,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (    t(3,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (    t(4,j,np2,ng2),j=1,4)
   write(72,*)
   write(72,*) '---------- REDUCED TMAT INNERTIAL (kcal/mol)'
   write(72,'(4(3x,f15.6))') (   tr(1,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (   tr(2,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (   tr(3,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (   tr(4,j,np2,ng2),j=1,4)
   write(72,*)
   write(72,*) '---------- TMAT INFINITY (kcal/mol)'
   write(72,'(4(3x,f15.6))') ( tinf(1,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') ( tinf(2,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') ( tinf(3,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') ( tinf(4,j,np2,ng2),j=1,4)
   write(72,*)
   write(72,*) '---------- REDUCED TMAT INFINITY (kcal/mol)'
   write(72,'(4(3x,f15.6))') (trinf(1,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (trinf(2,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (trinf(3,j,np2,ng2),j=1,4)
   write(72,'(4(3x,f15.6))') (trinf(4,j,np2,ng2),j=1,4)
   write(72,*)
   write(72,*) '---------- TOTAL TMAT (kcal/mol)'
   write(72,'(4(3x,f15.6))') ((t(1,j,np2,ng2) + tinf(1,j,np2,ng2)),j=1,4)
   write(72,'(4(3x,f15.6))') ((t(2,j,np2,ng2) + tinf(2,j,np2,ng2)),j=1,4)
   write(72,'(4(3x,f15.6))') ((t(3,j,np2,ng2) + tinf(3,j,np2,ng2)),j=1,4)
   write(72,'(4(3x,f15.6))') ((t(4,j,np2,ng2) + tinf(4,j,np2,ng2)),j=1,4)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate partial reorganization energies for PT and ET
   ! (ERPT and ERET) along with the cross coupling term (ERX)
   ! (Note that these quantities are supposed to be independent
   ! on the position of the quantum particle).
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   kpoint = npnts/2
   if (nkr.gt.1) then
      kgpoint = nkr/2
   else
      kgpoint = 1
   endif

   erpt = half*tr(2,2,kpoint,kgpoint)
   eret = half*tr(3,3,kpoint,kgpoint)
   erx  = half*tr(2,3,kpoint,kgpoint)

   write(6,'(/1x,''Partial reorganization energies (kcal/mol):'')')
   write(6,'(10x,''Lambda_PT:  '',f15.6)') erpt
   write(6,'(10x,''Lambda_ET:  '',f15.6)') eret
   write(6,'(10x,''Coupling:   '',f15.6)') erx

   write(72,'(/1x,''Partial reorganization energies (kcal/mol):'')')
   write(72,'(10x,''Lambda_PT:  '',f15.6)') erpt
   write(72,'(10x,''Lambda_ET:  '',f15.6)') eret
   write(72,'(10x,''Coupling:   '',f15.6)') erx

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! calculate the eigenvalues of the truncated reorganization
   ! energy matrix
   !
   !   erpt   erx
   !                  ---->     er1, er2  (er1<er2)
   !   erx   eret
   !
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ersum = eret + erpt
   erdif = eret - erpt
   ersq = sqrt(erdif*erdif + 4.d0*erx*erx)
   er1 = half*(ersum - ersq)
   er2 = half*(ersum + ersq)

   write(6,'(/1x,''Eigenvalues of the truncated reorganization energy matrix (kcal/mol):'')')
   write(6,'(10x,''Lambda_1:   '',f15.6)') er1
   write(6,'(10x,''Lambda_2:   '',f15.6)') er2

   write(72,'(/1x,''Eigenvalues of the truncated reorganization energy matrix (kcal/mol):'')')
   write(72,'(10x,''Lambda_1:   '',f15.6)') er1
   write(72,'(10x,''Lambda_2:   '',f15.6)') er2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! calculate the transformation angle
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tanth = half*(erdif - ersq)/erx
   sqtan = sqrt(1.d0 + tanth*tanth)
   sin_theta = tanth/sqtan
   cos_theta = 1.d0/sqtan

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! initialize transformation matrix for solvent coordinates:
   ! ztmat = S*C^{tr} so that (z1,z2) = [ztmat]*(zp,ze)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   sq1 = sqrt(2.d0*f0*er1)
   sq2 = sqrt(2.d0*f0*er2)
   ztmat(1,1) =  cos_theta/sq1
   ztmat(2,2) =  cos_theta/sq2
   ztmat(1,2) =  sin_theta/sq1
   ztmat(2,1) = -sin_theta/sq2

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! initialize transformation matrix for gradient:
   ! gtmat = S^{-1}*C^{tr} so that grad(z1,z2) = [gtmat]*grad(zp,ze)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   gtmat(1,1) =  sq1*cos_theta
   gtmat(2,2) =  sq2*cos_theta
   gtmat(1,2) =  sq1*sin_theta
   gtmat(2,1) = -sq2*sin_theta

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! initialize shifts for new solvent coordinates
   ! (d1,d2) = [S*C^{tr}]*(tr1a1b,tr1a2a)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   tr1a1b = tr(1,2,kpoint,kgpoint)
   tr1a2a = tr(1,3,kpoint,kgpoint)
   delta_1 =  cos_theta*tr1a1b/sq1 + sin_theta*tr1a2a/sq1
   delta_2 = -sin_theta*tr1a1b/sq2 + cos_theta*tr1a2a/sq2

   !===================================================================
   ! calculate the dipole moments of the diabatic charge distributions
   ! with the origin placed at the center of 1a charge distribution.
   ! (depends only on the proton coordinate)
   !===================================================================

   nhsol = iptsol(2)

   do kp=1,npnts
      xyzsol(1,nhsol) = rlist(k)*bohr2a
      xyzsol(2,nhsol) = 0.d0
      xyzsol(3,nhsol) = 0.d0
      call calculate_dipole_moment_diab(kp)
   enddo

   write(6,'(/1x,"Dipole moments of the diabatic charge distributions (proton at the middle of the PT interface)")')
   write(6,'( 1x,"   State  ","     d_x         d_y          d_z          |d|  ")')

   write(6,'( 1x,"    1a    ",4f12.6)') dipole_moment_diab_x(1,1,npnts/2), &
                                      & dipole_moment_diab_y(1,1,npnts/2), &
                                      & dipole_moment_diab_z(1,1,npnts/2), &
   &sqrt(dipole_moment_diab_x(1,1,npnts/2)**2 + dipole_moment_diab_y(1,1,npnts/2)**2 + dipole_moment_diab_z(1,1,npnts/2)**2)

   write(6,'( 1x,"    1b    ",4f12.6)') dipole_moment_diab_x(2,2,npnts/2), &
                                      & dipole_moment_diab_y(2,2,npnts/2), &
                                      & dipole_moment_diab_z(2,2,npnts/2), &
   &sqrt(dipole_moment_diab_x(2,2,npnts/2)**2 + dipole_moment_diab_y(2,2,npnts/2)**2 + dipole_moment_diab_z(2,2,npnts/2)**2)

   write(6,'( 1x,"    2a    ",4f12.6)') dipole_moment_diab_x(3,3,npnts/2), &
                                      & dipole_moment_diab_y(3,3,npnts/2), &
                                      & dipole_moment_diab_z(3,3,npnts/2), &
   &sqrt(dipole_moment_diab_x(3,3,npnts/2)**2 + dipole_moment_diab_y(3,3,npnts/2)**2 + dipole_moment_diab_z(3,3,npnts/2)**2)

   write(6,'( 1x,"    2b    ",4f12.6)') dipole_moment_diab_x(4,4,npnts/2), &
                                      & dipole_moment_diab_y(4,4,npnts/2), &
                                      & dipole_moment_diab_z(4,4,npnts/2), &
   &sqrt(dipole_moment_diab_x(4,4,npnts/2)**2 + dipole_moment_diab_y(4,4,npnts/2)**2 + dipole_moment_diab_z(4,4,npnts/2)**2)


   write(72,'(/1x,"Dipole moments of the diabatic charge distributions (proton at the middle of the PT interface)")')
   write(72,'( 1x,"   State  ","     d_x         d_y          d_z          |d|  ")')

   write(72,'( 1x,"    1a    ",4f12.6)') dipole_moment_diab_x(1,1,npnts/2), &
                                      & dipole_moment_diab_y(1,1,npnts/2), &
                                      & dipole_moment_diab_z(1,1,npnts/2), &
   &sqrt(dipole_moment_diab_x(1,1,npnts/2)**2 + dipole_moment_diab_y(1,1,npnts/2)**2 + dipole_moment_diab_z(1,1,npnts/2)**2)

   write(72,'( 1x,"    1b    ",4f12.6)') dipole_moment_diab_x(2,2,npnts/2), &
                                      & dipole_moment_diab_y(2,2,npnts/2), &
                                      & dipole_moment_diab_z(2,2,npnts/2), &
   &sqrt(dipole_moment_diab_x(2,2,npnts/2)**2 + dipole_moment_diab_y(2,2,npnts/2)**2 + dipole_moment_diab_z(2,2,npnts/2)**2)

   write(72,'( 1x,"    2a    ",4f12.6)') dipole_moment_diab_x(3,3,npnts/2), &
                                      & dipole_moment_diab_y(3,3,npnts/2), &
                                      & dipole_moment_diab_z(3,3,npnts/2), &
   &sqrt(dipole_moment_diab_x(3,3,npnts/2)**2 + dipole_moment_diab_y(3,3,npnts/2)**2 + dipole_moment_diab_z(3,3,npnts/2)**2)

   write(72,'( 1x,"    2b    ",4f12.6)') dipole_moment_diab_x(4,4,npnts/2), &
                                      & dipole_moment_diab_y(4,4,npnts/2), &
                                      & dipole_moment_diab_z(4,4,npnts/2), &
   &sqrt(dipole_moment_diab_x(4,4,npnts/2)**2 + dipole_moment_diab_y(4,4,npnts/2)**2 + dipole_moment_diab_z(4,4,npnts/2)**2)
   
   write(72,*)
   close(72)


!===================================================================
contains

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! calculates the dipole moments of the diabatic charge distributions
   ! (the origin is placed at the center of charge)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   subroutine calculate_dipole_moment_diab(kp)

      integer, intent(in) :: kp

      integer :: i, k, nhsol
      real(8) :: cocx, cocy, cocz, q_total
      real(8) :: dix, diy, diz
      real(8), dimension(natsol) :: xxx, yyy, zzz

      !-- move the origin to the center of charge distribution
      !   for the first (1a) diabatic state (with the proton excluded)
      
      cocx = 0.d0
      cocy = 0.d0
      cocz = 0.d0
      q_total = charge
      
      nhsol = iptsol(2)

      !do k=1,natsol
      !   if (k.eq.nhsol) cycle
      !   cocx = cocx + chrsol(1,k)*xyzsol(1,k)
      !   cocy = cocy + chrsol(1,k)*xyzsol(2,k)
      !   cocz = cocz + chrsol(1,k)*xyzsol(3,k)
      !enddo
      
      !if (q_total.ne.0.d0) then
      !   cocx = cocx/q_total
      !   cocy = cocy/q_total
      !   cocz = cocz/q_total
      !endif
      
      do k=1,natsol
         xxx(k) = xyzsol(1,k) - cocx
         yyy(k) = xyzsol(2,k) - cocy
         zzz(k) = xyzsol(3,k) - cocz
      enddo

      do i=1,nelst

         dix = 0.d0
         diy = 0.d0
         diz = 0.d0
         
         do k=1,natsol
            dix = dix + chrsol(i,k)*xxx(k)
            diy = diy + chrsol(i,k)*yyy(k)
            diz = diz + chrsol(i,k)*zzz(k)
         enddo

         dipole_moment_diab_x(i,i,kp) = dix
         dipole_moment_diab_y(i,i,kp) = diy
         dipole_moment_diab_z(i,i,kp) = diz

      enddo

   end subroutine calculate_dipole_moment_diab


!===================================================================
end subroutine setmat
