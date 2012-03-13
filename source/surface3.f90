subroutine surface3

!===================================================================C
!
!  Calculates the three-dimensional free energy surfaces
!  on rectangular grid in the space of solvent coordinates
!  ZP and ZE and gating coordinate (grid points).
!  Output data are written to the specified external
!  file in the format compatible with GNUPLOT plotting
!  program.
!
!  OPTIONS:
!
!  ADIAB - adiabatic electron/proton vibrational free energy surfaces
!
!  DIAB2 - ET diabatic free energy surfaces (two output files)
!
!  DIAB4 - diabatic free energy surfaces (four output files)
!
!  SCALE - scale solvent coordinates with the reorganization energies
!
!  TRANSFORM - use rotated and scaled frame in which the self energy
!              is diagonal (***NEW KEYWORD***)
!
!  ZP=ZP1/ZP2/NZP - grid along ZP: NZP points from ZP1 to ZP2
!
!  ZE=ZE1/ZE2/NZE - grid along ZE: NZE points from ZE1 to ZE2
!
!  NSTATES= - number of states to plot
!
!  OUTPUT= - name of the output files with surface data
!            (default filename "states.grid")
!
!  WEIGHTS - calculate EVB weights for NSTATES states
!
!  DKLIN= - name of the input file with a specification
!           of nonadiabatic couplings between adiabatic
!           states (only for ADIAB)
!
!  DKLOUT= - name of the output file with nonadiabatic couplings
!
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-03-13 21:58:25 $
!  $Revision: 5.5 $
!  $Log: not supported by cvs2svn $
!  Revision 5.4  2011/06/30 20:37:09  souda
!  allowing longer job directory names
!
!  Revision 5.3  2011/02/20 00:58:11  souda
!  Major additions/modifications:
!  (1) precalculation of the proton vibrational basis functions for METHOD=1
!  (2) Franck-Condon initial excitation added to DYNAMICS3
!  (3) addition of the module timers: module_timers.f90 (changes in Makefile!)
!
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!===================================================================C

   use pardim
   use keys
   use cst
   use timers
   use strings
   use solmat
   use control
   use quantum
   use feszz_3d
   use geogas, only: iptgas, xyzgas

   implicit none

   character(1024) :: options
   character( 160) :: fname
   character(  20) :: zdim
   character(   5) :: mode

   logical :: adiab, diab2, diab4, dkl
   logical :: scale, transform, weights, gatedim

   integer :: ikey, ielst, izp, ize, izp1, islash1, islash2, nzp
   integer :: ize1, nze, nzdim, instates, nstates, ioutput, lenf, ispa
   integer :: iweights, idklin, nlines, ndkl, idklout
   integer :: icount, i, j, k, np, kg, ndabf, ievb

   real(8) :: pscal, escal, zp1, zp2, ze1, ze2, zpstep, zestep
   real(8) :: rkg, zp, zp0, ze, ze0, z1, z2, dum1, dum2, dz1, dz2
   real(8) :: time0, time1, times

   integer, allocatable, dimension(:,:)   :: npair
   real(8), allocatable, dimension(:,:)   :: dkle, dklp

   real(8), allocatable, dimension(:)     :: fe
   real(8), allocatable, dimension(:,:)   :: z, enel, envib
   real(8), allocatable, dimension(:,:,:) :: psiel, psipr
   real(8), allocatable, dimension(:,:)   :: wght
   real(8), allocatable, dimension(:,:)   :: solv_coord_p, solv_coord_e

   !=== DEBUG variables =====================================
   !real(8) :: zinc, splus, sminus
   !real(8), allocatable, dimension(:,:) :: dkle_num, dklp_num
   !real(8), allocatable, dimension(:,:) :: zplus, zminus
   !=== END DEBUG variables =================================

   adiab   = .false.
   diab2   = .false.
   diab4   = .false.
   dkl     = .false.
   scale   = .false.
   weights = .false.
   transform = .false.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' SURF3(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in SURFACE3): You MUST specify options for SURF keyword ***"/)')
      stop
   else
      call getopt(keywrd,ikey+7,options)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Surface type (ADIAB, DIAB2, DIAB4)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' ADIAB').ne.0) then
      mode = 'ADIAB'
      adiab=.true.
      ielst = nelst
   elseif (index(options,' DIAB2').ne.0) then
      mode = 'DIAB2'
      diab2=.true.
      ielst = 2
   elseif (index(options,' DIAB4').ne.0) then
      mode = 'DIAB4'
      diab4=.true.
      ielst = 1
   else
      mode = 'ADIAB'
      adiab=.true.
      ielst = nelst
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Scaling of solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' SCALE').ne.0) then
      scale = .true.
      zdim = ' '
      pscal = 1.d0/erpt
      escal = 1.d0/eret
   else
      zdim = 'kcal/mol'
      pscal = 1.d0
      escal = 1.d0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Using transformed solvent coordinates (TRANSFORM)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (index(options,' TRANSFORM').ne.0) then
      transform = .true.
      zdim = '(kcal/mol)^(1/2)'
      if (scale) then
         write(*,'(/1x,"Warning in SURF3: SCALE keyword is incompatible with TRANSFORM and thus is ignored."/)')
         scale = .false.
         pscal = 1.d0
         escal = 1.d0
      endif
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Grid specification: solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp = index(options,' ZP=')
   ize = index(options,' ZE=')

   if (izp.ne.0) then

      izp1 = izp + 4
      islash1 = index(options(izp1:),'/')
      islash2 = index(options(izp1+islash1:),'/')
      zp1 = reada(options(izp1:izp1+islash1-2),1)
      zp2 = reada(options(izp1+islash1:izp1+islash1+islash2-1),1)
      nzp = reada(options,izp1+islash1+islash2)

   else

      write(*,'(/1x,"*** (in SURFACE3): You MUST specify ZP= option for SURF keyword ***"/)')
      stop

   endif

   if (ize.ne.0) then

      ize1 = ize + 4
      islash1 = index(options(ize1:),'/')
      islash2 = index(options(ize1+islash1:),'/')
      ze1 = reada(options(ize1:ize1+islash1-2),1)
      ze2 = reada(options(ize1+islash1:ize1+islash1+islash2-1),1)
      nze = reada(options,ize1+islash1+islash2)

   else

      write(*,'(/1x,"*** (in SURFACE3): You MUST specify ZE= option for SURF keyword ***"/)')
      stop

   endif

   if (.not.transform) then

      write(6,'(/1x,"Surface grid specification:",/,&
      &" along ZP: ",I3," points from ",F10.3," to ",F10.3,2X,A,/,&
      &" along ZE: ",I3," points from ",F10.3," to ",F10.3,2X,A)')&
      &  NZP,ZP1,ZP2,ZDIM,NZE,ZE1,ZE2,ZDIM

   else

      write(6,'(/1x,"Surface grid specification:",/,&
      &" along z1: ",I3," points from ",F10.3," to ",F10.3,2X,A,/,&
      &" along z2: ",I3," points from ",F10.3," to ",F10.3,2X,A)')&
      &  NZP,ZP1,ZP2,ZDIM,NZE,ZE1,ZE2,ZDIM

   endif

   if (scale) write(*,'(/1x,"The solvent coordinates ZP and ZE are scaled",/,&
              &          1x,"by the PT and ET partial reorganization",/,&
              &          1x,"energies (see above), respectively.")')

   if (transform) write(*,'(/1x,"The transformed solvent coordinates z1 and z2 are used:",/,&
                  &          1x,"the self energy is diagonal in this frame.")')

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Grid specification: gating coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (npntsg.gt.1) then
      gatedim = .true.
      write(6,'(/1x,"Gating grid specification:",/,&
      &" gating: ",i3," points from ",f10.3," to ",f10.3,2x,"A")')&
      &npntsg,glist(1)*bohr2a,glist(npntsg)*bohr2a
   else
      gatedim = .false.
      write(6,'(/1x,"Gating coordinate value:",f10.3,2x,"A")') abs(xyzgas(1,iptgas(3)) - xyzgas(1,iptgas(1)))
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of states to print out
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   instates = index(options,' NSTATES=')
   if (instates.ne.0) then
      nstates = reada(options,instates+9)
   else
      nstates = 10
   endif
   write(6,'(/1x,"Number of states to calculate: ",i3)') nstates

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output files for free energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' OUTPUT=')

   if (ioutput.eq.0) then

      fname = job(1:ljob)//'/surface3'
      lenf = ljob + 9
      write(6,'(/1x,"surface data are written to the file(s) <",a,">")') fname(1:lenf)

   else

      ispa = index(options(ioutput+8:),' ')
      fname = options(ioutput+8:ioutput+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      fname = job(1:ljob)//'/'//fname(1:lenf)
      lenf = lenf + ljob + 1
      write(6,'(/1x,"surface data are written to the file <",a,">")') fname(1:lenf)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Open output files
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (adiab) then

      open(1,file=fname(1:lenf)//'.dat',status='new')

   elseif (diab2) then

      open(1,file=fname(1:lenf)//'.1.dat',status='new')
      open(2,file=fname(1:lenf)//'.2.dat',status='new')

   elseif (diab4) then

      open(1,file=fname(1:lenf)//'.1a.dat',status='new')
      open(2,file=fname(1:lenf)//'.1b.dat',status='new')
      open(3,file=fname(1:lenf)//'.2a.dat',status='new')
      open(4,file=fname(1:lenf)//'.2b.dat',status='new')

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output file for EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   weights = index(options,' WEIGHTS').ne.0

   if (weights.and..not.diab4) then

      iweights = index(options,' WEIGHTS=')

      if (iweights.eq.0) then

         fname = job(1:ljob)//'/weights.grid'
         lenf = ljob + 13

      else

         ispa = index(options(iweights+9:),' ')
         fname = options(iweights+9:iweights+ispa+7)
         lenf = ispa - 1
         call locase(fname,lenf)
         fname = job(1:ljob)//'/'//fname(1:lenf)
         lenf = lenf + ljob + 1

      endif

      if (adiab) then
         open(21,file=fname(1:lenf),status='new')
      elseif (diab2) then
         open(21,file=fname(1:lenf)//'.1',status='new')
         open(22,file=fname(1:lenf)//'.2',status='new')
      endif

      write(6,'(/1x,"Results of the wavefunction analysis are written",/,&
      &" to the external file(s) <",a,">")') fname(1:lenf)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Parse options related to non-adiabatic couplings
   ! with respect to solvent coordinates.
   ! Note that this is meaningful only for adiabatic
   ! states (MODE='ADIAB')
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   idklin = index(options,' DKLIN=')

   if (adiab.and.idklin.ne.0) then

      dkl = .true.

      ispa = index(options(idklin+7:),' ')
      fname = options(idklin+7:idklin+ispa+5)
      lenf = ispa - 1
      call locase(fname,lenf)

      write(6,'(/1x,"Info about the nonadiabatic couplings",/,&
      &" with respect to solvent coordinates is read from",/,&
      &" the external file <",a,">")') fname(1:lenf)

      open(10,file=fname(1:lenf),status='old')

      read(10,*) nlines

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read pairs of states for which non-adiabatic couplings
      ! should be calculated
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (nlines.eq.0) then

         dkl = .false.
         close(10)

      else

         allocate(npair(nlines,2))

         ndkl = 0
         do i=1,nlines
            read(10,*) (npair(i,k),k=1,2)
            if (npair(i,1).lt.npair(i,2)) then
               ndkl = ndkl + 1
               npair(ndkl,1) = npair(i,1)
               npair(ndkl,2) = npair(i,2)
            endif
         enddo

         close(10)

         write(6,'(/1x,"The nonadiabatic couplings between the following adiabatic states",/,&
         &" will be calculated:",/,(1x,i2,"-",i2,","))') (npair(k,1),npair(k,2),k=1,ndkl)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Open output file for nonadiabatic couplings
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         idklout = index(options,' DKLOUT=')

         if (idklout.ne.0) then

            ispa = index(options(idklout+8:),' ')
            fname = options(idklout+8:idklout+ispa+6)
            lenf = ispa - 1
            call locase(fname,lenf)
            fname = job(1:ljob)//'/'//fname(1:lenf)
            lenf = lenf + ljob + 1

         else

            fname = job(1:ljob)//'/nacouplings.dat'
            lenf = ljob + 16

         endif

         write(6,'(/1x,"The nonadiabatic couplings on the grid",/,&
         &" will be written to the external file <",a,">")') fname(1:lenf)

         open(11,file=fname(1:lenf),status='new')

      endif

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate free energies on the grid of solvent coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (nzp.gt.1) then
      zpstep = (zp2-zp1)/(nzp-1)
   else
      zpstep = 0.d0
   endif

   if (nze.gt.1) then
      zestep = (ze2-ze1)/(nze-1)
   else
      zestep = 0.d0
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over the grid points
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   np = nzp*nze*npntsg
   write(6,'(/1x,"Calculating free energies at ",i6," grid points...")') np

   !-- generate two-dimensional array for solvent coordinates values

   allocate (solv_coord_p(nzp,nze))
   allocate (solv_coord_e(nzp,nze))

   do izp=1,nzp
      z1 = zp1 + (izp-1)*zpstep
      z1 = z1*pscal
      do ize=1,nze
         z2 = ze1 + (ize-1)*zestep
         z2 = z2*escal
         if (transform) then
            call z1z2_to_zpze(z1,z2,dum1,dum2)
            solv_coord_p(izp,ize) = dum1
            solv_coord_e(izp,ize) = dum2
         else
            solv_coord_p(izp,ize) = z1
            solv_coord_e(izp,ize) = z2
         endif
      enddo
   enddo

   time0 = second()

   ! Allocate arrays

   nzdim = nprst*ielst
   allocate (fe(nstates))
   allocate (z(nzdim,nzdim))
   !allocate (zplus(nzdim,nzdim), zminus(nzdim,nzdim))   !=== DEBUG
   allocate (psiel(ielst,npnts,ielst))
   allocate (psipr(ielst,nprst,npnts))
   allocate (enel(ielst,npnts))
   allocate (envib(ielst,nprst))
   if (weights) allocate (wght(nelst,nstates))
   if (dkl) then
      allocate (dklp(nstates,nstates),dkle(nstates,nstates))
      !allocate (dklp_num(nstates,nstates),dkle_num(nstates,nstates))  !=== DEBUG
   endif

   !-> Start of the loop over the gating grid points

   icount = 0
   write(*,*)
   write(*,'(49x,$)')

   loop_over_kg: do kg=1,npntsg

      rkg = glist(kg)*bohr2a

      loop_over_zp: do izp=1,nzp

         z1 = zp1 + (izp-1)*zpstep
         z1 = z1*pscal

         loop_over_ze: do ize=1,nze

            z2 = ze1 + (ize-1)*zestep
            z2 = z2*escal

            zp0 = solv_coord_p(izp,ize)
            ze0 = solv_coord_e(izp,ize)

            icount = icount + 1
            if (gatedim) then
               write(*,'(49("\b"),i10,":(",3f13.6,")",$)') icount,rkg,z1,z2
            else
               write(*,'(49("\b"),i10,":(",2f13.6,")",12x,$)') icount,z1,z2
            endif

            call feszz3(mode,1,kg,zp0,ze0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
            if (gatedim) then
               write(1,'(3f10.3,2x,40(1x,g16.9))') z1,z2,rkg,(fe(i),i=1,nstates)
            else
               write(1,'(2f10.3,2x,40(1x,g16.9))') z1,z2,(fe(i),i=1,nstates)     !, fe(2)-fe(1)  !===DEBUG
            endif

            if (weights) then
               call evbwei(mode,1,nstates,ndabf,z,ielst,psiel,psipr,wght)
               if (gatedim) then
                  write(21,'(3f10.3,2x,40g20.6)') z1,z2,rkg,((wght(ievb,i),ievb=1,4),i=1,nstates)
               else
                  write(21,'(2f10.3,2x,40g20.6)') z1,z2,((wght(ievb,i),ievb=1,4),i=1,nstates)
               endif
            endif

            if (dkl) then
               call coupzz3(mode,1,nstates,fe,nzdim,z,ielst,dklp,dkle)

               !=== DEBUG =========================================================
               !-- check numerical couplings
               !-------------------------------------------------------------------
               !zinc = 1.d-4
               !call feszz3(mode,1,kg,zp0+zinc,ze0,nstates,fe,nzdim,zplus,ndabf,ielst,enel,envib,psiel,psipr)
               !call feszz3(mode,1,kg,zp0-zinc,ze0,nstates,fe,nzdim,zminus,ndabf,ielst,enel,envib,psiel,psipr)
               !do i=1,nstates
               !   do j=i+1,nstates
               !      splus = 0.d0
               !      sminus = 0.d0
               !      do k=1,nzdim
               !         splus  = splus  + z(k,i)*zplus(k,j)
               !         sminus = sminus + z(k,i)*zminus(k,j)
               !      enddo
               !      dklp_num(i,j) = (splus - sminus)/(2.d0*zinc)
               !      dklp_num(j,i) = -dklp_num(i,j)
               !   enddo
               !enddo
               !call feszz3(mode,1,kg,zp0,ze0+zinc,nstates,fe,nzdim,zplus,ndabf,ielst,enel,envib,psiel,psipr)
               !call feszz3(mode,1,kg,zp0,ze0-zinc,nstates,fe,nzdim,zminus,ndabf,ielst,enel,envib,psiel,psipr)
               !do i=1,nstates
               !   do j=i+1,nstates
               !      splus = 0.d0
               !      sminus = 0.d0
               !      do k=1,nzdim
               !         splus  = splus  + z(k,i)*zplus(k,j)
               !         sminus = sminus + z(k,i)*zminus(k,j)
               !      enddo
               !      dkle_num(i,j) = (splus - sminus)/(2.d0*zinc)
               !      dkle_num(j,i) = -dkle_num(i,j)
               !   enddo
               !enddo
               !=== END DEBUG =====================================================

               if (transform) then
                  do i=1,nstates
                     do j=i+1,nstates
                        call gpge_to_g1g2(dklp(i,j),dkle(i,j),dz1,dz2)
                        dklp(i,j) =  dz1
                        dkle(i,j) =  dz2
                        dklp(j,i) = -dz1
                        dkle(j,i) = -dz2
                     enddo
                  enddo
               endif
               if (gatedim) then
                  write(11,'(2f10.3,40(1x,g16.9))') z1,z2,rkg,&
                  &(dklp(npair(i,1),npair(i,2)),dkle(npair(i,1),npair(i,2)),&
                  &sqrt(dklp(npair(i,1),npair(i,2))**2 + dkle(npair(i,1),npair(i,2))**2),i=1,ndkl)
               else
                  write(11,'(2f10.3,40(1x,g16.9))') z1,z2,&
                  &(dklp(npair(i,1),npair(i,2)),dkle(npair(i,1),npair(i,2)),&
                  &sqrt(dklp(npair(i,1),npair(i,2))**2 + dkle(npair(i,1),npair(i,2))**2),i=1,ndkl)
                  !&(dklp_num(npair(i,1),npair(i,2)),dkle_num(npair(i,1),npair(i,2)),&                        !===DEBUG
                  !&sqrt(dklp_num(npair(i,1),npair(i,2))**2 + dkle_num(npair(i,1),npair(i,2))**2),i=1,ndkl)   !===DEBUG
               endif
            endif

            if (diab2) then
               call feszz3(mode,2,kg,zp0,ze0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
               if (gatedim) then
                  write(2,'(3f10.3,2x,40g16.9)') z1,z2,rkg,(fe(i),i=1,nstates)
               else
                  write(2,'(2f10.3,2x,40g16.9)') z1,z2,(fe(i),i=1,nstates)
               endif
               if (weights) then
                  call evbwei(mode,2,nstates,ndabf,z,ielst,psiel,psipr,wght)
                  if (gatedim) then
                     write(22,'(3f10.3,2x,40g20.6)') z1,z2,rkg,((wght(ievb,i),ievb=1,4),i=1,nstates)
                  else
                     write(22,'(2f10.3,2x,40g20.6)') z1,z2,((wght(ievb,i),ievb=1,4),i=1,nstates)
                  endif
               endif
            endif

            if (diab4) then

               call feszz3(mode,2,kg,zp0,ze0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
               if (gatedim) then
                  write(2,'(3f10.3,2x,40g16.9)') z1,z2,rkg,(fe(i),i=1,nstates)
               else
                  write(2,'(2f10.3,2x,40g16.9)') z1,z2,(fe(i),i=1,nstates)
               endif

               call feszz3(mode,3,kg,zp0,ze0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
               if (gatedim) then
                  write(3,'(3f10.3,2x,40g16.9)') z1,z2,rkg,(fe(i),i=1,nstates)
               else
                  write(3,'(2f10.3,2x,40g16.9)') z1,z2,(fe(i),i=1,nstates)
               endif

               call feszz3(mode,4,kg,zp0,ze0,nstates,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
               if (gatedim) then
                  write(4,'(3f10.3,2x,40g16.9)') z1,z2,rkg,(fe(i),i=1,nstates)
               else
                  write(4,'(2f10.3,2x,40g16.9)') z1,z2,(fe(i),i=1,nstates)
               endif
            endif

         enddo loop_over_ze

         if (dkl) write(11,'(1x)')
         write(1,'(1x)')
         if (weights) write(21,'(1x)')
         if (diab2) then
            write(2,'(1x)')
            if (weights) write(22,'(1x)')
         endif
         if (diab4) then
            write(2,*)
            write(3,*)
            write(4,*)
         endif

      enddo loop_over_zp

   enddo loop_over_kg

   write(*,*)
   write(*,*)

   time1 = second()
   times = time1 - time0
   write(6,'(1x,"Done. Time elapsed   (sec): ",f12.3)') times
   write(6,'(1x,"      Time per point (sec): ",f12.3)') times/np

   if (dkl) close(11)
   close(1)
   if (weights) close(21)
   if (diab2) then
      close(2)
      if (weights) close(22)
   endif
   if (diab4) then
      close(2)
      close(3)
      close(4)
   endif

   ! Deallocate arrays

   deallocate (fe)
   deallocate (z)
   !deallocate (zplus, zminus)   !=== DEBUG
   deallocate (psiel,psipr)
   deallocate (enel,envib)
   deallocate (solv_coord_p, solv_coord_e)
   if (weights) deallocate (wght)
   if (dkl) then
      deallocate(dklp, dkle)
      !deallocate(dklp_num, dkle_num)  !=== DEBUG
      deallocate(npair)
   endif

   return

end subroutine surface3

