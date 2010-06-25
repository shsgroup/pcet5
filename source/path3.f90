subroutine path3
!===================================================================C
!  Calculates the free energy along the specified path
!  in the space of solvent coordinates ZP and ZE
!  and gating coordinate (defined on the grid).
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
!  P1=ZP1/ZE1 - first point (ZP1,ZE1) on the straight line path
!
!  P2=ZP2/ZE2 - second point (ZP2,ZE2) on the straight line path
!
!  NPATH=N - number of points along the straight line path
!           (N/2...P1...N...P2...N/2)
!
!  CURVE= name of the file with points along the (curvilinear) path
!
!  NSTATES= - number of states to plot
!
!  OUTPUT= - name of the output files with surface data
!            (default filename "states.path")
!
!  WEIGHTS - calculate EVB weights for NSTATES states
!
!  DKLIN= - name of the input file with a specification
!           of nonadiabatic couplings between adiabatic
!           states (only for ADIAB)
!
!  DKLOUT= - name of the output file with nonadiabatic couplings
!--------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  path3.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  path3.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.3  2008/04/11 00:07:20  souda
!  length of string OPTIONS increased to 1024
!  to accomodate more options (not critical)
!
!  Revision 1.2  2007/11/06 22:20:01  souda
!  new mmgen gas phase potential o-h o systems added
!
!  Revision 1.1.1.1  2004/01/13 20:08:17  souda
!  Initial PCET-4.0 Release
!
!
!===================================================================C
   use pardim
   use control
   use strings
   use keys
   use cst
   use quantum
   use feszz_3d

   implicit none

   character(1024) :: options
   character(  40) :: fname
   character(   5) :: mode
   logical         :: adiab, diab2, diab4, dkl, weights

   integer :: ikey, ielst, icurve, ip1, ip2, ip11, islash, ip21
   integer :: inpath, npath, mpath, npoints, i, ispa, lenf
   integer :: instates, nstates, ioutput, iweights
   integer :: idklin, ndkl, nlines, idklout
   integer :: iz, k, kg, ievb, ndabf, nzdim
   
   real(8) :: zp1, ze1, zp2, ze2, stepzp, stepze, zpi, zei
   real(8) :: time0, zp, ze, rr, time1, times
   real(8) :: second

   real(8), allocatable, dimension(:)     :: f
   real(8), allocatable, dimension(:,:)   :: z
   real(8), allocatable, dimension(:,:,:) :: psiel, psipr
   real(8), allocatable, dimension(:,:)   :: enel, envib
   real(8), allocatable, dimension(:,:)   :: wght

   integer, dimension(30,2)  :: npair
   real(8), dimension(30)    :: dkle, dklp
   real(8), dimension(2,500) :: point

   adiab   = .false.
   diab2   = .false.
   diab4   = .false.
   dkl     = .false.
   weights = .false.

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' PATH3(')

   if (ikey.eq.0) then
      write(*,'(/1x,''*** (in PATH): You MUST specify options for PATH keyword ***''/)')
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

   nzdim = ielst*nprst

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Path specification
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   icurve = index(options,' CURVE=')
   ip1 = index(options,' P1=')
   ip2 = index(options,' P2=')

   if (ip1.ne.0.and.ip2.ne.0) then

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! First point (ZP1,ZE1) on the straight line path
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ip11 = ip1+4
      islash = index(options(ip11:),'/')
      zp1 = reada(options(ip11:ip11+islash-2),1)
      ze1 = reada(options,ip11+islash)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Second point (ZP2,ZE2) on the straight line path
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ip21 = ip2+4
      islash = index(options(ip21:),'/')
      zp2 = reada(options(ip21:ip21+islash-2),1)
      ze2 = reada(options,ip21+islash)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Number of points between P1 and P2
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      inpath = index(options,' NPATH=')
      if (inpath.ne.0) then
         npath = reada(options,inpath+7)
      else
         npath = 100
      endif

      stepzp = (zp2-zp1)/(npath-1)
      stepze = (ze2-ze1)/(npath-1)

      mpath = npath/2

      npoints = 0
      do i=-mpath+1,npath+mpath
         npoints = npoints + 1
         zpi = zp1+(i-1)*stepzp
         zei = ze1+(i-1)*stepze
         point(1,npoints) = zpi
         point(2,npoints) = zei
      enddo

      write(6,'(/1x,''Straight line path through the points'')')
      write(6,'( 1x,''('',F7.3,'','',F7.3,'') and ('',F7.3,'','',F7.3,'')'')') zp1,ze1,zp2,ze2
      write(6,'( 1x,''Number of points in between: '',I3)') npath
      write(6,'( 1x,''Total number of points: '',I3)') npoints

   elseif (icurve.ne.0) then

      ispa = index(options(icurve+7:),' ')
      fname = options(icurve+7:icurve+ispa+5)
      lenf = ispa - 1
      call locase(fname,lenf)

      write(6,'(/1x,''Path information is read from the file <'',A,''>'')') fname(1:lenf)

      open(10,file=fname(1:lenf),status='old')

      npoints = 0
      10 read(10,*,end=11) zpi,zei
      npoints = npoints + 1
      point(1,npoints) = zpi
      point(2,npoints) = zei
      goto 10
      11 write(6,'( 1x,''Total number of points: '',i3)') npoints

   else

      write(*,'(/1x,''*** (in PATH3): '',/,&
      &''You MUST specify either P1= or CURVE= option for PATH3 keyword ***''/)')
      stop

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of states to print out
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   instates = index(options,' NSTATES=')
   if (instates.ne.0) then
      nstates = reada(options,instates+9)
   else
      nstates = 2
   endif

   write(6,'(/1x,''Number of states to calculate: '',i3)') nstates

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output files for free energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' OUTPUT=')

   if (ioutput.eq.0) then

      fname = job(1:ljob)//'/states.path3'
      lenf = ljob + 13

   else

      ispa = index(options(ioutput+8:),' ')
      fname = options(ioutput+8:ioutput+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      fname = job(1:ljob)//'/'//fname(1:lenf)
      lenf = lenf + ljob + 1

   endif

   write(6,'(/1x,''Path data are written to the file(s) <'',a,''>'')') fname(1:lenf)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Open output files for free energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (adiab) then

      open(1,file=fname(1:lenf),status='new')
      write(1,'("#----------------------------------------------")')
      write(1,'("# Adiabatic free energy surfaces along the path")')
      write(1,'("# Columns:")')
      write(1,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(1,'("#----------------------------------------------")')

   elseif (diab2) then

      open(1,file=fname(1:lenf)//'.1',status='new')
      write(1,'("#---------------------------------------------")')
      write(1,'("# Reactant free energy surfaces along the path")')
      write(1,'("# Columns:")')
      write(1,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(1,'("#---------------------------------------------")')
      open(2,file=fname(1:lenf)//'.2',status='new')
      write(2,'("#--------------------------------------------")')
      write(2,'("# Product free energy surfaces along the path")')
      write(2,'("# Columns:")')
      write(2,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(2,'("#--------------------------------------------")')

   elseif (diab4) then

      open(1,file=fname(1:lenf)//'.1a',status='new')
      write(1,'("#--------------------------------------------")')
      write(1,'("# 1a free energy surfaces along the path")')
      write(1,'("# Columns:")')
      write(1,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(1,'("#--------------------------------------------")')
      open(2,file=fname(1:lenf)//'.1b',status='new')
      write(2,'("#--------------------------------------------")')
      write(2,'("# 1b free energy surfaces along the path")')
      write(2,'("# Columns:")')
      write(2,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(2,'("#--------------------------------------------")')
      open(3,file=fname(1:lenf)//'.2a',status='new')
      write(3,'("#--------------------------------------------")')
      write(3,'("# 2a free energy surfaces along the path")')
      write(3,'("# Columns:")')
      write(3,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(3,'("#--------------------------------------------")')
      open(4,file=fname(1:lenf)//'.2b',status='new')
      write(4,'("#--------------------------------------------")')
      write(4,'("# 2b free energy surfaces along the path")')
      write(4,'("# Columns:")')
      write(4,'("# i Z_P(i) Z_E(i) R  U_1(i,R)  ...  U_N(i,R)")')
      write(4,'("#--------------------------------------------")')

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output file for EVB weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   weights = index(options,' WEIGHTS').ne.0

   if (weights.and..not.diab4) then

      iweights = index(options,' WEIGHTS=')

      if (iweights.eq.0) then

         fname = job(1:ljob)//'/weights.path3'
         lenf = ljob + 14

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

      write(6,'(/1x,''Results of the wavefunction analysis along the path are written'',/,&
      &'' to the external file(s) <'',a,''>'')') fname(1:lenf)

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

      write(6,'(/1x,''Information about the nonadiabatic couplings'',/,&
      &'' with respect to solvent coordinates is read from'',/,&
      &'' the external file <'',a,''>'')') fname(1:lenf)

      open(10,file=fname(1:lenf),status='old')

      read(10,*) nlines

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Read pairs of states for which non-adiabatic couplings
      ! should be calculated
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (nlines.eq.0) then

         dkl = .false.
         close(10)
         write(6,'(/1x,''The nonadiabatic couplings with respect to solvent coordinates'',/,&
                     &'' will not be calculated (NDKL=0)'')')

      else

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

         write(6,'(/1x,''The nonadiabatic couplings between the following adiabatic states'',/,&
                      &'' will be calculated:'',/,(1x,i2,''-'',i2,'',''))')&
                      &(npair(k,1),npair(k,2),k=1,ndkl)

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

            fname = job(1:ljob)//'/dklout'
            lenf = ljob + 7

         endif

         write(6,'(/1x,''The nonadiabatic couplings on the grid'',/,&
                     &'' will be written to the external file <'',a,''>'')') fname(1:lenf)

         open(11,file=fname(1:lenf),status='new')

      endif

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate free energies along the path
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   allocate (f(nstates))
   allocate (z(nzdim,nzdim))
   allocate (psiel(ielst,npnts,ielst))
   allocate (psipr(ielst,nprst,npnts))
   allocate (enel(ielst,npnts))
   allocate (envib(ielst,nprst))
   if (weights) allocate (wght(nelst,nstates))

   write(6,'(/1x,''Calculating free energies at '',I5,'' points along a path'')') npoints
   time0 = second()

   do iz=1,npoints

      zp = point(1,iz)
      ze = point(2,iz)

      do kg=1,npntsg

         rr = glist(kg)*bohr2a

         call feszz3(mode,1,kg,zp,ze,nstates,f,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
         write(1,'(i10,3f10.3,2x,20g15.9)') iz,zp,ze,rr,(f(i),i=1,nstates)

         if (weights) then
            call evbwei(mode,1,nstates,ndabf,z,ielst,psiel,psipr,wght)
            write(21,'(3f10.3,2x,40g20.6)') zp,ze,rr,((wght(ievb,i),ievb=1,4),i=1,nstates)
         endif

         if (diab2) then
            call feszz3(mode,2,kg,zp,ze,nstates,f,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
            write(2,'(i10,3f10.3,2x,20g15.9)') iz,zp,ze,rr,(f(i),i=1,nstates)
            if (weights) then
               call evbwei(mode,2,nstates,ndabf,z,ielst,psiel,psipr,wght)
               write(22,'(3f10.3,2x,40g20.6)') zp,ze,rr,((wght(ievb,i),ievb=1,4),i=1,nstates)
            endif
         endif

         if (diab4) then
            call feszz3(mode,2,kg,zp,ze,nstates,f,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
            write(2,'(i10,3f10.3,2x,20g15.9)') iz,zp,ze,rr,(f(i),i=1,nstates)
            call feszz3(mode,3,kg,zp,ze,nstates,f,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
            write(3,'(i10,3f10.3,2x,20g15.9)') iz,zp,ze,rr,(f(i),i=1,nstates)
            call feszz3(mode,4,kg,zp,ze,nstates,f,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)
            write(4,'(i10,3f10.3,2x,20g15.9)') iz,zp,ze,rr,(f(i),i=1,nstates)
         endif

         if (dkl) then
             call zcoup         !(ndabf,f,z,psiel,psipr,ndkl,npair,dklp,dkle)
             !write(11,'(2f10.3,30g15.9)') zp,ze, (dklp(i),dkle(i),dsqrt(dklp(i)**2+dkle(i)**2),i=1,ndkl)
         endif

      enddo

      if (npntsg.gt.1) then
         write(1,*)
         if (diab2) write(2,*)
         if (diab4) then
            write(2,*)
            write(3,*)
            write(4,*)
         endif
      endif

   enddo

   time1 = second()
   times = time1 - time0
   write(6,'(1x,''Done. Time elapsed (sec): '',f12.3)') times

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

   deallocate (f,z,psiel,psipr)
   deallocate (enel,envib)
   if (weights) deallocate (wght)

   return

end subroutine path3

