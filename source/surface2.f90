subroutine surface2

!===================================================================C
!
!  Calculates the two-dimensional free energy surfaces
!  on rectangular grid in the space of solvent coordinates
!  ZP and ZE. Gating coordinate is quantized.
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
!              is diagonal (***NEW KEYWORD***)!
!
!  ZP=ZP1/ZP2/NZP - grid along ZP: NZP points from ZP1 to ZP2
!
!  ZE=ZE1/ZE2/NZE - grid along ZE: NZE points from ZE1 to ZE2
!
!  NSTATES= - number of states to plot
!            (in case of MGQUANT=3 it is a number of electron-proton
!             vibronic states, so the total number of states will be
!             NSTATES*NGAST)
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
!  $Date: 2011-06-30 20:37:09 $
!  $Revision: 5.4 $
!  $Log: not supported by cvs2svn $
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
   use cst
   use timers
   use keys
   use strings
   use solmat
   use control
   use quantum
   use feszz_2d

   implicit none

   character(1024) :: options
   character( 160) :: fname
   character(  20) :: zdim
   character(   5) :: mode

   logical :: adiab, diab2, diab4, dkl
   logical :: scale, transform, weights

   integer :: ikey, ielst, izp, ize, izp1, islash1, islash2, nzp, ize1, nze
   integer :: instates, nstates_tot, nf, nstates, ngatst, ioutput, lenf, ispa
   integer :: iweights, idklin, nlines, ndkl, nzdim, nz, npsiga
   integer :: i, ndabf, ievb

   real(8) :: pscal, escal, zp1, zp2, ze1, ze2, zpstep, zestep, z1, z2, dum1, dum2
   real(8) :: zp, zp0, ze, ze0, timep1, timep2, time0, time1, times

   integer, dimension(30,2) :: npair
   real(8), dimension(30)   :: dkle, dklp

   real(8), allocatable, dimension(:)       :: f, z, psiga
   real(8), allocatable, dimension(:,:,:,:) :: psiel, psipr
   real(8), allocatable, dimension(:,:,:)   :: enel, envib, envibg
   real(8), allocatable, dimension(:,:)     :: wght
   real(8), allocatable, dimension(:,:)   :: solv_coord_p, solv_coord_e

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Check whether gating quantization is specified
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (.not.gquant) then
      write(*,'(/1x,"*** (in SURFACE2): You MUST specify QUANTUM for GATING keyword ***"/)')
      stop
   endif

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
   ikey = index(keywrd,' SURF2(')

   if (ikey.eq.0) then
      write(*,'(/1x,"*** (in SURFACE2): You MUST specify options for SURF2 keyword ***"/)')
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
         write(*,'(/1x,"Warning in SURF2: SCALE keyword is incompatible with TRANSFORM and thus is ignored."/)')
         scale = .false.
         pscal = 1.d0
         escal = 1.d0
      endif
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Grid specification
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp = index(options,' ZP=')
   ize = index(options,' ZE=')

   if (izp.ne.0) then

      izp1 = izp+4
      islash1 = index(options(izp1:),'/')
      islash2 = index(options(izp1+islash1:),'/')
      zp1 = reada(options(izp1:izp1+islash1-2),1)
      zp2 = reada(options(izp1+islash1:izp1+islash1+islash2-1),1)
      nzp = reada(options,izp1+islash1+islash2)

   else

      write(*,'(/1x,"*** (in SURFACE2): You MUST specify ZP= option for SURF keyword ***"/)')
      stop

   endif

   if (ize.ne.0) then

      ize1 = ize+4
      islash1 = index(options(ize1:),'/')
      islash2 = index(options(ize1+islash1:),'/')
      ze1 = reada(options(ize1:ize1+islash1-2),1)
      ze2 = reada(options(ize1+islash1:ize1+islash1+islash2-1),1)
      nze = reada(options,ize1+islash1+islash2)

   else

      write(*,'(/1x,"*** (in SURFACE2): You MUST specify ZE= option for SURF keyword ***"/)')
      stop

   endif

   if (.not.transform) then

      write(6,'(/1x,"Surface grid specification:",/,&
      &" along ZP: ",I3," points from ",F7.3," to ",F7.3,2X,A,/,&
      &" along ZE: ",I3," points from ",F7.3," to ",F7.3,2X,A)')&
      &  NZP,ZP1,ZP2,ZDIM,NZE,ZE1,ZE2,ZDIM

   else

      write(6,'(/1x,"Surface grid specification:",/,&
      &" along z1: ",I3," points from ",F7.3," to ",F7.3,2X,A,/,&
      &" along z2: ",I3," points from ",F7.3," to ",F7.3,2X,A)')&
      &  NZP,ZP1,ZP2,ZDIM,NZE,ZE1,ZE2,ZDIM

   endif

   if (scale) write(*,'(/1x,"The solvent coordinates ZP and ZE are scaled",/,&
             &1X,"by the PT and ET partial reorganization",/,&
             &1X,"energies (see above), respectively.")')

   if (transform) write(*,'(/1x,"The transformed solvent coordinates z1 and z2 are used:",/,&
                     &          1x,"the self energy is diagonal in this frame.")')


   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Number of states to print out
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   instates = index(options,' NSTATES=')
   if (instates.ne.0) then
      nstates_tot = reada(options,instates+9)
   else
      nstates_tot = 10
   endif
   write(6,'(/1x,"Number of states to calculate: ",i3)') nstates_tot

   nf = nstates_tot

   if (mgquant.eq.3) then

      nstates = (nstates_tot-1)/ngast + 1
      ngatst = mod(nstates_tot-1,ngast) + 1
      if (nstates.gt.1) then
         write(6,'( 1x,"(",i4," gating states for each of the first",&
                          &i4," electron-proton vibronic states and",&
                          &i4," gating states for the",i4,a3," vibronic state)")')&
                          &ngast,nstates-1,ngatst,nstates,th(nstates)
      else
         write(6,'( 1x,"(",i4," gating states for the",i4,a3," vibronic state)")')&
         &ngatst,nstates,th(nstates)
      endif

   else

      nstates = nstates_tot

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Output files for free energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ioutput = index(options,' OUTPUT=')

   if (ioutput.eq.0) then

      fname = job(1:ljob)//'/states.grid'
      lenf = ljob + 12
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

      open(1,file=fname(1:lenf),status='new')

   elseif (diab2) then

      open(1,file=fname(1:lenf)//'.1',status='new')
      open(2,file=fname(1:lenf)//'.2',status='new')

   elseif (diab4) then

      open(1,file=fname(1:lenf)//'.1a',status='new')
      open(2,file=fname(1:lenf)//'.1b',status='new')
      open(3,file=fname(1:lenf)//'.2a',status='new')
      open(4,file=fname(1:lenf)//'.2b',status='new')

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

      write(6,'(/1x,"The nonadiabatic couplings with respect to solvent coordinates",/,&
      &"in case of quantum gating are not implemented yet")')

      dkl = .false.

      ! dkl = .true.

      ! ispa = index(options(idklin+7:),' ')
      ! fname = options(idklin+7:idklin+ispa+5)
      ! lenf = ispa - 1
      ! call locase(fname,lenf)

      ! write(6,'(/1x,"info about the nonadiabatic couplings",/,&
      ! &" with respect to solvent coordinates is read from",/,&
      ! &" the external file <",a,">")') fname(1:lenf)

      ! open(10,file=fname(1:lenf),status='old')

      ! read(10,*) nlines

      !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !! Read pairs of states for which non-adiabatic couplings
      !! should be calculated
      !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      ! if (nlines.eq.0) then
      !
      !    dkl = .false.
      !    close(10)
      !
      ! else
      !
      !    ndkl = 0
      !    do i=1,nlines
      !       read(10,*) (npair(i,k),k=1,2)
      !       if (npair(i,1).lt.npair(i,2)) then
      !          ndkl = ndkl + 1
      !          npair(ndkl,1) = npair(i,1)
      !          npair(ndkl,2) = npair(i,2)
      !       endif
      !    enddo
      !
      !    close(10)
      !
      !    write(6,'(/1x,"The nonadiabatic couplings between the following states",/,&
      !    &" will be calculated:",/,(1X,I2,"-",I2,","))') (npair(k,1),npair(k,2),k=1,ndkl)
      !
      !    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !    ! Open output file for nonadiabatic couplings
      !    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !
      !    idklout = index(options,' DKLOUT=')
      !
      !    if (idklout.ne.0) then
      !
      !       ispa = index(options(idklout+8:),' ')
      !       fname = options(idklout+8:idklout+ispa+6)
      !       lenf = ispa - 1
      !       call locase(fname,lenf)
      !       fname = job(1:ljob)//'/'//fname(1:lenf)
      !       lenf = lenf + ljob + 1
      !
      !    else
      !
      !       fname = job(1:ljob)//'/dklout'
      !       lenf = ljob + 7
      !
      !    endif
      !
      !    write(6,'(/1x,"The nonadiabatic couplings on the grid",/,&
      !    &" will be written to the external file <",a,">")') fname(1:lenf)
      !
      !    open(11,file=fname(1:lenf),status='new')
      !
      ! endif

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate free energies on the grid
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

   ! Allocate local arrays

   if (mgquant.eq.1.or.mgquant.eq.2) then
      nzdim = ielst*nprst*ngast
      nz = nzdim*nzdim
   elseif (mgquant.eq.3) then
      nzdim = ielst*nprst
      nz = nzdim*nzdim*npntsg
   endif

   allocate (z(nz))
   allocate (f(nf))

   if (weights) allocate (wght(nelst,nstates_tot))

   npsiga = ielst*nprst*ngast*npntsg

   allocate (psiel(ielst,npnts,npntsg,ielst))
   allocate (psipr(ielst,nprst,npnts,npntsg))
   allocate (psiga(npsiga))
   allocate (enel(ielst,npnts,npntsg))
   allocate (envib(ielst,nprst,npntsg))
   allocate (envibg(ielst,nprst,ngast))

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loops over the grid points
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

   write(6,'(/1x,"Calculating free energies at ",i6," grid points..."/)') nzp*nze
   time0 = second()

   !write(*,*)
   !icount = 0

   do izp=1,nzp

      z1 = zp1 + (izp-1)*zpstep
      z1 = z1*pscal

      do ize=1,nze

         z2 = ze1 + (ize-1)*zestep
         z2 = z2*escal

         zp0 = solv_coord_p(izp,ize)
         ze0 = solv_coord_e(izp,ize)

         !icount = icount + 1
         !write(*,'(i10,2f20.6)') icount,zp0,ze0

         timep1 = second()
         call feszz2(mode,1,zp0,ze0,nstates,nf,f,nz,z,ndabf,&
                    &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
         timep2 = second()
         write(1,'(2f10.3,2x,20g15.9)') z1,z2,(f(i),i=1,nstates_tot)
         write(*,'(f20.3,2x,2f10.3,2x,10g15.9)') timep2-timep1,z1,z2,(f(i),i=1,nstates_tot)

         if (weights) then
            call evbweig(mode,1,nstates_tot,ndabf,nz,z,ielst,psiel,psipr,npsiga,psiga,wght)
            write(21,'(2f10.3,2x,40g20.6)') z1,z2,((wght(ievb,i),ievb=1,4),i=1,nstates_tot)
         endif

         if (diab2) then

            timep1 = second()
            call feszz2(mode,2,zp0,ze0,nstates,nf,f,nz,z,ndabf,&
                       &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
            timep2 = second()
            write(2,'(2f10.3,2x,20g15.9)') z1,z2,(f(i),i=1,nstates_tot)
            write(*,'(f20.3,2x,2f10.3,2x,20g15.9)') timep2-timep1,z1,z2,(f(i),i=1,nstates_tot)

            if (weights) then
               call evbweig(mode,2,nstates_tot,ndabf,nz,z,ielst,psiel,psipr,npsiga,psiga,wght)
               write(22,'(2f10.3,2x,40g20.6)') z1,z2,((wght(ievb,i),ievb=1,4),i=1,nstates_tot)
            endif

         endif

         if (diab4) then
            call feszz2(mode,3,zp0,ze0,nstates,nf,f,nz,z,ndabf,&
                       &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
            write(3,'(2f10.3,2x,20g15.9)') z1,z2,(f(i),i=1,nstates_tot)
            call feszz2(mode,4,zp0,ze0,nstates,nf,f,nz,z,ndabf,&
                       &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
            write(4,'(2f10.3,2x,20g15.9)') z1,z2,(f(i),i=1,nstates_tot)
         endif

         if (dkl) then
            call zcoup   !(ndabf1,f1,z,psiel,psipr,ndkl,npair,dklp,dkle)
            !write(11,'(2f10.3,30g15.9)') z1,z2,(dklp(i),dkle(i),dsqrt(dklp(i)**2+dkle(i)**2),i=1,ndkl)
         endif

      enddo

      if (dkl) write(11,'(1x)')
      write(1,'(1x)')
      if (weights) write(21,'(1x)')
      if (diab2) then
         write(2,'(1x)')
         if (weights) write(22,'(1x)')
      endif
      if (diab4) then
         write(2,'(1x)')
         write(3,'(1x)')
         write(4,'(1x)')
      endif
       write(*,*)

   enddo

   time1 = second()
   times = time1 - time0
   write(6,'(1x,"Done. Time elapsed (sec): ",f12.3)') times

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

   deallocate (f,z,psiel,psipr,psiga)
   deallocate (enel,envib,envibg)
   deallocate (solv_coord_p, solv_coord_e)
   if (weights) deallocate (wght)

   return

end subroutine surface2

