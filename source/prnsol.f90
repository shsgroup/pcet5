subroutine prnsol

!===================================================================C
!
!  Print out solvated energy profiles along the proton coordinate
!  at specified values of solvent coordinates ZP and ZE
!  (if PTSOL keyword is specified)
!
!  June 16, 2000: Alexander Soudackov, University of Notre Dame
!
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!
!===================================================================C

   use pardim
   use cst
   use keys
   use strings
   use quantum
      
   implicit none
   character(1024) :: options
   character(40)   :: fname
   integer         :: ikey, izp, ize, iptsol, ispa, lenf, k, kg, i
   real*8          :: zp, ze
   real*8, dimension(nelst)       :: hd4, hd2, haa
   real*8, dimension(nelst,nelst) :: cu

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract options
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ikey = index(keywrd,' PTSOL(')

   if (ikey.eq.0) then

      write(*,'(/1X,''*** (in PRNSOL): You MUST specify PTSOL keyword with options ***''/)')
      stop

   else

      call getopt(keywrd,ikey+7,options)

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Extract the values of ZP and ZE coordinates
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   izp = index(options,' ZP=')
   ize = index(options,' ZE=')

   if (izp.eq.0.or.ize.eq.0) then
      write(*,'(/1X,''*** (in PRNSOL): You MUST specify ZP= and ZE= options ***''/)')
      stop
   else
      zp = reada(options,izp+4)
      ze = reada(options,ize+4)
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Open the external file for writing the profiles
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   iptsol = index(options,' OUTPUT=')

   IF (IPTSOL.EQ.0) THEN

      WRITE(6,'(/1X,''Solvated profiles along H coordinate '',&
      &''for zp='',f7.3,'' and ze='',f7.3,/,1x,&
      &''are written to the external file <'',a,''>'')')&
      &zp,ze,job(1:ljob)//'/'//'ptsol.dat'
      write(6,'(1x,''zp='',f7.3,'' ze='',f7.3)') zp,ze
      open(1,file=job(1:ljob)//'/'//'ptsol.dat',status='new')

   else

      ispa = index(options(iptsol+8:),' ')
      fname = options(iptsol+8:iptsol+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1X,''Solvated profiles along H coordinate '',&
      &''for ZP='',F7.3,'' and ZE='',F7.3,/,1X,&
      &''are written to the external file <'',a,''>'')')&
      &zp,ze,job(1:ljob)//'/'//fname(1:lenf)
      open(1,file=job(1:ljob)//'/'//fname(1:lenf),status='new')

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Caculate the energies of the diabatic and adiabatic states
   ! on the grid along the proton coordinate at ZP and ZE
   !
   ! The data structure is as follows:
   !
   ! col 1:       RLIST(k)  - proton coordinate at k-th grid point (A)
   ! col 2:       GLIST(k)  - gating coordinate at k-th grid point (A)
   ! col 3-6:     HD4(i,k)  - diabatic energies    (DIAB4) (kcal/mol)
   ! col 7-10:    HD2(i,k)  - ET diabatic energies (DIAB2) (kcal/mol)
   ! col 11-14:   HAA(i,k)  - adiabatic energies   (ADIAB) (kcal/mol)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(1,'("#==================================================================")')
   write(1,'("# col 1:     RLIST(k)  - proton coordinate at k-th grid point (A)")')
   write(1,'("# col 2:     GLIST(k)  - gating coordinate at k-th grid point (A)")')
   write(1,'("# col 3-6:   HD4(i,k)  - diabatic energies    (DIAB4) (kcal/mol)")')
   write(1,'("# col 7-10:  HD2(i,k)  - ET diabatic energies (DIAB2) (kcal/mol)")')
   write(1,'("# col 11-14: HAA(i,k)  - adiabatic energies   (ADIAB) (kcal/mol)")')
   write(1,'("#==================================================================")')

   do kg=1,npntsg
      do k=1,npnts
         call usol('DIAB4',k,kg,zp,ze,hd4,cu)
         call usol('DIAB2',k,kg,zp,ze,hd2,cu)
         call usol('ADIAB',k,kg,zp,ze,haa,cu)
         write(1,'(20g20.10)') rlist(k)*bohr2a,glist(kg)*bohr2a,&
         &(hd4(i),i=1,nelst),(hd2(i),i=1,nelst),(haa(i),i=1,nelst)
      enddo
      write(1,*)
   enddo

   close(1)
   return

end subroutine prnsol

