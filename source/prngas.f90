subroutine prngas

!===================================================================C
!
!  Print out the gas-phase and electronically solvated
!  energy profiles along the proton coordinate
!  (if PTGAS keyword is specified)
!
!  June 9, 2000: Alexander Soudackov, University of Notre Dame
!
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!===================================================================C

   use pardim
   use cst
   use keys
   use strings
   use quantum

   implicit none
   character(40) :: fname
   integer       :: iptgas, ispa, lenf, k, kg, i
   real*8, dimension(nelst) :: hd, ha12, ha, hdes, ha12es, haes

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Open the external file
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   iptgas = index(keywrd,' PTGAS=')

   if (iptgas.eq.0) then

      write(6,'(/1X,''Gas phase energy surfaces are written to the external file <'',A,''>'')')&
      &job(1:ljob)//'/'//'ptgas'
      open(1,file=job(1:ljob)//'/'//'ptgas.dat',status='new')

   else

      ispa = index(keywrd(iptgas+7:),' ')
      fname = keywrd(iptgas+7:iptgas+ispa+6)
      lenf = ispa - 1
      call locase(fname,lenf)
      write(6,'(/1X,''Gas phase energy surfaces are written to the external file <'',A,''>'')')&
      &job(1:ljob)//'/'//fname(1:lenf)
      open(1,file=job(1:ljob)//'/'//fname(1:lenf),status='new')

   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Caculate the energies of the diabatic and adiabatic states
   ! on the grid along the proton coordinate for the gas phase
   ! and including the electronic solvation.
   !
   ! The data structure is as follows:
   !
   ! col 1:     RLIST(k)    - proton coordinate at k-th grid point
   ! col 2:     GLIST(k)    - gating coordinate at k-th grid point
   ! col 3-6:   HD(i,k)     - gas-phase diabatic energies
   ! col 7-10:  HA12(i,k)   - gas-phase ET-diabatic/PT-adiabatic energies
   ! col 11-14: HA(i,k)     - gas-phase adiabatic energies
   ! col 15-18: HDES(i,k)   - electronically solvated diabatic energies
   ! col 19-22: HA12ES(i,k) - electronically solvated ET-diabatic/PT-adiabatic energies
   ! col 23-26: HAES(i,k)   - electronically solvated adiabatic energies
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   write(1,'("#====================================================================================")')
   write(1,'("# col 1:     RLIST(k)    - proton coordinate at k-th grid point")')
   write(1,'("# col 2:     GLIST(k)    - gating coordinate at k-th grid point")')
   write(1,'("# col 3-6:   HD(i,k)     - gas-phase diabatic energies")')
   write(1,'("# col 7-10:  HA12(i,k)   - gas-phase ET-diabatic/PT-adiabatic energies")')
   write(1,'("# col 11-14: HA(i,k)     - gas-phase adiabatic energies")')
   write(1,'("# col 15-18: HDES(i,k)   - electronically solvated diabatic energies")')
   write(1,'("# col 19-22: HA12ES(i,k) - electronically solvated ET-diabatic/PT-adiabatic energies")')
   write(1,'("# col 23-26: HAES(i,k)   - electronically solvated adiabatic energies")')
   write(1,'("#====================================================================================")')

   do kg=1,npntsg
      do k=1,npnts
         call ugas(k, kg, hd, ha12, ha, hdes, ha12es, haes)
         write(1,'(2f10.6,2x,30g20.10)') rlist(k)*bohr2a,glist(kg)*bohr2a,&
        &(hd(i),i=1,nelst),(ha12(i),i=1,nelst),(ha(i),i=1,nelst),&
	&(hdes(i),i=1,nelst),(ha12es(i),i=1,nelst),(haes(i),i=1,nelst)
      enddo
      write(1,*)
   enddo

   close(1)

   return

end subroutine prngas
