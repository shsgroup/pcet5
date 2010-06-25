subroutine evbwei(mode,iset,nstates,ndabf,z,ielst,psiel,psipr,wght)
!============================================================================!
!  Calculates the EVB weights of the wavefunction
!
!  MODE - type of the state (ADIAB/DIAB2/DIAB4)
!  ISET - set of states (1 for ADIAB, 1/2 for DIAB2)
!  NSTATES - number of states to analyze (only one for METHOD=2/3)
!  NDABF - dimension of the total Hamiltonian
!  Z - total wavefunction coefficients (eigenvectors)
!  PSIEL - electronic wavefunctions
!  PSIPR - proton vibrational wavefunctions
!----------------------------------------------------------------------------!
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  evbwei.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  evbwei.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2003/12/19 16:49:42  souda
!  Initial PCET-4.0 Release
!
!
!============================================================================!
   use pardim
   use control
   use quantum

   implicit none

   character(5), intent(in)  :: mode
   integer,      intent(in)  :: iset,nstates,ndabf,ielst
   real(8),      intent(in)  :: z(ndabf,ndabf),&
                                psiel(ielst,npnts,ielst),&
                                psipr(ielst,nprst,npnts)
   real(8),      intent(out) :: wght(nelst,nstates)

   integer :: ielst0, istate, i, j, ishift, mu, k, kk
   real(8) :: zji, toti, ck, c2k, protwf

   real(8), dimension(ndabf) :: vect2

   wght = 0.d0  !array operation

   if (mode.eq.'ADIAB') then
      ielst0 = nelst
   elseif (mode.eq.'DIAB2') then
      ielst0 = 2
   elseif (mode.eq.'DIAB4') then
      ielst0 = 1
      do istate=1,nstates
         wght(iset,istate) = 1.d0
      enddo
      return
   endif

   if (ielst0.ne.ielst) then
      write(6,*) ' wrong input for evbwei: ielst=',ielst
      stop
   endif

   if (method.eq.1) then

      do istate=1,nstates

         vect2 = 0.d0

         do j=1,ndabf
            zji = z(j,istate)
            vect2(j) = zji*zji
         enddo

         do i=1,ielst

            ishift = i + 2*iset - 2
            toti = 0.d0

            do mu=1,nprst
               kk = (i-1)*nprst + mu
               toti = toti + vect2(kk)
            enddo

            wght(ishift,istate) = toti*100.d0

         enddo

      enddo

   else

      do i=1,ielst

         ishift = i + 2*iset - 2
         toti = 0.d0

         do k=1,npnts

            ck = psiel(1,k,i)
            c2k = ck*ck
            protwf = psipr(1,1,k)
            toti = toti + protwf*c2k*protwf

         enddo

         wght(ishift,1) = toti*100.d0

      enddo

   endif

   return

end subroutine evbwei
