subroutine evbweig(mode,iset,nstates,ndabf,nz,z,&
                   ielst,psiel,psipr,npsiga,psiga,wght)
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
!  PSIGA - gating vibrational wavefunctions
!-------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  evbweig.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  evbweig.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2003/12/19 16:51:27  souda
!  Initial PCET-4.0 Release
!
!
!============================================================================!
   use pardim
   use control
   use quantum

   implicit none

   character(5), intent(in)  :: mode
   integer,      intent(in)  :: iset,nstates,ndabf,nz,ielst,npsiga
   real(8),      intent(in)  :: z(nz),&
                                psiel(ielst,npnts,npntsg,ielst),&
                                psipr(ielst,nprst,npnts,npntsg),&
                                psiga(npsiga)
   real(8),      intent(out) :: wght(nelst,nstates)

   integer :: ielst_, istate, j, i, ishift, mu, m, kk
   integer :: nvibst, ivibst, igast, kg
   real(8) :: zji, toti, zji2

   real(8), allocatable, dimension(:,:)   :: ztemp2
   real(8), allocatable, dimension(:,:,:) :: ztemp3, psiga3_tmp
   real(8), allocatable, dimension(:)     :: vect, vect2


   wght = 0.d0

   if (mode.eq.'ADIAB') then
      ielst_ = nelst
   elseif (mode.eq.'DIAB2') then
      ielst_ = 2
   elseif (mode.eq.'DIAB4') then
      ielst_ = 1
      do istate=1,nstates
         wght(iset,istate) = 1.d0
      enddo
      return
   endif

   if (method.eq.1) then

      allocate (vect(ndabf),vect2(ndabf))

      if (mgquant.eq.1.or.mgquant.eq.2) then

         allocate (ztemp2(ndabf,ndabf))
         ztemp2 = reshape(z,(/ndabf,ndabf/))

      elseif (mgquant.eq.3) then

         nvibst = (nstates-1)/ngast + 1
         allocate (ztemp3(ndabf,ndabf,npntsg))
         ztemp3 = reshape(z,(/ndabf,ndabf,npntsg/))
         allocate (psiga3_tmp(nvibst,ngast,npntsg))
         psiga3_tmp = reshape(psiga,(/nvibst,ngast,npntsg/))

      endif

      do istate=1,nstates

         vect = 0.d0
         vect2 = 0.d0

         if (mgquant.eq.1.or.mgquant.eq.2) then

            do j=1,ndabf
               zji = ztemp2(j,istate)
               vect(j) = zji
               vect2(j) = zji*zji
            enddo

            do i=1,ielst

               ishift = i + 2*iset - 2
               toti = 0.d0

               do mu=1,nprst
                  do m=1,ngast
                     kk = (i-1)*nprst*ngast + (mu-1)*nprst + m
                     toti = toti + vect2(kk)
                  enddo
               enddo

               wght(ishift,istate) = toti*100.d0

            enddo

         elseif (mgquant.eq.3) then

            ivibst = (istate-1)/ngast + 1
            igast = mod(istate-1,ngast) + 1

            do j=1,ndabf
               zji2 = 0.d0
               do kg=1,npntsg
                  zji = ztemp3(j,ivibst,kg)
                  zji2 = zji2 + psiga3_tmp(ivibst,igast,kg)*&
                                zji*zji*&
                                psiga3_tmp(ivibst,igast,kg)
               enddo
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

         endif

      enddo

      deallocate (vect,vect2)
      if (mgquant.eq.1.or.mgquant.eq.2) then
         deallocate (ztemp2)
      elseif (mgquant.eq.3) then
         deallocate (ztemp3,psiga3_tmp)
      endif

   else

      write(6,*) ' in evbweig: method=2/3 not done yet'

      !do i=1,ielst
      !
      !  ishift = i + 2*iset - 2
      !  toti = 0.d0
      !
      !  do k=1,npnts
      !
      !     ck = psiel(1,k,i)
      !     c2k = ck*ck
      !     protwf = psipr(1,1,k)
      !     toti = toti + protwf*c2k*protwf
      !
      !  enddo
      !
      !  wght(ishift,1) = toti*100.d0
      !
      !enddo

   endif

   return

end subroutine evbweig
