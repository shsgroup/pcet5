subroutine weight2(nout,zp,ze,mode,iset,istate)

!====================================================================
!
!  Prints out the weights of the wavefunction
!
!  ZP,ZE - values of solvent coordinates
!  NOUT - unit for the output file
!  MODE - type of the state (ADIAB/DIAB2/DIAB4)
!  ISET - set of states (1 for ADIAB, 1/2 for DIAB2)
!  ISTATE - state to analyze
!
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-04-13 23:49:48 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:37  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!====================================================================

   use pardim
   use control
   use quantum
   use feszz_2d, only: feszz2
   use sorting

   implicit none

   integer,      intent(in) :: nout, iset, istate
   character(5), intent(in) :: mode
   real(8),      intent(in) :: zp, ze

   character(2), dimension(4) :: evb = (/'1a','1b','2a','2b'/)
   logical :: adiab, diab2, diab4
   integer :: ielst, nzdim, nz, nstates, nfe, npsiga, ndabf
   integer :: i, ishift, j, k, kk, ievb, ievbshift, ivib, ivibg, mu, m
   integer :: ivibst, igast, kg, jndex
   real(8) :: zji, sumc, contrib, toti, zji2, zji_sq, psiga3_sq
   real(8) :: ck, c2k, protwf, gatewf

   real(8), allocatable, dimension(:)       :: fe, z, psiga
   real(8), allocatable, dimension(:,:)     :: ztemp1
   real(8), allocatable, dimension(:,:,:)   :: ztemp3, psiga3
   real(8), allocatable, dimension(:,:,:,:) :: psiel, psipr, psiga4
   real(8), allocatable, dimension(:,:,:)   :: enel, envib, envibg
   real(8), allocatable, dimension(:)       :: vect2
   integer, allocatable, dimension(:)       :: ivect2


   adiab = .false.
   diab2 = .false.
   diab4 = .false.

   if (mode.eq.'ADIAB') then
      adiab=.true.
      ielst = nelst
   elseif (mode.eq.'DIAB2') then
      diab2=.true.
      ielst = 2
   elseif (mode.eq.'DIAB4') then
      diab4=.true.
      ielst = 1
   endif

   write(nout,'(/1x,60("-"))')
   write(nout,'(1x,"WAVEFUNCTION ANALYSIS (quantum gating mode)")')
   write(nout,'(1x,"mode=",a5,3x,"iset=",i1,3x,"istate=",i2,3x,"zp=",f8.3,3x,"ze=",f8.3)')&
   &mode,iset,istate,zp,ze
   write(nout,'(1x,60("-"))')

   if (diab4) then
      write(nout,'(/1x,"the information requested is meaningless..."/)')
      return
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the wavefunction
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (mgquant.eq.1.or.mgquant.eq.2) then
      nzdim = ielst*nprst*ngast
      nz = nzdim*nzdim
      nstates = istate
   elseif (mgquant.eq.3) then
      nzdim = ielst*nprst
      nz = nzdim*nzdim*npntsg
      nstates = (istate-1)/ngast + 1
   endif
   allocate (z(nz))

   nfe = istate
   npsiga = ielst*nprst*ngast*npntsg

   allocate (fe(nfe))
   allocate (psiel(ielst,npnts,npntsg,ielst))
   allocate (psipr(ielst,nprst,npnts,npntsg))
   allocate (psiga(npsiga))
   allocate (enel(ielst,npnts,npntsg))
   allocate (envib(ielst,nprst,npntsg))
   allocate (envibg(ielst,nprst,ngast))

   call feszz2(mode,iset,zp,ze,nstates,nfe,fe,nz,z,ndabf,&
              &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (method.eq.1) then

      write(nout,'(1x,"Diabatic electronic basis (Method=1)")')
      write(nout,'(1x,79("-"))')

      write(nout,'(/1x,"Major contributions:"/)')
      write(nout,'( t5, "EVB-index",&
                  & t20,"proton vib. st.",&
                  & t40,"gating vib. st.",&
                  & t60,"contribution(%)",/,&
                  & 1x,79("-"))')

      if (mgquant.eq.1.or.mgquant.eq.2) then

         allocate  (vect2(ndabf))
         allocate (ivect2(ndabf))
         vect2 = 0.d0
         ivect2 = 0

         allocate (ztemp1(ndabf,ndabf))
         ztemp1 = reshape(z,(/ndabf,ndabf/))

         do j=1,ndabf
            zji = ztemp1(j,istate)
            vect2(j) = zji*zji
         enddo

         call bubbli(ndabf,vect2,ivect2)

         sumc = 0.d0
         k = 0

         do while(sumc.le.0.95d0)

            k = k + 1
            jndex = ivect2(k)
            contrib = vect2(jndex)
            ievb = (jndex-1)/(nprst*ngast) + 1
            ievbshift = ievb + 2*iset - 2
            ivib  = (jndex - (ievb-1)*nprst*ngast - 1)/ngast + 1
            ivibg = jndex - (ievb-1)*nprst*ngast - (ivib-1)*ngast

            sumc = sumc + contrib

            write(nout,'(1x,t7,a2,t20,i3,t40,i3,t62,f8.3)')&
            &evb(ievbshift),ivib,ivibg,contrib*100.d0

         enddo
         write(nout,'(1x,60("-"))')

         write(nout,'(/1x,"Total EVB composition (summed up over the vibrational states):"/)')

         do i=1,ielst

            ishift = i + 2*iset - 2
            toti = 0.d0
            do mu=1,nprst
               do m=1,ngast
                  kk = (i-1)*nprst*ngast + (mu-1)*ngast + m
                  toti = toti + vect2(kk)
               enddo
            enddo

            write(nout,'(1x,"evb state (",a2,"): ",f8.3," %")') evb(ishift),toti*100.d0

         enddo

         deallocate  (vect2)
         deallocate (ivect2)
         deallocate (ztemp1)

      elseif (mgquant.eq.3) then

         allocate  (vect2(ndabf))
         allocate (ivect2(ndabf))
         vect2 = 0.d0
         ivect2 = 0

         allocate (psiga3(nstates,ngast,npntsg))
         psiga3 = reshape(psiga,(/nstates,ngast,npntsg/))

         allocate (ztemp3(ndabf,ndabf,npntsg))
         ztemp3 = reshape(z,(/ndabf,ndabf,npntsg/))

         ! weights are averaged over the corresponding gating wavefunction

         ivibst = (istate-1)/ngast + 1
         igast = mod(istate-1,ngast) + 1

         do j=1,ndabf
            zji2 = 0.d0
            do kg=1,npntsg
               zji = ztemp3(j,istate,kg)
               zji_sq = zji*zji
               psiga3_sq = psiga3(ivibst,igast,kg)*psiga3(ivibst,igast,kg)
               zji2 = zji2 + psiga3_sq*zji_sq
            enddo
            vect2(j) = zji*zji
         enddo

         call bubbli(ndabf,vect2,ivect2)

         sumc = 0.d0
         k = 0

         do while(sumc.le.0.95d0)

            k = k + 1
            jndex = ivect2(k)
            contrib = vect2(jndex)
            ievb = (jndex-1)/nprst + 1
            ievbshift = ievb + 2*iset - 2
            ivib  = jndex - (ievb-1)*nprst

            sumc = sumc + contrib

            write(nout,'(1x,t7,a2,t20,i3,t42,f8.3)') evb(ievbshift),ivib,contrib*100.d0

         enddo
         write(nout,'(1x,60("-"))')

         write(nout,'(/1x,"Total EVB composition (summed up over the vibrational states):"/)')

         do i=1,ielst

            ishift = i + 2*iset - 2
            toti = 0.d0
            do mu=1,nprst
                  kk = (i-1)*nprst + mu
                  toti = toti + vect2(kk)
            enddo

            write(nout,'(1x,"EVB state (",a2,"): ",f8.3," %")') evb(ishift),toti*100.d0

         enddo

         deallocate  (vect2)
         deallocate (ivect2)
         deallocate (ztemp3,psiga3)

      endif

   !------------------------------------------------------------------
   else

      write(nout,'(1X,"Adiabatic electronic basis (Method=2/3)")')
      write(nout,'(1X,"Weights are averaged over the ground state vibrational wavefunctions")')
      write(nout,'(1X,60("-"))')

      write(nout,'(/1X,"EVB weights for the ground state:"/)')
      write(nout,'( t5,"EVB-index",t20,"contribution(%)",/,1x,40("-"))')

      if (mgquant.eq.1.or.mgquant.eq.2) then
         allocate (psiga4(ielst,nprst,ngast,npntsg))
         psiga4 = reshape(psiga,(/ielst,nprst,ngast,npntsg/))
      elseif (mgquant.eq.3) then
         allocate (psiga3(nstates,ngast,npntsg))
         psiga3 = reshape(psiga,(/nstates,ngast,npntsg/))
      endif

      do i=1,ielst

         ishift = i + 2*iset - 2
         toti = 0.d0

         do k=1,npnts
            do kg=1,npntsg

               ck = psiel(1,k,kg,i)
               c2k = ck*ck
               protwf = psipr(1,1,k,kg)

               if (mgquant.eq.1.or.mgquant.eq.2) then
                  gatewf = psiga4(1,1,1,kg)
               elseif (mgquant.eq.3) then
                  gatewf = psiga3(1,1,kg)
               endif

               toti = toti + gatewf*protwf*c2k*protwf*gatewf

            enddo
         enddo

         write(nout,'(1x,"EVB state (",a2,"): ",f8.3," %")') evb(ishift),toti*100.d0

      enddo

      if (mgquant.eq.1.or.mgquant.eq.2) then
         deallocate (psiga4)
      elseif (mgquant.eq.3) then
         deallocate (psiga3)
      endif

   endif

   deallocate (fe,z,psiel,psipr,psiga)
   deallocate (enel,envib,envibg)

   return

end subroutine weight2
