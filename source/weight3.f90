subroutine weight3(nout,kg,zp,ze,mode,iset,istate)

!====================================================================
!
!  Prints out the weights of the wavefunction
!
!  KG     - value (grid point) of gating coordinate
!  ZP,ZE  - values of solvent coordinates
!  NOUT   - unit for the output file
!  MODE   - type of the state (ADIAB/DIAB2/DIAB4)
!  ISET   - set of states (1 for ADIAB, 1/2 for DIAB2)
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
   use cst
   use control
   use quantum
   use feszz_3d, only: feszz3
   use sorting

   implicit none

   integer,      intent(in) :: nout, kg, iset, istate
   character(5), intent(in) :: mode
   real(8),      intent(in) :: zp, ze

   character(2), dimension(4) :: evb = (/'1a','1b','2a','2b'/)
   logical :: adiab, diab2, diab4

   integer :: ielst, nzdim, ndabf
   integer :: i, ishift, j, mu, k, kk, ievb, ievbshift, ivib, jndex
   real(8) :: rr, zji, sumc, contrib, toti, ck, c2k, protwf

   real(8), allocatable, dimension(:)     :: fe
   real(8), allocatable, dimension(:,:)   :: z
   real(8), allocatable, dimension(:,:,:) :: psiel, psipr
   real(8), allocatable, dimension(:,:)   :: enel, envib
   real(8), allocatable, dimension(:)     :: vect2
   integer, allocatable, dimension(:)     :: ivect2


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

   rr = glist(kg)*bohr2a

   write(nout,'(/1x,60("-"))')
   write(nout,'(1x,"WAVEFUNCTION ANALYSIS at R=",f6.3)') rr
   write(nout,'(1x,"mode=",a5,3x,"iset=",i1,3x,"istate=",i2,3x,"zp=",f8.3,3x,"ze=",f8.3)')&
   &mode,iset,istate,zp,ze
   write(nout,'(1x,60("-"))')

   if (diab4) then
      write(nout,'(/1x,"The information requested is meaningless..."/)')
      return
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the wavefunction
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   nzdim = ielst*nprst
   allocate (fe(istate))
   allocate (z(nzdim,nzdim))
   allocate (psiel(ielst,npnts,ielst))
   allocate (psipr(ielst,nprst,npnts))
   allocate (enel(ielst,npnts))
   allocate (envib(ielst,nprst))

   call feszz3(mode,iset,kg,zp,ze,istate,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the weights
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (method.eq.1) then

      WRITE(NOUT,'(1X,"Diabatic electronic basis (Method=1)")')
      WRITE(NOUT,'(1X,60("-"))')

      WRITE(NOUT,'(/1X,"Major contributions:"/)')
      write(nout,'( t5,"EVB-index",t20,"vibr. state",t40,"contribution(%)",/,1x,60("-"))')

      allocate  (vect2(ndabf))
      allocate (ivect2(ndabf))
      vect2 = 0.d0
      ivect2 = 0

      do j=1,ndabf
         zji = z(j,istate)
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
         ivib = jndex - (ievb-1)*nprst
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

   else

      write(nout,'(1x,"Adiabatic electronic basis (Method=2/3)")')
      write(nout,'(1x,"Weights are averaged over the ground state vibrational wavefunction")')
      write(nout,'(1x,60("-"))')

      write(nout,'(/1x,"EVB weights for the ground state:"/)')
      write(nout,'( t5,"EVB-index",t20,"contribution(%)",/,1x,40("-"))')

      do i=1,ielst

         ishift = i + 2*iset - 2
         toti = 0.d0

         do k=1,npnts

            ck = psiel(1,k,i)
            c2k = ck*ck
            protwf = psipr(1,1,k)
            toti = toti + protwf*c2k*protwf

         enddo

         write(nout,'(1x,"EVB state (",a2,"): ",f8.3," %")') evb(ishift),toti*100.d0

      enddo

   endif

   deallocate (fe,z,psiel,psipr)
   deallocate (enel,envib)

   return

end subroutine weight3


