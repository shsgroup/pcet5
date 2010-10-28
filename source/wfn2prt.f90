subroutine wfn2prt(nout,noutg,zp,ze,mode,iset,nvibst,ngatst)

!=====================================================================
!
!    Writes out the proton  vibrational wavefunctions
!
!    NOUT   - file unit for proton vibrational functions
!    NOUTG  - file unit for gating vibrational functions
!    ZP,ZE  - values of solvent coordinates
!    MODE   - type of the state (ADIAB/DIAB2/DIAB4)
!    ISET   - set of states (1 for ADIAB, 1/2 for DIAB2)
!    NVIBST - number of proton vibrational states to print
!    NGATST - number of gating vibrational states to print
!
!    Description of the output file NOUT (proton):
!
!    Columns:
!  -------------------------------------------------------------------
!  proton    electr.    vibr.      vibrational
!  coord.    energy     energy     wavefunction
!  -------------------------------------------------------------------
!  RLIST(K) GLIST(KG)
!  (ENEL(I,K,KG),(ENVIB(I,N,KG),PSIPR(I,N,K,KG),N=1,NVIBST),I=1,IELST)
!  -------------------------------------------------------------------
!
!    Description of the output file NOUTG (gating):
!
!    Columns:
!  -------------------------------------------------------------------
!  gating    proton vibr.    gating vibr.      vibrational
!  coord.    energy          energy            wavefunction
!  -------------------------------------------------------------------
!  if MGQUANT = 1 or 2
!
!  GLIST(KG)
!  ((ENVIB(I,N,KG),PSIGA(I,N,L,KG),L=1,NGATST),N=1,NVIBST),I=1,IELST)
!
!  if MGQUANT = 3
!
!  GLIST(KG) ((PSIGA(1,L),L=1,NGATST),I=1,IELST)
!  -------------------------------------------------------------------
!
!    IELST - number of electronic basis functions
!
!----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:37 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!=====================================================================

   use pardim
   use control
   use cst
   use quantum
   use feszz_2d, only: feszz2
   use feszz_3d, only: feszz3

   implicit none
   character(5), intent(in) :: mode
   integer,      intent(in) :: nout,noutg,iset,nvibst,ngatst
   real(8),      intent(in) :: zp, ze

   logical adiab,diab2,diab4
   integer :: ielst, nstates, nzdim, nz, nfe, npsiga, ndabf, ndabf3, nz3dim
   integer :: k, kg, i, n, l

   real(8), allocatable, dimension(:)       :: fe, z, psiga
   real(8), allocatable, dimension(:,:,:,:) :: psiel, psipr, psiga4
   real(8), allocatable, dimension(:,:,:)   :: enel, envib, envibg, psiga3

   real(8), allocatable, dimension(:)     :: f3
   real(8), allocatable, dimension(:,:)   :: z3, enel3, envib3
   real(8), allocatable, dimension(:,:,:) :: psiel3, psipr3


   if (mode.eq.'ADIAB') then
      adiab = .true.
      ielst = nelst
   elseif (mode.eq.'DIAB2') then
      diab2 = .true.
      ielst = 2
   elseif (mode.eq.'DIAB4') then
      diab4 = .true.
      ielst = 1
   endif

   !~~~~~~~~~~~~~~~~~~~~~
   ! Calculate the states
   !~~~~~~~~~~~~~~~~~~~~~

   nstates = 1
   if (mgquant.eq.1.or.mgquant.eq.2) then
      nzdim = ielst*nprst*ngast
      nz = nzdim*nzdim
      nfe = nstates
   elseif (mgquant.eq.3) then
      nzdim = ielst*nprst
      nz = nzdim*nzdim*npntsg
      nfe = nstates*ngatst
   endif
   allocate (z(nz))

   npsiga = ielst*nprst*ngast*npntsg
   allocate (fe(nfe))
   allocate (psiel(ielst,npnts,npntsg,ielst), psipr(ielst,nprst,npnts,npntsg), psiga(npsiga))
   allocate (enel(ielst,npnts,npntsg), envib(ielst,nprst,npntsg), envibg(ielst,nprst,ngast))

   call feszz2(mode,iset,zp,ze,nstates,nfe,fe,nz,z,ndabf,&
              &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)

   !~~~~~~~~~~~~~~~~~~~~~
   ! Print out to NOUT
   !~~~~~~~~~~~~~~~~~~~~~
   do kg=1,npntsg
      do k=1,npnts
         write(nout,'(1x,50g15.6)') rlist(k)*bohr2a,glist(kg)*bohr2a,&
        &(enel(i,k,kg),(envib(i,n,kg),envib(i,n,kg)+psipr(i,n,k,kg)*100.d0,n=1,nvibst),i=1,ielst)
      enddo
      write(nout,*)
   enddo

   !~~~~~~~~~~~~~~~~~~~~~
   ! Print out to NOUTG
   !~~~~~~~~~~~~~~~~~~~~~
   if (mgquant.eq.1.or.mgquant.eq.2) then

      allocate (psiga4(ielst,nprst,ngast,npntsg))
      psiga4 = reshape(psiga,(/ielst,nprst,ngast,npntsg/))
      do kg=1,npntsg
         write(noutg,'(1x,50g15.6)') glist(kg)*bohr2a,&
       &((envib(i,n,kg),n=1,nvibst),i=1,ielst),&
      &(((envibg(i,n,l),envibg(i,n,l)+psiga4(i,n,l,kg)*50.d0,&
         &l=1,ngatst),n=1,nvibst),i=1,ielst)
      enddo
      deallocate (psiga4)

   elseif (mgquant.eq.3) then

      allocate (psiga3(nstates,ngast,npntsg))
      psiga3 = reshape(psiga,(/nstates,ngast,npntsg/))

      nz3dim = nprst*ielst
      allocate (f3(nstates))
      allocate (z3(nz3dim,nz3dim))
      allocate (psiel3(ielst,npnts,ielst),psipr3(ielst,nprst,npnts))
      allocate (enel3(ielst,npnts),envib3(ielst,nprst))

      do kg=1,npntsg
         call feszz3(mode,iset,kg,zp,ze,nstates,f3,nz3dim,z3,ndabf3,&
                    &ielst,enel3,envib3,psiel3,psipr3)
         write(noutg,'(1x,50g15.6)') glist(kg)*bohr2a,&
        &(f3(i),(fe(i+l-1),fe(i+l-1)+psiga3(i,l,kg)*50.d0,l=1,ngatst),i=1,nstates)
      enddo

      deallocate (psiga3)
      deallocate (f3,z3)
      deallocate (psiel3,psipr3)
      deallocate (enel3,envib3)

   endif

   deallocate (fe,z)
   deallocate (psiel,psipr,psiga)
   deallocate (enel,envib,envibg)

   return

end subroutine wfn2prt
