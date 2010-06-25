subroutine wfn3prt(nout,zp,ze,mode,iset,nvibst)
!===================================================================
!  Writes out the proton vibrational wavefunctions
!  as two-dimensional functions of proton and gating coordinates
!
!  NOUT - file unit
!  ZP,ZE - values of solvent coordinates
!  MODE - type of the state (ADIAB/DIAB2/DIAB4)
!  ISET - set of states (1 for ADIAB, 1/2 for DIAB2)
!  NVIBST - number of vibrational states to print
!
!  Description of the output file:
!
!  Columns:
!-------------------------------------------------------------------
!  proton   gating     electr.    vibr.      vibrational
!  coord.   coord.     energy     energy     wavefunction
!-------------------------------------------------------------------
!  RLIST(K) GLIST(KG) (ENEL(I,K),(ENVIB(I,N),PSIPR(I,N,K),N=1,NVIBST),I=1,IELST
!-------------------------------------------------------------------
!
!  IELST - number of electronic basis functions
!
!--------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:37
!  4.1
!  Exp
!  wfn3prt.f90,v 4.1 2010/06/25 20:02:37 souda Exp
!  wfn3prt.f90,v
!  Revision 4.1  2010/06/25 20:02:37  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:13:52  souda
!  Initial PCET-4.0 Release
!
!
!======================================================================
   use pardim
   use cst
   use quantum
   use feszz_3d, only: feszz3

   implicit none
   character(5), intent(in) :: mode
   integer,      intent(in) :: nout, iset, nvibst
   real(8),      intent(in) :: zp, ze

   logical :: adiab, diab2, diab4
   integer :: ielst, nzdim, kg, k, i, n, ndabf

   real(8), allocatable, dimension(:)     :: fe
   real(8), allocatable, dimension(:,:)   :: z
   real(8), allocatable, dimension(:,:)   :: enel, envib
   real(8), allocatable, dimension(:,:,:) :: psiel, psipr

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

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Loop over grid points along the gating coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   nzdim = ielst*nprst
   allocate (fe(1), z(nzdim,nzdim))
   allocate (psiel(ielst,npnts,ielst), psipr(ielst,nprst,npnts))
   allocate (enel(ielst,npnts), envib(ielst,nprst))

   do kg=1,npntsg

      ! calculate the states
      call feszz3(mode,iset,kg,zp,ze,1,fe,nzdim,z,ndabf,ielst,enel,envib,psiel,psipr)

      do k=1,npnts
         write(nout,'(1x,40g15.6)') rlist(k)*bohr2a,glist(kg)*bohr2a,&
         &(enel(i,k),(envib(i,n),envib(i,n)+psipr(i,n,k)*20.d0,n=1,nvibst),i=1,ielst)
      enddo
      write(nout,*)

   enddo

   deallocate (fe,z)
   deallocate (psiel,psipr)
   deallocate (enel,envib)

   return

end subroutine wfn3prt
