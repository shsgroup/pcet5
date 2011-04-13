module feszz_2d

!======================================================================     
!  2D free energy surface (quantization for proton and gating modes)
!======================================================================     
!
!  $Author: souda $
!  $Date: 2011-04-13 23:49:48 $
!  $Revision: 5.4 $
!  $Log: not supported by cvs2svn $
!  Revision 5.3  2011/02/20 00:58:11  souda
!  Major additions/modifications:
!  (1) precalculation of the proton vibrational basis functions for METHOD=1
!  (2) Franck-Condon initial excitation added to DYNAMICS3
!  (3) addition of the module timers: module_timers.f90 (changes in Makefile!)
!
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!======================================================================     

   use pardim
   use cst
   use timers
   use control
   use gasmat
   use solmat
   use quantum
   use eispack, only: rs
   use lapack_wrappers
   use turbomole_wrappers

   implicit none
   private
   
   public :: feszz2
   public :: evbweig
   
contains

   subroutine feszz2(mode,iset,zp,ze,nstates,nfe,fe,nz,z,ndabf,&
                    &ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
   !======================================================================
   !
   !     Calculates free energy as a function of two solvent variables
   !     for proton-coupled electron transfer reactions.
   !     The continuum model is utilized. The electron and
   !     the proton are treated quantum-mechanically.
   !     The gating mode quantum states are adiabatic with respect
   !     to the solvent coordinates and double adiabatic with
   !     respect to the electron-proton vibronic states.
   !
   !     MODE - calculation mode:
   !         ='ADIAB' - adiabatic electron/proton/gating vibrational
   !                    free energies;
   !         ='DIAB2' - ET diabatic free energy surfaces;
   !         ='DIAB4' - diabatic (EVB) free energy surfaces.
   !
   !     ISET - set of states to calculate (only one set is available
   !            for ADIAB, two sets for DIAB2 and four sets for DIAB4)
   !                                                                          
   !     ZP, ZE are values of the effective medium coordinates (kcal/mol)     
   !
   !     NSTATES - number of states to calculate
   !               (in case of MGQUANT=3 it is a number of electron-proton
   !                vibronic states, so the total number of states will be
   !                NSTATES*NGAST)
   !
   !     FE(I) - free energy of vibronic (including gating) states (kcal/mol)
   !        I - vibronic quantum number
   !
   !     Z(I,J) - eigenvectors of the total vibronic
   !              (including gating mode) Hamiltonian
   !
   !     NZDIM - dimension of the total Hamiltonian matrix
   !             (dimension of the basis) (input)
   !
   !     NDABF - dimension of the total Hamiltonian matrix
   !             (dimension of the basis) (output)
   !
   !     IELST - number of electronic states included in calculation
   !
   !     ENEL(I,K,KG) - electronic energies (kcal/mol)
   !          I - electronic quantum number
   !            K - the point number along the proton coordinate
   !              KG - the point number along the gating coordinate
   !
   !     ENVIB(I,N,KG) - proton vibrational levels (kcal/mol)
   !           I - electronic quantum number
   !             N - vibrational quantum number
   !               KG - the point number along the gating coordinate
   !
   !     ENVIBG(I,N,M) - gating vibrational levels (kcal/mol)
   !            I - electronic quantum number
   !              N - proton vibrational quantum number
   !                M - gating mode quantum number
   !
   !     PSIEL(I,K,KG,L) - electronic wave function (coefficients)
   !           I - quantum number
   !             K - the point number along the proton coordinate
   !               KG - the point number along the gating coordinate
   !                  L - index of the EVB basis state
   !
   !     PSIPR(I,N,K,KG) - protonic wave function (coefficients)
   !           I - electronic quantum number
   !             N - protonic quantum number
   !               K - the point number along the proton coordinate
   !                 KG - the point number along the gating coordinate
   !
   !c     PSIGA(I,N,M,KG) - gating vibrational wave function (coefficients)
   !c           I - electronic quantum number
   !c             N - protonic quantum number
   !c               M - gating mode quantum number
   !c                 KG - the point index of the gating mode grid
   !
   !     PSIGA(*) - gating vibrational wave function (coefficients)
   !           one-dimensional array with a structure depending on
   !           MGQANT:
   !           1 or 2 - ((((KG=1,NPNTSG),M=1,NGAST),N=1,NPRST),I=1,IELST)
   !                        I  - electronic quantum number
   !                        N  - protonic quantum number
   !                        M  - gating mode quantum number
   !                        KG - the point index of the gating mode grid
   !
   !           3 - (((KG=1,NPNTSG),M=1,NGAST),I=1,NSTATES)
   !                  I  - vibronic quantum number
   !                  M  - gating mode quantum number
   !                  KG - the point index of the gating mode grid
   !
   !======================================================================     

   character(5), intent(in) :: mode
   integer,      intent(in) :: iset, nstates, nfe, nz, ielst, npsiga
   real(8),      intent(in) :: zp, ze

   integer, intent(out)                                      :: ndabf
   real(8), intent(out), dimension(nfe)                      :: fe
   real(8), intent(out), dimension(nz)                       :: z
   real(8), intent(out), dimension(ielst,npnts,npntsg,ielst) :: psiel
   real(8), intent(out), dimension(ielst,nprst,npnts,npntsg) :: psipr
   real(8), intent(out), dimension(npsiga)                   :: psiga
   real(8), intent(out), dimension(ielst,npnts,npntsg)       :: enel
   real(8), intent(out), dimension(ielst,nprst,npntsg)       :: envib
   real(8), intent(out), dimension(ielst,nprst,ngast)        :: envibg

   ! Local saved variables
   ! the following is only for maintaining the phases
   ! of the wavefunctions along the solvent coordinates.
   ! It is not necessary unless we want to calculate
   ! non-adiabatic coupling vector for solvent dynamics
   !
   !integer, save                          :: numr    = 0
   !real*8, save, dimension(maxsta,maxsta) :: zprev2  = 0.d0
   !real*8, save, dimension(nelst,nelst)   :: zeprev2 = 0.d0
   !real*8, save, dimension(maxpnt,maxpnt) :: zpprev2 = 0.d0

   ! Local variables

   character(5) :: elmode
   logical :: adiab, diab2, diab4

   integer :: ielst0, i, j, n, m, l, lg, k, kg, kg1, kg2, kk, ll
   integer :: ierr, iel, ipr, iga, jel, jpr, jga, nzdim, nzdim2
   integer :: ielshift, jelshift, izstart, izend, istate, kz
   real(8) :: gm, p, termg, termd, termd1, termd2, dterm, gterm
   real(8) :: hfij, dtij, gtij, psipr_ovl, psiga4_ovl

   integer, allocatable, dimension(:)       :: istel, istpr, istga
   real(8), allocatable, dimension(:,:)     :: hf, zzprev2
   real(8), allocatable, dimension(:)       :: wfe(:), fv1(:), fv2(:)
   real(8), allocatable, dimension(:,:,:,:) :: dij, gij
   real(8), allocatable, dimension(:,:,:,:) :: dimunu, gimunu
   real(8), allocatable, dimension(:,:)     :: envibg3, zsquare, fe3
   real(8), allocatable, dimension(:,:,:)   :: psiga3
   real(8), allocatable, dimension(:,:,:,:) :: psiga4

   ! increment number of calls
   !numr = numr + 1

   adiab = mode.eq.'ADIAB'
   diab2 = mode.eq.'DIAB2'
   diab4 = mode.eq.'DIAB4'

   ! Reduced mass for gating coordinate (A.U.) 
   gm = dm*am/(dm+am)

   ! Number of electronic basis states

   if (adiab) then
      ielst0 = nelst
   elseif(diab2) then
      ielst0 = 2
   elseif (diab4) then
      ielst0 = 1
   else
      ielst0 = nelst
   endif

   if (ielst0.ne.ielst) then
      write(6,*) ' Wrong input for feszz2: ielst=',ielst
      stop
   endif

   ! Calculation of the electronic states
   ! U     - solvated adiabatic states
   ! CU    - solvated adiabatic eigenvectors (coefficients)

   !enel  = 0.d0
   !psiel = 0.d0
   !envib = 0.d0
   !psipr = 0.d0
   !psiga = 0.d0

   if (method.eq.1) then

      elmode = 'DIAB4'

   elseif (method.eq.2.or.method.eq.3) then

      if (adiab) elmode = 'ADIAB'
      if (diab2) elmode = 'DIAB2'
      if (diab4) elmode = 'DIAB4'

   else

      elmode = 'DIAB4'

   endif

   !----------------------------------------------------------------------
   ! calculate electronic (enel, psiel) and protonic (envib, psipr) states
   !----------------------------------------------------------------------
   call electron_proton_states

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   ! building the d_ij and g_ij matrix elements as functions   
   ! of the proton and gating coordinates                      
   ! (needed only for method=2)                                
   !                                                           
   ! dij(i,j,k,kg), k,kg are points along the proton and       
   !                gating coordinates                         
   ! gij(i,j,k,kg), k,kg are points along the proton and       
   !                gating coordinates                         
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

   if (method.eq.2.and..not.diab4) then

      allocate (dij(ielst,ielst,npnts,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('dij')

      allocate (gij(ielst,ielst,npnts,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('gij')

      call calcdij2
      call calcgij2

   endif

   !======================================================================
   ! building of the electronic/proton/gating vibrational basis
   ! and calculating the free energies of states
   !======================================================================

   !--------------------------------------------------------------
   ! mgquant = 1/2
   !--------------------------------------------------------------
   if (mgquant.eq.1.or.mgquant.eq.2) then

      ! the total hamiltonian is built in the basis of the products
      ! of electronic, proton vibrational and gating vibrational
      ! states:
      !
      ! phi(i,mu,m) = psi(i) x chi(i,mu) x phi(i,mu,m)
      !
      ! mgquant=1 all the non-adiabatic terms originating from
      !           the kinetic energy operator for gating motion
      !           are included.
      !
      ! mgquant=2 all the non-adiabatic terms originating from
      !           the kinetic energy operator for gating motion
      !           are ignored (the electron-proton non-adiabatic
      !           terms are treated according to the method)

      nzdim = ielst*nprst*ngast
      allocate (istel(nzdim),istpr(nzdim),istga(nzdim),stat=ierr)
      if (ierr.ne.0) call alloc_error('istel,istpr,istga')
      istel = 0
      istpr = 0
      istga = 0

      ! ordering of the basis states
      !
      ! istel(k) - the quantum number of the electronic state
      ! istpr(k) - the quantum number of the protonic state
      ! istga(k) - the quantum number of the gating state
 
      ndabf = 0
      do i=1,ielst
         do n=1,nprst
            do m=1,ngast
               ndabf = ndabf + 1
               istel(ndabf) = i
               istpr(ndabf) = n
               istga(ndabf) = m
            enddo
         enddo
      enddo

      ! building the d_imn(r) and g_imn(r) matrix elements
      ! as functions of the gating coordinate
      ! (needed only if mqquant=1)
      !
      ! dij(i,j,k,kg), k,kg are points along the proton and
      !                gating coordinates
      ! gij(i,j,k,kg), k,kg are points along the proton and
      !                gating coordinates

      if (mgquant.eq.1) then
         allocate (dimunu(ielst,nprst,nprst,npntsg),&
                   gimunu(ielst,nprst,nprst,npntsg),stat=ierr)
         if (ierr.ne.0) call alloc_error('dimunu,gimunu')
         dimunu = 0.d0
         gimunu = 0.d0
         call calcdimn
         call calcgimn
      endif

      ! calculation of the gating vibrational states for each
      ! basis electron-proton state

      allocate (psiga4(ielst,nprst,ngast,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('psiga4')
      psiga4 = 0.d0

      call gating_states1
      psiga = reshape(psiga4,(/npsiga/),(/0.d0/))

      if (method.eq.3.and.mgquant.eq.2) then

         ! no need to build a hamiltonian matrix,
         ! it is diagonal

         do i=1,nstates
           iel = istel(i)
           ipr = istpr(i)
           iga = istga(i)
           fe(i) = envibg(iel,ipr,iga)
         enddo

         ! deallocate arrays and return

         ! clean_memory
         if (allocated(istel  )) deallocate(istel  )
         if (allocated(istpr  )) deallocate(istpr  )
         if (allocated(istga  )) deallocate(istga  )
         if (allocated(dij    )) deallocate(dij    )
         if (allocated(gij    )) deallocate(gij    )
         if (allocated(dimunu )) deallocate(dimunu )
         if (allocated(gimunu )) deallocate(gimunu )
         if (allocated(psiga4 )) deallocate(psiga4 )
         return

      endif

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! building of the total hamiltonian matrix hf(i,j)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      allocate (hf(ndabf,ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('hf')
      hf = 0.d0

      do i=1,ndabf

         iel = istel(i)
         ielshift = iel + 2*iset - 2
         ipr = istpr(i)
         iga = istga(i)

         do j=i,ndabf

            jel = istel(j)
            jelshift = jel + 2*iset - 2
            jpr = istpr(j)
            jga = istga(j)

            hfij = 0.d0

            if (iel.eq.jel.and.ipr.eq.jpr.and.iga.eq.jga) then
               hfij = hfij + envibg(iel,ipr,iga)*cal2au
            endif

            if (iel.eq.jel.and.mgquant.eq.1) then

               ! terms with derivatives of proton states
               ! with respect to the gating coordinate
               ! (proton-gating non-adiabatic terms)

               termg = 0.d0
               do kg=1,npntsg
                  psiga4_ovl = psiga4(iel,ipr,iga,kg)*psiga4(iel,jpr,jga,kg)
                  termg = termg + psiga4_ovl*gimunu(iel,ipr,jpr,kg)
               enddo

               termd = 0.d0
               do lg=1,ngast

                  termd1 = 0.d0
                  do kg=1,npntsg
                     psiga4_ovl = psiga4(iel,ipr,iga,kg)*psiga4(iel,jpr,lg,kg)
                     termd1 = termd1 + psiga4_ovl*dimunu(iel,ipr,jpr,kg)
                  enddo

                  termd2 = 0.d0
                  do kg1=1,npntsg
                     do kg2=1,npntsg
                        psiga4_ovl = psiga4(iel,jpr,lg,kg1)*psiga4(iel,jpr,jga,kg2)
                        termd2 = termd2 + psiga4_ovl*dgx(kg1,kg2)
                     enddo
                  enddo

                  termd = termd + termd1*termd2

               enddo
               termd = termd*2.d0

               hfij = hfij - hbar*hbar/(2.d0*gm)*(termg + termd)

               ! terms with derivatives of electronic states
               ! with respect to the gating coordinate
               ! (electron-gating non-adiabatic terms)

               ! ignore them... too small and complicated

            endif

            ! electron-proton non-adiabatic terms

            if (method.eq.1) then

               if (iel.ne.jel) then
                  do kg=1,npntsg
                     do k=1,npnts
                        psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,k,kg)
                        psiga4_ovl = psiga4(iel,ipr,iga,kg)*psiga4(jel,jpr,jga,kg)
                        hfij = hfij + cal2au*h0(ielshift,jelshift,k,kg)*psipr_ovl*psiga4_ovl
                     enddo
                  enddo
               endif

            elseif (method.eq.2) then

               dterm = 0.d0
               if (iel.ne.jel) then
                  do kg=1,npntsg
                     do k=1,npnts
                        do l=1,npnts
                           psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,l,kg)
                           psiga4_ovl = psiga4(iel,ipr,iga,kg)*psiga4(jel,jpr,jga,kg)
                           dterm = dterm + psipr_ovl*psiga4_ovl*dij(iel,jel,k,kg)*dx(l,k)
                        enddo
                     enddo
                  enddo
                  dtij = -hbar*hbar*dterm/pm
                  hfij = hfij + dtij
               endif

               gterm = 0.d0
               do kg=1,npntsg
                  do k=1,npnts
                     psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,k,kg)
                     psiga4_ovl = psiga4(iel,ipr,iga,kg)*psiga4(jel,jpr,jga,kg)
                     gterm = gterm + psipr_ovl*psiga4_ovl*gij(iel,jel,k,kg)
                  enddo
               enddo

               gtij = -hbar*hbar*gterm/(2.d0*pm)
               hfij = hfij + gtij

            endif

            hf(i,j) = hfij
            if (i.ne.j) hf(j,i) = hfij

         enddo
      enddo

      deallocate (istel,istpr,istga)
      deallocate (psiga4)
      if (mgquant.eq.1) deallocate (dimunu,gimunu)
      if (method.eq.2.and..not.diab4) deallocate (dij,gij)

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! diagonalization of the hamiltonian matrix hf(i,j)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      allocate (zsquare(ndabf,ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('zsquare')
      zsquare = 0.d0

      allocate (fv1(ndabf),fv2(ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('fv1,fv2')
      fv1 = 0.d0
      fv2 = 0.d0

      allocate (wfe(ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('wfe')
      wfe = 0.d0

      !call rs(ndabf,ndabf,hf,wfe,ndabf,zsquare,fv1,fv2,ierr)
      !call dspevx_wrapper(ndabf,hf,wfe,z,ierr)
      !call dspevd_wrapper(ndabf,hf,wfe,z,ierr)
      !call rdiag_wrapper(ndabf,hf,wfe,z,ierr)
      call dsyevr_wrapper(ndabf,hf,wfe,z,ierr)

      if (ierr.ne.0) write(*,*) "Diagonalization error in feszz2, IERR =", ierr

      do i=1,nstates
         fe(i) = wfe(i)*au2cal
      enddo

      deallocate (hf,wfe)
      deallocate (fv1,fv2)

      ! check phases

      !if (numr.gt.1) then
      !   do kk=1,ndabf
      !      p = dot_product(zsquare(1:ndabf,kk),zprev2(1:ndabf,kk))
      !      if (p.lt.0d0) then
      !         do ll=1,ndabf
      !            zsquare(ll,kk) = -zsquare(ll,kk)
      !         enddo
      !      endif
      !   enddo
      !endif
      !zprev2(1:ndabf,1:ndabf) = zsquare

      nzdim2 = ndabf*ndabf
      z = reshape(zsquare,(/nzdim2/),(/0.d0/))
      deallocate (zsquare)

   !--------------------------------------------------------------
   ! mgquant = 3
   !--------------------------------------------------------------
   elseif (mgquant.eq.3) then

      ! the total hamiltonian is built in the basis of the products
      ! of electronic and proton vibrational states.
      ! gating states are calculated on top of them
      !
      ! phi(i,mu,m) = [psi(i) x chi(i,mu)]_k(zp,ze) x phi(k,m)

      nzdim = ielst*nprst
      allocate (istel(nzdim),istpr(nzdim),stat=ierr)
      if (ierr.ne.0) call alloc_error('istel,istpr')
      istel = 0
      istpr = 0

      ! ordering of the basis states
      !
      ! istel(k) - the quantum number of the electronic state
      ! istpr(k) - the quantum number of the protonic state

      ndabf  = 0
      do i=1,ielst
         do n=1,nprst
            ndabf = ndabf + 1
            istel(ndabf) = i
            istpr(ndabf) = n
         enddo
      enddo

      allocate (hf(ndabf,ndabf),wfe(ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('hf,wfe')
      hf = 0.d0
      wfe = 0.d0
      allocate (fv1(ndabf),fv2(ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('fv1,fv2')
      fv1 = 0.d0
      fv2 = 0.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! now, build hamiltonian matrix for electron-proton states
      ! at each gating grid point
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      allocate (fe3(ndabf,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('fe3')
      fe3 = 0.d0

      allocate (zsquare(ndabf,ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('zsquare')
      zsquare = 0.d0

      allocate (zzprev2(ndabf,ndabf),stat=ierr)
      if (ierr.ne.0) call alloc_error('zzprev2')
      zzprev2 = 0.d0

      do kg=1,npntsg

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! in double adiabatic approximation (method=3): all d- and g-terms
         ! are zero and the total hamiltonian matrix is diagonal.
         ! also if mode='diab4' then there is no off-diagonal matrix elements
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (method.eq.3.or.diab4) then

            do i=1,nstates
              iel = istel(i)
              ipr = istpr(i)
              fe3(i,kg) = envib(iel,ipr,kg)
            enddo

            do kz=1,ndabf
               zsquare(kz,kz) = 1.d0
            enddo

            izstart = (kg-1)*ndabf*ndabf + 1
            izend   = kg*ndabf*ndabf
            z(izstart:izend) = reshape(zsquare,(/ndabf*ndabf/))
            cycle

         endif

         do i=1,ndabf

            iel = istel(i)
            ielshift = iel + 2*iset - 2
            ipr = istpr(i)

            do j=i,ndabf

               jel = istel(j)
               jelshift = jel + 2*iset - 2
               jpr = istpr(j)

               hfij = 0.d0

               if (iel.eq.jel.and.ipr.eq.jpr) then
                  hfij = hfij + envib(iel,ipr,kg)*cal2au
               endif

               if (method.eq.1) then

                  if (iel.ne.jel) then
                     do k=1,npnts
                        psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,k,kg)
                        hfij = hfij + h0(ielshift,jelshift,k,kg)*cal2au*psipr_ovl
                     enddo
                  endif

               else

                  dterm = 0.d0
                  if (iel.ne.jel) then
                     do k=1,npnts
                        do l=1,npnts
                           psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,l,kg)
                           dterm = dterm + psipr_ovl*dij(iel,jel,k,kg)*dx(l,k)
                        enddo
                     enddo
                     dtij = -hbar**2*dterm/pm
                     hfij = hfij + dtij
                  endif

                  gterm = 0.d0
                  do k=1,npnts
                     psipr_ovl = psipr(iel,ipr,k,kg)*psipr(jel,jpr,k,kg)
                     gterm = gterm + psipr_ovl*gij(iel,jel,k,kg)
                  enddo

                  gtij = -hbar**2*gterm/(2.d0*pm)
                  hfij = hfij + gtij

               endif

               hf(i,j) = hfij
               if (i.ne.j) hf(j,i) = hfij

            enddo
         enddo

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! diagonalization of the hamiltonian matrix hf(i,j)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !call rs(ndabf,ndabf,hf,wfe,ndabf,zsquare,fv1,fv2,ierr)
         !call dspevx_wrapper(ndabf,hf,wfe,z,ierr)
         !call dspevd_wrapper(ndabf,hf,wfe,z,ierr)
         !call rdiag_wrapper(ndabf,hf,wfe,z,ierr)
         call dsyevr_wrapper(ndabf,hf,wfe,z,ierr)

         if (ierr.ne.0) write(*,*) "Diagonalization error in feszz2, IERR =", ierr

         ! check phases: maintain the phase along the gating grid

         !if (kg.eq.1) then
         !   if (numr.gt.1) then
         !      do kk=1,ndabf
         !         p = dot_product(zsquare(1:ndabf,kk),zprev2(1:ndabf,kk))
         !         if (p.lt.0d0) then
         !            do ll=1,ndabf
         !               zsquare(ll,kk) = -zsquare(ll,kk)
         !            enddo
         !         endif
         !      enddo
         !   endif
         !   zprev2(1:ndabf,1:ndabf) = zsquare
         !else

         if (kg.gt.1) then
            do kk=1,ndabf
               p = dot_product(zsquare(1:ndabf,kk),zzprev2(1:ndabf,kk))
               if (p.lt.0d0) then
                  do ll=1,ndabf
                     zsquare(ll,kk) = -zsquare(ll,kk)
                  enddo
               endif
            enddo
         endif
         zzprev2 = zsquare

         ! fill the array for the free energy on the gating greed
         do i=1,ndabf
            fe3(i,kg) = wfe(i)*au2cal
         enddo

         ! append the matrix of coefficients to the output array
         izstart = (kg-1)*ndabf*ndabf + 1
         izend   = kg*ndabf*ndabf
         z(izstart:izend) = reshape(zsquare,(/ndabf*ndabf/))

      enddo

      deallocate (zsquare,stat=ierr)
      if (ierr.ne.0) call dealloc_error('zsquare')

      deallocate (zzprev2,stat=ierr)
      if (ierr.ne.0) call dealloc_error('zzprev2')

      deallocate (hf,wfe,stat=ierr)
      if (ierr.ne.0) call dealloc_error('hf,wfe')

      deallocate (fv1,fv2,stat=ierr)
      if (ierr.ne.0) call dealloc_error('fv1,fv2')

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! calculation of the gating vibrational states for each
      ! electron-proton vibronic state
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      allocate (psiga3(nstates,ngast,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('psiga3')

      allocate (envibg3(nstates,ngast),stat=ierr)
      if (ierr.ne.0) call alloc_error('envibg3')

      call gating_states3

      istate = 0
      do i=1,nstates
         do m=1,ngast
            istate = istate + 1
            if (istate.le.nfe) fe(istate) = envibg3(i,m)
         enddo
      enddo

      psiga = reshape(psiga3,(/npsiga/),(/0.d0/))

      deallocate (envibg3,stat=ierr)
      if (ierr.ne.0) call dealloc_error('envibg3')

      deallocate (psiga3,stat=ierr)
      if (ierr.ne.0) call dealloc_error('psiga3')

      deallocate (fe3,stat=ierr)
      if (ierr.ne.0) call dealloc_error('fe3')

      deallocate (istel,istpr,stat=ierr)
      if (ierr.ne.0) call dealloc_error('istel,istpr')

   else

      write(6,*) 'wrong mgquant value passed to feszz2:',mgquant
      write(6,*) '(from feszz2...)'
      ! clean_memory
      if (allocated(istel  )) deallocate(istel  )
      if (allocated(istpr  )) deallocate(istpr  )
      if (allocated(istga  )) deallocate(istga  )
      if (allocated(hf     )) deallocate(hf     )
      if (allocated(zzprev2)) deallocate(zzprev2)
      if (allocated(wfe    )) deallocate(wfe    )
      if (allocated(fv1    )) deallocate(fv1    )
      if (allocated(fv2    )) deallocate(fv2    )
      if (allocated(dij    )) deallocate(dij    )
      if (allocated(gij    )) deallocate(gij    )
      if (allocated(dimunu )) deallocate(dimunu )
      if (allocated(gimunu )) deallocate(gimunu )
      if (allocated(envibg3)) deallocate(envibg3)
      if (allocated(psiga3 )) deallocate(psiga3 )
      if (allocated(psiga4 )) deallocate(psiga4 )
      if (allocated(zsquare)) deallocate(zsquare)
      if (allocated(fe3    )) deallocate(fe3    )
      stop

   endif

   !write(*,*) 'leaving feszz2... any memory leaks?'
   !if (allocated(istel))   write(*,*) '*****> istel is still allocated'
   !if (allocated(istpr))   write(*,*) '*****> istpr is still allocated'
   !if (allocated(istga))   write(*,*) '*****> istga is still allocated'
   !if (allocated(u))       write(*,*) '*****> u is still allocated'
   !if (allocated(cu))      write(*,*) '*****> cu is still allocated'
   !if (allocated(vpr))     write(*,*) '*****> vpr is still allocated'
   !if (allocated(w))       write(*,*) '*****> w is still allocated'
   !if (allocated(ham))     write(*,*) '*****> ham is still allocated'
   !if (allocated(cp))      write(*,*) '*****> cp is still allocated'
   !if (allocated(vprg))    write(*,*) '*****> vprg is still allocated'
   !if (allocated(wg))      write(*,*) '*****> wg is still allocated'
   !if (allocated(hamg))    write(*,*) '*****> hamg is still allocated'
   !if (allocated(cpg))     write(*,*) '*****> cpg is still allocated'
   !if (allocated(ceprev))  write(*,*) '*****> ceprev is still allocated'
   !if (allocated(cpprev))  write(*,*) '*****> cpprev is still allocated'
   !if (allocated(hf))      write(*,*) '*****> hf is still allocated'
   !if (allocated(zzprev2)) write(*,*) '*****> zzprev2 is still allocated'
   !if (allocated(wfe))     write(*,*) '*****> wfe is still allocated'
   !if (allocated(fv1))     write(*,*) '*****> fv1 is still allocated'
   !if (allocated(fv2))     write(*,*) '*****> fv2 is still allocated'
   !if (allocated(dij))     write(*,*) '*****> dij is still allocated'
   !if (allocated(gij))     write(*,*) '*****> gij is still allocated'
   !if (allocated(dimunu))  write(*,*) '*****> dimunu is still allocated'
   !if (allocated(gimunu))  write(*,*) '*****> gimunu is still allocated'
   !if (allocated(psiga3))  write(*,*) '*****> psiga3 is still allocated'
   !if (allocated(psiga4))  write(*,*) '*****> psiga4 is still allocated'
   !if (allocated(zsquare)) write(*,*) '*****> zsquare is still allocated'
   !if (allocated(fe3))     write(*,*) '*****> fe3 is still allocated'

   ! clean_memory                                
   if (allocated(istel  )) deallocate(istel  )   
   if (allocated(istpr  )) deallocate(istpr  )   
   if (allocated(istga  )) deallocate(istga  )   
   if (allocated(hf     )) deallocate(hf     )   
   if (allocated(wfe    )) deallocate(wfe    )   
   if (allocated(fv1    )) deallocate(fv1    )   
   if (allocated(fv2    )) deallocate(fv2    )   
   if (allocated(dij    )) deallocate(dij    )   
   if (allocated(gij    )) deallocate(gij    )   
   if (allocated(dimunu )) deallocate(dimunu )   
   if (allocated(gimunu )) deallocate(gimunu )   
   if (allocated(envibg3)) deallocate(envibg3)   
   if (allocated(psiga3 )) deallocate(psiga3 )   
   if (allocated(psiga4 )) deallocate(psiga4 )   
   if (allocated(zsquare)) deallocate(zsquare)   
   if (allocated(fe3    )) deallocate(fe3    )   
   return

   contains

   !============================================================================!
   subroutine electron_proton_states
   !============================================================================!
   ! calculates electronic and protonic states on 2d grid
   !============================================================================!

      ! local variables
      integer :: ierr, kg, k, i, ishift, kk, ll, l, n
      real*8  :: p

      real*8,  allocatable :: u(:), cu(:,:)
      real*8,  allocatable :: vpr(:),  w(:),  ham(:,:),  cp(:,:)
      real*8,  allocatable :: ceprev(:,:), cpprev(:,:)

      ! allocate some work arrays

      allocate (u(nelst),cu(nelst,nelst),stat=ierr)
      if (ierr.ne.0) call alloc_error('u,cu')
      u = 0.d0
      cu = 0.d0

      allocate (ceprev(nelst,nelst),stat=ierr)
      if (ierr.ne.0) call alloc_error('ceprev')
      ceprev = 0.d0

      allocate (cpprev(npnts,npnts),stat=ierr)
      if (ierr.ne.0) call alloc_error('cpprev')
      cpprev = 0.d0

      allocate (vpr(npnts),w(npnts),ham(npnts,npnts),cp(npnts,npnts),stat=ierr)
      if (ierr.ne.0) call alloc_error('vpr,w,ham,cp')
      vpr = 0.d0
      w = 0.d0
      ham = 0.d0
      cp = 0.d0

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! loop over the grid points along the gating coordinate
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do kg=1,npntsg

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! loop over the grid points along the proton coordinate
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do k=1,npnts

            ! electronic profiles and wavefunctions
            call usol(elmode,k,kg,zp,ze,u,cu)

            if (adiab) then

               ! all four states are needed (ielst must be 4)

               do i=1,ielst
                  enel(i,k,kg) =  u(i)
               enddo

            elseif (diab2) then

               ! only two states are needed (ielst must be 2):
               ! 1a,1b if iset=1
               ! 2a,2b if iset=2

               do i=1,ielst
                  ishift = i + 2*iset - 2
                  enel(i,k,kg) = u(ishift)
               enddo

            elseif (diab4) then

               ! only one state is needed:

               enel(1,k,kg) = u(iset)

            endif

            ! check phases of the wavefunctions

            !if (k.eq.1.and.kg.eq.1) then
            !
            !   if (numr.gt.1) then
            !      do kk=1,nelst
            !         p = dot_product(cu(1:nelst,kk),zeprev2(1:nelst,kk))
            !         if (p.lt.0d0) then
            !            do ll=1,nelst
            !               cu(ll,kk) = -cu(ll,kk)
            !            enddo
            !         endif
            !      enddo
            !   endif
            !   zeprev2 = cu
            !
            !else

            if (k.gt.1.or.kg.gt.1) then

               do kk=1,nelst
                  p = dot_product(cu(1:nelst,kk),ceprev(1:nelst,kk))
                  if (p.lt.0d0) then
                     do ll=1,nelst
                        cu(ll,kk) = -cu(ll,kk)
                     enddo
                  endif
               enddo

            endif

            ceprev = cu

            ! assign the electronic wavefunction

            if (adiab) then

               do l=1,ielst
                  do i=1,ielst
                     psiel(i,k,kg,l) = cu(l,i)
                  enddo
               enddo

            elseif (diab2) then

               do l=1,ielst
                  do i=1,ielst
                     ishift = i + 2*iset - 2
                     psiel(i,k,kg,l) = cu(l,ishift)
                  enddo
               enddo

            elseif (diab4) then

               psiel(1,k,kg,1) = 1.d0

            endif

         enddo

         !===============================================================
         ! calculation of the protonic vibrational states (basis states)
         ! method = 1 - in  diabatic electronic potentials
         !        = 2 - in adiabatic electronic potentials
         !        = 3 - in adiabatic electronic potentials
         !===============================================================

         do i=1,ielst

            vpr = 0.d0
            ham = 0.d0
            cp = 0.d0
            w = 0.d0

            do k=1,npnts
               vpr(k) = enel(i,k,kg)*cal2au
            enddo
            call calcham(npnts,vpr,hke,ham)
            call calcwf(npnts,ham,cp,w)

            ! check phases

            !if (i.eq.1.and.kg.eq.1) then
            !
            !   if (numr.gt.1) then
            !      do kk=1,nprst
            !         p = dot_product(cp(1:npnts,kk),zpprev2(1:npnts,kk))
            !         if (p.lt.0d0) then
            !            do ll=1,npnts
            !               cp(ll,kk) = -cp(ll,kk)
            !            enddo
            !         endif
            !      enddo
            !   endif
            !   zpprev2(1:npnts,1:npnts) = cp
            !
            !else

            if (i.gt.1.or.kg.gt.1) then

               do kk=1,nprst
                  p = dot_product(cp(1:npnts,kk),cpprev(1:npnts,kk))
                  if (p.lt.0d0) then
                     do ll=1,npnts
                        cp(ll,kk) = -cp(ll,kk)
                     enddo
                  endif
               enddo

            endif

            cpprev = cp

            do n=1,nprst
               envib(i,n,kg) = w(n)*au2cal
               do l=1,npnts
                  psipr(i,n,l,kg) = cp(l,n)
               enddo
            enddo

         enddo   ! loop over proton grid

      enddo  ! loop over gating grid

      deallocate (u,cu,stat=ierr)
      if (ierr.ne.0) call dealloc_error('u,cu')

      deallocate (ceprev,cpprev,stat=ierr)
      if (ierr.ne.0) call dealloc_error('ceprev,cpprev')

      deallocate (vpr,ham,cp,w,stat=ierr)
      if (ierr.ne.0) call dealloc_error('vpr,ham,cp,w')

      return

   end subroutine electron_proton_states

   !============================================================================!
   subroutine gating_states1
   !============================================================================!
   ! calculation of the gating vibrational states for each
   ! basis electron-proton state
   !============================================================================!
      integer :: ierr, i, n, kg, kk, ll, m, l
      real*8  :: p

      real*8,  allocatable :: vprg(:), wg(:), hamg(:,:), cpg(:,:)
      real*8,  allocatable :: cpprev(:,:)

      allocate (vprg(npntsg),wg(npntsg),&
                hamg(npntsg,npntsg),cpg(npntsg,npntsg),stat=ierr)
      if (ierr.ne.0) call alloc_error('vprg,wg,hamg,cpg')
      vprg = 0.d0
      wg   = 0.d0
      hamg = 0.d0
      cpg  = 0.d0

      allocate (cpprev(npntsg,npntsg),stat=ierr)                              
      if (ierr.ne.0) call alloc_error('cpprev')                               
      cpprev = 0.d0                                                           

      psiga4 = 0.d0

      do i=1,ielst
         do n=1,nprst                                                         

            hamg = 0.d0                                                       
            cpg = 0.d0                                                        
            wg = 0.d0                                                         

            do kg=1,npntsg                                                    
               vprg(kg) = envib(i,n,kg)*cal2au                                
            enddo                                                             
            call calcham(npntsg,vprg,hgke,hamg)                               
            call calcwf(npntsg,hamg,cpg,wg)                                   

            ! check phases (check this very carefully!!!)                     

            !if (i.eq.1.and.n.eq.1) then                      
            !
            !   if (numr.gt.1) then                                            
            !      do kk=1,ngast                                               
            !         p = dot_product(cpg(1:npntsg,kk),zpprev2(1:npntsg,kk))   
            !         if (p.lt.0d0) then                                       
            !            do ll=1,npntsg                                        
            !               cpg(ll,kk) = -cpg(ll,kk)                           
            !            enddo                                                 
            !         endif                                                    
            !      enddo                                                       
            !   endif                                                          
            !   zpprev2(1:npntsg,1:npntsg) = cpg                               
            !
            !else                                                              

            if (i.gt.1.or.n.gt.1) then

               do kk=1,ngast                                                  
                  p = dot_product(cpg(1:npntsg,kk),cpprev(1:npntsg,kk))       
                  if (p.lt.0d0) then                                          
                     do ll=1,npntsg                                           
                        cpg(ll,kk) = -cpg(ll,kk)                              
                     enddo                                                    
                  endif                                                       
               enddo                                                          

            endif                                                             

            cpprev = cpg                                                      

            do m=1,ngast                                                      
               envibg(i,n,m) = wg(m)*au2cal                                   
               do l=1,npntsg                                                  
                  psiga4(i,n,m,l) = cpg(l,m)                                  
               enddo                                                          
            enddo                                                             

         enddo                                                                
      enddo                                                                   

      deallocate (vprg,hamg,cpg,cpprev,wg,stat=ierr)                          
      if (ierr.ne.0) call dealloc_error('vprg,hamg,cpg,cpprev,wg')                   

      return

   end subroutine gating_states1

   !============================================================================!
   subroutine gating_states3
   !============================================================================!
   ! calculation of the gating vibrational states for each
   ! electron-proton state
   !============================================================================!
      integer :: ierr, i, kg, kk, ll, m
      real*8  :: p

      real*8,  allocatable :: vprg(:), wg(:), hamg(:,:), cpg(:,:)
      real*8,  allocatable :: cpprev(:,:)

      envibg3 = 0.d0
      psiga3  = 0.d0

      allocate (vprg(npntsg),wg(npntsg),hamg(npntsg,npntsg),cpg(npntsg,npntsg),stat=ierr)   
      if (ierr.ne.0) call alloc_error('vprg,wg,hamg,cpg')                                   
      vprg = 0.d0                                                                           
      wg   = 0.d0                                                                           
      hamg = 0.d0                                                                           
      cpg  = 0.d0                                                                           

      allocate (cpprev(npntsg,npntsg),stat=ierr)                                            
      if (ierr.ne.0) call alloc_error('cpprev')                                             
      cpprev = 0.d0                                                                         

      do i=1,nstates

         hamg = 0.d0
         cpg  = 0.d0
         wg   = 0.d0

         do kg=1,npntsg
            vprg(kg) = fe3(i,kg)*cal2au
         enddo
         call calcham(npntsg,vprg,hgke,hamg)
         call calcwf(npntsg,hamg,cpg,wg)

         ! check phases

         !if (i.eq.1) then
         !
         !   if (numr.gt.1) then
         !      do kk=1,ngast
         !         p = dot_product(cpg(1:npntsg,kk),zpprev2(1:npntsg,kk))
         !         if (p.lt.0d0) then
         !            do ll=1,npntsg
         !               cpg(ll,kk) = -cpg(ll,kk)
         !            enddo
         !         endif
         !      enddo
         !   endif
         !   zpprev2(1:npntsg,1:npntsg) = cpg
         !
         !else

         if (i.gt.1) then

            do kk=1,ngast
               p = dot_product(cpg(1:npntsg,kk),cpprev(1:npntsg,kk))
               if (p.lt.0d0) then
                  do ll=1,npntsg
                     cpg(ll,kk) = -cpg(ll,kk)
                  enddo
               endif
            enddo

         endif

         cpprev = cpg

         do m=1,ngast
            envibg3(i,m) = wg(m)*au2cal
            do kg=1,npntsg
               psiga3(i,m,kg) = cpg(kg,m)
            enddo
         enddo

      enddo

      deallocate (vprg,hamg,cpg,cpprev,wg,stat=ierr)
      if (ierr.ne.0) call dealloc_error('vprg,hamg,cpg,cpprev,wg')

      return

   end subroutine gating_states3

   !============================================================================!
   subroutine calcdij2
   !====================================================================
   ! Calculates an array DIJ(I,J,K,KG) representing the                    
   ! the matrix elements <Psi_i|Nabla_p|Psi_j>                             
   ! depending on the proton coordinate (point K in RLIST)                 
   ! and the gating coordinate (point KG in GLIST)                         
   !                                                                       
   ! The two options are implemented:                                      
   ! 1) direct numerical calculation of <Psi_i|Nabla_p|Psi_j>              
   ! 2) calculation in a line with the exact expressions                   
   !====================================================================   
      integer :: option, kg, k, k1, k2, i, j, l, l_shift, ll, ll_shift
      real*8  :: step, d, psi, dpsi, tll

      dij = 0.d0
      option = 2

      select case(option)

         case(1)
         !~~~~~~~~~~~~~~~~~~~~~~~
         ! option (1) - numerical
         !~~~~~~~~~~~~~~~~~~~~~~~

         step = rlist(2) - rlist(1)

         do kg=1,npntsg
            do k=2,npnts-1

               k1 = k - 1
               k2 = k + 1

               do i=1,ielst
                  do j=i,ielst

                     d = 0.d0

                     if (i.ne.j) then

                        do l=1,ielst
                           psi = psiel(i,k,kg,l)
                           dpsi = psiel(j,k2,kg,l)-psiel(j,k1,kg,l)
                           d = d + psi*dpsi
                        enddo

                        d = d/(2.d0*step)

                     endif

                     dij(i,j,k,kg) = d
                     if (i.ne.j) dij(j,i,k,kg) = -d

                  enddo
               enddo

            enddo
         enddo

         case(2)
         !~~~~~~~~~~~~~~~~~~~~~~~~
         ! option (2) - analytical
         !~~~~~~~~~~~~~~~~~~~~~~~~

         do kg=1,npntsg
            do k=1,npnts

               do i=1,ielst
                  do j=i,ielst

                     d = 0.d0

                     if (i.ne.j) then

                        do l=1,ielst
                           l_shift = l + 2*iset - 2
                           do ll=1,ielst
                              ll_shift = ll + 2*iset - 2
                              if (l.eq.ll) then
                                 tll = dh0(l_shift,l_shift,k,kg) - dtinf(l_shift,k,kg)/2.d0
                              else
                                 tll = dh0(l_shift,ll_shift,k,kg)
                              endif
                              d = d + psiel(i,k,kg,l)*psiel(j,k,kg,ll)*tll
                           enddo
                        enddo

                        d = d/(enel(j,k,kg)-enel(i,k,kg))/cal2au

                     endif

                     dij(i,j,k,kg) = d
                     if (i.ne.j) dij(j,i,k,kg) = -d

                  enddo

               enddo

            enddo
         enddo

      end select

      return

   end subroutine calcdij2

   !============================================================================!
   subroutine calcgij2
   !========================================================================
   !     Calculates an array GIJ(I,J,K,KG) representing the
   !     the matrix elements <Psi_i|Nabla_p^2|Psi_j>
   !     depending on the proton coordinate (point K in PRLIST)
   !     and the gating distance (point KG in GLIST)
   !
   !     The two options are implemented:
   !     1) direct numerical calculation of <Psi_i|Nabla_p^2|Psi_j>
   !     2) calculation in a line with the exact expressions
   !========================================================================
      integer :: option, kg, k, k1, k2, k3, i, j, l, kk, ll, kkshift, llshift
      real*8  :: step, g, psi, d2psi, ds1, t1kl, tr1l, tr1k
      real*8  :: g1, ejj, tkl, tkl1, tkl2, g2, g3, gik3
      real*8, dimension(2) :: ym

      gij = 0.d0
      option = 2

      select case(option)

         case(1)
         !---------------------------------------------------------------
         ! option (1): numerical
         !---------------------------------------------------------------

         step = rlist(2) - rlist(1)

         do kg=1,npntsg
            do k=2,npnts-1
               k1 = k - 1
               k2 = k + 1
               do i=1,ielst
                  do j=1,ielst
                     g = 0.d0
                     do l=1,ielst
                        psi = psiel(i,k,kg,l)
                        d2psi = psiel(j,k2,kg,l)-2.d0*psiel(j,k,kg,l)+psiel(j,k1,kg,l)
                        g = g + psi*d2psi
                     enddo
                     g = g/(step*step)
                     gij(i,j,k,kg) = g
                  enddo
               enddo
            enddo
         enddo

         case(2)
         !---------------------------------------------------------------
         ! option (2): analytical
         !---------------------------------------------------------------

         do kg=1,npntsg
            do k=1,npnts

               ! derivatives of the self-energy term

               ym(1) = zp*cal2au
               ym(2) = ze*cal2au
               ds1 = 0.d0
               do kk=1,2
                  do ll=1,2
                     t1kl = t1(kk,ll,k,kg)/cal2au
                     tr1l = tr(1,ll+1,k,kg)*cal2au
                     tr1k = tr(1,kk+1,k,kg)*cal2au
                     ds1 = ds1 + dtr(kk+1,k,kg)*t1kl*(ym(ll)+tr1l) + (ym(kk)+tr1k)*t1kl*dtr(ll+1,k,kg)
                  enddo
               enddo
               ds1 = 0.5d0*ds1

               !~~~~~~~~~~~~~~~~~~~~~~
               ! analitical expression
               !~~~~~~~~~~~~~~~~~~~~~~

               do i=1,ielst
                  do j=1,ielst

                     g = 0.d0

                     if (i.eq.j) then

                        do kk=1,ielst
                           if (kk.ne.i) then
                              g = g - dij(i,kk,k,kg)*dij(i,kk,k,kg)
                           endif
                        enddo

                     else

                        g1 = 0.d0
                        ejj = 0.d0

                        do kk=1,ielst
                           kkshift = kk + 2*iset - 2
                           do ll=1,ielst
                              llshift = ll + 2*iset - 2
                              if (kk.eq.ll) then
                                 tkl1 =  dh0(kkshift,kkshift,k,kg) -  dtinf(kkshift,k,kg)/2.d0
                                 tkl2 = d2h0(kkshift,kkshift,k,kg) - d2tinf(kkshift,k,kg)/2.d0
                              else
                                 tkl1 =  dh0(kkshift,llshift,k,kg)
                                 tkl2 = d2h0(kkshift,llshift,k,kg)
                              endif
                              g1 = g1 + psiel(i,k,kg,kk)*psiel(j,k,kg,ll)*tkl2
                              ejj = ejj + psiel(j,k,kg,kk)*psiel(j,k,kg,ll)*tkl1
                           enddo
                        enddo

                        g1 = g1/(enel(j,k,kg)-enel(i,k,kg))/cal2au
                        ejj = ejj + ds1 - dt(1,k,kg)/2.d0
                        ejj = ejj/(enel(j,k,kg)-enel(i,k,kg))/cal2au
                        g3 = ejj*dij(i,j,k,kg)

                        g2 = 0.d0

                        do k3=1,ielst
                           if (k3.ne.j) then
                              gik3 = 0.d0
                              do kk=1,ielst
                                 kkshift = kk + 2*iset - 2
                                 do ll=1,ielst
                                    llshift = ll + 2*iset - 2
                                    if (kk.eq.ll) then
                                       tkl = dh0(kkshift,kkshift,k,kg) - dtinf(kkshift,k,kg)/2.d0
                                    else
                                       tkl = dh0(kkshift,llshift,k,kg)
                                    endif
                                    gik3 = gik3 + psiel(i,k,kg,kk)*psiel(k3,k,kg,ll)*tkl
                                 enddo
                              enddo
                              if (i.eq.k3) gik3=gik3+ds1-dt(1,k,kg)/2.d0
                              g2 = g2 + gik3*dij(k3,j,k,kg)/(enel(j,k,kg)-enel(i,k,kg))/cal2au
                           endif
                        enddo

                        g = g1 + 2.d0*g2 - 2.d0*g3

                     endif

                     gij(i,j,k,kg) = g

                  enddo
               enddo
            enddo
         enddo

      end select

      return

   end subroutine calcgij2

   !============================================================================!
   subroutine calcdimn
   !========================================================================
   ! Calculates an array DIMUNU(I,MU,NU,KG) representing the
   ! the matrix elements <Chi_imu|Nabla_R|Chi_inu>
   ! depending on the gating coordinate (point KG in GLIST)
   !
   ! The two options are implemented:
   ! 1) direct numerical calculation
   ! 2) calculation in a line with the exact expressions
   !========================================================================
      integer :: option, i, mu, nu, kg, kg1, kg2, k
      integer :: ishift
      real*8  :: step, r2step, d, chimu, dchinu
      real*8  :: chinu, dhii

      dimunu = 0.d0
      option = 2

      select case(option)

         case(1)
         !~~~~~~~~~~~~~~~~~~~~~~~
         ! option (1) - numerical
         !~~~~~~~~~~~~~~~~~~~~~~~

         step = glist(2) - glist(1)
         r2step = 0.5d0/step

         do i=1,ielst

            do mu=1,nprst-1
               do nu=mu+1,nprst

                  do kg=2,npntsg-1

                     kg1 = kg - 1
                     kg2 = kg + 1

                     d = 0.d0
                     do k=1,npnts
                        chimu = psipr(i,mu,k,kg)
                        dchinu = psipr(i,nu,k,kg2) - psipr(i,nu,k,kg1)
                        d = d + chimu*dchinu
                     enddo
                     d = d*r2step

                     dimunu(i,mu,nu,kg) =  d
                     dimunu(i,nu,mu,kg) = -d

                  enddo

               enddo
            enddo

         enddo

         case(2)
         !~~~~~~~~~~~~~~~~~~~~~~~~
         ! option (2) - analytical
         !~~~~~~~~~~~~~~~~~~~~~~~~

         do i=1,ielst

            ishift = i + 2*iset - 2

            do mu=1,nprst-1
               do nu=mu+1,nprst

                  do kg=1,npntsg

                     d = 0.d0
                     do k=1,npnts
                        chimu = psipr(i,mu,k,kg)
                        chinu = psipr(i,nu,k,kg)
                        dhii  = dgh0(ishift,ishift,k,kg)
                        d = d + chimu*dhii*chinu
                     enddo
                     d = d/(envib(i,nu,kg)-envib(i,mu,kg))/cal2au

                     dimunu(i,mu,nu,kg) =  d
                     dimunu(i,nu,mu,kg) = -d

                  enddo

               enddo
            enddo

         enddo

      end select

      return

   end subroutine calcdimn

   !============================================================================!
   subroutine calcgimn
   !========================================================================
   !     Calculates an array GIJ(I,MU,NU,KG) representing the
   !     the matrix elements <Chi_imu|Nabla_R^2|Chi_inu>
   !     depending on the gating coordinate (point KG in GLIST)
   !
   !     The two options are implemented:
   !     1) direct numerical calculation
   !     2) calculation in a line with the exact expressions
   !========================================================================
      integer :: option, i, mu, nu, kg, kg1, kg2, k
      real*8  :: step, rstepq, g, chimu, d2chinu, d

      option = 1

      select case(option)

         case(1)
         !~~~~~~~~~~~~~~~~~~~~~~~
         ! option (1) - numerical
         !~~~~~~~~~~~~~~~~~~~~~~~

         step  = glist(2) - glist(1)
         rstepq = 1.d0/(step*step)

         do i=1,ielst

            do mu=1,nprst
               do nu=1,nprst

                  do kg=2,npntsg-1

                     kg1 = kg - 1
                     kg2 = kg + 1

                     g = 0.d0
                     do k=1,npnts
                        chimu = psipr(i,mu,k,kg)
                        d2chinu = psipr(i,nu,k,kg2) - 2.d0*psipr(i,nu,k,kg) + psipr(i,nu,k,kg1)
                        d = d + chimu*d2chinu
                     enddo
                     g = g*rstepq

                     gimunu(i,mu,nu,kg) =  g

                  enddo

               enddo
            enddo

         enddo

         case(2)
         !~~~~~~~~~~~~~~~~~~~~~~~~
         ! option (2) - analytical
         !~~~~~~~~~~~~~~~~~~~~~~~~

         write(*,*) 'sorry, analytical gimunu is not coded yet'
         write(*,*) 'anyway, it must be slower than numerical...'
         write(*,*) 'please, change option in calcgimn.'
         stop

      end select

      return

   end subroutine calcgimn

   !============================================================================!
   subroutine alloc_error(string_)
      character(*), intent(in) :: string_
      write(*,*) '======> allocation error: '//string_
      stop
   end subroutine alloc_error
      
   !============================================================================!
   subroutine dealloc_error(string_)
      character(*), intent(in) :: string_
      write(*,*) '======> deallocation error: '//string_
      stop
   end subroutine dealloc_error

   !============================================================================!

   end subroutine feszz2

   !============================================================================!
   subroutine evbweig(mode,iset,nstates,ndabf,nz,z,&
                     &ielst,psiel,psipr,npsiga,psiga,wght)
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
   !============================================================================!
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

end module feszz_2d
