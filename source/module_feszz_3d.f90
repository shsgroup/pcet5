module feszz_3d
!======================================================================
!  3D free energy surface
!======================================================================
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_feszz_3d.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_feszz_3d.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 19:58:33  souda
!  Initial PCET-4.0 Release
!
!======================================================================
   use pardim
   use cst
   use control
   use gasmat
   use solmat
   use quantum
   use eispack, only: rs

   implicit none
   private
   
   public :: feszz3
   public :: evbwei
   
   contains

   !======================================================================
   subroutine feszz3(mode,iset,kg,zp,ze,nstates,fe,nz,z,ndabf,&
                    &ielst,enel,envib,psiel,psipr)
   !======================================================================
   !
   !     Calculates free energy as a function of two solvent variables
   !     and a gating coordinate
   !     for proton-coupled electron transfer reactions.
   !     The continuum model is utilized. The electron and
   !     the proton are treated quantum-mechanically.
   !
   !     MODE - calculation mode:
   !         ='ADIAB' - adiabatic electron/proton vibrational
   !                    free energies;
   !         ='DIAB2' - ET diabatic free energy surfaces;
   !         ='DIAB4' - diabatic (EVB) free energy surfaces.
   !
   !     ISET - set of states to calculate (only one set is available
   !            for ADIAB, two sets for DIAB2 and four sets for DIAB4)
   !
   !     KG - grid point along the gating coordinate
   !
   !     ZP, ZE are values of the effective medium coordinates (kcal/mol)
   !
   !     NSTATES - number of states to calculate
   !
   !     FE(I) - free energy of states (kcal/mol)
   !
   !     NZ - dimension of the Z matrix
   !
   !     Z(I,J) - eigenvectors of the total Hamiltonian
   !
   !     NDABF - dimension of the total Hamiltonian matrix
   !
   !     IELST - number of electronic basis states
   !
   !     ENEL(I,K) - electronic energies (kcal/mol)
   !          I - electronic quantum number
   !            K - the point number along the proton coordinate
   !
   !     ENVIB(I,N) - proton vibrational levels (kcal/mol)
   !           I - electronic quantum number
   !             N - vibrational quantum number
   !
   !     PSIEL(I,K,L) - electronic wave function (coefficients)
   !           I - quantum number
   !             K - the point number along the proton coordinate
   !               L - index of the EVB basis state
   !
   !     PSIPR(I,N,K) - protonic wave function (coefficients)
   !           I - electronic quantum number
   !             N - protonic quantum number
   !               K - the point number along the proton coordinate
   !
   !======================================================================
   character(5), intent(in) :: mode
   integer,      intent(in) :: iset, kg, nstates, nz, ielst
   real(8),      intent(in) :: zp, ze

   integer, intent(out)                               :: ndabf
   real(8), intent(out), dimension(nstates)           :: fe
   real(8), intent(out), dimension(nz,nz)             :: z
   real(8), intent(out), dimension(ielst,npnts,ielst) :: psiel
   real(8), intent(out), dimension(ielst,nprst,npnts) :: psipr
   real(8), intent(out), dimension(ielst,npnts)       :: enel
   real(8), intent(out), dimension(ielst,nprst)       :: envib

   logical :: adiab, diab2, diab4
   character(5) :: elmode

   integer :: nzdim, ielst0, i, j, n, k, l, kk, ll, ierr
   integer :: iel, ipr, ielshift, jel, jpr, jelshift
   real(8) :: p, hfij, psipr_ovl, dterm, dtij, gterm, gtij

   integer, allocatable, dimension(:)     :: istel, istpr
   real(8), allocatable, dimension(:,:,:) :: dij, gij
   real(8), allocatable, dimension(:,:)   :: hf
   real(8), allocatable, dimension(:)     :: wfe, fv1, fv2

   ! the following is only for maintaining the phases
   ! of the wavefunctions along the solvent coordinates.
   ! It is not necessary unless we want to calculate
   ! non-adiabatic coupling vector for solvent dynamics
   !
   !integer, save                           :: numr   = 0
   !real*8,  save, dimension(maxsta,maxsta) :: zprev  = 0.d0
   !real*8,  save, dimension(nelst,nelst)   :: zeprev = 0.d0
   !real*8,  save, dimension(maxpnt,maxpnt) :: zpprev = 0.d0

   nzdim = nprst*ielst
   fe = 0.d0

   !numr = numr + 1

   adiab = mode.eq.'ADIAB'
   diab2 = mode.eq.'DIAB2'
   diab4 = mode.eq.'DIAB4'

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! number of electronic basis states
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
      write(6,*) ' wrong input for feszz3: ielst=',ielst
      stop
   endif

   !======================================================================
   !     calculation of the electronic states
   !======================================================================

   if (method.eq.1) then
      elmode = 'DIAB4'
   elseif (method.eq.2.or.method.eq.3) then
      if (adiab) elmode = 'ADIAB'
      if (diab2) elmode = 'DIAB2'
      if (diab4) elmode = 'DIAB4'
   else
      elmode = 'DIAB4'
   endif

   call electronic_states

   !======================================================================
   !     calculation of the protonic vibrational states (basis states)
   !     method = 1 - in  diabatic electronic potentials
   !            = 2 - in adiabatic electronic potentials
   !            = 3 - in adiabatic electronic potentials
   !======================================================================
   call proton_states

   !======================================================================
   !     building of the electronic/proton vibrational basis
   !======================================================================

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !     ordering of the basis states
   !
   !     ndabf - the size of the basis
   !     istel(k) - the quantum number of the electronic state
   !     istpr(k) - the quantum number of the protonic state
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   allocate (istel(nzdim),istpr(nzdim))

   ndabf = 0
   do i=1,ielst
      do n=1,nprst
         ndabf = ndabf + 1
         istel(ndabf) = i
         istpr(ndabf) = n
      enddo
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !     building the d_ij and g_ij matrix elements as functions
   !     of the proton coordinate (needed only for method=2)
   !
   !     dij(i,j,k), k is a point along the proton coordinate
   !     gij(i,j,k), k is a point along the proton coordinate
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   if (method.eq.2.and..not.diab4) then
      allocate (dij(ielst,ielst,npnts),gij(ielst,ielst,npnts))
      dij = 0.d0
      gij = 0.d0
      call calcdij3(2)
      call calcgij3(2)
   endif

   !======================================================================
   !     building of the total hamiltonian matrix hf(i,j)
   !======================================================================

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !     in double adiabatic approximation (method=3): all d- and g-terms
   !     are zero and the total hamiltonian matrix is diagonal.
   !     also if mode='diab4' then there is no offdiagonal matrix elements
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   if (method.eq.3.or.diab4) then
      do i=1,nstates
         iel = istel(i)
         ipr = istpr(i)
         fe(i) = envib(iel,ipr)
      enddo
      ! clean_memory
      if (allocated(dij   )) deallocate(dij   )
      if (allocated(gij   )) deallocate(gij   )
      if (allocated(istel )) deallocate(istel )
      if (allocated(istpr )) deallocate(istpr )
      return
   endif

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   !     otherwise we have to build the total hamiltonian matrix and
   !     diagonalize it
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   allocate (hf(ndabf,ndabf),wfe(ndabf))
   allocate (fv1(ndabf),fv2(ndabf))

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
            hfij = hfij + envib(iel,ipr)*cal2au
         endif

         if (method.eq.1) then

            if (iel.ne.jel) then
               do k=1,npnts
                  psipr_ovl = psipr(iel,ipr,k)*psipr(jel,jpr,k)
                  hfij = hfij + psipr_ovl*h0(ielshift,jelshift,k,kg)*cal2au
               enddo
            endif

         else

            dterm = 0.d0
            if (iel.ne.jel) then
               do k=1,npnts
                  do l=1,npnts
                     psipr_ovl = psipr(iel,ipr,k)*psipr(jel,jpr,l)
                     dterm = dterm + psipr_ovl*dij(iel,jel,k)*dx(l,k)
                  enddo
               enddo
               dtij = -hbar*hbar*dterm/pm
               hfij = hfij + dtij
            endif

            gterm = 0.d0
            do k=1,npnts
               psipr_ovl = psipr(iel,ipr,k)*psipr(jel,jpr,k)
               gterm = gterm + psipr_ovl*gij(iel,jel,k)
            enddo

            gtij = -hbar*hbar*gterm/(2.d0*pm)
            hfij = hfij + gtij

         endif

         hf(i,j) = hfij
         if (i.ne.j) hf(j,i) = hfij

      enddo
   enddo

   if (method.eq.2.and..not.diab4) then
      if (allocated(dij)) deallocate (dij)
      if (allocated(gij)) deallocate (gij)
   endif
   deallocate (istel,istpr)

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! diagonalization of the hamiltonian matrix hf(i,j)
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   fv1 = 0.d0
   fv2 = 0.d0
   z = 0.d0
   call rs(ndabf,ndabf,hf,wfe,ndabf,z,fv1,fv2,ierr)

   do i=1,nstates
      fe(i) = wfe(i)*au2cal
   enddo

   deallocate (hf,wfe)
   deallocate (fv1,fv2)

   ! check phases between the calls

   !if (numr.gt.1) then
   !   do kk=1,ndabf
   !      p = dot_product(z(1:ndabf,kk),zprev(1:ndabf,kk))
   !      if (p.lt.0d0) then
   !         do ll=1,ndabf
   !            z(ll,kk) = -z(ll,kk)
   !         enddo
   !      endif
   !   enddo
   !endif
   !zprev(1:ndabf,1:ndabf) = z(1:ndabf,1:ndabf)

   return

   contains

   !============================================================================!
   subroutine electronic_states
   !======================================================================
   !  calculation of the electronic states
   !======================================================================
      integer :: k, i, ishift, kk, ll, l
      real*8  :: p
      real*8,  allocatable :: u(:), cu(:,:), ceprev(:,:)

      enel  = 0.d0
      psiel = 0.d0

      allocate (u(nelst),cu(nelst,nelst),ceprev(nelst,nelst))

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !     loop over the grid points along the proton coordinate:
      !     calculate electronic states
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do k=1,npnts

         ! electronic profiles and wavefunctions
         call usol(elmode,k,kg,zp,ze,u,cu)

         if (adiab) then

            ! all four states are needed

            do i=1,nelst
               enel(i,k) =  u(i)
            enddo

         elseif (diab2) then

            ! only two states are needed:
            ! 1a,1b if iset=1
            ! 2a,2b if iset=2

            do i=1,2
               ishift = i + 2*iset - 2
               enel(i,k) = u(ishift)
            enddo

         elseif (diab4) then

            ! only one state is needed:
            enel(1,k) = u(iset)

         endif

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! check phases of the wavefunctions
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !if (k.eq.1) then
         !
         !   if (numr.gt.1) then
         !      do kk=1,nelst
         !         p = dot_product(cu(1:nelst,kk),zeprev(1:nelst,kk))
         !         if (p.lt.0d0) then
         !            do ll=1,nelst
         !               cu(ll,kk) = -cu(ll,kk)
         !            enddo
         !         endif
         !      enddo
         !   endif
         !   zeprev = cu

         if (k.gt.1) then

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

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! assign the electronic wavefunction
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (adiab) then

            do l=1,nelst
               do i=1,nelst
                  psiel(i,k,l) = cu(l,i)
               enddo
            enddo

         elseif(diab2) then

            do l=1,2
               do i=1,2
                  ishift = i + 2*iset - 2
                  psiel(i,k,l) = cu(l,ishift)
               enddo
            enddo

         elseif(diab4) then

            psiel(1,k,1) = 1.d0

         endif

      enddo

      deallocate (u,cu,ceprev)

   end subroutine electronic_states

   !============================================================================!
   subroutine proton_states
   !======================================================================
   !     calculation of the protonic vibrational states (basis states)       
   !     method = 1 - in  diabatic electronic potentials                     
   !            = 2 - in adiabatic electronic potentials                     
   !            = 3 - in adiabatic electronic potentials                     
   !======================================================================
      integer :: i, k, kk, ll, n, l
      real*8  :: p

      real*8, allocatable, dimension(:)   :: vpr, w
      real*8, allocatable, dimension(:,:) :: ham, cp, cpprev

      allocate (vpr(npnts),w(npnts))
      allocate (ham(npnts,npnts),cp(npnts,npnts),cpprev(npnts,npnts))

      psipr = 0.d0
      envib = 0.d0

      do i=1,ielst

         ham = 0.d0
         cp = 0.d0
         w = 0.d0
         do k=1,npnts
            vpr(k) = enel(i,k)*cal2au
         enddo
         call calcham(npnts,vpr,hke,ham)
         call calcwf(npnts,ham,cp,w)

         ! check phases

         !if (i.eq.1) then
         !
         !   if (numr.gt.1) then
         !      do kk=1,nprst
         !         p = dot_product(cp(1:npnts,kk),zpprev(1:npnts,kk))
         !         if (p.lt.0d0) then
         !            do ll=1,npnts
         !               cp(ll,kk) = -cp(ll,kk)
         !            enddo
         !         endif
         !      enddo
         !   endif
         !   zpprev(1:npnts,1:npnts) = cp

         if (i.gt.1) then

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
            envib(i,n) = w(n)*au2cal
            do l=1,npnts
               psipr(i,n,l) = cp(l,n)
            enddo
         enddo

      enddo

      deallocate (vpr,w,ham)
      deallocate (cp,cpprev)

   end subroutine proton_states

   !============================================================================!
   subroutine calcdij3(option)
   !========================================================================
   !     Calculates an array DIJ(I,J,K) representing the
   !     the matrix elements <Psi_i|Nabla_p|Psi_j>
   !     depending on the proton coordinate (point K in RLIST)
   !     and AT THE FIXED gating coordinate (point KG in GLIST)
   !
   !     The two options are implemented:
   !     1) direct numerical calculation of <Psi_i|Nabla_p|Psi_j>
   !     2) calculation in a line with the exact expressions
   !========================================================================
      integer, intent(in) :: option
      integer :: k, k1, k2, i, j, l, l_shift, ll, ll_shift
      real*8  :: step, d, psi, dpsi, tll

      dij = 0.d0
      !option = 2

      select case(option)

         case(1)
         !~~~~~~~~~~~~~~~~~~~~~~~
         ! option (1) - numerical
         !~~~~~~~~~~~~~~~~~~~~~~~

         step = rlist(2) - rlist(1)

         do k=2,npnts-1

            k1 = k - 1
            k2 = k + 1

            do i=1,ielst
               do j=i,ielst

                  d = 0.d0

                  if (i.ne.j) then

                     do l=1,ielst
                        psi = psiel(i,k,l)
                        dpsi = psiel(j,k2,l)-psiel(j,k1,l)
                        d = d + psi*dpsi
                     enddo

                     d = d/(2.d0*step)

                  endif

                  dij(i,j,k) = d
                  if (i.ne.j) dij(j,i,k) = -d

               enddo
            enddo

         enddo

         case(2)
         !~~~~~~~~~~~~~~~~~~~~~~~~
         ! option (2) - analytical
         !~~~~~~~~~~~~~~~~~~~~~~~~

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
                           d = d + psiel(i,k,l)*psiel(j,k,ll)*tll
                        enddo
                     enddo

                     d = d/(enel(j,k)-enel(i,k))/cal2au

                  endif

                  dij(i,j,k) = d
                  if (i.ne.j) dij(j,i,k) = -d

               enddo

            enddo

         enddo

      end select

      return

   end subroutine calcdij3

   !============================================================================!
   subroutine calcgij3(option)
   !========================================================================
   !     Calculates an array GIJ(I,J,K) representing the
   !     the matrix elements <Psi_i|Nabla_p^2|Psi_j>
   !     depending on the proton coordinate (point K in PRLIST)
   !     and at the fixed gating distance (point KG in GLIST)
   !
   !     The two options are implemented:
   !     1) direct numerical calculation of <Psi_i|Nabla_p^2|Psi_j>
   !     2) calculation in a line with the exact expressions
   !========================================================================
      integer, intent(in) :: option
      integer :: i, j, k, k1, k2, k3, l, kk, kkshift, ll, llshift
      real*8  :: step, g, psi, d2psi, ds1, g1, g2, g3, gik3, ejj
      real*8  :: tkl, tkl1, tkl2, t1kl, tr1l, tr1k

      real*8, dimension(2) :: ym(2)

      gij = 0.d0
      !option = 2

      select case(option)

         case(1)
         !---------------------------------------------------------------
         !     option (1): numerical
         !---------------------------------------------------------------

         step = rlist(2) - rlist(1)

         do k=2,npnts-1
            k1 = k - 1
            k2 = k + 1
            do i=1,ielst
               do j=1,ielst
                  g = 0.d0
                  do l=1,ielst
                     psi = psiel(i,k,l)
                     d2psi = psiel(j,k2,l)-2.d0*psiel(j,k,l) + psiel(j,k1,l)
                     g = g + psi*d2psi
                  enddo
                  g = g/(step*step)
                  gij(i,j,k) = g
               enddo
            enddo
         enddo

         case(2)
         !---------------------------------------------------------------
         !     option (2): analytical
         !---------------------------------------------------------------

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
                           g = g - dij(i,kk,k)*dij(i,kk,k)
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
                              tkl1 = dh0(kkshift,kkshift,k,kg) - dtinf(kkshift,k,kg)/2.d0
                              tkl2 = d2h0(kkshift,kkshift,k,kg) - d2tinf(kkshift,k,kg)/2.d0
                           else
                              tkl1 = dh0(kkshift,llshift,k,kg)
                              tkl2 = d2h0(kkshift,llshift,k,kg)
                           endif
                           g1 = g1 + psiel(i,k,kk)*psiel(j,k,ll)*tkl2
                           ejj = ejj + psiel(j,k,kk)*psiel(j,k,ll)*tkl1
                        enddo
                     enddo

                     g1 = g1/(enel(j,k)-enel(i,k))/cal2au
                     ejj = ejj + ds1 - dt(1,k,kg)/2.d0
                     ejj = ejj/(enel(j,k)-enel(i,k))/cal2au
                     g3 = ejj*dij(i,j,k)

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
                                 gik3=gik3+psiel(i,k,kk)*psiel(k3,k,ll)*tkl
                              enddo
                           enddo
                           if (i.eq.k3) gik3 = gik3 + ds1 - dt(1,k,kg)/2.d0
                           g2 = g2 + gik3*dij(k3,j,k)/(enel(j,k)-enel(i,k))/cal2au
                        endif
                     enddo

                     g = g1 + 2.d0*g2 - 2.d0*g3

                  endif

                  gij(i,j,k) = g

               enddo
            enddo
         enddo

      end select

      return

   end subroutine calcgij3

   !============================================================================!
   end subroutine feszz3

   !============================================================================!
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
   !============================================================================!
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
   !============================================================================!

end module feszz_3d
