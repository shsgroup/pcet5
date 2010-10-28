module elcm

!---------------------------------------------------------------------
! Contains the routines related to the ellipsoidal cavity model
!---------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!---------------------------------------------------------------------

   use pardim
   use cst
   use parsol
   use geosol

   !---------------------------------------------------------------------
   implicit none
   private

   !---------------------------------------------------------------------
   public :: tmat, dtm, dtinfm, d2tinfm

   !---------------------------------------------------------------------
   contains

   !---------------------------------------------------------------------
   subroutine tmat(t,tinf,tr,trinf)
   !======================================================================C
   !  Calculates inertial and non-inertial reorganization
   !  energy matrices at the current configuration of the solute.
   !
   !  The reduced matrices (TR,TRINF) are defined as
   !
   !     t'_{11} = t_{11}
   !     t'_{1i} = t_{1i} - t_{11}                    (i .ne. 1)
   !     t'_{ij} = t_{ij} + t_{11} - t_{i1} - t_{j1}  (i,j .ne. 1)
   !
   !  and
   !
   !     t_{ij} = T_{ii,jj} = -INT{\phi_{ii} \rho_{jj}}
   !
   !  The simple electrostatic ellipsoidal model is utilized
   !  for calculation of polarization potentials \phi_{ii}
   !-----------------------------------------------------------------------
      implicit none
      real*8, intent(out), dimension(4,4) :: t, tinf, tr, trinf

      integer :: i, j, l
      real*8  :: tij0, tij8, xl, xlal, xmul
      real*8  :: phi0j, phi8j, t00, t008

      do i=1,4
         do j=1,4
            tr   (i,j) = zero
            trinf(i,j) = zero
            t    (i,j) = zero
            tinf (i,j) = zero
         enddo
      enddo

      do i=1,4
         do j=1,i

            tij0 = zero
            tij8 = zero

            do l=1,natsol

               xl = xyzsol(1,l)
               xlal = one
               xmul = two*xl/r
               call polpot(j,eps0,xlal,xmul,phi0j)
               call polpot(j,eps8,xlal,xmul,phi8j)
               tij0 = tij0 - chrsol(i,l)*(phi0j - phi8j)
               tij8 = tij8 - chrsol(i,l)*phi8j

            enddo

            t(i,j) = tij0
            if (i.ne.j) t(j,i) = tij0
            tinf(i,j) = tij8
            if (i.ne.j) tinf(j,i) = tij8

         enddo
      enddo

      ! (HYD) 19 July 2001: new option to take care of solvation matrices
      ! (HYD) now taking care of solvation matrix options
      ! all manipulations done using one triangle of matrices

      ! SYMT option

      if (symt) then

          ! first two line to ensure property
          ! amongst matrix elements respected

          t(2,3) = t(1,2) + t(1,3) - t(2,2)
          t(1,4) = t(1,2) + t(1,3) - t(1,1)
          tinf(2,3) = tinf(1,2) + tinf(1,3) - tinf(2,2)
          tinf(1,4) = tinf(1,2) + tinf(1,3) - tinf(1,1)

          ! now symmetrizing matrix for a symmetric system

          t(4,4) = t(1,1)
          t(3,3) = t(2,2)
          t(3,4) = t(1,2)
          t(2,4) = t(1,3)
          tinf(4,4) = tinf(1,1)
          tinf(3,3) = tinf(2,2)
          tinf(3,4) = tinf(1,2)
          tinf(2,4) = tinf(1,3)

       endif

      ! SYMPT option

      if (sympt) then
         t(2,2) = t(1,1)
         tinf(2,2) = tinf(1,1)
      endif

      ! SYMET option

      if (symet) then
         t(3,3) = t(1,1)
         tinf(3,3) = tinf(1,1)
      endif

      ! (HYD) now taking care of removing 2b dependence
      !       using a set of reduced densities

      if (reddens) then
         t(1,4) = t(1,2) + t(1,3) - t(1,1)
         t(2,4) = t(2,2) + t(2,3) - t(1,2)
         t(3,4) = t(2,3) + t(3,3) - t(1,3)
         t(4,4) = t(1,1) + t(2,2) + t(3,3) + 2.d0*(t(2,3)-t(1,2)-t(1,3))
         tinf(1,4) = tinf(1,2) + tinf(1,3) - tinf(1,1)
         tinf(2,4) = tinf(2,2) + tinf(2,3) - tinf(1,2)
         tinf(3,4) = tinf(2,3) + tinf(3,3) - tinf(1,3)
         tinf(4,4) = tinf(1,1) + tinf(2,2) + tinf(3,3) + 2.d0*(tinf(2,3)-tinf(1,2)-tinf(1,3))
      endif

      ! (HYD) end solvation matrix options

      !(HYD) now symmetrizing relative to diagonal
      !      must always be done else whole approach not validated
      !      nosymd should be used to check only to see how valid approach is

      if (.not.nosymd) then
         do i=1,4
            if ((i+1).lt.4) then
               do j=i+1,4
                  t(j,i) = t(i,j)
                  tinf(j,i) = tinf(i,j)
               enddo
            endif
         enddo
      endif

      ! (end HYD section)

      ! Construction of the reduced matrices [t']

      t00  = t(1,1)
      t008 = tinf(1,1)
      tr(1,1) = t00
      trinf(1,1) = t008

      do i=2,4
         tr(i,1) = t(i,1) - t00
         tr(1,i) = tr(i,1)
         trinf(i,1) = tinf(i,1) - t008
         trinf(1,i) = trinf(i,1)
      enddo

      do i=2,4
         do j=2,4
            tr(i,j) = t(i,j) + t00 - t(i,1) - t(j,1)
            trinf(i,j) = tinf(i,j) + t008 - tinf(i,1) - tinf(j,1)
         enddo
      enddo

      ! conversion a.u. --> kcal/mol

      do i=1,4
         do j=1,4
            t(i,j) = t(i,j)*au2cal
            tinf(i,j) = tinf(i,j)*au2cal
            tr(i,j) = tr(i,j)*au2cal
            trinf(i,j) = trinf(i,j)*au2cal
         enddo
      enddo

      return

   end subroutine tmat

   !---------------------------------------------------------------------
   subroutine dtm(dt)
   !======================================================================C
   !     Calculates the derivatives (au) of the inertial
   !     reorganization energy matrix elements at the current
   !     configuration of the solute
   !
   !     DT(I) = DT(1,I)
   !
   !     The simple electrostatic ellipsoidal model is utilized
   !     for calculation of polarization potentials \phi_{ii}
   !
   !======================================================================C
      implicit none
      real*8, intent(out), dimension(4) :: dt

      logical :: prot
      integer :: i, l
      real*8  :: dti, xcl, xlal, xmul, dphi0, dphi8

      dt = zero

      do i=1,4

         dti = zero

         do l=1,natsol

            if (l.eq.iptsol(2)) then
               prot=.true.
            else
               prot=.false.
            endif

            xcl = xyzsol(1,l)
            xlal = one
            xmul = two*xcl/r
            call dpolpot(prot,i,eps0,r,xlal,xmul,dphi0)
            call dpolpot(prot,i,eps8,r,xlal,xmul,dphi8)
            dti = dti - chrsol(1,l)*(dphi0-dphi8)

         enddo

         dt(i) = dti

      enddo

      return

   end subroutine dtm

   !---------------------------------------------------------------------
   subroutine dtinfm(dtinfk)
   !=====================================================================
   !  Calculates the derivatives (au) of the non-inertial reorganization
   !  energy matrix elements at the current configuration of the solute.
   !
   !  The simple electrostatic ellipsoidal model is utilized
   !  for calculation of polarization potentials \phi_{ii}
   !
   !=====================================================================
      implicit none
      real*8, intent(out), dimension(4) :: dtinfk

      logical :: proton
      integer :: i, l
      real*8  :: dtii8, xcl, xlal, xmul, dphi8


      dtinfk = zero

      do i=1,4

         dtii8 = zero

         do l=1,natsol

            if (l.eq.iptsol(2)) then
               proton=.true.
            else
               proton=.false.
            endif

            xcl = xyzsol(1,l)
            xlal = one
            xmul = two*xcl/r
            call dpolpot(proton,i,eps8,r,xlal,xmul,dphi8)
            dtii8 = dtii8 - chrsol(i,l)*dphi8

         enddo

         dtinfk(i) = dtii8

      enddo

      return

   end subroutine dtinfm

   !---------------------------------------------------------------------
   subroutine d2tinfm(d2tinf)
   !=====================================================================
   !  Calculates the second derivatives (au) of the non-inertial
   !  reorganization energy matrix elements at the current
   !  configuration of the solute.
   !
   !  The simple electrostatic ellipsoidal model is utilized
   !  for calculation of polarization potentials \phi_{ii}
   !
   !=====================================================================
      implicit none
      real*8, intent(out), dimension(4) :: d2tinf

      logical :: prot
      integer :: i, l
      real*8  :: d2tii8, xcl, xlal, xmul, d2phi8

      d2tinf = zero

      do i=1,4

         d2tii8 = zero

         do l=1,natsol

            if (l.eq.iptsol(2)) then
               prot=.true.
            else
               prot=.false.
            endif

            xcl = xyzsol(1,l)
            xlal = one
            xmul = two*xcl/r
            call d2polpot(prot,i,eps8,r,xlal,xmul,d2phi8)
            d2tii8 = d2tii8 - chrsol(i,l)*d2phi8

         enddo

         d2tinf(i) = d2tii8

      enddo

      return

   end subroutine d2tinfm

   !---------------------------------------------------------------------
   subroutine bcoef(n,j,eps,bn)
   !=====================================================================
   !  Calculates the expansion coefficient B(N)
   !  for the polarization potential at the point
   !  (XLA,XMU) inside the ellipsoidal cavity
   !  with the charge distribution (J) from the surrounding
   !  medium with dielectric constant EPS
   !=====================================================================
      implicit none
      integer, intent(in)  :: n, j
      real*8,  intent(in)  :: eps
      real*8,  intent(out) :: bn
      
      integer :: k
      real*8  :: pn, pn1, qn, qn1, xj, xa, xck, xk, pnk, rau

      pn  = plgndr(n,  l0)
      pn1 = plgndr(n-1,l0)
      qn  = qlgndr(n,  l0)
      qn1 = qlgndr(n-1,l0)

      if (n.eq.0) then
         xj = 1.d0
      else
         xj = one - (l0 - pn1/pn)/(l0 - qn1/qn)/eps
      endif

      xa = zero

      do k=1,natsol
         xck = xyzsol(1,k)
         xk = two*xck/r
         pnk = plgndr(n,xk)
         xa = xa + chrsol(j,k)*pnk
      enddo

      rau = r*a2bohr
      bn = (two/rau)*(one/eps-one)*(two*n+one)*xa*qn/xj/pn

      return

   end subroutine bcoef

   !---------------------------------------------------------------------
   subroutine dbcoef(n,j,eps,dbn)
   !=====================================================================
   !  Calculates the derivative of the expansion coefficient B(N)
   !  for the polarization potential at point
   !  (XLA,XMU) inside the ellipsoidal cavity
   !  with charge distribution (J) from the surrounding
   !  medium with dielectric constant EPS
   !  (with respect to the proton ellipsoidal coordinate)
   !=====================================================================
      implicit none
      
      integer, intent(in)  :: n, j
      real*8,  intent(in)  :: eps
      real*8,  intent(out) :: dbn

      integer :: numh
      real*8  :: pn, pn1, qn, qn1, xj, xch, xmuh
      real*8  :: pnh, pn1h, dpn, rau, chrjh

      pn  = plgndr(n,  l0)
      pn1 = plgndr(n-1,l0)
      qn  = qlgndr(n,  l0)
      qn1 = qlgndr(n-1,l0)

      if (n.eq.0) then
         xj = 1.d0
      else
         xj = one - (l0 - pn1/pn)/(l0 - qn1/qn)/eps
      endif

      numh = iptsol(2)
      xch = xyzsol(1,numh)

      xmuh = two*xch/r
      pnh = plgndr(n,xmuh)
      pn1h = plgndr(n-1,xmuh)
      dpn = n*(xmuh*pnh-pn1h)/(xmuh*xmuh-1.d0)

      rau = r*a2bohr
      chrjh = chrsol(j,numh)
      dbn = (two/rau)*(one/eps-one)*(two*n+one)*chrjh*dpn*qn/xj/pn

      return

   end subroutine dbcoef

   !---------------------------------------------------------------------
   subroutine d2bcoef(n,j,eps,d2bn)
   !======================================================================C
   !     Calculates the second derivative of the expansion
   !     coefficient B(N) for the polarization potential at point
   !     (XLA,XMU) inside the ellipsoidal cavity
   !     with charge distribution (J) from the surrounding
   !     medium with dielectric constant EPS
   !     (with respect to the proton ellipsoidal coordinate)
   !======================================================================C
      implicit none

      integer, intent(in)  :: n, j
      real*8,  intent(in)  :: eps
      real*8,  intent(out) :: d2bn

      integer :: numh
      real*8  :: pn, pn1, qn, qn1, xj, xch, xmuh
      real*8  :: pnh, pn1h, pn2h, d2pn, rau, chrjh

      pn  = plgndr(n,  l0)
      pn1 = plgndr(n-1,l0)
      qn  = qlgndr(n,  l0)
      qn1 = qlgndr(n-1,l0)

      if (n.eq.0) then
         xj = 1.d0
      else
         xj = one - (l0 - pn1/pn)/(l0 - qn1/qn)/eps
      endif

      numh = iptsol(2)
      xch = xyzsol(1,numh)

      xmuh = two*xch/r
      pnh = plgndr(n,xmuh)
      pn1h = plgndr(n-1,xmuh)
      pn2h = plgndr(n-2,xmuh)
      d2pn = (n-1)*pn2h+(3-2*n)*xmuh*pn1h+(n-1)*xmuh**2*pnh-pnh
      d2pn = d2pn/(xmuh**2.d0-1.d0)**2.d0

      rau = r*a2bohr
      chrjh = chrsol(j,numh)
      d2bn = (two/rau)*(one/eps-one)*n*(two*n+one)*chrjh*d2pn*qn/xj/pn

      return

   end subroutine d2bcoef

   !---------------------------------------------------------------------
   subroutine polpot(j,eps,xla,xmu,phi)
   !======================================================================C
   !     Calculates polarization potential PHI (a.u.) at point
   !     (XLA,XMU) inside the ellipsoidal cavity
   !     with charge distribution (J) from the surrounding
   !     medium with dielectric constant EPS
   !======================================================================C
      implicit none

      integer, intent(in)  :: j
      real*8,  intent(in)  :: eps, xla, xmu
      real*8,  intent(out) :: phi

      real*8,  parameter :: acc = 1.d-9
      integer, parameter :: nmax = 200

      integer :: n
      real*8  :: phi0, bn, pnla, pnmu

      phi = zero
      phi0 = -999.d0
      n = -1
      do while (dabs(phi-phi0).ge.acc.and.n.lt.nmax)
         n = n + 1
         if (mod(n,2).ne.0) phi0 = phi
         call bcoef(n,j,eps,bn)
         pnla = plgndr(n,xla)
         pnmu = plgndr(n,xmu)
         phi = phi + bn*pnla*pnmu
      enddo

      if (n.ge.nmax) then
         write(*,'(/1x,''no convergence in polpot''/)')
         stop 'in polpot...'
      endif

      return

   end subroutine polpot

   !---------------------------------------------------------------------
   subroutine dpolpot(proton,j,eps,r,xla,xmu,dphi)
   !======================================================================C
   !     Calculates the derivative of the polarization potential PHI
   !     at (XLA,XMU) inside the ellipsoidal cavity
   !     with charge distribution (J) from the surrounding
   !     medium with dielectric constant EPS
   !     (with respest to the proton position)
   !======================================================================C
      implicit none

      logical, intent(in)  :: proton
      integer, intent(in)  :: j
      real*8,  intent(in)  :: eps, r, xla, xmu
      real*8,  intent(out) :: dphi

      real*8,  parameter :: acc = 1.d-9
      integer, parameter :: nmax = 50
      
      integer :: n
      real*8  :: dphi0, dbn, pnla, pnmu, bn, pnmu1, rau

      dphi = zero
      dphi0 = -999.d0
      n = -1

      do while (dabs(dphi-dphi0).ge.acc.and.n.lt.nmax)

         n = n + 1
         if (mod(n,2).ne.0) dphi0 = dphi
         call dbcoef(n,j,eps,dbn)
         pnla = plgndr(n,xla)
         pnmu = plgndr(n,xmu)
         if (.not.proton) then
            dphi = dphi + dbn*pnla*pnmu
         else
            call bcoef(n,j,eps,bn)
            pnmu1 = plgndr(n-1,xmu)
            dphi = dphi + dbn*pnla*pnmu + n*bn*(xmu*pnmu-pnmu1)/(xmu**2-1)
         endif

      enddo

      rau = r*a2bohr
      dphi = 2.d0*dphi/rau

      if (n.ge.nmax) then
         write(*,'(/1x,''no convergence in dpolpot''/)')
         stop 'in dpolpot...'
      endif

      return

   end subroutine dpolpot

   !---------------------------------------------------------------------
   subroutine d2polpot(prot,j,eps,r,xla,xmu,d2phi)
   !=====================================================================
   !  Calculates the second derivative of the polarization
   !  potential PHI at (XLA,XMU) inside the ellipsoidal cavity
   !  with charge distribution (J) from the surrounding
   !  medium with dielectric constant EPS
   !  (with respest to the proton position)
   !=====================================================================
      implicit none
      logical, intent(in)  :: prot
      integer, intent(in)  :: j
      real*8,  intent(in)  :: eps, r, xla, xmu
      real*8,  intent(out) :: d2phi

      real*8,  parameter :: acc = 1.d-9
      integer, parameter :: nmax = 50

      integer :: n
      real*8  :: d2phi0, d2bn, pnla, pnmu
      real*8  :: dbn, bn, pnmu1, pnmu2, dpn, d2pn, rau

      d2phi = zero
      d2phi0 = -999.d0
      n = -1

      do while (dabs(d2phi-d2phi0).ge.acc.and.n.lt.nmax)

         n = n + 1
         if (mod(n,2).ne.0) d2phi0 = d2phi
         call d2bcoef(n,j,eps,d2bn)
         pnla = plgndr(n,xla)
         pnmu = plgndr(n,xmu)
         if (.not.prot) then
            d2phi = d2phi + d2bn*pnla*pnmu
         else
            call dbcoef(n,j,eps,dbn)
            call  bcoef(n,j,eps,bn)
            pnmu1 = plgndr(n-1,xmu)
            pnmu2 = plgndr(n-2,xmu)
            dpn = n*(xmu*pnmu-pnmu1)/(xmu*xmu-1.d0)
            d2pn = (n-1)*pnmu2+(3-2*n)*xmu*pnmu1+(n-1)*xmu*xmu*pnmu-pnmu
            d2pn = d2pn/(xmu*xmu-1.d0)**2.d0
            d2phi = d2phi + d2bn*pnla*pnmu + 2.d0*dbn*dpn + bn*d2pn
         endif

      enddo

      rau = r*a2bohr
      d2phi = 4.d0*d2phi/(rau*rau)

      if (n.ge.nmax) then
         write(*,'(/1x,''no convergence in d2polpot''/)')
         stop 'in d2polpot: no convergence...'
      endif

      return

   end subroutine d2polpot

   !---------------------------------------------------------------------
   real*8 function plgndr(l,x)
   !---------------------------------------------------------------------
   ! (C) Copr. 1986-92 Numerical Recipes Software $$$.
   !---------------------------------------------------------------------

      implicit none
      integer, intent(in)  :: l
      real*8,  intent(in)  :: x

      integer :: ll
      real*8  :: pll, pmm, pmmp1

      pmm = 1.d0

      if (l.le.0) then

         plgndr = pmm

      else

         pmmp1 = x*pmm

         if (l.eq.1) then
            plgndr = pmmp1
         else
            do ll=2,l
               pll = (x*dble(2*ll-1)*pmmp1-dble(ll-1)*pmm)/dble(ll)
               pmm = pmmp1
               pmmp1 = pll
            enddo
            plgndr = pll
         endif

      endif

      return

   end function plgndr

   !---------------------------------------------------------------------
   real*8 function qlgndr(l,x)
   !---------------------------------------------------------------------
   !  (C) Copr. 1986-92 Numerical Recipes Software $$$.
   !---------------------------------------------------------------------
      implicit none
      integer, intent(in)  :: l
      real*8,  intent(in)  :: x

      integer :: ll
      real*8 ::  qll, qmm, qmmp1

      qmm = 0.5d0*dlog((x+1.d0)/(x-1.d0))

      if (l.le.0) then

         qlgndr = qmm

      else

         qmmp1 = -1.d0 + x*qmm

         if (l.eq.1) then
            qlgndr = qmmp1
         else
            do ll=2,l
               qll = (x*dble(2*ll-1)*qmmp1-dble(ll-1)*qmm)/dble(ll)
               qmm = qmmp1
               qmmp1 = qll
            enddo
            qlgndr = qll
         endif

      endif

      return

   end function qlgndr
   !---------------------------------------------------------------------

end module elcm
