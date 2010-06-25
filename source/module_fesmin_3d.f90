module fesmin_3d
!=======================================================================
!     Routines related to the finfing on minima
!     on 3-dimensional free-energy surfaces
!-----------------------------------------------------------------------
!
!     souda
!     2010/06/25 20:02:35
!     4.1
!     Exp
!     module_fesmin_3d.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!     module_fesmin_3d.f90,v
!     Revision 4.1  2010/06/25 20:02:35  souda
!     Release 4.1
!
!     Revision 1.6  2007/11/06 22:20:00  souda
!     new mmgen gas phase potential o-h o systems added
!
!     Revision 1.5  2007/03/12 23:08:03  souda
!     Modifications related to using LBFGS minimization method.
!
!     Revision 1.4  2004/06/21 16:39:48  souda
!     some beautification...
!
!     Revision 1.3  2004/06/14 15:58:24  souda
!     decrease amount of out put in verbose mode in rate*
!
!     Revision 1.2  2004/05/15 03:32:45  souda
!     Added Borgis-Hynes rate routine
!
!     Revision 1.1.1.1  2003/12/19 17:16:58  souda
!     Initial PCET-4.0 Release
!
!
!=======================================================================
   use pardim
   use cst
   use eispack
   use control
   use quantum
   
   implicit none
   private

   character(5), save :: mode
   integer, save      :: iset, ist

   !-- maximum number of iterations in the minimization
   integer, public, save :: maxit=100

   !-- Newton-Raphson minimization parameters
   real(8), public, save :: slim=1.d-1   ! the maximum length of the Newton-Raphson step
   real(8), public, save :: acc=1.d-6    ! minimization accuracy (gradient norm criterion)

   !-- LBFGS minimization parameters
   integer, public, save :: iprint=-1      ! output level
   integer, public, save :: mcorr=5        ! number of corrections used in the limited memory matrix
   real(8), public, save :: factr=1.0d+6   ! tolerance for the function (factor to machine precision)
   real(8), public, save :: pgtol=1.0d-6   ! tolerance for components of the projected gradient
   
   public :: fes3min
   public :: fes3gmin

   contains

   !=======================================================================
   subroutine fes3min(mode_,iset_,istate_,zp0,zps,ze0,zes,r0,&
                      ierr,zpmin,zemin,kgmin,rmin,femin,grad,hess)
   !=======================================================================
   !     Finds minima on the spesified three-dimensional
   !     free energy surface using canonical Newton-Raphson
   !     minimization algorithm for solvent coordinates and
   !     simple minimization on the grid for the gating coordinate.
   !
   !     MODE - (ADIAB, DIAB2 or DIAB4) - type of the surface
   !
   !     ISET - a set of the surfaces (1 for ADIAB, 1/2 for DIAB2,
   !            1/2/3/4 for DIAB4)
   !
   !     ISTATE - quantum number of the desired state
   !
   !     ZP0,ZE0 - initial values of solvent coordinates (kcal/mol)
   !
   !     ZPS,ZES - scaling factors for the solvent coordinates
   !
   !     R0 - initial value of gating coordinate (A)
   !
   !     IERR - error code (0 - normal)
   !
   !     ZPMIN,ZEMIN - solvent coordinates at the minimum
   !
   !     KGMIN,RMIN - grid point index and the gating coordinate
   !                  at the minimum
   !
   !     FEMIN - free energy at the minimum (kcal/mol)
   !
   !     GRAD - derivatives (gradient) at the minimum
   !
   !     HESS - second derivatives (Hessian) at the minimum
   !
   !=======================================================================

      character(5), intent(in) :: mode_
      integer,      intent(in) :: iset_, istate_
      real*8,       intent(in) :: zp0, zps, ze0, zes, r0

      integer, intent(out) :: ierr, kgmin
      real*8,  intent(out) :: zpmin, zemin, rmin, femin
      real*8,  intent(out), dimension(3)   :: grad
      real*8,  intent(out), dimension(3,3) :: hess

      integer, parameter :: nrdim = 2
      logical :: done

      integer :: i, n, kg0, kg, igdir, ichdir, iter, nf, ng, nfc
      real*8  :: delgmin, delg, ginc, ginc2, fx
      real*8  :: dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin
      real*8  :: fright, fmiddl, fleft

      real*8, dimension(nrdim)       :: x0, x, xs, g, s, w, ev
      real*8, dimension(nrdim)       :: gleft, gright
      real*8, dimension(nrdim,nrdim) :: f, d

      n = nrdim
      mode = mode_
      iset = iset_
      ist  = istate_

      ! find a gating grid index corresponding to the initial
      ! value of the gating coordinate

      delgmin = 999999.d0
      kg0 = 1
      do kg=1,npntsg
         delg = dabs(glist(kg)*bohr2a - r0)
         if (delg.le.delgmin) then
             delgmin = delg
             kg0 = kg
         endif
      enddo

      ! initial direction of the search along the gating grid

      igdir = 1
      done = .false.
      kg = kg0
      kgmin = kg0
      x(1) = zp0*zps
      x(2) = ze0*zes
      xs(1) = zps
      xs(2) = zes
      femin = 999999.d0
      ichdir = 0
      ginc = bohr2a*(glist(2)-glist(1))
      ginc2 = ginc*ginc

      do while (.not.done)

         kg = kg + igdir
         write(*,'(1x,"--> gating grid point ",i3,": ",f8.3)') kg,glist(kg)*bohr2a

         if (kg.lt.1.or.kg.gt.npntsg) then
            write(*,'(/1x,"warning: minimization failed...")')
            write(*,'( 1x,"search along the gating coordinate leads out of bounds...")')
            write(*,'( 1x,"consider larger gating interval...")')
            stop
         endif

         x0(1) = x(1)
         x0(2) = x(2)

         ierr = 0

         if (iminim.eq.1) then

            !-- Newton-Raphson optimization
            call newton3(n,kg,x0,x,xs,g,f,d,ev,w,s,fx,slim,acc,maxit,ierr,iter,nf,ng,nfc)

	 elseif (iminim.eq.2) then

            !-- LBFGS optimization
            call lbfgs3(n,kg,mcorr,x0,x,xs,g,fx,iprint,factr,pgtol,maxit,ierr,iter,nf)

	 endif

         if (ierr.ne.0) then
            write(6,'(/1x,''warning: minimization failed...'')')
            write(6,'(1x,10(''-''),''input'',10(''-''))')
            write(6,'(1x,''mode = '',a5)') mode_
            write(6,'(1x,''iset = '',i1)') iset_
            write(6,'(1x,''ist  = '',i2)') istate_
            write(6,'(1x,''zp0  = '',f10.3)') zp0
            write(6,'(1x,''zps  = '',f10.3)') zps
            write(6,'(1x,''ze0  = '',f10.3)') ze0
            write(6,'(1x,''zes  = '',f10.3)') zes
            write(6,'(1x,''kg,r = '',i5,f10.3)') kg,glist(kg)*bohr2a
            write(6,'(1x,''slim = '',g20.10)') slim
            write(6,'(1x,''acc  = '',g20.10)') acc
            write(6,'(1x,10(''-''),''output'',10(''-''))')
            write(6,'(1x,''ierr = '',i2)') ierr
            write(6,'(1x,''iter = '',i3)') iter
            write(6,'(1x,''nf   = '',i3)') nf
            write(6,'(1x,''ng   = '',i3)') ng
            write(6,'(1x,''nfc  = '',i3)') nfc
            write(6,'(1x,''fx   = '',f10.3)') fx
            write(6,'(1x,''x    = '',2f15.6)') (x(i),i=1,2)
            write(6,'(1x,''g    = '',2f15.6)') (g(i),i=1,2)
            write(6,'(1x,''ev   = '',2f15.6)') (ev(i),i=1,2)
            write(6,'(1x,10(''-''),''output'',10(''-'')/)')
         endif

         if (fx.lt.femin) then

            zpmin = x(1)/zps
            zemin = x(2)/zes
            kgmin = kg
            rmin = bohr2a*glist(kg)
            femin = fx
            dzpmin = zps*g(1)
            dzemin = zes*g(2)

            d2zpmin = zps*zps*f(1,1)
            d2zemin = zes*zes*f(2,2)
            d2zpzemin = zps*zes*f(1,2)

         elseif (kg.ne.kgmin) then

            igdir = -igdir
            ichdir = ichdir + 1

         endif

         done = ichdir.ge.2

      enddo

      grad = 0.d0
      hess = 0.d0

      ! gradient at the minimum

      grad(1) = dzpmin
      grad(2) = dzemin
      grad(3) = 0.d0

      ! second derivative with respect to the solvent coordinates

      hess(1,1) = d2zpmin
      hess(2,2) = d2zemin
      hess(1,2) = d2zpzemin
      hess(2,1) = hess(1,2)

      ! second derivative with respect to the gating coordinate
      ! (on the grid, very approximate)

      if (kg-2.lt.1.or.kg+2.gt.npntsg) then
         write(6,'(/1x,"warning: 2-nd derivative with respect")')
         write(6,'( 1x,"to the gating coordinate can not be evaluated (set to zero).")')
         write(6,'( 1x,"(too close to the boundaries)")')
         hess(3,3) = 0.d0
      else
         call fes3(n,kgmin+2,x,xs,fright)
         call fes3(n,kgmin  ,x,xs,fmiddl)
         call fes3(n,kgmin-2,x,xs,fleft)
         hess(3,3) = 0.25d0*(fright + fleft - 2.d0*fmiddl)/ginc2
      endif

      ! second mixed derivative with respect to the gating coordinate
      ! and solvent coordinates
      ! (on the grid, very approximate)

      if (kg-1.lt.1.or.kg+1.gt.npntsg) then
         write(6,'(/1x,"warning: 2-nd mixed derivative with respect")')
         write(6,'( 1x,"can not be evaluated (set to zero).")')
         write(6,'( 1x,"(too close to the boundaries)")')
      else
         call d1fes3(n,kgmin+1,x,xs,gleft)
         call d1fes3(n,kgmin-1,x,xs,gright)
         hess(1,3) = half*(gright(1) - gleft(1))/ginc
	 hess(3,1) = hess(1,3)
         hess(2,3) = half*(gright(2) - gleft(2))/ginc
	 hess(3,2) = hess(2,3)
      endif
 
      return

   end subroutine fes3min

   !=======================================================================
   subroutine fes3gmin(mode_,iset_,istate_,kg,zp0,zps,ze0,zes,&
                       ierr,zpmin,zemin,femin,&
                       dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)
   !=======================================================================
   !     Finds constrained minima on the spesified three-dimensional
   !     free energy surface using canonical Newton-Raphson
   !     minimization algorithm.
   !
   !     MODE - (ADIAB, DIAB2 or DIAB4) - type of the surface
   !
   !     ISET - a set of the surfaces (1 for ADIAB, 1/2 for DIAB2,
   !            1/2/3/4 for DIAB4)
   !
   !     ISTATE - quantum number of the desired state
   !
   !     KG - grid point along the gating coordinate (constraint)
   !
   !     ZP0,ZE0 - initial values of solvent coordinates (kcal/mol)
   !
   !     ZPS,ZES - scaling factors for the solvent coordinates
   !
   !     IERR - error code (0 - normal)
   !
   !     ZPMIN,ZEMIN - solvent coordinates at the minimum
   !
   !     FEMIN - free energy at the minimum (kcal/mol)
   !
   !     DZPMIN,DZEMIN - derivatives at the minimum
   !
   !     D2ZPMIN,D2ZEMIN,D2ZPZEMIN - second derivatives at the minimum
   !=======================================================================

      character(5), intent(in) :: mode_
      integer,      intent(in) :: iset_, istate_, kg
      real*8,       intent(in) :: zp0, zps, ze0, zes
      
      integer, intent(out) :: ierr
      real*8,  intent(out) :: zpmin, zemin, femin
      real*8,  intent(out) :: dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin

      integer, parameter :: nrdim = 2
      integer :: i, n, iter, nf, ng, nfc
      real*8  :: fx

      real*8, dimension(nrdim)       :: x0, x, xs, g, s, w, ev
      real*8, dimension(nrdim,nrdim) :: f, d

      n = nrdim
      mode = mode_
      iset = iset_
      ist  = istate_

      x0(1) = zp0*zps
      x0(2) = ze0*zes
      xs(1) = zps
      xs(2) = zes

      ierr = 0

      if (iminim.eq.1) then

         !-- Newton-Rafson optimization
         call newton3(n,kg,x0,x,xs,g,f,d,ev,w,s,fx,slim,acc,maxit,ierr,iter,nf,ng,nfc)

      elseif (iminim.eq.2) then

         !-- LBFGS optimization
         call lbfgs3(n,kg,mcorr,x0,x,xs,g,fx,iprint,factr,pgtol,maxit,ierr,iter,nf)

      endif

      zpmin = x(1)/zps
      zemin = x(2)/zes
      femin = fx
      dzpmin = zps*g(1)
      dzemin = zes*g(2)

      if (iminim.eq.1) then

         d2zpmin = zps*zps*f(1,1)
         d2zemin = zes*zes*f(2,2)
         d2zpzemin = zps*zes*f(1,2)

      elseif (iminim.eq.1) then

         call d2fes3(n,kg,x,xs,f)
         d2zpmin   = f(1,1)
         d2zemin   = f(2,2)
         d2zpzemin = f(1,2)
                                                                                                                                                      
      endif


      if (ierr.ne.0) then

         write(6,'(/1x,''warning: minimization failed...'')')
         write(6,'(1x,10(''-''),''input'',10(''-''))')
         write(6,'(1x,''mode = '',a5)') mode_
         write(6,'(1x,''iset = '',i1)') iset_
         write(6,'(1x,''ist  = '',i2)') istate_
         write(6,'(1x,''zp0  = '',f10.3)') zp0
         if(iminim.eq.1) write(6,'(1x,''zps  = '',f10.3)') zps
         write(6,'(1x,''ze0  = '',f10.3)') ze0
         if(iminim.eq.1) write(6,'(1x,''zes  = '',f10.3)') zes

         if (iminim.eq.1) then
            write(6,'(1x,''slim = '',g20.10)') slim
            write(6,'(1x,''acc  = '',g20.10)') acc
         elseif (iminim.eq.2) then
            write(6,'(1x,"factr = ",g20.10)') factr
            write(6,'(1x,"pgtol = ",g20.10)') pgtol
	 endif

         write(6,'(1x,10(''-''),''output'',10(''-''))')
         write(6,'(1x,''ierr = '',i2)') ierr

         if (iminim.eq.2) then
            !-- error information from LBFGS routine
         endif

         write(6,'(1x,''iter = '',i3)') iter
         write(6,'(1x,''nf   = '',i3)') nf
         if(iminim.eq.1) write(6,'(1x,''ng   = '',i3)') ng
         if(iminim.eq.1) write(6,'(1x,''nfc  = '',i3)') nfc
         write(6,'(1x,''fx   = '',f10.3)') fx
         write(6,'(1x,''x    = '',2f15.6)') (x(i),i=1,2)
         write(6,'(1x,''g    = '',2f15.6)') (g(i),i=1,2)
         if(iminim.eq.1) write(6,'(1x,''ev   = '',2f15.6)') (ev(i),i=1,2)
         write(6,'(1x,10(''-''),''output'',10(''-'')/)')

      endif

      return

   end subroutine fes3gmin

   !=======================================================================
   subroutine newton3&
      (n,kg,x0,x,xs,g,f,d,ev,w,s,fx,slim,eps,lim,istat,iter,nf,ng,nfc)
 
      integer, intent(in)  :: n, kg, lim
      real*8,  intent(in)  :: slim, eps
      real*8,  intent(in)  :: x0(n), xs(n)
      integer, intent(out) :: iter, nf, ng, nfc, istat
      real*8,  intent(out) :: fx
      real*8,  intent(out) :: x(n), g(n), f(n,n), s(n), d(n,n), w(n), ev(n)
      
      logical :: conv
      integer :: ierr, i, j
      real*8  :: ss, sg, gnorm

      x = x0
      s = 0.d0

      iter = 0
      nf = 0
      ng = 0
      nfc = 0
      istat = 1

      call fes3(n,kg,x,xs,fx)
      nf = nf + 1
      call d1fes3(n,kg,x,xs,g)
      ng = ng + 1
      call d2fes3(n,kg,x,xs,f)
      nfc = nfc + 1

      d = f

      call treql(n,d,ev,w,ierr)
      if (ierr.ne.0) then
         istat = -3
         return
      endif

      ! minimization loop
      do

         call invers(n,f,w,ierr)
         if (ierr.ne.0) then
            istat = -2
            exit
         endif

         do i=1,n
            s(i)=0
            do j=1,n
               s(i) = s(i) - f(i,j)*g(j)
            enddo
         enddo

         ss = slen(n,s)
         sg = slen(n,g)

         if (ss.gt.slim) then
            do i=1,n
               s(i)=s(i)/ss*slim
            enddo
         endif

         x = x + s

         iter = iter + 1
         call fes3(n,kg,x,xs,fx)
         nf = nf + 1
         call d1fes3(n,kg,x,xs,g)
         ng = ng + 1
         call d2fes3(n,kg,x,xs,f)
         nfc = nfc + 1

         d = f

         call treql(n,d,ev,w,ierr)
         if (ierr.ne.0) then
            istat = -3
            exit
         endif

         conv = .true.

         do i=1,n
            if (dabs(s(i)).gt.eps) then
               istat=1
               conv=.false.
            endif
         enddo

         gnorm = slen(n,g)
         write(*,'(1x,i4,2x,8f15.6)') iter,(x(i)/xs(i),i=1,n),fx,gnorm,(ev(i),i=1,n)

         if (iter.gt.lim) istat=-1

         if (conv) then

            istat = 0
            exit

         else

            if (istat.lt.0) exit

         endif

      enddo

      return

   end subroutine newton3


   !===================================================================!
   subroutine lbfgs3(n,kg,m,x0,x,xs,g,fx,iprint,factr,pgtol,maxit,istat,iter,nf)

      integer, intent(in)                  :: n, kg, m, maxit, iprint
      real*8,  intent(in)                  :: factr, pgtol
      real*8,  intent(in),  dimension(n)   :: x0, xs
      integer, intent(out)                 :: iter, nf, istat
      real*8,  intent(out)                 :: fx
      real*8,  intent(out), dimension(n)   :: x, g

      !-- local variables

      integer, parameter :: mmax=17     ! maximum number of limited memory corrections
 
      character(len=60)       :: task, csave
      logical, dimension(4)   :: lsave
      integer, dimension(n)   :: nbd
      integer, dimension(3*n) :: iwa
      integer, dimension(44)  :: isave

      real(kind=8) :: f
      real(kind=8), dimension(n)  :: l, u
      real(kind=8), dimension(29) :: dsave
      real(kind=8), dimension(2*m*n+4*n+12*m*m+12*m) :: wa

      real(kind=8) :: t1, t2
      integer      :: i

      istat = 0

      !-- We wish to have output at every iteration.
      ! iprint = 1

      !-- tolerances in the stopping criteria.
      !
      !   factr = 1.d+12 for low accuracy;
      !           1.d+7  for moderate accuracy; 
      !           1.d+1  for extremely high accuracy.
      !
      !   pgtol >= 0
      !   The iteration will stop when max{|proj g_i|, i=1,...,n} <= pgtol
      !   where pg_i is the ith component of the projected gradient.
      !   The user can suppress this termination test by setting pgtol=0.

      !-- m - number of corrections used in the limited memory matrix.
      !   It is not altered by the routine.  Values of m < 3  are
      !   not recommended, and large values of m can result in excessive
      !   computing time. The range  3 <= m <= 20 is recommended. 

      if (m.lt.3.or.m.gt.20) then
         write(*,*)
         write(*,*) 'Number of corrections used in the limited memory matrix (parameter m)'
         write(*,*) 'is outside of the recommended range (3 <= m <= 20). Optimization aborted...'
         return
      endif
 
      !-- no bounds

      do i=1,n
         nbd(i)=0
         l(i)=0.0d0
         u(i)=0.0d0
      enddo

      !-- We now define the starting point.

      do i=1,n
         x(i) = x0(i)
      enddo
 
      !-- start the iteration by initializing task
 
      task = 'START'

      !------- the beginning of the loop ----------
 
 111  continue

      optimization_loop: do

         !-- This is the call to the L-BFGS-B code.
         call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)
 
         if (task(1:2).eq.'FG') then

            !-- the minimization routine has returned to request the
            !   function f and gradient g values at the current x.

            !-- Compute function value f
            call fes3(n,kg,x,xs,f)

            !-- Compute gradient vector g
            call d1fes3(n,kg,x,xs,g)

         elseif (task(1:5) .eq. 'NEW_X') then

            !-- the minimization routine has returned with a new iterate,
            !   and we have opted to continue the iteration.

            iter = isave(30)

            !-- Terminate if the total number of f and g evaluations exceeds 10,000
            if (isave(34).ge.10000) then
               task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS 10,000'
               istat = 1
            endif

            !-- Terminate if the number of iterations exceeds maxit
            if (isave(30).ge.maxit) then
               task='STOP: TOTAL NO. of iterations EXCEEDS maxit'
               istat = 2
            endif

            !-- Terminate if  |proj g|/(1+|f|) < 1.0d-10, where 
            !   "proj g" denoted the projected gradient

            if(dsave(13).le.1.d-10*(1.0d0 + abs(f))) &
            & task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'

            !-- print the following information at each iteration:
            !
            ! - the current iteration number, isave(30),
            ! - the total number of f and g evaluations, isave(34),
            ! - the value of the objective function f,
            ! - the norm of the projected gradient,  dsve(13)

            write(*,'(2(a,i5,4x),a,d12.5,4x,a,d12.5)') &
            & 'Iteration',isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)

         else

            exit

         endif

      enddo optimization_loop

      !-- If the run is to be terminated, we print also the information
      !   contained in task as well as the final value of x.

      fx = f
      nf = isave(34)

      if (task(1:4) .eq. 'STOP') then

         write (*,*) task  
         write (*,*) 'Final X='
         write (*,'((1x,6(1x,f15.6)))') (x(i),i = 1,n)

      elseif (task(1:4) .eq. 'CONV') then
 
         istat = 0

      elseif (task(1:4) .eq. 'ABNO') then
 
         istat = -1
         write (*,*) 'Abnormal termination: termination conditions not satisfied'
         write (*,*) 'Final X='
         write (*,'((1x,6(1x,f15.6)))') (x(i),i = 1,n)
         write (*,*) 'Final function:',fx
         write (*,*) 'Final gradient:',g

      elseif (task(1:4) .eq. 'ERROR') then
 
         istat = -2
         write (*,*) 'Error termination: check input parameters'
         write (*,*) 'current X='
         write (*,'((1x,6(1x,f15.6)))') (x(i),i = 1,n)
         write (*,*) 'Current function:',fx
         write (*,*) 'Current gradient:',g

      endif

   end subroutine lbfgs3


   !=======================================================================
   subroutine fes3(n,kg,x,xs,e)
   !=======================================================================
   !     Calculates free energy of the particular state
   !     N    - number of variables (almost always two: ZP and ZE)
   !     KG   - grid point along the gating coordinate
   !     X    - vector of variables
   !     XS   - vector of scaling factors for variables
   !     E    - free energy in kcal/mol
   !=======================================================================
      use feszz_3d, only: feszz3
      integer, intent(in) :: n, kg
      real*8,  intent(in), dimension(n) :: x, xs
      real*8,  intent(out) :: e

      integer :: ielst, nzdim, ndabf
      real*8  :: zp, ze

      real*8, allocatable :: fe(:),z(:,:),psipr(:,:,:)
      real*8, allocatable :: psiel(:,:,:)
      real*8, allocatable :: enel(:,:),envib(:,:)

      zp = x(1)/xs(1)
      ze = x(2)/xs(2)

      if (mode.eq.'ADIAB') then
         ielst = nelst
      elseif (mode.eq.'DIAB2') then
         ielst = 2
      elseif (mode.eq.'DIAB4') then
         ielst = 1
      else
         write(6,'(//1x,"in fes3: mode is invalid: ",a5)') mode
         stop
      endif

      nzdim = ielst*nprst
      allocate (fe(ist),z(nzdim,nzdim),psipr(ielst,nprst,npnts))
      allocate (psiel(ielst,npnts,ielst))
      allocate (enel(ielst,npnts),envib(ielst,nprst))

      call feszz3(mode,iset,kg,zp,ze,ist,fe,nzdim,z,ndabf,&
                  ielst,enel,envib,psiel,psipr)
      e = fe(ist)

      deallocate (fe,z,psipr)
      deallocate (psiel)
      deallocate (enel,envib)

      return

   end subroutine fes3

   !=======================================================================
   subroutine d1fes3(n,kg,x,xs,g)

      integer, intent(in) :: n, kg
      real*8,  intent(in),  dimension(n) :: x, xs
      real*8,  intent(out), dimension(n) :: g

      real*8, parameter :: xinc = 1.d-4
      integer :: i
      real*8, dimension(n) :: xt
      
      real*8 :: ep2, ep1

      do i=1,n
         xt(i) = x(i)
      enddo

      do i=1,n
         xt(i) = x(i) + xinc*xs(i)/2.d0
         call fes3(n,kg,xt,xs,ep2)
         xt(i) = x(i) - xinc*xs(i)/2.d0
         call fes3(n,kg,xt,xs,ep1)
         xt(i) = x(i)
         g(i) = (ep2-ep1)/xinc/xs(i)
      enddo

      return

   end subroutine d1fes3

   !=======================================================================
   subroutine d2fes3(n,kg,x,xs,f)

      integer, intent(in) :: n, kg
      real*8,  intent(in),  dimension(n)   :: x, xs
      real*8,  intent(out), dimension(n,n) :: f

      real*8, parameter :: xinc = 1.d-4, two = 2.d0
      real*8, dimension(n) :: xt, grad, gnext
      integer :: i, j

      call d1fes3(n,kg,x,xs,grad)

      xt = x

      do i=1,n
         xt(i) = x(i) + xinc*xs(i)
         call d1fes3(n,kg,xt,xs,gnext)
         xt(i) = x(i)
         do j=1,n
            f(i,j)= (gnext(j)-grad(j))/xinc/xs(j)
         enddo
      enddo

      do i=1,n
         do j=1,i-1
            f(i,j) = (f(i,j)+f(j,i))/two
            f(j,i) = f(i,j)
         enddo
      enddo

      return

   end subroutine d2fes3

end module fesmin_3d
