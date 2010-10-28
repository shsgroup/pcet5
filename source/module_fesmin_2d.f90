module fesmin_2d

!-------------------------------------------------------------------
!  minimization on a 2D free energy surface
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!-------------------------------------------------------------------

   use pardim
   use cst
   use eispack
   use control
   use quantum
   
   implicit none
   private

   character(5), save :: mode
   integer,      save :: iset, ist

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

   public :: fes2min

   contains

   !=======================================================================
   subroutine fes2min(mode_,iset_,istate_,zp0,zps,ze0,zes,&
                      ierr,zpmin,zemin,femin,&
                      dzpmin,dzemin,d2zpmin,d2zemin,d2zpzemin)
   !=======================================================================
   !     Finds minima on the spesified two-dimensional
   !     free energy surface using canonical Newton-Raphson
   !     minimization algorithm or LBFGS algorithm
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
   !     IERR - error code (0 - normal)
   !
   !     ZPMIN,ZEMIN - solvent coordinates at the minimum
   !
   !     FEMIN - free energy at the minimum (kcal/mol)
   !
   !     DZPMIN,DZEMIN - derivatives at the minimum
   !
   !     D2ZPMIN,D2ZEMIN,D2ZPZEMIN - second derivatives at the minimum
   !
   !=======================================================================
      implicit none

      character(5), intent(in) :: mode_
      integer,      intent(in) :: iset_, istate_
      real*8,       intent(in) :: zp0, zps, ze0, zes
      
      integer, intent(out) :: ierr
      real*8,  intent(out) :: zpmin, zemin, femin
      real*8,  intent(out) :: dzpmin, dzemin, d2zpmin, d2zemin, d2zpzemin

      integer, parameter :: nrdim2 = 2

      integer :: i, n, iter, nf, ng, nfc
      real*8  :: fx, factr, pgtol

      real*8, dimension(nrdim2)        :: x0, x, xs, g, s, w, ev
      real*8, dimension(nrdim2,nrdim2) :: f, d

      n = nrdim2
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
         call newton2(n,x0,x,xs,g,f,d,ev,w,s,fx,slim,acc,maxit,ierr,iter,nf,ng,nfc)

      elseif (iminim.eq.2) then

         !-- LBFGS optimization
         call lbfgs2(n,mcorr,x0,x,xs,g,fx,iprint,factr,pgtol,maxit,ierr,iter,nf)

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

         call d2fes2(n,x,xs,f)
         d2zpmin   = f(1,1)
         d2zemin   = f(2,2)
         d2zpzemin = f(1,2)

      endif

      if (ierr.ne.0) then

         write(6,'(/1x,"warning: minimization failed...")')
         write(6,'(1x,10("-"),"input",10("-"))')
         write(6,'(1x,"mode = ",a5)') mode_
         write(6,'(1x,"iset = ",i1)') iset_
         write(6,'(1x,"ist  = ",i2)') istate_
         write(6,'(1x,"zp0  = ",f10.3)') zp0
         if(iminim.eq.1) write(6,'(1x,"zps  = ",f10.3)') zps
         write(6,'(1x,"ze0  = ",f10.3)') ze0
         if(iminim.eq.1) write(6,'(1x,"zes  = ",f10.3)') zes

         if (iminim.eq.1) then
	    write(6,'(1x,"slim = ",g20.10)') slim
            write(6,'(1x,"acc  = ",g20.10)') acc
         elseif (iminim.eq.2) then
	    write(6,'(1x,"factr = ",g20.10)') factr
            write(6,'(1x,"pgtol = ",g20.10)') pgtol
	 endif

         write(6,'(1x,10("-"),"output",10("-"))')
         write(6,'(1x,"ierr = ",i2)') ierr

         if (iminim.eq.2) then
            !-- error information from LBFGS routine
	 endif

         write(6,'(1x,"iter = ",i3)') iter
         write(6,'(1x,"nf   = ",i3)') nf
         if(iminim.eq.1) write(6,'(1x,"ng   = ",i3)') ng
         if(iminim.eq.1) write(6,'(1x,"nfc  = ",i3)') nfc
         write(6,'(1x,"fx   = ",f10.3)') fx
         write(6,'(1x,"x    = ",2f15.6)') (x(i),i=1,2)
         write(6,'(1x,"g    = ",2f15.6)') (g(i),i=1,2)
         if(iminim.eq.1) write(6,'(1x,"ev   = ",2f15.6)') (ev(i),i=1,2)
         write(6,'(1x,10("-"),"output",10("-")/)')

      endif

      return

   end subroutine fes2min

   !===================================================================!
   subroutine newton2 &
      (n,x0,x,xs,g,f,d,ev,w,s,fx,slim,eps,lim,istat,iter,nf,ng,nfc)

      integer, intent(in)                  :: n, lim
      real*8,  intent(in)                  :: slim, eps
      real*8,  intent(in),  dimension(n)   :: x0, xs
      integer, intent(out)                 :: iter, nf, ng, nfc, istat
      real*8,  intent(out)                 :: fx
      real*8,  intent(out), dimension(n)   :: x, g, s, w, ev
      real*8,  intent(out), dimension(n,n) :: f, d
      
      logical :: conv
      integer :: ierr, i, j
      real*8  :: ss, sg, gnorm

      x = x0
      s = 0.d0

      iter=0
      nf=0
      ng=0
      nfc=0
      istat=1

      call fes2(n,x,xs,fx)
      nf = nf + 1
      call d1fes2(n,x,xs,g)
      ng = ng + 1
      call d2fes2(n,x,xs,f)
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
            s(i) = 0.d0
            do j=1,n
               s(i) = s(i) - f(i,j)*g(j)
            enddo
         enddo

         ss = slen(n,s)
         sg = slen(n,g)

         if (ss.gt.slim) then
            do i=1,n
               s(i) = s(i)/ss*slim
            enddo
         endif

         x = x + s

         iter = iter + 1
         call fes2(n,x,xs,fx)
         nf = nf + 1
         call d1fes2(n,x,xs,g)
         ng = ng + 1
         call d2fes2(n,x,xs,f)
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
               istat = 1
               conv = .false.
            endif
         enddo

         gnorm = slen(n,g)
         !write(*,'(1x,i4,2x,8f15.6)') iter,(x(i)/xs(i),i=1,n),fx,gnorm,(ev(i),i=1,n)

         if (iter.gt.lim) istat = -1

         if (conv) then

            istat = 0
            exit

         else

            if (istat.lt.0) exit

         endif

      enddo

      return

   end subroutine newton2

   !===================================================================!
   subroutine lbfgs2(n,m,x0,x,xs,g,fx,iprint,factr,pgtol,maxit,istat,iter,nf)

      integer, intent(in)                  :: n, m, maxit, iprint
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
            call fes2(n,x,xs,f)

            !-- Compute gradient vector g
            call d1fes2(n,x,xs,g)

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

            write(*,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') &
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

   end subroutine lbfgs2

   !===================================================================!
   subroutine fes2(n,x,xs,e)
   !===================================================================!
   !     Calculates free energy of the particular state
   !     N    - number of variables (almost always two: ZP and ZE)
   !     X    - vector of variables
   !     XS   - vector of scaling factors for variables
   !     E    - free energy in kcal/mol
   !===================================================================!
      use feszz_2d, only: feszz2
      integer, intent(in) :: n
      real*8,  intent(in), dimension(n) :: x, xs
      real*8,  intent(out) :: e

      integer :: ielst, nzdim, nz, nfe, nstates, ierr, npsiga, ndabf
      real*8  :: zp, ze

      real*8, allocatable :: fe(:),z(:),psipr(:,:,:,:),psiga(:)
      real*8, allocatable :: psiel(:,:,:,:)
      real*8, allocatable :: enel(:,:,:),envib(:,:,:),envibg(:,:,:)

      zp = x(1)/xs(1)
      ze = x(2)/xs(2)

      if (mode.eq.'ADIAB') then
         ielst = nelst
      elseif (mode.eq.'DIAB2') then
         ielst = 2
      elseif (mode.eq.'DIAB4') then
         ielst = 1
      else
         ielst = nelst
      endif

      if (mgquant.eq.1.or.mgquant.eq.2) then
         nzdim = ielst*nprst*ngast
         nz = nzdim*nzdim
         nstates = ist
         nfe = nstates
      elseif (mgquant.eq.3) then
         nzdim = ielst*nprst
         nz = nzdim*nzdim*npntsg
         nstates = (ist-1)/ngast + 1
         nfe = ist
      endif

      npsiga = ielst*nprst*ngast*npntsg

      allocate (z(nz),&
                fe(nfe),&
                psipr(ielst,nprst,npnts,npntsg),&
                psiga(npsiga),&
                psiel(ielst,npnts,npntsg,ielst),&
                enel(ielst,npnts,npntsg),&
                envib(ielst,nprst,npntsg),&
                envibg(ielst,nprst,ngast),stat=ierr)

      if (ierr.ne.0) then
         write(*,*) '#### allocation error in fes2'
         stop
      endif

      ndabf = 0
      call feszz2(mode,iset,zp,ze,nstates,nfe,fe,nz,z,ndabf,&
                  ielst,enel,envib,envibg,psiel,psipr,npsiga,psiga)
      e = fe(ist)

      deallocate (z,fe,psipr,psiga,psiel,enel,envib,envibg,stat=ierr)
      if (ierr.ne.0) then
         write(*,*) '#### deallocation error in fes2'
         stop
      endif

      return

   end subroutine fes2

   !===================================================================!
   subroutine d1fes2(n,x,xs,g)

      integer, intent(in) :: n
      real*8,  intent(in),  dimension(n) :: x, xs
      real*8,  intent(out), dimension(n) :: g

      real*8, parameter :: xinc = 1.d-4
      integer :: i
      real*8, dimension(n) :: xt
      
      real*8 :: ep2, ep1

      xt = x

      do i=1,n
         xt(i) = x(i) + xinc*xs(i)/2.d0
         call fes2(n,xt,xs,ep2)
         xt(i) = x(i) - xinc*xs(i)/2.d0
         call fes2(n,xt,xs,ep1)
         xt(i) = x(i)
         g(i) = (ep2-ep1)/xinc/xs(i)
      enddo

      return

   end subroutine d1fes2

   !===================================================================!
   subroutine d2fes2(n,x,xs,f)

      integer, intent(in)                     :: n
      real*8,  intent(in),     dimension(n)   :: x, xs
      real*8,  intent(out),    dimension(n,n) :: f

      real*8, parameter :: xinc = 1.d-4, two = 2.d0
      real*8, dimension(n) :: xt, grad, gnext
      integer :: i, j

      call d1fes2(n,x,xs,grad)
      
      xt = x

      do i=1,n
         xt(i) = x(i) + xinc*xs(i)
         call d1fes2(n,xt,xs,gnext)
         xt(i) = x(i)
         do j=1,n
            f(i,j) = (gnext(j)-grad(j))/xinc/xs(j)
         enddo
      enddo

      do i=1,n
         do j=1,i-1
            f(i,j) = (f(i,j)+f(j,i))/two
            f(j,i) = f(i,j)
         enddo
      enddo

      return

   end subroutine d2fes2

end module fesmin_2d
