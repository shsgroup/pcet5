module quantum
!===============================================================
!  Quantities for 1D-schroedinger equation (a.u.)
!---------------------------------------------------------------
!    Quantities for proton (a.u.)
!---------------------------------------------------------------
!    PM       - mass of the quantum particle
!    ALIM     - left integration limit
!    BLIM     - right integration limit
!    RLIST    - grid points along the quantum proton coordinate
!    NPNTS    - number of grid points along the quantum proton
!               coordinate
!    NPNTSSOL - number of grid points along the quantum proton
!               coordinate for solvation calculations
!    NPRST    - number of proton vibrational states to include
!    HKE      - kinetic energy matrix in the grid representation
!               (for unit-mass particle).
!    DX       - matrix of the derivative (moment) in the grid
!               representation.
!---------------------------------------------------------------
!    Quantities for gating coordinate (a.u.)
!---------------------------------------------------------------
!    DM       - mass of the proton donor
!    AM       - mass of the proton acceptor
!    AGLIM    - left grid limit
!    BGLIM    - right grid limit
!    GLIST    - grid points along the gating coordinate
!    HGKE     - kinetic energy matrix in the grid representation
!               (for unit-mass particle).
!    DGX      - matrix of the first derivative in the grid
!               representation
!    NPNTSG   - number of grid points along the gating
!               coordinate
!    NPNTSSOLG- number of grid points along the gating
!               coordinate for solvation calculations
!    NGAST    - number of gating vibrational states
!               per vibronic state
!---------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!===============================================================

   use cst
   use control
   use eispack

   implicit none
   public
   save

   ! proton quantities
   integer                             :: npnts, npntssol, nprst
   real*8                              :: pm, alim, blim, rhmin
   real*8, allocatable, dimension(:)   :: rlist
   real*8, allocatable, dimension(:,:) :: hke, dx

   ! gating quantities
   integer                             :: npntsg, npntssolg, ngast
   real*8                              :: dm, am, aglim, bglim
   real*8, allocatable, dimension(:)   :: glist
   real*8, allocatable, dimension(:,:) :: hgke, dgx

contains
!===============================================================

   subroutine alloc_pquant(ndim,deriv)

      integer, intent(in) :: ndim
      logical, intent(in) :: deriv

      integer :: k
      real*8  :: step, x

      if (ndim.le.0) then
         write(*,*) "ERROR in alloc_pquant: Wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(rlist) ) then

         allocate (rlist(ndim+1))
         step = a2bohr*(blim-alim)/(ndim-1)
         x = alim*a2bohr
         do k=1,ndim+1
            rlist(k) = x
            x = x + step
         enddo

      else

         write(*,*) "ERROR in alloc_rquant: RLIST is already allocated..."
         stop

      endif

      if ( .not.allocated(hke) ) then
         allocate (hke(ndim,ndim))
         hke = zero
      else
         write(*,*) "ERROR in alloc_pquant: HKE is already allocated..."
         stop
      endif

      if (deriv) then

         if ( .not.allocated(dx) ) then

            allocate (dx(ndim,ndim))
            dx = zero

         else

            write(*,*) "ERROR in alloc_pquant: DX is already allocated..."
            stop

         endif

      endif

   end subroutine alloc_pquant

   subroutine dealloc_pquant
      if (allocated(rlist)) deallocate (rlist)
      if (allocated(hke) )  deallocate (hke)
      if (allocated(dx)  )  deallocate (dx)
   end subroutine dealloc_pquant

   subroutine alloc_gquant(ngdim,derivg)
   
      integer, intent(in) :: ngdim
      logical, intent(in) :: derivg
      
      integer :: k
      real*8  :: x, stepg

      if (ngdim.le.0) then
         write(*,*) "ERROR in alloc_gquant: Wrong dimension passed..."
         stop
      endif

      if ( .not.allocated(glist) ) then

         allocate (glist(ngdim+1))

         if (ngdim.gt.1) then

            stepg = a2bohr*(bglim-aglim)/(ngdim-1)
            x = aglim*a2bohr
            do k=1,ngdim+1
               glist(k) = x
               x = x + stepg
            enddo

         else

            glist(1) = 0.d0

         endif

      else

         write(*,*) "ERROR in alloc_gquant: GLIST is already allocated..."
         stop

      endif

      if (gquant) then

         if ( .not.allocated(hgke) ) then                                      
            allocate (hgke(ngdim,ngdim))                                       
            hgke = zero                                                        
         else                                                                  
            write(*,*) "ERROR in alloc_gquant: HGKE is already allocated..."   
            stop                                                               
         endif                                                                 

         if (derivg) then                                                      

            if ( .not.allocated(dgx) ) then                                    

               allocate (dgx(ngdim,ngdim))                                     
               dgx = zero                                                      

            else                                                               

               write(*,*) "ERROR in alloc_gquant: DGX is already allocated..." 
               stop                                                            

            endif                                                              

         endif                                                                 

      endif

   end subroutine alloc_gquant

   subroutine dealloc_gquant
      if (allocated(glist)) deallocate (glist)
      if (allocated(hgke) ) deallocate (hgke)
      if (allocated(dgx)  ) deallocate (dgx)
   end subroutine dealloc_gquant

   pure subroutine gridke(npnts_,xint_,hke_)
   !-----------------------------------------------------------------------
   ! This subroutine calculates a kinetic energy matrix on the grid
   ! (npnts_ SHOULD BE AN INTEGER POWER OF TWO !!!)
   !-----------------------------------------------------------------------

      integer, intent(in) :: npnts_
      real*8,  intent(in) :: xint_
      real*8, intent(out), dimension(npnts_,npnts_) :: hke_

      integer :: nn,i,j,k,kk
      real*8  :: xintf, delk
      real*8, allocatable :: pmom(:),delf(:)

      allocate (pmom(npnts_))
      allocate (delf(2*npnts_))

      ! delk=two*pi/xint_  this was not accounted for additional grid point needed
      !                    by FFT routine

      xintf = xint_ + xint_/(npnts_-1)
      delk = two*pi/xintf

      nn = npnts_/2
      do i=1,nn
         pmom(i) = (i-1)*delk*hbar
         pmom(i+nn) = (i-1-nn)*delk*hbar
      enddo

      do i=1,npnts_

         do j=1,npnts_
            if (j.eq.i) then
               delf(2*j-1) = one
	       delf(2*j)   = zero
            else
               delf(2*j-1) = zero
	       delf(2*j)   = zero
            endif
         enddo

         !---- do fourier transform from the position
         !     to the momentum representation

         call four1(delf,npnts_,1)

         do k=1,2*npnts_
	    kk = k/2 + mod(k,2)
            delf(k)=delf(k)*pmom(kk)**2/two
         enddo

         !---- do inverse fourier transform from the momentum
         !     representation to the postion representation

         call four1(delf,npnts_,-1)

         do j=1,npnts_
            hke_(j,i) = delf(2*j-1)/dble(npnts_)
         enddo

      enddo

      deallocate (pmom,delf)

      return
   end subroutine gridke

   pure subroutine griddx(npnts_,xint_,dx_)
   !-----------------------------------------------------------------------
   ! This subroutine calculates a matrix of the first derivative
   ! (npnts_ SHOULD BE AN INTEGER POWER OF TWO !!!)
   !-----------------------------------------------------------------------

      integer, intent(in) :: npnts_
      real*8,  intent(in) :: xint_
      real*8,  intent(out), dimension(npnts_,npnts_) :: dx_

      integer :: nn, i, j, k
      real*8  :: xintf, delk, delf_re, delf_im
      real*8, allocatable :: delf(:), pmom(:)

      allocate (pmom(npnts_))
      allocate (delf(2*npnts_))

      ! delk=two*pi/xint_  this was not accounted for additional grid point needed
      !                    by FFT routine

      xintf = xint_ + xint_/(npnts_-1)
      delk = two*pi/xintf

      nn = npnts_/2
      do i=1,nn
         pmom(i)=(i-1)*delk*hbar
         pmom(i+nn)=(i-1-nn)*delk*hbar
      enddo

      do i=1,npnts_

         do j=1,npnts_
            if (j.eq.i) then
               delf(2*j-1) = one
               delf(2*j)   = zero
            else
               delf(2*j-1) = zero
               delf(2*j)   = zero
            endif
         enddo

         !---- do fourier transform from the position
         !     to the momentum representation

         call four1(delf,npnts_,1)

         do k=1,npnts_
	    delf_re = delf(2*k-1)
	    delf_im = delf(2*k)
	    delf(2*k-1) = -delf_im*pmom(k)/hbar
            delf(2*k)   =  delf_re*pmom(k)/hbar
         enddo

         !---- do inverse fourier transform from the momentum
         !     representation to the postion representation

         call four1(delf,npnts_,-1)

         do j=1,npnts_
            dx_(j,i) = delf(2*j-1)/dble(npnts_)
         enddo

      enddo

      deallocate (pmom,delf)

      return
   end subroutine griddx

   pure subroutine four1(data,nn,isign)
   !---1D FFT routine (Numerical Recipes)

      integer, intent(in) :: nn, isign
      real*8,  intent(inout), dimension(2*nn) :: data

      integer :: n, j, i, m, mmax, istep
      real*8  :: wr, wi, wpr, wpi, wtemp, theta, tempr, tempi

      n = 2*nn
      j = 1

      do i=1,n,2
         if (j.gt.i) then
            tempr = data(j)
            tempi = data(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         endif
         m = n/2
         do while ((m.ge.2).and.(j.gt.m))
            j = j - m
            m = m/2
         enddo
         j = j + m
      enddo

      mmax = 2
  
      do while (n.gt.mmax)
         istep = 2*mmax
         theta = 6.28318530717959d0/(isign*mmax)
         wpr = -2.d0*dsin(0.5d0*theta)**2
         wpi = dsin(theta)
         wr = 1.d0
         wi = 0.d0
         do m=1,mmax,2
            do i=m,n,istep
               j = i + mmax
               tempr = wr*data(j) - wi*data(j+1)
               tempi = wr*data(j+1) + wi*data(j)
               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            enddo
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         enddo
         mmax = istep
      enddo

      return
 
   end subroutine four1

   !============================================================================!
   pure subroutine calcham(npnts_,v_,hke_,ham_)
   !============================================================================!
   ! Initializes the Hamiltonian matrix on the grid
   !============================================================================!

      ! input/output variables
      integer, intent(in)                              :: npnts_
      real*8,  intent(in),    dimension(npnts_)        :: v_
      real*8,  intent(in),    dimension(npnts_,npnts_) :: hke_
      real*8,  intent(inout), dimension(npnts_,npnts_) :: ham_

      ! local variables
      integer :: i

      ! add potential and kinetic energy matrices to form
      ! the Hamiltonian matrix ham_

      ham_ = 0.d0
      do i=1,npnts_
         ham_(i,i) = v_(i)
      enddo

      ham_ = ham_ + hke_

      return

   end subroutine calcham

   !============================================================================!
   pure subroutine calcwf(n_,h_,c_,w_)
   !============================================================================!
   !  Calculates the wavefunctions and eigen values
   !  of the Hamiltonian matrix given on the grid
   !============================================================================!

      ! input
      integer, intent(in)                      :: n_
      real*8,  intent(inout), dimension(n_,n_) :: h_
      real*8,  intent(inout), dimension(n_,n_) :: c_
      real*8,  intent(inout), dimension(n_)    :: w_
      
      integer :: ierr
      real*8, dimension(n_) :: work1, work2

      c_ = 0.d0
      w_ = 0.d0
      work1 = 0.d0
      work2 = 0.d0
      ierr = 0

      !** diagonalize Hamiltonian matrix **!
      call rs(n_,n_,h_,w_,n_,c_,work1,work2,ierr)
      
      return

   end subroutine calcwf

end module quantum
