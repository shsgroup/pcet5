module random_generators

!--------------------------------------------------------------------
!  Random number generators
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!--------------------------------------------------------------------

   implicit none
   private

   real(8), parameter :: a1=3.949846138d0
   real(8), parameter :: a3=0.252408784d0
   real(8), parameter :: a5=0.076542912d0
   real(8), parameter :: a7=0.008355968d0
   real(8), parameter :: a9=0.029899776d0

   !-- parameters from ran2nr routine

   integer, parameter :: im1=2147483563
   integer, parameter :: im2=2147483399
   integer, parameter :: imm1=im1-1
   integer, parameter :: ia1=40014
   integer, parameter :: ia2=40692
   integer, parameter :: iq1=53668
   integer, parameter :: iq2=52774
   integer, parameter :: ir1=12211
   integer, parameter :: ir2=3791
   integer, parameter :: ntab=32
   integer, parameter :: ndiv=1+imm1/ntab
   real(4), parameter :: am=1.0/im1
   real(4), parameter :: eps=1.2e-7
   real(4), parameter :: rnmx=1.-eps

   !-- initial random seed
   integer, public :: iseed

   public :: set_random_seed
   public :: ran2nr
   public :: gaussdist
   public :: gaussdist_boxmuller

contains

   !--------------------------------------------------------------------
   subroutine set_random_seed(iseed_)
      integer, intent(in) :: iseed_
      iseed = iseed_
   end subroutine set_random_seed

   !--------------------------------------------------------------------
   subroutine gaussdist(v,idum)

      implicit none
      integer, intent(inout) :: idum
      real(8), intent(out)   :: v
      integer :: irantype
      real(4) :: r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, rsum
      real(8) :: randum, rrr, rr2

      !-- use different random number generator depending on value of irantype
      !   irantype=1 => ran2 from numerical recipes (from soo young)
      !   irantype=2 => dl_poly random number generator (from alexander)

      irantype = 1

      if (irantype.eq.1) then

         call ran2nr(r1,idum)
         call ran2nr(r2,idum)
         call ran2nr(r3,idum)
         call ran2nr(r4,idum)
         call ran2nr(r5,idum)
         call ran2nr(r6,idum)
         call ran2nr(r7,idum)
         call ran2nr(r8,idum)
         call ran2nr(r9,idum)
         call ran2nr(r10,idum)
         call ran2nr(r11,idum)
         call ran2nr(r12,idum)

      else

         randum = duni()
         r1 = duni()
         r2 = duni()
         r3 = duni()
         r4 = duni()
         r5 = duni()
         r6 = duni()
         r7 = duni()
         r8 = duni()
         r9 = duni()
         r10 = duni()
         r11 = duni()
         r12 = duni()

      endif

      rsum = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11 + r12
      rrr = (real(rsum,8) - 6.0d0)/4.0d0
      rr2 = rrr*rrr
      v = rrr*(a1 + rr2*(a3 + rr2*(a5 + rr2*(a7 + rr2*a9))))

   end subroutine gaussdist
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   subroutine gaussdist_boxmuller(v, idum)

      implicit none

      integer, intent(inout) :: idum
      real(8), intent(out)   :: v

      logical, save :: lnew=.false.
      real(8), save :: store
      integer :: irantype
      real(4) :: wr1, wr2
      real(8) :: r, sqln, w1, w2, randum

      irantype=1

      if (lnew) then

         v = store
         lnew = .false.

      else

         !-- repeat until (r<1.0 and r<>0.0)
         10 continue

         if (irantype.eq.1) then
         
            call ran2nr(wr1,idum)  
            call ran2nr(wr2,idum)  
            
         else

            randum = duni()
            wr1 = duni()
            wr2 = duni()

         endif

         w1 = 2.d0*real(wr1,8) - 1.d0
         w2 = 2.d0*real(wr2,8) - 1.d0
         
         r = w1*w1 + w2*w2
         
         if ((r.GE.1.d0) .OR. (r.EQ.0.d0)) goto 10
         !-- end of repeat until loop
         
         sqln = sqrt(-2.d0*log(r)/r)
         v = w1*sqln
         store = w2*sqln
         lnew = .true.

      endif

   end subroutine gaussdist_boxmuller
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! numerical recipes: the art of scientific  computing (2nd ed.) 
   ! written by W. H. Press et. al.                                
   !--------------------------------------------------------------------
   ! Long period (lt.2E18) random generator of L'Ecuyer with Bays-Durham
   ! shuffle and added safe guards: Returns a uniform random deviate
   ! between 0.0 and 1.0 (exclusive of the endpoint values). Call with 
   ! idum a negative integer to initialize; thereafter, do not alter idum
   ! between successive deviates in a sequence. RNMX should approximate
   ! the largest floating value that is less than 1.
   ! [computation takes about 1.5 times longer than ran1nr(idum)]
   !--------------------------------------------------------------------
   subroutine ran2nr(ran2,idum)

      !-- based on Soo Young's code

      implicit none

      integer, intent(inout) :: idum
      real(4), intent(out)   :: ran2
      
      integer :: j, k
      integer, dimension(ntab), save :: iv=0
      integer, save :: iy=0, idum2=123456789

      !-- when initializing
      !if (iy.eq.0) write(*,*) 'when iy=0 in ran2nr, idum=',idum

      if (idum.le.0) then
         !write(*,*) 'restart the sequence of random number in ran2nr'
         !write(*,*) 'in ran2nr, idum=',idum
         !write(*,*)
         idum = max(-idum,1)     ! prevent idum=0
         idum2 = idum
         !-- load the shuffle table (after 8 warm-ups)         
         do j=ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum - k*iq1) - k*ir1
            if(idum.lt.0) idum = idum + im1
            if(j.le.ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif

      !-- compute idum=mod(ia1*idum,im1) & idum2=mod(ia2*idum2,im2)
      !   without overflows by Scharge's method

      k = idum/iq1
      idum = ia1*(idum - k*iq1) - k*ir1
      if(idum.lt.0) idum = idum + im1
      k = idum2/iq2
      idum2 = ia2*(idum2 - k*iq2) - k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      
      j = 1 + iy/ndiv              ! in the range 1:ntab
      !-- idum is shuffled, idum & idum2 are combined to generate output      
      iy = iv(j) - idum2
      iv(j) = idum
      if(iy.lt.1) iy = iy + imm1
      ran2 = min(am*iy,rnmx)     ! because users do not expect 
                                 ! endpoint values

   end subroutine ran2nr
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !     
   !     dl_poly random number generator based on the universal
   !     random number generator of marsaglia, zaman and tsang
   !     (stats and prob. lett. 8 (1990) 35-39.) it must be
   !     called once to initialise parameters u,c,cd,cm
   !     
   !     copyright daresbury laboratory 1992
   !     author -  w.smith         july 1992
   !     
   !--------------------------------------------------------------------
   function duni() result(random_number)

      logical, save :: new=.true.
      real(8) :: random_number

      integer :: i, j, k, l, m, ii, jj
      real(4) :: s, t

      integer, save :: ir, jr
      real(4), save :: u(97), c, cd, cm, uni

      if (new) then

         !--  initial values of i,j,k must be in range 1 to 178 (not all 1)
         !    initial value of l must be in range 0 to 168.

         i = 12
         j = 34
         k = 56
         l = 78

         ir = 97
         jr = 33

         new = .false.

         do ii=1,97
            s = 0.0
            t = 0.5
            do jj=1,24
               m = mod(mod(i*j,179)*k,179)
               i = j
               j = k
               k = m
               l = mod(53*l+1,169)
               if(mod(l*m,64).ge.32) s = s + t
               t = 0.5*t
            enddo
            u(ii) = s
         enddo

         c  =   362436.0/16777216.0
         cd =  7654321.0/16777216.0
         cm = 16777213.0/16777216.0

      else

         !-- calculate random number

         uni = u(ir) - u(jr)
         if(uni.lt.0.0) uni = uni + 1.0
         u(ir) = uni
         ir = ir - 1
         if(ir.eq.0) ir = 97
         jr = jr - 1
         if(jr.eq.0) jr = 97
         c = c - cd
         if(c.lt.0.0) c = c + cm
         uni = uni - c
         if(uni.lt.0.0) uni = uni + 1.0
         random_number = real(uni,8)

      endif

   end function duni
   !--------------------------------------------------------------------
 

end module random_generators
