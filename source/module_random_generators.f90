module random_generators

!--------------------------------------------------------------------
!  Random number generators
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-03-28 20:57:28 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!--------------------------------------------------------------------

   implicit none
   private

   !-- type of the random number generator to use in gaussdist routines
   !   =1 - ran2nr() routine from numerical recipes
   !   =2 - duni() routine from DL_POLY
   integer :: irantype=1

   !-- parameters for gaussdist

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

   !-- principal random seed
   !   dynamic, changed during each call to ran2nr
   !   if negative the seeds are reinitialized
   integer, public :: iseed

   !-- ran2nr seeds
   integer :: iy, idum2
   integer, dimension(ntab) :: iv

   !-- Random seeds for duni: i_seed, j_seed, k_seed, l_seed
   !   initial values of i_seed, j_seed, k_seed must be
   !   in the range from 1 to 178 (not all 1);
   !   initial value of l_seed must be in the range from 0 to 168.

   integer, public, save :: i_seed=12
   integer, public, save :: j_seed=34
   integer, public, save :: k_seed=56
   integer, public, save :: l_seed=78

   !-- Additional seeds for duni (not necessary to change)

   integer, save :: ir=97
   integer, save :: jr=33

   real(kind=4), dimension(97), save :: u

   real(kind=4), save :: c  =   362436.0/16777216.0
   real(kind=4), save :: cd =  7654321.0/16777216.0
   real(kind=4), save :: cm = 16777213.0/16777216.0
   real(kind=4), save :: uni

   interface set_duni_random_seeds
      module procedure set_duni_random_seeds_input
      module procedure set_duni_random_seeds_clock
   end interface

   interface set_random_seed
      module procedure set_random_seed_input
      module procedure set_random_seed_clock
      module procedure set_random_seed_pbsid
   end interface

   public :: set_duni_random_seeds, initialize_duni
   public :: set_random_seed
   public :: reset_random_seed
   public :: save_random_seeds
   public :: restore_random_seeds
   public :: ran2nr, duni
   public :: gaussdist
   public :: gaussdist_boxmuller

contains

   !--------------------------------------------------------------------
   subroutine set_random_seed_input(iseed_)
      integer, intent(in) :: iseed_
      iseed = iseed_
   end subroutine set_random_seed_input

   !---------------------------------------------------------------------
   subroutine set_random_seed_clock()
      integer, dimension(8) :: current_values
      call date_and_time(values=current_values)
      iseed = -(current_values(8) + 1)
   end subroutine set_random_seed_clock

   !---------------------------------------------------------------------
   subroutine set_random_seed_pbsid(pbsvar)

      character(len=*), intent(in) :: pbsvar
      character(len=240) :: varstring
      integer :: ierr, getenvqq, intvar, idot

      ierr = getenvqq(pbsvar,varstring)

      if (ierr.lt.0) then
         write(*,'(1x," *** Error code while getting the value of the environment variable ",a," : ",i4)') trim(pbsvar), ierr
         stop
      elseif (ierr.eq.0) then
         write(*,'(1x," *** (Random seed generation) The environment variable ",a," is not defined")') trim(pbsvar)
         stop
      endif

      !-- extract actual ID by cutting off the hostname
      idot = index(varstring,".")
      varstring = varstring(1:idot-1)

      !-- convert to integer
      read(varstring,'(i)') intvar

      !-- use three last digits
      iseed = -mod(intvar,1000)

   end subroutine set_random_seed_pbsid

   !--------------------------------------------------------------------
   subroutine reset_random_seed
      if (iseed.gt.0) iseed = -iseed
   end subroutine reset_random_seed

   !--------------------------------------------------------------------
   subroutine save_random_seeds(ichannel)
      !-- writes the seeds to the external checkpoint binary file
      !   (file should be opened as new)
      integer, intent(in) :: ichannel
      integer :: k
      write(ichannel) i_seed, j_seed, k_seed, l_seed, (u(k),k=1,97), uni
      write(ichannel) iseed, iy, idum2, (iv(k),k=1,ntab)
   end subroutine save_random_seeds

   !--------------------------------------------------------------------
   subroutine restore_random_seeds(ichannel)
      !-- reads the seeds from the external checkpoint binary file
      !   (file should exist and be opened as old)
      integer, intent(in) :: ichannel
      integer :: k
      read(ichannel) i_seed, j_seed, k_seed, l_seed, (u(k),k=1,97), uni
      read(ichannel) iseed, iy, idum2, (iv(k),k=1,ntab)
   end subroutine restore_random_seeds

   !--------------------------------------------------------------------
   function gaussdist() result(v)

      implicit none
      real(8) :: v
      real(4) :: r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, rsum
      real(8) :: randum, rrr, rr2

      v = 0.d0

      !-- use different random number generator depending on value of irantype
      !   irantype=1 => ran2nr from numerical recipes (from soo young)
      !   irantype=2 => dl_poly random number generator (from alexander)

      if (irantype.eq.1) then

          r1 = ran2nr()
          r2 = ran2nr()
          r3 = ran2nr()
          r4 = ran2nr()
          r5 = ran2nr()
          r6 = ran2nr()
          r7 = ran2nr()
          r8 = ran2nr()
          r9 = ran2nr()
         r10 = ran2nr()
         r11 = ran2nr()
         r12 = ran2nr()

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

   end function gaussdist
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   function gaussdist_boxmuller() result(v)

      implicit none

      real(8) :: v

      logical, save :: lnew=.false.
      real(8), save :: store
      real(4) :: wr1, wr2
      real(8) :: r, sqln, w1, w2, randum

      v = 0.d0

      if (lnew) then

         v = store
         lnew = .false.

      else

         !-- repeat until (r<1.0 and r<>0.0)
         10 continue

         if (irantype.eq.1) then
         
            wr1 = ran2nr()
            wr2 = ran2nr()
            
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

   end function gaussdist_boxmuller
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! numerical recipes: the art of scientific  computing (2nd ed.) 
   ! written by W. H. Press et. al.                                
   !--------------------------------------------------------------------
   ! Long period (lt.2E18) random generator of L'Ecuyer with Bays-Durham
   ! shuffle and added safe guards: Returns a uniform random deviate
   ! between 0.0 and 1.0 (exclusive of the endpoint values). Call with 
   ! iseed a negative integer to initialize; thereafter, do not alter iseed
   ! between successive deviates in a sequence. RNMX should approximate
   ! the largest floating value that is less than 1.
   ! [computation takes about 1.5 times longer than ran1nr]
   !--------------------------------------------------------------------
   function ran2nr() result(ran2)

      !-- based on Soo Young's code

      implicit none

      real(4) :: ran2
      
      integer :: j, k
      !integer, dimension(ntab), save :: iv=0
      !integer, save :: iy=0, idum2=123456789

      !-- when initializing
      !if (iy.eq.0) write(*,*) 'when iy=0 in ran2nr, iseed=',iseed

      if (iseed.le.0) then
         !write(*,*) 'restart the sequence of random number in ran2nr'
         !write(*,*) 'in ran2nr, iseed=',iseed
         !write(*,*)
         iseed = max(-iseed,1)     ! prevent iseed=0
         idum2 = iseed
         !-- load the shuffle table (after 8 warm-ups)         
         do j=ntab+8,1,-1
            k = iseed/iq1
            iseed = ia1*(iseed - k*iq1) - k*ir1
            if(iseed.lt.0) iseed = iseed + im1
            if(j.le.ntab) iv(j) = iseed
         enddo
         iy = iv(1)
      endif

      !-- compute iseed=mod(ia1*iseed,im1) & idum2=mod(ia2*idum2,im2)
      !   without overflows by Scharge's method

      k = iseed/iq1
      iseed = ia1*(iseed - k*iq1) - k*ir1
      if(iseed.lt.0) iseed = iseed + im1
      k = idum2/iq2
      idum2 = ia2*(idum2 - k*iq2) - k*ir2
      if(idum2.lt.0) idum2 = idum2 + im2
      
      j = 1 + iy/ndiv              ! in the range 1:ntab
      !-- iseed is shuffled, iseed & idum2 are combined to generate output      
      iy = iv(j) - idum2
      iv(j) = iseed
      if(iy.lt.1) iy = iy + imm1
      ran2 = min(am*iy,rnmx)     ! because users do not expect 
                                 ! endpoint values

   end function ran2nr
   !--------------------------------------------------------------------

   !---------------------------------------------------------------------
   subroutine set_duni_random_seeds_input(i_, j_, k_, l_)

      integer, intent(in) :: i_, j_, k_, l_
      integer :: i_tmp, j_tmp, k_tmp, l_tmp

      i_tmp = mod(abs(i_),178)
      j_tmp = mod(abs(j_),178)
      k_tmp = mod(abs(k_),178)
      l_tmp = mod(abs(l_),168)

      if (i_tmp.ne.0) then
         i_seed = i_tmp
      else
         i_seed = 178
      endif

      if (j_tmp.ne.0) then
         j_seed = j_tmp
      else
         j_seed = 178
      endif

      if (k_tmp.ne.0) then
         k_seed = k_tmp
      else
         k_seed = 178
      endif

      if (i_seed.eq.1.and.j_seed.eq.1.and.k_seed.eq.1) k_seed = 111

   end subroutine set_duni_random_seeds_input

   !---------------------------------------------------------------------
   subroutine set_duni_random_seeds_clock()

      integer, dimension(8) :: current_values
      
      call date_and_time(values=current_values)
      
      i_seed = 3*current_values(6) + 1
      j_seed = 3*current_values(7) + 1
      k_seed = mod(current_values(8),178) + 1
      l_seed = 7*current_values(5)

      if (i_seed.eq.1.and.j_seed.eq.1.and.k_seed.eq.1) k_seed = 111

   end subroutine set_duni_random_seeds_clock

   !---------------------------------------------------------------------
   subroutine initialize_duni()

      integer :: ii, jj, m_seed
      real(kind=4) :: s, t

      do ii=1,97

          s = 0.0
          t = 0.5

          do jj=1,24

             m_seed = mod(mod(i_seed*j_seed,179)*k_seed,179)
             i_seed = j_seed
             j_seed = k_seed
             k_seed = m_seed
             l_seed = mod(53*l_seed+1,169)
             if (mod(l_seed*m_seed,64).ge.32) s = s + t
             t = 0.5*t

          enddo

          u(ii) = s

      enddo

   end subroutine initialize_duni

   !*********************************************************************
   !     
   !     dl_poly random number generator based on the universal
   !     random number generator of marsaglia, zaman and tsang
   !     (stats and prob. lett. 9 (1990) 35-39.) it must be
   !     called once to initialise parameters u,c,cd,cm
   !     
   !     copyright daresbury laboratory 1992
   !     aurhor -  w.smith         july 1992
   !     
   !*********************************************************************
   function duni() result(random_number)

      real(8) :: random_number

      !-- calculate random number

      uni = u(ir) - u(jr)
      if (uni.lt.0.d0) uni = uni + 1.d0
      u(ir) = uni
      ir = ir - 1
      if (ir.eq.0) ir = 97
      jr = jr - 1
      if (jr.eq.0) jr = 97
      c = c - cd
      if (c.lt.0.d0) c = c + cm
      uni = uni - c
      if (uni.lt.0.d0) uni = uni + 1.d0
      random_number = dble(uni)

   end function duni


end module random_generators
