module rates_et2

!===============================================================================
!   Rate constants for two-state ET model
!-------------------------------------------------------------------------------

   use cst
   use parsol
   use polyroots

   implicit none

contains

   !----------------------------------------------------------------------------
   !  Kramers-Grote-Hynes rate for Debye-1 relaxation
   !----------------------------------------------------------------------------
   function kgh_rate_debye1(temp,fr,fb,eb) result(rate)
      real(kind=8) :: temp, fr, fb, eb
      real(kind=8) :: rate
      rate = 1.d12*(0.5d0/pi)*sqrt(fr*fb)/(f0*taul)*exp(-eb/kb/temp)
   end function kgh_rate_debye1

   !----------------------------------------------------------------------------
   !  Kramers-Grote-Hynes rate for Debye-2 relaxation
   !----------------------------------------------------------------------------
   function kgh_rate_debye2(temp,fr,fb,eb) result(rate)
      real(kind=8) :: temp, fr, fb, eb
      real(kind=8) :: rate
      real(kind=8) :: omegab
      omegab = (fb - gamma)*taualpha - etax + sqrt(((fb-gamma)*taualpha - etax)**2 + 4.d0*fb*etax*taualpha)
      omegab = omegab/(2.d0*etax*taualpha)
      rate = 1.d12*(0.5d0/pi)*sqrt(fr/fb)*omegab*exp(-eb/kb/temp)
   end function kgh_rate_debye2

   !----------------------------------------------------------------------------
   !  Kramers-Grote-Hynes rate for Onodera-1 relaxation
   !----------------------------------------------------------------------------
   function kgh_rate_onodera1(temp,fr,fb,eb) result(rate)
      real(kind=8) :: temp, fr, fb, eb
      real(kind=8) :: rate
      real(kind=8) :: factor1, factor2
      factor1 = -0.5d0*etax/sqrt(effmass1*fb) + sqrt(1.d0 + etax*etax/4.d0/effmass1/fb)
      factor2 = sqrt(fr/effmass1)
      rate = 1.d12*(0.5d0/pi)*factor1*factor2*exp(-eb/kb/temp)
   end function kgh_rate_onodera1

   !----------------------------------------------------------------------------
   !  Kramers-Grote-Hynes rate for Onodera-2 relaxation
   !----------------------------------------------------------------------------
   function kgh_rate_onodera2(temp,fr,fb,eb) result(rate)
      real(kind=8) :: temp, fr, fb, eb
      real(kind=8) :: rate
      integer :: nroot, i
      real(kind=8) :: a, b, c, d, omegab
      complex(kind=8), dimension(3) :: root

      !-- Solve characteristic equation
      a = effmass1*taualpha
      b = effmass1 + etax*taualpha
      c = etax + gamma*taualpha - fb*taualpha
      d = -fb
      call cubic(a, b, c, d, root, nroot)

      omegab = 0.d0
      do i=1,nroot
         if (real(root(i)).gt.0.and.abs(dimag(root(i))).lt.1.d-12) then
            omegab = dreal(root(i))
            exit
         endif
      enddo

      rate = 1.d12*(0.5d0/pi)*sqrt(fr/fb)*omegab*exp(-eb/kb/temp)

   end function kgh_rate_onodera2

end module rates_et2
