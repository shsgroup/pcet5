!=======================================================================
!  Parameters for Runge-Kutta methods
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-26 21:06:21 $
!  $Revision: 5.1 $
!  $Log: not supported by cvs2svn $
!
!=======================================================================
module rk_parameters

   implicit none
   public

   real(8), parameter, dimension(4) :: c_rk4=(/ 0.0d0, &
                                             &  0.5d0, &
                                             &  0.5d0, &
                                             &  1.d0   /)

   real(8), parameter, dimension(4) :: b_rk4=(/ 1.d0/6.d0, &
                                             &  1.d0/3.d0, &
                                             &  1.d0/3.d0, &
                                             &  1.d0/6.d0  /)

   real(8), parameter, dimension(4,4) :: a_rk4 = reshape((/ 0.0d0, 0.0d0, 0.0d0, 0.0d0,  &
                                                         &  0.5d0, 0.0d0, 0.0d0, 0.0d0,  &
                                                         &  0.0d0, 0.5d0, 0.0d0, 0.0d0,  &
                                                         &  0.0d0, 0.0d0, 1.0d0, 0.0d0   /), &
                                                         &  shape=(/4,4/), &
                                                         &  order=(/2,1/)  )

end module rk_parameters
