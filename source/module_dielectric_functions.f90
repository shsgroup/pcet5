module dielectric_functions

   !---------------------------------------------------------------------
   ! Contains the routines for frequency dependent dielectric functions
   ! The following models are implemented:
   ! - Debye dielectric function
   ! - combination of Debye functions
   ! - Onodera refinement of the Debye function
   ! - Onodera model with two relaxation times
   !---------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2011-02-25 19:11:25 $
   !  $Revision: 5.2 $
   !  $Log: not supported by cvs2svn $
   !  Revision 5.1  2010/10/26 21:06:20  souda
   !  new routines/modules
   !
   !
   !---------------------------------------------------------------------

   use cst
   use parsol

   !---------------------------------------------------------------------
   implicit none
   private

   complex(kind=8), parameter :: ii=(zero,one)

   public :: debye_epsilon
   public :: debye_f
   public :: onodera_epsilon
   public :: onodera_f

contains

   !---------------------------------------------------------------------
   ! complex-valued Debye dielectric function
   !---------------------------------------------------------------------
   function debye_epsilon(omega) result(epsilon)
      real(kind=8), intent(in) :: omega    ! frequency (1/ps)
      complex(kind=8) :: epsilon
      epsilon = eps8_dyn + (eps0_dyn-eps8_dyn)/(one - ii*omega*taud)
   end function debye_epsilon

   !---------------------------------------------------------------------
   ! complex-valued Debye f(omega) = 1/epsilon(omega)
   !---------------------------------------------------------------------
   function debye_f(omega) result(f)
      real(kind=8), intent(in) :: omega    ! frequency (1/ps)
      complex(kind=8) :: f
      f = f0*(one - ii*omega*taul)
   end function debye_f

   !---------------------------------------------------------------------
   ! complex-valued Onodera dielectric function
   !---------------------------------------------------------------------
   function onodera_epsilon(omega) result(epsilon)
      real(kind=8), intent(in) :: omega    ! frequency (1/ps)
      complex(kind=8) :: epsilon
      epsilon = eps8_dyn + (eps0_dyn-eps8_dyn)/((one - ii*omega*tau0)*(one - ii*omega*taud))
   end function onodera_epsilon

   !---------------------------------------------------------------------
   ! complex-valued Onodera f(omega) = 1/epsilon(omega)
   !---------------------------------------------------------------------
   function onodera_f(omega) result(f)
      real(kind=8), intent(in) :: omega    ! frequency (1/ps)
      complex(kind=8) :: f
      complex(kind=8) :: eps
      eps = onodera_epsilon(omega)
      f = four*pi*eps*eps8_dyn/(eps-eps8_dyn)
   end function onodera_f

end module dielectric_functions
