module laser

   !---------------------------------------------------------------------
   !
   ! Contains the data and routines for pump laser characteristics
   !
   ! PUMP_SHAPE - spectral shape of the pump laser pulse
   !              0 - rectangular shape
   !              1 - Lorentzian shape
   !              2 - Gaussian shape
   !
   ! PUMP_ENERGY - energy at the center of the laser pulse (eV)
   !
   ! PUMP_WIDTH - full width at half maximum (FWHM) (eV)
   !
   !---------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2010-11-03 06:22:16 $
   !  $Revision: 5.1 $
   !  $Log: not supported by cvs2svn $
   !
   !---------------------------------------------------------------------

   use cst
   use solmat, only: dipole_moment_diab_x, &
                   & dipole_moment_diab_y, &
                   & dipole_moment_diab_z

   !---------------------------------------------------------------------
   implicit none
   private

   integer :: pump_shape
   real(8) :: pump_energy
   real(8) :: pump_width

   public :: set_laser
   public :: spectral_shape

contains

   !---------------------------------------------------------------------
   ! initialization routine
   !---------------------------------------------------------------------
   subroutine set_laser(shape,energy,width)
      integer, intent(in) :: shape
      real(8), intent(in) :: energy, width
      pump_shape = shape
      pump_energy = energy
      pump_width = width
   end subroutine set_laser

   !---------------------------------------------------------------------
   ! spectral lineshape of the pump laser pulse
   !---------------------------------------------------------------------

   function spectral_shape(energy) result(value)

      real(8), intent(in) :: energy     ! energy in eV
      real(8) :: value

      real(8) :: sigma, sigma_sq, prefactor

      value = 0.d0

      select case(pump_shape)

         case(0)

            if (energy .lt. pump_energy - pump_width/2.d0 .or. &
               &energy .gt. pump_energy + pump_width/2.d0) then
               value = 0.d0
            else
               value = 1.d0/pump_width
            endif

         case(1)

            value = (1.d0/pi)*(pump_width/2.d0)/((energy-pump_energy)**2 + pump_width**2/4.d0)

         case(2)

            sigma = 0.5d0*pump_width/sqrt(2.d0*log(2.d0))
            sigma_sq = sigma*sigma
            prefactor = 1.d0/sqrt(2.d0*pi)/sigma
            value = prefactor*exp(-(energy-pump_energy)**2/(2.d0*sigma_sq))

      end select

   end function spectral_shape

end module laser
