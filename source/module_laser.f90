module laser

   !---------------------------------------------------------------------
   !
   !  Contains the data and routines for pump laser characteristics
   !
   !  PUMP_SHAPE - spectral shape of the pump laser pulse
   !               0 - rectangular shape
   !               1 - Lorentzian shape
   !               2 - Gaussian shape
   !
   !  PUMP_ENERGY - energy at the center of the laser pulse (eV)
   !
   !  PUMP_WIDTH - full width at half maximum (FWHM) (eV)
   !
   !---------------------------------------------------------------------
   !
   !  $Author: souda $
   !  $Date: 2010-11-10 21:14:21 $
   !  $Revision: 5.3 $
   !  $Log: not supported by cvs2svn $
   !  Revision 5.2  2010/11/04 22:43:08  souda
   !  Next iteration... and two additional Makefiles for building the code with debug options.
   !
   !  Revision 5.1  2010/11/03 06:22:16  souda
   !  New module for laser pump characteristics
   !
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
   public :: lorentzian, gaussian
   public :: print_laser_spectrum

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

      !real(8) :: sigma, sigma_sq, prefactor

      value = 0.d0

      select case(pump_shape)

         case(0)

            if (energy .lt. pump_energy - pump_width/2.d0 .or. &
               &energy .gt. pump_energy + pump_width/2.d0) then
               value = 0.d0
            else
               value = 1.d0/pump_width
            endif

         case(1)   !- Lorentzian lineshape

            value = lorentzian(energy,pump_energy,pump_width)
            !value = (1.d0/pi)*(pump_width/2.d0)/((energy-pump_energy)**2 + pump_width**2/4.d0)

         case(2)

            value = gaussian(energy,pump_energy,pump_width)

            !sigma = 0.5d0*pump_width/sqrt(2.d0*log(2.d0))
            !sigma_sq = sigma*sigma
            !prefactor = 1.d0/sqrt(2.d0*pi)/sigma
            !value = prefactor*exp(-(energy-pump_energy)**2/(2.d0*sigma_sq))

      end select

   end function spectral_shape

   !---------------------------------------------------------------------
   ! Lorentzian function (normalized)
   !---------------------------------------------------------------------
   function lorentzian(x,x0,w) result(value)
      real(8), intent(in) :: x, x0, w
      real(8) :: value
      real(8) :: dx
      dx = x - x0
      value = (1.d0/pi)*(w/2.d0)/(dx*dx + w*w/4.d0)
   end function lorentzian

   !---------------------------------------------------------------------
   ! Gaussian function (normalized)
   !---------------------------------------------------------------------
   function gaussian(x,x0,w) result(value)
      real(8), intent(in) :: x, x0, w
      real(8) :: value
      real(8) :: dx, dx_sq, sigma, sigma_sq, prefactor
      dx = x - x0
      dx_sq = dx*dx
      sigma = 0.5d0*w/sqrt(2.d0*log(2.d0))
      sigma_sq = sigma*sigma
      prefactor = 1.d0/sqrt(2.d0*pi)/sigma
      value = prefactor*exp(-dx_sq/(2.d0*sigma_sq))
   end function gaussian


   !------------------------------------------------------------------------
   ! output spectral lineshape of the pump laser pulse to the external file
   !------------------------------------------------------------------------
   subroutine print_laser_spectrum(ichannel)
      integer, intent(in) :: ichannel

      integer, parameter :: npoints=200
      integer :: k
      real(8) :: energy, exc_en_min, exc_en_max, range, den, p

      !-- set the excitation energy range

      exc_en_min = pump_energy - 10.d0*pump_width
      if (exc_en_min.lt.0.d0) exc_en_min = 0.d0
      exc_en_max = pump_energy + 10.d0*pump_width
      range = exc_en_max - exc_en_min
      den = range/(npoints-1)

      !-- print the spectrum

      write(ichannel,'("#",71("="))')
      write(ichannel,'("#","   Pump laser power spectrum")')
      write(ichannel,'("#",71("-"))')
      write(ichannel,'("#",t10,"Energy (eV)",t35,"Normalized intensity")')
      write(ichannel,'("#",71("-"))')
      
      do k=1,npoints

         energy = exc_en_min + (k-1)*den
         p = spectral_shape(energy)

         write(ichannel,'(f20.6,t35,f15.6)') energy, p

      enddo

      write(ichannel,'("#",71("="))')

   end subroutine print_laser_spectrum



end module laser
