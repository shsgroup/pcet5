module cst
!===================================================================
!     Constants and conversion factors
!-------------------------------------------------------------------
!  HBAR   - Planck constant in atomic units
!  KB     - Boltzmann constant in kcal/(mol*K)
!  DALTON - Atomic Mass Unit (in electron mass units)
!  HMASS  - Mass of proton (in electron mass units)
!  DMASS  - Mass of deiteron (in electron mass units)
!  BOHR2A - Bohr to Angstrom conversion factor
!  A2BOHR - Angstrom to Bohr conversion factor
!  AU2CAL - atomic units to kcal/mol conversion factor
!  CAL2AU - kcal/mol to atomic units conversion factor
!  EV2CM  - electron volts to inverse centimeters
!  CM2EV  - inverse centimeters to electron volts
!  AU2EV  - atomic units to electron volts
!  EV2AU  - electron volts to atomic units
!  HZ2CM  - Hertz to inverse centimeters
!  CM2HZ  - inverse centimeters to Hertz
!-------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:35
!  4.1
!  Exp
!  module_cst.f90,v 4.1 2010/06/25 20:02:35 souda Exp
!  module_cst.f90,v
!  Revision 4.1  2010/06/25 20:02:35  souda
!  Release 4.1
!
!  Revision 1.2  2004/05/15 03:32:45  souda
!  Added Borgis-Hynes rate routine
!
!  Revision 1.1.1.1  2004/01/13 19:37:29  souda
!  Initial PCET-4.0 Release
!
!
!===================================================================
   implicit none
   public
   save

   real(8), parameter :: zero  = 0.0d0,&
                       & half  = 0.5d0,&
                       & one   = 1.0d0,&
                       & two   = 2.0d0,&
                       & three = 3.0d0,&
                       & four  = 4.0d0,&
                       & five  = 5.0d0

   real(8) :: pi, kb, hbar, hbarps, dalton, hmass, dmass
   real(8) :: bohr2a, a2bohr, au2cal, cal2au, ev2cm, cm2ev
   real(8) :: au2ev, ev2au, hz2cm, cm2hz, e2, ev2cal, pico

   public :: init_cst

  !===================================================================
   contains

   subroutine init_cst
   !======================================================================C
   !     Initialization
   !======================================================================C

      ! Constants

      ! Pi
      pi = four*atan(one)

      ! Planck constant (atomic units)
      hbar  = one

      ! Planck constant (kcal/mol)*picosecond
      hbarps = 4.7685388d-2/pi

      ! Boltzmann constant in kcal/(mol*K)
      kb = 1.9872159d-3

      ! Proton and Deiteron masses (electron mass units)
      hmass = 1836.1527d0
      dmass = 3672.3054d0

      ! Dalton (electron mass units)
      dalton = 1822.8880d0

      ! Conversion factors
      bohr2a =    0.529167d0
      a2bohr =    one/bohr2a
      au2cal =  627.5095d0
      cal2au =    one/au2cal
      ev2cm  = 8065.54477d0
      cm2ev  =    one/ev2cm
      au2ev  =   27.2113834d0
      ev2au  =    one/au2ev
      hz2cm  =    0.333564d-10
      cm2hz  =    one/hz2cm
      e2     =   14.39982283d0
      ev2cal =   23.04512014d0
      pico   =    1.d12

   end subroutine init_cst

end module cst

