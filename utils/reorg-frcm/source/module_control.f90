module control

!===============================================================
!     Control parameters
!---------------------------------------------------------------
!     METHOD  - method of calculation of vibronic states
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     1 - using the DIABATIC approach: diagonalization
!         of the total Hamiltonian in the basis of the
!         products of DIABATIC electronic states (EVB)
!         and vibrational states (in the DIABATIC
!         EVB potentials).
!
!     2 - using the ADIABATIC approach: diagonalization
!         of the total Hamiltonian in the basis of the
!         products of ADIABATIC electronic states
!         and vibrational states (in the ADIABATIC
!         potentials). Note that this approach involves
!         the calculation of non-adiabatic "d" and "g"
!         derivative terms.
!
!     3 - using the DOUBLE ADIABATIC approach: the
!         non-adiabatic "d" and "g" derivative terms
!         are assumed to be zero.
!
!     GQUANT  - quantum (TRUE) vs. classical (FALSE)
!               description of the gating mode
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     MGQUANT  - method of incorporation of gating
!                vibrational states
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     1 - gating states are calculated for each basis
!         vibronic (electron-proton) state. Then the full
!         Hamiltonian is constructed including all
!         additional terms arising from non-adiabatic
!         interactions between proton and gating motion.
!         Final states are obtained by diagonalization
!         of the full Hamiltonian.
!
!     2 - gating states are calculated for each basis
!         vibronic (electron-proton) state. Non-adiabatic
!         interaction between vibronic and gating
!         states are ignored.
!         Final states are obtained by diagonalization
!         of the full Hamiltonian.
!
!     3 - gating states are calculated for each final
!         vibronic (electron-proton) state. Non-adiabatic
!         interaction between vibronic and gating
!         states are ignored.
!         Final states are obtained by diagonalization
!         of the full Hamiltonian.
!
!     IGAS  - type of the gas phase potential
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     1 - MM1:   MM EVB potential 1 (for Nocera-like systems)
!     2 - WATER: MM EVB potential for water cluster (Voth-Schmidt)
!     3 - MM10:  MM EVB potential 10 (for other systems)
!     4 - LEPS:  LEPS EVB potential
!
!     ISOLV  - type of the solvation model
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     1 - ELLIPSE: simple ellipsoidal model
!     2 - FRCM:    Frequency Resolved Cavity Model
!
!     IMINIM - minimization method for 2D and 3D surfaces
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     1 - Newton-Rafson
!     2 - Limited mempry BFGS method (LBFGS)
!
!     CHARGE - total charge of the solute
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!-------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!==================================================================

   implicit none
   public
   save

   integer :: method, mgquant, igas, isolv, iminim
   real*8, allocatable, dimension(:) :: charge
   logical :: gquant
   !===============================================================

contains

   subroutine alloc_charge(n_)
      integer, intent(in) :: n_
      allocate (charge(n_))
   end subroutine alloc_charge

   subroutine dealloc_charge
      if (allocated(charge)) deallocate (charge)
   end subroutine dealloc_charge

end module control
