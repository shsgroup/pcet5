module parsol
!=======================================================================
!     Solvent and cavity parameters for solvation calculations
!-----------------------------------------------------------------------
!     EPS0    - static dielectric constant
!     EPS8    - optical dielectric constant
!     KAPPA   - factor for VdW radii
!     DELTA   - the width of the layer between two cavities
!     A       - major semiaxis of the ellipsoidal cavity
!     B       - minor semiaxis of the ellipsoidal cavity
!     R       - interfocal distance of the ellipsoidal cavity
!     L0      - R/2A
!     SYMT    - SYMMETRIZATION OF T MATRICES
!     REDDENS - reduced density is used
!     SYMPT   - symmetrize T-matrices for PT only
!     SYMET   - symmetrize T-matrices for ET only
!     NOSYMD  - ???
!-----------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_parsol.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_parsol.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:04:40  souda
!  Initial PCET-4.0 Release
!
!
!=======================================================================
   implicit none
   public
   save

   logical :: symt, reddens, sympt, symet, nosymd
   real(8) :: eps0, eps8, kappa, delta, a, b, r, l0

   !=======================================================================

end module parsol
