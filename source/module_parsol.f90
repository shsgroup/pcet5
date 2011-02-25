module parsol

!=======================================================================
!     Solvent and cavity parameters for solvation calculations
!-----------------------------------------------------------------------
!
!     Static model
!     ------------
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
!
!
!     Debye dynamical model
!     -----------------------------
!     f0      - inverse Pekar factor = 4*\pi*eps0*eps8/(eps0-eps8)
!     taud    - Debye relaxation time (ps)
!     taul    - longitudinal relaxation time (ps)
!
!     Debye-Onodera dynamical model
!     -----------------------------
!     f0      - inverse Pekar factor = 4*\pi*eps0*eps8/(eps0-eps8)
!     taud    - Debye relaxation time (ps)
!     taul    - longitudinal relaxation time (ps)
!     tau0    - Onodera time (ps)
!     tau0l   - Onodera longitudinal time (ps)
!
!     Debye-Onodera dynamical model with two relaxation periods
!     ---------------------------------------------------------
!
!     eps1    - additional dielectric constant
!     tau1    - first relaxation time (ps)
!     tau2    - second relaxation time (ps)
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-02-25 19:11:25 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!=======================================================================

   use cst

   implicit none
   public
   save

   logical :: symt, reddens, sympt, symet, nosymd
   real(8) :: eps0, eps8, kappa, delta, a, b, r, l0

   real(8) :: eps0_dyn, eps8_dyn, f0, taud, taul, tau0, tau0l, effmass
   real(8) :: eps1, tau1, tau2

contains

   !=======================================================================
   subroutine set_debye_model_parameters()
      taul  = eps8_dyn*taud/eps0_dyn
      effmass = 0.d0
   end subroutine set_debye_model_parameters

   !=======================================================================
   subroutine set_debye2_model_parameters()
      !-- (not implemented yet)
      taul  = eps8_dyn*taud/eps0_dyn
      effmass = 0.d0
   end subroutine set_debye2_model_parameters

   !=======================================================================
   subroutine set_onodera_model_parameters()
      taul  = eps8_dyn*taud/eps0_dyn
      tau0l = eps8_dyn*tau0/eps0_dyn
      effmass = f0*tau0*taul
   end subroutine set_onodera_model_parameters

   !=======================================================================
   subroutine set_onodera2_model_parameters()
      f0 = four*pi*eps0_dyn*eps8_dyn/(eps0_dyn-eps8_dyn)
      !-- to be added...
      !eps1 = 
      !tau1 =
      !tau2 =
      !...
   end subroutine set_onodera2_model_parameters

end module parsol
