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
!  $Date: 2012-03-13 22:03:04 $
!  $Revision: 5.4 $
!  $Log: not supported by cvs2svn $
!  Revision 5.3  2011/02/25 19:11:25  souda
!  Now using a separate set of dielectric constant for solvent dynamics.
!
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

   real(8) :: eps0_dyn, eps8_dyn, f0, taud, taul, tau0, tau0l, effmass1, effmass2
   real(8) :: eps1_dyn, tau1, tau2
   real(8) :: taualpha, gamma, etax, etay

contains

   !=======================================================================
   subroutine set_debye_model_parameters()
      taul  = eps8_dyn*taud/eps0_dyn
      !f0 = four*pi*eps0_dyn*eps8_dyn/(eps0_dyn - eps8_dyn)
      effmass1 = 0.d0
      effmass2 = 0.d0
   end subroutine set_debye_model_parameters

   !=======================================================================
   subroutine set_debye2_model_parameters()

      real(8) :: alpha, alpha2, alpha3, gamma1

      !f0 = four*pi*eps0_dyn*eps8_dyn/(eps0_dyn - eps8_dyn)

      alpha = (eps1_dyn - eps8_dyn)*tau1 + (eps0_dyn - eps1_dyn)*tau2
      alpha2 = alpha*alpha

      etax = 4.d0*pi*eps8_dyn*eps8_dyn*tau1*tau2/alpha

      taualpha = alpha/(eps0_dyn - eps8_dyn)

      gamma1 = (eps0_dyn - eps1_dyn)*(eps1_dyn - eps8_dyn)*(tau1 - tau2)*(tau1 - tau2)
      gamma = 4.d0*pi*eps8_dyn*eps8_dyn*gamma1/(eps0_dyn - eps8_dyn)/alpha2

      etay = gamma*taualpha

      effmass1 = 0.d0
      effmass2 = 0.d0

   end subroutine set_debye2_model_parameters

   !=======================================================================
   subroutine set_onodera_model_parameters()
      taul  = eps8_dyn*taud/eps0_dyn
      tau0l = eps8_dyn*tau0/eps0_dyn
      !f0 = four*pi*eps0_dyn*eps8_dyn/(eps0_dyn - eps8_dyn)
      if (tau0.gt.0) then
         effmass1 = f0*tau0*taul
         effmass2 = f0*tau0*taul
      endif
   end subroutine set_onodera_model_parameters

   !=======================================================================
   subroutine set_onodera2_model_parameters()

      real(8) :: eta_up, alpha, alpha2, alpha3, gamma1, gamma2

      !f0 = four*pi*eps0_dyn*eps8_dyn/(eps0_dyn - eps8_dyn)

      alpha = (eps1_dyn - eps8_dyn)*tau1 + (eps0_dyn - eps1_dyn)*tau2                                                                            
      alpha2 = alpha*alpha                                                                                                   
      alpha3 = alpha*alpha2                                                                                                  

      if (tau0.gt.0) then
         effmass1 = 4.d0*pi*eps8_dyn*eps8_dyn*tau0*tau1*tau2/alpha
         effmass2 = effmass1
      endif

      eta_up = (eps0_dyn - eps1_dyn)*(tau0 + tau1)*tau2*tau2 + (eps1_dyn - eps8_dyn)*(tau0 + tau2)*tau1*tau1
      etax = 4.d0*pi*eps8_dyn*eps8_dyn*eta_up/alpha2

      taualpha = alpha/(eps0_dyn - eps8_dyn)

      gamma1 = (eps0_dyn - eps1_dyn)*(eps1_dyn - eps8_dyn)*(tau1 - tau2)*(tau1 - tau2)
      gamma2 = taualpha - tau0
      gamma = 4.d0*pi*eps8_dyn*eps8_dyn*gamma1*gamma2/alpha3

      etay = gamma*taualpha

   end subroutine set_onodera2_model_parameters

end module parsol
