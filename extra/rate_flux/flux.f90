module flux

!==============================================================
!  computes the reactive flux (integrand)
!  Exp[-chi*t^2/2 + p*(Cos[t]-1) + i*(Sin[t]+theta*t]
!==============================================================
   use constants
   use rate_pars

   implicit none
   public

   ! temperature dependent flux parameters
   real(8), save :: beta, zeta, chi, p, q, theta

   contains

   !----------------------------------------------------------------------------
   subroutine set_flux_pars(temp_)
   
   ! Calculates the flux parameters at given temperature
   
      real(8), intent(in) :: temp_     ! temperature (K)

      beta = au2kcal/kb/temp_                                     ! atomic units
      zeta = 1.d0/tanh(half*beta*freq)
      chi = two*lambda_z/(beta*freq*freq)
      q = lambda_r + lambda_a - two*zeta*sqrt(lambda_r*lambda_a)
      q = q/freq
      p = zeta*(lambda_r+lambda_a) - two*sqrt(lambda_r*lambda_a)
      p = p/freq
      theta = (deltag + lambda_z)/freq

   end subroutine set_flux_pars

   !----------------------------------------------------------------------------
   function flux_real(tau_) result(flux_real_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: flux_real_

      ! local variables
      real(8) :: tau2, exp_factor, phase

      tau2 = tau_*tau_
      exp_factor = exp( -half*chi*tau2 + p*(cos(tau_)-one) )
      phase = q*sin(tau_) + theta*tau_
      flux_real_ = exp_factor*cos(phase)

   end function flux_real

   !----------------------------------------------------------------------------
   function flux_imag(tau_) result(flux_imag_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: flux_imag_

      ! local variables
      real(8) :: tau2, exp_factor, phase

      tau2 = tau_*tau_
      exp_factor = exp( -half*chi*tau2 + p*(cos(tau_)-one) )
      phase = q*sin(tau_) + theta*tau_
      flux_imag_ = exp_factor*sin(phase)

   end function flux_imag
   !----------------------------------------------------------------------------
   function flux_real_transformed(t_) result(flux_real_transformed_)

      ! input parameters      
      real(8), intent(in) :: t_       ! transformed scaled time (dimensionless)

      ! output
      real(8) :: flux_real_transformed_

      ! local variables
      real(8) :: x1, x2, t2, flux_re1, flux_re2

      x1 = (one-t_)/t_
      x2 = (t_-one)/t_
      t2 = t_*t_
      flux_re1 = flux_real(x1)
      flux_re2 = flux_real(x2)
      flux_real_transformed_ = (flux_re1 + flux_re2)/t2

   end function flux_real_transformed
   !----------------------------------------------------------------------------
   function damp_factor_full(tau_) result(damp_factor_full_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: damp_factor_full_
      
      ! local
      real(8) :: tau2

      tau2 = tau_*tau_
      damp_factor_full_ = exp( -half*chi*tau2 + p*(cos(tau_)-one) )

   end function damp_factor_full

   !----------------------------------------------------------------------------
   function damp_factor(tau_) result(damp_factor_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: damp_factor_
      
      ! local
      real(8) :: tau2

      tau2 = tau_*tau_
      damp_factor_ = exp( -half*chi*tau2 )

   end function damp_factor

   !----------------------------------------------------------------------------
   function phase_argument(tau_) result(phase_argument_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: phase_argument_

      phase_argument_ = q*sin(tau_) + theta*tau_

   end function phase_argument

   !----------------------------------------------------------------------------
   function cos_phase_argument(tau_) result(cos_phase_argument_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: cos_phase_argument_

      cos_phase_argument_ = cos(q*sin(tau_)+theta*tau_)

   end function cos_phase_argument
   !----------------------------------------------------------------------------
   function fct_dqh(t_) result(fct_dqh_)
   
   ! flux function fct for DQH64 integration
   !
   ! Exp[-t^2]*fct[t]
   
      ! input
      real(8), intent(in)  :: t_

      ! output
      real(8) :: fct_dqh_

      ! local variables
      real(8) :: tau

      tau = sqrt(two/chi)*t_
      fct_dqh_ = exp(p*(cos(tau)-one))*cos_phase_argument(tau)

   end function fct_dqh
      
end module flux
