module rate_flux
!-------------------------------------------------------------------------------
! Tests various numerical integration methods to calculate
! the flux integral in the rate expression
!
! Integral to calculate is:
!
! Int[-inf,inf] Exp[-x*t^2 + p*(Cos[t]-1)]*Cos[q*Sin[t]+theta*t]
!
! Methods used:
!
! (1) DQSF - standard Gaussian quadrature for the whole function
!
! (2) DQAGS - adaptive algorithm QAGS (from QUADPACK) for the whole function
!
! (3) DQH64 (from SSPLIB) - for integrals Exp[-t^2]*F[x] - 64 points quadrature
!
! (5) Monte-Carlo integration with importance sampling (Doll)
!
! (6) High-temperature analytical expression
!
! (7) Low-temperature analytical expression
!
!-------------------------------------------------------------------------------
   use cst
   use monte_carlo
   use ssplib

   implicit none
   private

   ! free energy surface parameters (all in atomic units)
   
   real(8), save :: homega     ! R-mode frequency in atomic units (hbar*omega)
   real(8), save :: mass       ! R-mode reduced mass in atomic units
   real(8), save :: alpha      ! Coupling parameter in 1/Bohr
   real(8), save :: dr         ! Delta R in a.u.
   real(8), save :: lambda_a   ! coupling reorganization energy
   real(8), save :: lambda_r   ! R-mode reorganization energy
   real(8), save :: lambda_z   ! solvent reorganization energy
   real(8), save :: deltag     ! reaction free energy (energy bias) in a.u.
   real(8), save :: v0_sq      ! squared constant coupling in a.u.

   ! temperature dependent flux parameters
   real(8), save :: beta, zeta, chi, p, q, theta

   ! rate prefactors
   real(8) :: prefactor, rate_exact_prefactor
   
   public :: set_flux_pars
   !public :: rate_exact_dqsf, rate_exact_dqh, rate_exact_qags, rate_exact_mc, rate_high, rate_low
   public :: rate_exact, rate_high, rate_low

   interface rate_exact
      module procedure rate_exact_dqsf, &
                     & rate_exact_dqh,  &
                     & rate_exact_qags, &
                     & rate_exact_mc
   end interface

   contains

   !----------------------------------------------------------------------------
   subroutine set_flux_pars(homega_,mass_,lambda_r_,alpha_,lambda_z_,deltag_,v0_sq_,temp_)
   
      real(8), intent(in) :: homega_, mass_, lambda_z_, lambda_r_, alpha_, deltag_, v0_sq_
      real(8), intent(in) :: temp_     ! temperature (K)

      ! free energy surface parameters

      homega = homega_*cm2ev*ev2au
      mass = mass_*dalton
      alpha = alpha_*bohr2a

      lambda_z = lambda_z_*cal2au
      lambda_a = half*alpha*alpha/mass
                                                !lambda_r = half*mass*homega*homega*dr*dr
      lambda_r = lambda_r_*cal2au

      deltag = deltag_*cal2au
      v0_sq = v0_sq_*cal2au*cal2au

      ! temperature independent rate prefactor (V^2/hbar) in a.u./sec
      prefactor = pico*v0_sq/(hbarps*cal2au)

      ! flux parameters at given temperature
      beta = au2cal/kb/temp_                                     ! atomic units
      zeta = one/tanh(half*beta*homega)
      chi = two*lambda_z/(beta*homega*homega)
      q = lambda_r + lambda_a - two*zeta*sqrt(lambda_r*lambda_a)
      q = q/homega
      p = zeta*(lambda_r+lambda_a) - two*sqrt(lambda_r*lambda_a)
      p = p/homega
      theta = (deltag + lambda_z)/homega

      ! temperature dependent rate prefactor (V^2/hbar) in 1/sec
      rate_exact_prefactor = (prefactor/homega)*dexp(two*lambda_a*zeta/homega)

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
   function trig_factor(tau_) result(trig_factor_)

      ! input parameters      
      real(8), intent(in) :: tau_       ! frequency scaled time (dimensionless)

      ! output
      real(8) :: trig_factor_
      
      ! local
      real(8) :: tau2

      tau2 = tau_*tau_
      trig_factor_ = exp(p*(cos(tau_)-one))*cos_phase_argument(tau_)

   end function trig_factor

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
      fct_dqh_ = sqrt(two/chi)*trig_factor(tau)

   end function fct_dqh

   !----------------------------------------------------------------------------
   ! exact rates
   !----------------------------------------------------------------------------
   subroutine rate_exact_dqsf(ndim_, error_, result_)
   
      ! input
      integer, intent(in) :: ndim_         ! number of integration points
      
      ! output
      real(8), intent(out) :: error_       ! estimated error of integration
      real(8), intent(out) :: result_      ! final rate (1/sec)

      ! local variables
      real(8) :: step, tau
      integer :: i
      real(8), allocatable, dimension(:) :: integrand, integral

      step = one/(ndim_-1)
      allocate (integrand(ndim_),integral(ndim_))

      integrand(1) = zero
      tau = step
      do i=2,ndim_
         integrand(i) = flux_real_transformed(tau)
	 tau = tau + step
      enddo

      call dqsf(step,integrand,integral,ndim_)
      deallocate (integrand,integral)
      error_ = abs(integral(ndim_)-integral(ndim_-1))
      result_ = rate_exact_prefactor*integral(ndim_)

      return

   end subroutine rate_exact_dqsf

   !----------------------------------------------------------------------------
   subroutine rate_exact_dqh(result_)

      real(8), intent(out) :: result_

      real(8) :: integral

      call dqh64(fct_dqh,integral)
      result_ = rate_exact_prefactor*integral
      
      return
      
   end subroutine rate_exact_dqh

   !----------------------------------------------------------------------------
   subroutine rate_exact_qags(epsabs_ , epsrel_,&
                            & limit_, abserr_, neval_, ier_, alist_, blist_,&
                            & rlist_ , elist_, iord_, last_, result_)

      ! input
      real(8), intent(in) :: epsabs_, epsrel_
      integer, intent(in) :: limit_

      ! output
      real(8), intent(out) :: result_, abserr_
      integer, intent(out) :: neval_, ier_, last_
      real(8), dimension(limit_), intent(out) :: alist_, blist_, rlist_, elist_
      integer, dimension(limit_), intent(out) :: iord_

      ! local variables
      real(8) :: integral, a, b

      a = zero
      b = one

      call dqagse(flux_real_transformed,a,b,epsabs_,epsrel_,limit_,&
		& integral,abserr_,neval_,ier_,alist_,&
		& blist_,rlist_,elist_,iord_,last_)
      result_ = rate_exact_prefactor*integral

      if (ier_ .ne. 0) then

         write(*,*)
         write(*,*) "Abnormal termination of DQAGS integration routine: ier=",ier_
         write(*,*) "Requested accuracy is not achieved"

         select case(ier_)

            case(1)
               write(*,*) "Maximum number of subdivisions has been achieved"

            case(2)
               write(*,*) "Occurrence of round-off error is detected"

            case(3)
               write(*,*) "Extremely bad integrand behaviour occurs at some points"

            case(4)
               write(*,*) "The algorithm does not converge"

            case(5)
               write(*,*) "The interval is probably divergent or slowly convergent"

            case(6)
               write(*,*) "Input is invalid: epsabs < 0 and epsrel < max(50*rel.mach.acc.,0.5d-28)"

         end select

      endif

      return

   end subroutine rate_exact_qags

   !----------------------------------------------------------------------------
   subroutine rate_exact_mc(nsteps_, maxstep_, nsteps_accepted_, tau_min_, tau_max_, result_)

      ! input
      real(8), intent(in) :: maxstep_
      integer, intent(in) :: nsteps_
      
      ! output
      integer, intent(out) :: nsteps_accepted_
      real(8), intent(out) :: tau_min_, tau_max_, result_

      ! local variables
      real(8) :: integral, distr_norm

      nsteps_accepted_ = 0
      tau_min_ = zero
      tau_max_ = zero
      result_ = zero
      distr_norm = sqrt(two*pi/chi)

      call mc_integrate1(damp_factor, distr_norm, trig_factor,&
                       & nsteps_, maxstep_, nsteps_accepted_,&
                       & tau_min_, tau_max_, integral)

      if (integral <= zero) then
         write(*,*)
         write(*,*) "Monte Carlo integration did not converge..."
         result_ = zero
      else
         result_ = rate_exact_prefactor*integral
      endif

      return

   end subroutine rate_exact_mc

   !----------------------------------------------------------------------------
   ! asymptotic rates
   !----------------------------------------------------------------------------
   subroutine rate_high(result_)
   
      real(8), intent(out) :: result_

      ! local variables
      real(8) :: lambda_total, sq_lambda, upper_term, sq_denom
      real(8) :: exp_term

      lambda_total = lambda_z + lambda_r + lambda_a
      sq_lambda = sqrt(lambda_a*lambda_r)
      upper_term = deltag+lambda_z-four*sq_lambda/(beta*homega)
      !--- sq_denom = lambda_total + beta*homega*sq_lambda*upper_term/lambda_total ---(old version of 2005 ????)
      sq_denom = lambda_total
      exp_term = upper_term*upper_term/(four*lambda_total)
      result_ = prefactor*exp(four*lambda_a/beta/homega/homega)*sqrt(pi*beta/sq_denom)*exp(-beta*exp_term)

      return
      
   end subroutine rate_high

   !----------------------------------------------------------------------------
   subroutine rate_low(result_)
   
      real(8), intent(out) :: result_

      ! local variables
      real(8) :: exp_term1, exp_term2, exp_term3, sq_diff

      exp_term1 = two*lambda_a/homega
      sq_diff = sqrt(lambda_a) - sqrt(lambda_r)
      exp_term2 = sq_diff*sq_diff/homega
      exp_term3 = (deltag+lambda_z)*(deltag+lambda_z)/(four*lambda_z)
      result_ = prefactor*exp(exp_term1)*exp(-exp_term2)*sqrt(pi*beta/lambda_z)*exp(-beta*exp_term3)

      return
      
   end subroutine rate_low
   
end module rate_flux
