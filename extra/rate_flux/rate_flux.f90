program rate_flux

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

   use constants
   use rate_pars
   use flux   
   use monte_carlo
   use ssplib

   implicit none

   ! parameters
   character(9), parameter :: rate_file = "rates.dat"   ! rates vs temperature
   real(8) :: pico = 1.d12

   ! variables
   real(8) :: prefactor, rate_exact_prefactor, temp, temp_low, temp_high, temp_step
   real(8) :: tau
   integer :: i, ntemp, itemp
   real(8) :: second, start_time, end_time, elapsed_time

   ! DQSF variables
   real(8) :: rate_exact_dqsf
   real(8) :: step_dqsf, error_dqsf
   integer :: ndim_dqsf
   real(8), allocatable, dimension(:) :: integrand_dqsf, integral_dqsf

   ! DQAGS variables
   real(8) :: rate_exact_dqags
   real(8) :: a_qags, b_qags, epsabs_qags, epsrel_qags
   integer :: limit_qags, last_qags
   real(8) :: result_qags, abserr_qags
   integer :: neval_qags, ier_qags
   real(8), allocatable, dimension(:) :: alist_qags, blist_qags, rlist_qags, elist_qags
   integer, allocatable, dimension(:) :: iord_qags
   real(8) :: rate_exact_qags

   ! DQH64 variables
   real(8) :: integral_dqh, rate_exact_dqh

   ! Monte Carlo variables
   real(8) :: integral_mc
   real(8) :: step_mc, tau_max, tau_min
   real(8) :: rate_exact_mc
   integer :: nsteps_mc, nsteps_accepted
   
   ! other variables
   real(8) :: lambda_total, sq_lambda, upper_term, sq_denom
   real(8) :: exp_term, exp_term1, exp_term2, exp_term3, sq_diff
   real(8) :: rate_high_t, rate_low_t

   ! read and set parameters
   call read_pars

   ! temperature independent rate prefactor (V^2/hbar) in (kcal/mol)/sec
   prefactor = pico*v*v/hbarps

   ! numerical integration parameters

   write(*,*)
   write(*,'("Numerical integration parameters")')

   write(*,*)
   write(*,'("Number of points in [0,1] for DQSF: ",$)')
   read(*,*) ndim_dqsf

   write(*,*)
   write(*,'("Absolute accuracy requested in DQAGS: ",$)')
   read(*,*) epsabs_qags
   write(*,'("Relative accuracy requested in DQAGS: ",$)')
   read(*,*) epsrel_qags
   write(*,'("Upperbound on the number of subintervals in DQAGS: ",$)')
   read(*,*) limit_qags
   if (limit_qags <= 1) limit_qags = 100
   a_qags = zero
   b_qags = one

   write(*,*)
   write(*,'("Monte Carlo maximum step: ",$)')
   read(*,*) step_mc
   write(*,'("Total number of Monte Carlo steps: ",$)')
   read(*,*) nsteps_mc

   ! temperature range

   write(*,*)
   write(*,'("Temperature range - from: ",$)')
   read(*,*) temp_low
   write(*,'("                      to: ",$)')
   read(*,*) temp_high
   write(*,'("        number of points: ",$)')
   read(*,*) ntemp

   temp_step = dabs(temp_high - temp_low)/(ntemp-1)

   ! open output file
   open(1,file=rate_file)
   write (1,'("#---------------------------------------------------------------")')
   write (1,'("# Log10[Rate(1/s)] vs. temperature")')
   write (1,'("#")')
   write (1,'("# columns:")')
   write (1,'("#")')
   write (1,'("# T(K)  100/T  DQSF  DQH64  QAGS  Monte-Carlo  High-T  Low-T")') 
   write (1,'("#---------------------------------------------------------------")')

   ! Loop over temperature

   temp = temp_low - temp_step

   do itemp=1,ntemp

      temp = temp + temp_step
      call set_flux_pars(temp)

      write(*,*)
      write(*,*) "*** Temperature",temp," K ***************************"
      write(*,*) "    kT: ",au2kcal/beta

      !------------!
      ! Exact rate !
      !------------!

      rate_exact_prefactor = (prefactor/freq)*dexp(two*lambda_a*zeta/freq)

      ! numerical calculation of the flux integral
      
      !---(1)--- DQSF (from SSPLIB)

      step_dqsf = one/(ndim_dqsf-1)
      allocate (integrand_dqsf(ndim_dqsf),integral_dqsf(ndim_dqsf))

      integrand_dqsf(1) = zero
      tau = zero + step_dqsf
      do i=2,ndim_dqsf
         integrand_dqsf(i) = flux_real_transformed(tau)
	 tau = tau + step_dqsf
      enddo

      start_time = second()
      call dqsf(step_dqsf,integrand_dqsf,integral_dqsf,ndim_dqsf)
      end_time = second()
      elapsed_time = end_time - start_time
      error_dqsf = abs(integral_dqsf(ndim_dqsf)-integral_dqsf(ndim_dqsf-1))
      rate_exact_dqsf = rate_exact_prefactor*integral_dqsf(ndim_dqsf)

      write(*,*)
      write(*,*) "--- DQSF integration"
      write(*,*) "    Integral: ",integral_dqsf(ndim_dqsf)
      write(*,*) "    Error   : ",error_dqsf
      write(*,*) "    Rate(1/s): ",rate_exact_dqsf
      write(*,*) "Elapsed time: ",elapsed_time

      deallocate (integrand_dqsf,integral_dqsf)

      !---(2)--- DQH64 (from SSPLIB)
      
      start_time = second()
      call dqh64(fct_dqh,integral_dqh)
      integral_dqh = integral_dqh*sqrt(two/chi)
      end_time = second()
      elapsed_time = end_time - start_time
      rate_exact_dqh = rate_exact_prefactor*integral_dqh

      write(*,*)
      write(*,*) "--- DQH64 integration"
      write(*,*) "    Integral: ",integral_dqh
      write(*,*) "    Rate(1/s): ",rate_exact_dqh
      write(*,*) "Elapsed time: ",elapsed_time

      !---(3)--- DQAGS (from QUADPACK)

      allocate (alist_qags(limit_qags),blist_qags(limit_qags),rlist_qags(limit_qags),elist_qags(limit_qags))
      allocate (iord_qags(limit_qags))
      start_time = second()
      call dqagse(flux_real_transformed,&
                & a_qags,b_qags,epsabs_qags,epsrel_qags,limit_qags,&
		& result_qags,abserr_qags,neval_qags,ier_qags,alist_qags,&
		& blist_qags,rlist_qags,elist_qags,iord_qags,last_qags)
      end_time = second()
      elapsed_time = end_time - start_time
      rate_exact_qags = rate_exact_prefactor*result_qags

      write(*,*)
      write(*,*) "--- DQAGS integration"
      write(*,*) "    Integral: ",result_qags
      write(*,*) "    Absolute error: ",abserr_qags
      write(*,*) "    ierror indicator: ",ier_qags
      write(*,*) "    Number of function evaluations: ",neval_qags
      write(*,*) "    Number of subintervals used: ",last_qags
      write(*,*) "    Rate(1/s): ",rate_exact_qags
      write(*,*) "Elapsed time: ",elapsed_time

      deallocate (alist_qags,blist_qags,rlist_qags,elist_qags)
      deallocate (iord_qags)

      !---(4)--- Simple Monte Carlo

      start_time = second()
      call mc_integrate1(nsteps_mc,step_mc,nsteps_accepted,tau_min,tau_max,integral_mc)
      end_time = second()
      elapsed_time = end_time - start_time

      rate_exact_mc = rate_exact_prefactor*integral_mc
      if (rate_exact_mc <= zero) rate_exact_mc = one

      write(*,*)
      write(*,*) "--- Monte Carlo integration"
      write(*,*) "    Total number of steps: ",nsteps_mc
      write(*,*) "    Number of accepted steps: ",nsteps_accepted
      write(*,*) "    Sampling range: ",tau_min,tau_max
      write(*,*) "    Integral: ",integral_mc
      write(*,*) "    Rate(1/s): ",rate_exact_mc
      write(*,*) "Elapsed time: ",elapsed_time

      !------------------!
      ! Approximate rate !
      !------------------!

      !---(5)--- High-temperature analytic expression

      lambda_total = lambda_z + lambda_r + lambda_a
      sq_lambda = sqrt(lambda_a*lambda_r)
      upper_term = deltag+lambda_z-four*sq_lambda/(beta*freq)
      sq_denom = lambda_total + beta*freq*sq_lambda*upper_term/lambda_total
      exp_term = upper_term*upper_term*(two-lambda_z/lambda_total)/(four*lambda_total)
      rate_high_t = prefactor*exp(four*lambda_a/beta/freq/freq)*sqrt(pi*beta/sq_denom)*exp(-beta*exp_term)
      write(*,*)
      write(*,*) "--- High-temperature limit"
      write(*,*) "    Rate(1/s): ",rate_high_t

      !---(6)--- Low-temperature analytic expression
      
      exp_term1 = two*lambda_a/freq
      sq_diff = sqrt(lambda_a) - sqrt(lambda_r)
      exp_term2 = sq_diff*sq_diff/freq
      exp_term3 = (deltag+lambda_z)*(deltag+lambda_z)/(four*lambda_z)
      rate_low_t = prefactor*exp(exp_term1)*exp(-exp_term2)*sqrt(pi*beta/lambda_z)*exp(-beta*exp_term3)
      write(*,*)
      write(*,*) "--- Low-temperature limit"
      write(*,*) "    Rate(1/s): ",rate_low_t

      !------- write to output file
      write(1,'(2f8.3,2x,7g12.6)') temp, 100.d0/temp,&
                            &   log10(rate_exact_dqsf),&
			    &   log10(rate_exact_dqh),&
			    &   log10(rate_exact_qags),&
			    &   log10(abs(rate_exact_mc)),&
			    &   log10(rate_high_t),&
			    &   log10(rate_low_t)

   enddo  ! end temperature loop

   close(1)
   
end program rate_flux