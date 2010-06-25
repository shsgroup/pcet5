module monte_carlo

   use constants
   use rate_pars
   use flux
   
   implicit none
   private
   
   public :: mc_integrate1
   
   contains
   
   !--------------------------------------------------------------------
   subroutine mc_integrate1(n_, step_, n_accepted_, tau_min_, tau_max_, result_)
   
      implicit none

      ! input
      integer, intent(in) :: n_
      real(8), intent(in) :: step_
      
      ! output
      integer, intent(out) :: n_accepted_
      real(8), intent(out) :: tau_min_, tau_max_, result_
      
      ! local variables
      integer :: i
      real(8) :: rnumber1, rnumber2
      real(8) :: tau0, taui, delta
      real(8) :: energy_mc_0, energy_mc, function_mc, function_mc_0
      logical :: accept
      
      result_ = zero
      n_accepted_ = 0
      tau0 = zero
      energy_mc_0 = half*chi*tau0*tau0
      function_mc_0 = exp(p*(cos(tau0)-one))*cos_phase_argument(tau0)

      tau_max_ = -99999.d0
      tau_min_ =  99999.d0

      ! initialize pseudo random number generator
      call random_seed

      ! start Metropolis random walk
      do i=1,n_

         call random_number(rnumber1)
         taui = tau0 + (two*rnumber1-one)*step_

         energy_mc = half*chi*taui*taui
         function_mc = exp(p*(cos(taui)-one))*cos_phase_argument(taui)
         delta = energy_mc - energy_mc_0

         accept = .false.
         if (delta < zero) then
            accept = .true.
         else
            call random_number(rnumber2)
            if (exp(-delta) > rnumber2) accept = .true.
         endif

         if (accept) then
            n_accepted_ = n_accepted_ + 1
            if (taui > tau_max_) tau_max_ = taui
            if (taui < tau_min_) tau_min_ = taui
            tau0 = taui
            energy_mc_0 = energy_mc
            result_ = result_ + function_mc
         endif

      enddo

      result_ = sqrt(two*pi/chi)*result_/n_accepted_
      
      return
      
   end subroutine mc_integrate1
   
end module monte_carlo
