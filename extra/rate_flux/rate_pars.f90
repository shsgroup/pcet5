module rate_pars

!-----------------------------------------------------------------------
! Contains subroutines and functions relevant to the calculation
! or initialization of the parameters used in the rate calculation
!-----------------------------------------------------------------------
   use constants

   implicit none
   public

   ! parameters
   
   real(8), save :: freq       ! frequency in atomic units (hbar*omega)
   real(8), save :: mass       ! mass in atomic units
   real(8), save :: alpha      ! Coupling parameter in 1/Bohr
   real(8), save :: dr         ! Delta R in a.u.
   real(8), save :: lambda_a   ! coupling reorganization energy
   real(8), save :: lambda_r   ! R-mode reorganization energy
   real(8), save :: lambda_z   ! solvent reorganization energy
   real(8), save :: deltag     ! reaction free energy (energy bias) in a.u.
   real(8), save :: v          ! constant coupling in a.u.

   contains
   
   subroutine read_pars
      
      real(8) :: freq_inp, mass_inp, lambda_z_inp, dr_inp, alpha_inp, deltag_inp, v_inp

      write(*,*)
      write(*,*) "====> Get ready to input parameters (Hit enter when ready)"
      read(*,*)
      write(*,*)
      write(*,'("R-mode frequency (1/cm)                 : ",$)')
      read(*,*) freq_inp
      write(*,'("R-mode reduced mass (Daltons)           : ",$)')
      read(*,*) mass_inp
      write(*,'("R-mode shift, Delta R (A)               : ",$)')
      read(*,*) dr_inp
      write(*,'("Coupling parameter alpha (1/A)          : ",$)')
      read(*,*) alpha_inp
      write(*,'("Solvent reorganization energy (kcal/mol): ",$)')
      read(*,*) lambda_z_inp
      write(*,'("Reaction free energy (kcal/mol)         : ",$)')
      read(*,*) deltag_inp
      write(*,'("Constant coupling (kcal/mol)            : ",$)')
      read(*,*) v_inp

      ! conversion to atomic units and additional calculations

      freq = freq_inp/ev2cm/au2ev
      mass = mass_inp*dalton
      dr = dr_inp/bohr2a
      alpha = alpha_inp*bohr2a

      lambda_z = lambda_z_inp/au2kcal
      lambda_a = half*alpha*alpha/mass
      lambda_r = half*mass*freq*freq*dr*dr

      deltag = deltag_inp/au2kcal
      v = v_inp/au2kcal

   end subroutine read_pars

end module rate_pars