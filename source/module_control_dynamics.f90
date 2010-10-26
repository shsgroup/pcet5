module control_dynamics
!=======================================================================
!     Control parameters for solvent dynamics
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-26 21:06:20 $
!  $Revision: 5.1 $
!  $Log: not supported by cvs2svn $
!
!=======================================================================
   implicit none
   public
   save

   character(len=10) :: solvent_model
   character(len=10) :: interpolation="LINEAR"
   logical :: mdqt=.false.     ! flag for MDQT dynamics
   integer :: initial_state=1  ! initial state
   integer :: initial_set=1    ! initial set of states
   integer :: nstates          ! number of states in MDQT
   integer :: nsteps=100       ! number of steps
   integer :: ndump=1          ! dump trajectory every ndump steps
   integer :: ntraj=1          ! number of trajectories
   integer :: nqsteps          ! number of timesteps in TDSE
   real(8) :: tstep            ! time step for solvent dynamics (ps)
   real(8) :: qtstep           ! timestep for TDSE (ps)
   real(8) :: temp             ! temperature (K)
!=======================================================================

contains

   !-------------------------------------
   ! set MDQT variables (TDSE)
   !-------------------------------------
   subroutine set_tdse_timestep(nqsteps_,tstep_)
      integer, intent(in) :: nqsteps_
      real(8), intent(in) :: tstep_
      nqsteps = nqsteps_
      qtstep = tstep_/real(nqsteps_)
   end subroutine set_tdse_timestep

end module control_dynamics
