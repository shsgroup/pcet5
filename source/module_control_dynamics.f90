module control_dynamics
!=======================================================================
!     Control parameters for solvent dynamics
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2012-04-06 22:38:46 $
!  $Revision: 5.4 $
!  $Log: not supported by cvs2svn $
!  Revision 5.3  2011/03/01 23:54:03  souda
!  Variable timestep for quantum propagation implemented (thanks to Sharon) - that fixes the problems with the conservation of the norm of the time-dependent wavefunction.
!
!  Revision 5.2  2010/11/04 22:43:08  souda
!  Next iteration... and two additional Makefiles for building the code with debug options.
!
!  Revision 5.1  2010/10/26 21:06:20  souda
!  new routines/modules
!
!
!=======================================================================
   implicit none
   public
   save

   character(len=10) :: solvent_model
   character(len=10) :: interpolation="LINEAR"

   logical :: mdqt=.false.        ! flag for MDQT dynamics
   logical :: phase_corr=.false.  ! flag for MDQT dynamics with phase correction (Shenvi-Subotnik)
   logical :: decoherence=.false. ! flag for MDQT dynamics with decoherence (AFSSH) (Landry-Shenvi-Subotnik)
   integer :: nsteps=100          ! number of steps
   integer :: ndump=1             ! dump trajectory every ndump steps
   integer :: ntraj=1             ! number of trajectories
   integer :: nqsteps             ! number of timesteps in TDSE
   integer :: maxnqsteps          ! maximum number of timesteps in TDSE
   real(8) :: tstep               ! time step for solvent dynamics (ps)
   real(8) :: qtstep              ! timestep for TDSE (ps)
   real(8) :: temp                ! temperature (K)
   real(8) :: dzeta=1.d0          ! coefficient for decoherence rate (A-FSSH, Eq.32)
!=======================================================================

contains

   !-------------------------------------
   ! set MDQT variables (TDSE)
   !-------------------------------------
   subroutine set_tdse_timestep(nqsteps_,maxnqsteps_,tstep_)
      integer, intent(in) :: nqsteps_, maxnqsteps_
      real(8), intent(in) :: tstep_
      nqsteps = nqsteps_
      maxnqsteps = maxnqsteps_
      qtstep = tstep_/real(nqsteps_)
   end subroutine set_tdse_timestep

end module control_dynamics
