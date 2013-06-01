subroutine banner(version)

!-------------------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:35 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!-------------------------------------------------------------------------------

   implicit none
   character(len=*), intent(in) :: version

   !-- Prints banner

   write(6,'(//)')
   write(6,'(1x,"          2b, or not 2b: that is the question                                     ")')
   write(6,'(/)')
   write(6,'(1x,"                        William Shakespear, Hamlet (Act 3, Scene 1)               ")')
   write(6,'(/)')
   write(6,'(/)')
   write(6,'(1x,"___/\\\\\\\\\\\\\__________/\\\\\\\\\__/\\\\\\\\\\\\\\\__/\\\\\\\\\\\\\\\_ ",a)') trim(version)
   write(6,'(1x," __\/\\\/////////\\\_____/\\\////////__\/\\\///////////__\///////\\\/////__       ")')
   write(6,'(1x,"  __\/\\\_______\/\\\___/\\\/___________\/\\\___________________\/\\\_______      ")')
   write(6,'(1x,"   __\/\\\\\\\\\\\\\/___/\\\_____________\/\\\\\\\\\\\___________\/\\\_______     ")')
   write(6,'(1x,"    __\/\\\/////////____\/\\\_____________\/\\\///////____________\/\\\_______    ")')
   write(6,'(1x,"     __\/\\\_____________\//\\\____________\/\\\___________________\/\\\_______   ")')
   write(6,'(1x,"      __\/\\\______________\///\\\__________\/\\\___________________\/\\\_______  ")')
   write(6,'(1x,"       __\/\\\________________\////\\\\\\\\\_\/\\\\\\\\\\\\\\\_______\/\\\_______ ")')
   write(6,'(1x,"__________\///____________________\/////////__\///////////////________\///________")')
   write(6,'(/)')
   write(6,'(1x,"   Version ",a)') trim(version)
   write(6,'(/)')
   write(6,'(1x,"   Hammes-Schiffer Group                                                          ")')
   write(6,'(1x,"   Department of Chemistry                                                        ")')
   write(6,'(1x,"   University of Illinois at Urbana-Champaign                                     ")')
   write(6,'(1x,"   Urbana, IL 61801                                                               ")')
   write(6,'(1x,"__________________________________________________________________________________")')
   write(6,'(/)')
   write(6,'(1x,"   The program calculates free energies, nonadiabatic rates and solvent dynamics  ")')
   write(6,'(1x,"   for general proton-coupled electron transfer (PCET) reaction in polar solvents.")')
   write(6,'(1x,"   The solvent is described in dielectric continuum approximation, the solute is  ")')
   write(6,'(1x,"   described in terms of the four-state EVB model, the proton and electrons are   ")')
   write(6,'(1x,"   treated quantum mechanically.                                                  ")')
   write(6,'(/)')
   write(6,'(1x,"   REFERENCES                                                                     ")')
   write(6,'(/)')
   write(6,'(1x,"   1. A. V. Soudackov and S. Hammes-Schiffer.                                     ")')
   write(6,'(1x,"      Multistate Continuum Theory for Multiple Charge Transfer                    ")')
   write(6,'(1x,"      Reactions in Solution.                                                      ")')
   write(6,'(1x,"      J. Chem. Phys. 111 (1999) 4672-4687.                                        ")')
   write(6,'(/)')
   write(6,'(1x,"   2. A. V. Soudackov and S. Hammes-Schiffer.                                     ")')
   write(6,'(1x,"      Theoretical Study of Photoinduced Proton-Coupled Electron Transfer          ")')
   write(6,'(1x,"      through Asymmetric Salt Bridges.                                            ")')
   write(6,'(1x,"      J. Am. Chem. Soc. 121 (1999) 10598-10607.                                   ")')
   write(6,'(/)')
   write(6,'(1x,"__________________________________________________________________________________")')
   write(6,'(/)')

end subroutine banner
