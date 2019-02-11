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
   write(6,'(1x,"        __\///____________________\/////////__\///////////////________\///________")')
   write(6,'(/)')
   write(6,'(1x,"   Version ",a)') trim(version)
   write(6,'(/)')
   write(6,'(1x,"   Hammes-Schiffer Group                                                          ")')
   write(6,'(1x,"   Department of Chemistry                                                        ")')
   write(6,'(1x,"   Yale University                                                                ")')
   write(6,'(1x,"   New Haven, CT 06520                                                            ")')
   write(6,'(1x,"__________________________________________________________________________________")')
   write(6,'(/)')
   write(6,'(1x,"   The program calculates solvent reorganization energy matrices for a set of     ")')
   write(6,'(1x,"   fixed charge distributions given as collections of partial charges.            ")')
   write(6,'(1x,"   The solvent is described in dielectric continuum approximation within the      ")')
   write(6,'(1x,"   Frequency Resolved Cavity Model (FRCM) or Ellipsoidal Cavity Model (ELCM).     ")')
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
