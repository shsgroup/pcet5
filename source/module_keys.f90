module keys
!===============================================================
!  Keywords, titles, job name, and date
!-----------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_keys.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_keys.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.3  2007/11/06 22:20:00  souda
!  new mmgen gas phase potential o-h o systems added
!
!  Revision 1.2  2004/01/23 19:10:00  souda
!  Longer JOB string
!
!  Revision 1.1.1.1  2004/01/13 20:03:06  souda
!  Initial PCET-4.0 Release
!
!
!
!===============================================================
   implicit none
   public
   save

   character(1024) :: keywrd
   character(32)   :: strdat
   character(80)   :: title, title2
   character(80)   :: job
   integer         :: ljob

   !===============================================================

end module keys

