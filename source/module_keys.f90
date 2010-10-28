module keys

!===============================================================
!  Keywords, titles, job name, and date
!---------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
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

end module keys

