module keys

!===============================================================
!  Keywords, titles, job name, and date
!---------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2011-03-28 20:57:59 $
!  $Revision: 5.3 $
!  $Log: not supported by cvs2svn $
!  Revision 5.2  2010/10/28 21:29:36  souda
!  First (working and hopefully bug-free) source of PCET 5.x
!
!
!===============================================================

   implicit none
   public
   save

   character(1024) :: keywrd
   character(32)   :: strdat
   character(160)  :: title, title2
   character(160)  :: job
   integer         :: ljob

end module keys

