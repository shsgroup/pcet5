subroutine getmin(in,nrmin,nsmin,rzpi,rzei,rzpj,rzej)
!===================================================================C
!  Reads Minima of the states for the rate calculation
!===================================================================C
!  Parameters:
!  IN      - File containing location of the minima
!  NRMIN   - Number of Precurssor minima to read in
!  NSMIN   - Number of successor minima to read in
!  RZPI    - Read in zp value of precurssor
!  RZEI    - Read in ze value of precurssor
!  RZPJ    - Read in zp value of successor
!  RZEJ    - Read in ze value of successor
!-------------------------------------------------------------------
!  Input file should have the following format
!  Card 1: Title line
!  Card 2: Comment line (ie presurssor or successor)
!  Card 3 -> n: Read in state zp and ze values
!            (n is nrmin or npmin)
!  Card n+1: Comment line (ie presurssor or successor)
!  Card n+2 -> m: Read in state zp and ze values
!            (m is nrmin or npmin)
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:35 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!===================================================================C
   implicit none

   integer, intent(in) :: in, nrmin, nsmin
   real(8), intent(out), dimension(nrmin) :: rzpi, rzei
   real(8), intent(out), dimension(nsmin) :: rzpj, rzej

   character(80) :: tit, kom
   integer       :: i

   rzpi = 0.d0
   rzei = 0.d0
   rzpj = 0.d0
   rzej = 0.d0

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! read title and comments cards
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   read(in,'(a)') tit

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! read in precurssor states minima
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   read(in,'(a)') kom
   do i=1,nrmin
      read(in,*) rzpi(i),rzei(i)
   enddo

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! read in successor states minima
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   read(in,'(a)') kom
   do i=1,nsmin
      read(in,*) rzpj(i),rzej(i)
   enddo

   return

end subroutine getmin
