subroutine initmat

!===================================================================C
!  Initializes the matrices on the grid
!--------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:35 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!===================================================================C

   use pardim
   use control
   use gasmat
   use solmat
   use quantum

   implicit none
   logical :: deriv, derivg

   !-- do we need to allocate the arrays for derivatives?

   deriv  = method.eq.2.and.(igas.eq.1.or.igas.eq.4.or.igas.eq.5).and.isolv.eq.1
   derivg = gquant.and.mgquant.eq.1.and.(igas.eq.4.or.igas.eq.5)

   !-- allocate arrays for matrices on the grid
   !   by calling corresponding module procedures

   call alloc_gasmat(npnts,npntsg,deriv,derivg)
   call alloc_solmat(npnts,npntsg,deriv)
   call alloc_pquant(npnts,deriv)
   call alloc_gquant(npntsg,derivg)

   return

end subroutine initmat

