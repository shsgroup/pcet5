module elmnts

!=======================================================================
!     Chemical element symbols
!-----------------------------------------------------------------------
!
!  Chemical elements symbols (H - No) including special symbols:
!  XX ( 99) - dummy atom
!  De (103) - general electron donor
!  Ae (104) - general electron acceptor
!  Dp (105) - general proton donor
!  Ap (106) - general proton acceptor
!  Ps (107) - general pseudo-atom
!  XX (108) - dummy atom
!
!-----------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!=======================================================================

   public
   save

   character(4), dimension(108) :: ELSYM =                             &
   (/ '  H ',                                                 '  He',  &
      '  Li', '  Be', '  B ', '  C ', '  N ', '  O ', '  F ', '  Ne',  &
      '  Na', '  Mg', '  Al', '  Si', '  P ', '  S ', '  Cl', '  Ar',  &
      '  K ', '  Ca', '  Sc', '  Ti', '  V ', '  Cr', '  Mn',          &
                      '  Fe', '  Co', '  Ni', '  Cu', '  Zn',          &
                      '  Ga', '  Ge', '  As', '  Se', '  Br', '  Kr',  &
      '  Rb', '  Sr', '  Y ', '  Zr', '  Nb', '  Mo', '  Tc',          &
                      '  Ru', '  Rh', '  Pd', '  Ag', '  Cd',          &
                      '  In', '  Sn', '  Sb', '  Te', '  I ', '  Xe',  &
      '  Cs', '  Ba',                                                  &
              '  La', '  Ce', '  Pr', '  Nd', '  Pm', '  Sm', '  Eu',  &
              '  Gd', '  Tb', '  Dy', '  Ho', '  Er', '  Tm', '  Yb',  &
                      '  Lu', '  Hf', '  Ta', '  W ', '  Re',          &
                      '  Os', '  Ir', '  Pt', '  Au', '  Hg',          &
                      '  Tl', '  Pb', '  Bi', '  Po', '  At', '  Rn',  &
      '  Fr', '  Ra',                                                  &
              '  Ac', '  Th', '  Pa', '  U ', '  Np', '  Pu', '  Am',  &
              '  Cm', '  Bk', '  Cf', '  XX', '  Fm', '  Md', '  No',  &
                      '  De', '  Ae', '  Dp', '  Ap', '  Ps', '    '     /)

end module elmnts
