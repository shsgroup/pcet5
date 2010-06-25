module msevb_water
!---------------------------------------------------------------------
!   MSEVB routines for water clusters (Helene Decornez)
!-------------------------------------------------------------------
!
!  souda
!  2010/06/25 20:02:36
!  4.1
!  Exp
!  module_msevb_water.f90,v 4.1 2010/06/25 20:02:36 souda Exp
!  module_msevb_water.f90,v
!  Revision 4.1  2010/06/25 20:02:36  souda
!  Release 4.1
!
!  Revision 1.1.1.1  2004/01/13 20:03:30  souda
!  Initial PCET-4.0 Release
!
!
!---------------------------------------------------------------------
   use pardim
   use elmnts
   use geogas

   !---------------------------------------------------------------------
   implicit none
   private
   save

   integer, PARAMETER :: NATMAX=31
   integer, PARAMETER :: MAXVBSTATE=15

   ! Parameters for VET with a water chain (Henderson and Cave)
   real*8, parameter :: somea=0.3d0*627.5d0
   real*8, parameter :: somebeta=1.4d0

   ! For TIP3P Water from:
   ! Dang and Pettitt, JPC, 1987, 91,3349-3354

   ! Intermolecular contributions
   ! units for Esq are kcal*A/mol
   ! units for Asq are kcal*A^12/mol
   ! units for Csq are kcal*A^6/mol
   real*8, parameter :: Asq  =  580000.0d0
   real*8, parameter :: Csq  = -525.0d0
   real*8, parameter :: Esqr =  332.1022d0

   ! Intramolecular contributions
   real*8, parameter :: knsOWHW    = 1059.162d0
   real*8, parameter :: knsHWOWHW  = 68.087d0
   real*8, parameter :: distHWOW   = 0.96d0
   real*8, parameter :: anglHWOWHW = 104.5d0*3.14159265359d0/180.d0

   ! Now parameters for water hydronium interactions from the EVB potential
   ! parameters in Voth paper

   real*8, parameter :: msDOH   =  266.3d0
   real*8, parameter :: AOH     =  1.285d0
   real*8, parameter :: R0OH    =  0.98d0
   real*8, parameter :: Kalpha  =  73.27d0
   real*8, parameter :: alpha0  =  116.d0*3.14159265359d0/180.d0
   real*8, parameter :: epsOP   =  0.1535936196592d0
   real*8, parameter :: sgmOP   =  (3.164d0+3.1506d0)/2.d0
   real*8, parameter :: capB    =  2.591d0
   real*8, parameter :: smlb    =  3.50d0
   real*8, parameter :: D0OO    =  2.50d0
   real*8, parameter :: alpha   =  15.0d0
   real*8, parameter :: smlr0OO =  1.90d0
   real*8, parameter :: beta    =  4.50d0
   real*8, parameter :: capR0OO =  3.14d0
   real*8, parameter :: DOO     =  2.875d0
   real*8, parameter :: knsP    =  0.27d0
   real*8, parameter :: knsK    =  11.5d0
   real*8, parameter :: gamma   =  1.85d0
   real*8, parameter :: Vijkns  = -32.925d0

   ! from waterpar.inc: COMMON/water001/NAT
   integer :: NAT

   ! from evbato.inc
   integer     :: NVB,VBMOL(MAXVBSTATE,MAXVBSTATE,NATMAX)               ! COMMON/EVBATO001/
   character*4 :: VBATO(MAXVBSTATE,MAXVBSTATE,NATMAX)                   ! COMMON/EVBATO003/
   character*4 :: INPSYM(NATMAX)
   character*4 :: SYMET(2)
   logical     :: DIAG(MAXVBSTATE,MAXVBSTATE)                           ! COMMON/EVBATO004/
   logical     :: ETEVB,EQUIL                                           ! COMMON/EVBATO008/
   real*8      :: QVB(MAXVBSTATE,MAXVBSTATE,NATMAX),DIST(NATMAX,NATMAX) ! COMMON/EVBATO002/
   real*8      :: XDON,YDON,ZDON,XACP,YACP,ZACP                         ! COMMON/EVBATO006/
   real*8      :: QDON,QACP                                             ! COMMON/EVBATO009/

   real*8,  public, dimension(NATMAX) :: X0EVB, Y0EVB, Z0EVB            ! COMMON/EVBATO000/
   real*8,  public :: DIST_ET=0.d0, DIST_PT=0.d0, DELTAG=0.d0           ! COMMON/EVBATO005/
   integer, public :: Cl, NATET, VEPTMETH=1                             ! COMMON/EVBATO010/

   ! from evbham.inc
   real*8 :: VBHAM(2,2)
   real*8 :: DHDX(MAXVBSTATE,MAXVBSTATE,NATMAX)
   real*8 :: DHDY(MAXVBSTATE,MAXVBSTATE,NATMAX)
   real*8 :: DHDZ(MAXVBSTATE,MAXVBSTATE,NATMAX)
   real*8 :: cvbmat(maxvbstate,maxvbstate)
   real*8 :: evbmat(maxvbstate)

   !---------------------------------------------------------------------
   public :: initevb, msevb


   !---------------------------------------------------------------------
   contains

   !***********************************************************************
   !=======================================================================
   ! INIT_EVB.f
   !=======================================================================
   ! This file is called in the beginning to build all the valence bond
   ! states to be used by the potential and determines all the values
   ! necessary for both the diagonal and the offdiagonal states
   ! The diagonal states will have an H3O+ as a central molecule
   ! and the offdiagonal states will have an h5o2+ as a central
   ! molecule
   !=======================================================================
   subroutine INITEVB

      implicit none

      ! local variables
      CHARACTER SYMB*4
      INTEGER IA,JA,KA,LA,STATE,STT, IC,IAS
      INTEGER I,J,K,VBSTATEI,VBSTATEJ,WATER
      INTEGER TMPIA,TMPJA,TMPKA,TMPLA
      INTEGER NPD,NPA
      REAL*8  MINDIST,TMPDIST
      LOGICAL VB(NATMAX),HYDR(NATMAX),USED(NATMAX)

      !~~~~~~~~~~~ AVS modification ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! initialize number of atoms and cartesian coordinates
      ! from the common blocks of the combined interface:
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      NPD = IPTGAS(1)
      NPA = IPTGAS(3)
      DIST_PT = DABS(XYZGAS(1,NPA)-XYZGAS(1,NPD))

      nat = NATGAS

      DO IA=1,NAT

         x0evb(IA) = XYZGAS(1,IA)
         y0evb(IA) = XYZGAS(2,IA)
         z0evb(IA) = XYZGAS(3,IA)

         IAS = 1
         SYMB = ELSYM(LABGAS(IA))
         DO IC=1,4
            IF (SYMB(IC:IC).EQ.' ') THEN
               IAS = IAS + 1
            ELSE
               EXIT
            ENDIF
         ENDDO

         inpsym(IA) = SYMB(IAS:)

      ENDDO

      !========================================================================
      !  1. FIND PIVOT H3O+ MOLECULE (ONE WITH CLOSEST 3 HYDROGEN ATOMS)
      !  2. DEFINE THIS OXYGEN AS 'OP' OF STATE 1 WITH 3 'HP' HYDROGENS
      !  3. DEFINE CLOSEST OXYGEN TO ONE IN 2 AS 'OP' OF STATE 2, DO ALL STATES
      !     (THERE WILL BE THE SAME NUMBER OF STATES AS THERE ARE OXYGEN ATOMS)
      !  4. ASSIGN ALL OXYGENS IN EACH STATE THAT IS NOT AN 'OP' TO BE AN 'OW'
      !  5. FOR EACH STATE:
      !     START WITH 'OP'
      !=========================================================================

      VBSTATEI=0
      WATER=2
      NVB=0
      TMPDIST=0.

      !=========================================================================
      ! CALCULATE ALL INTERATOMIC DISTANCES
      !=========================================================================
      DO IA=1,NAT
         VB(IA)=.FALSE.
         USED(IA)=.FALSE.
         HYDR(IA)=.FALSE.
         DO JA=1,NAT
            DIST(IA,JA)=0.
            DIST(IA,JA)=DST(IA,JA)
         ENDDO
      ENDDO

      !=========================================================================
      !  1. FIND PIVOT H3O+ MOLECULE
      !     THIS WILL BE THE OXYGEN ATOM WITH THE CLOSEST 3 HYDROGENS
      !=========================================================================
      MINDIST=100000.
      VBSTATEI=1
      NVB=1
      DO IA=1,NAT
         IF (INPSYM(IA).EQ.'O') THEN
            DO JA=1,NAT
              IF (INPSYM(JA).EQ.'H') THEN
                DO KA=1,NAT
                  IF (INPSYM(KA).EQ.'H'.AND.KA.NE.JA) THEN
                    DO LA=1,NAT
                      IF (INPSYM(LA).EQ.'H'.AND.(LA.NE.JA.AND.LA.NE.KA)) THEN
                        TMPDIST=DIST(IA,JA)+DIST(IA,KA)+DIST(IA,LA)
                        IF (TMPDIST.LT.MINDIST) THEN
                          MINDIST=TMPDIST
                          TMPIA=IA
                          TMPJA=JA
                          TMPKA=KA
                          TMPLA=LA
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
         ENDIF
      ENDDO
      USED(TMPIA)=.TRUE.
      VBMOL(VBSTATEI,VBSTATEI,TMPIA)=1
      VBMOL(VBSTATEI,VBSTATEI,TMPJA)=1
      VBMOL(VBSTATEI,VBSTATEI,TMPKA)=1
      VBMOL(VBSTATEI,VBSTATEI,TMPLA)=1
      VBATO(VBSTATEI,VBSTATEI,TMPIA)='OP'
      VBATO(VBSTATEI,VBSTATEI,TMPJA)='HP'
      VBATO(VBSTATEI,VBSTATEI,TMPKA)='HP'
      VBATO(VBSTATEI,VBSTATEI,TMPLA)='HP'
      QVB(VBSTATEI,VBSTATEI,TMPIA)=-0.5d0
      QVB(VBSTATEI,VBSTATEI,TMPJA)=0.5d0
      QVB(VBSTATEI,VBSTATEI,TMPKA)=0.5d0
      QVB(VBSTATEI,VBSTATEI,TMPLA)=0.5d0
      VB(TMPIA)=.TRUE.
      VB(TMPJA)=.TRUE.
      VB(TMPKA)=.TRUE.
      VB(TMPLA)=.TRUE.
      HYDR(TMPJA)=.TRUE.
      HYDR(TMPKA)=.TRUE.
      HYDR(TMPLA)=.TRUE.

      DO IA=1,NAT
         IF (INPSYM(IA).EQ.'O'.AND..NOT.VB(IA)) THEN
            NVB=NVB+1
            VBMOL(VBSTATEI,VBSTATEI,IA)=WATER
            VBATO(VBSTATEI,VBSTATEI,IA)='OW'
            QVB(VBSTATEI,VBSTATEI,IA)=-0.834d0
            VB(IA)=.TRUE.
            MINDIST=10000.
            DO JA=1,NAT
               IF (INPSYM(JA).EQ.'H'.AND..NOT.VB(JA)) THEN
               DO KA=1,NAT
                  IF (INPSYM(KA).EQ.'H'.AND.(.NOT.VB(KA).AND.KA.NE.JA)) THEN
                     TMPDIST=DIST(IA,JA)+DIST(IA,KA)
                     IF (TMPDIST.LT.MINDIST) THEN
                       MINDIST=TMPDIST
                       TMPJA=JA
                       TMPKA=KA
                     ENDIF
                   ENDIF
                 ENDDO
               ENDIF
             ENDDO
             VB(TMPJA)=.TRUE.
             VB(TMPKA)=.TRUE.
             VBMOL(VBSTATEI,VBSTATEI,TMPJA)=WATER
             VBMOL(VBSTATEI,VBSTATEI,TMPKA)=WATER
             VBATO(VBSTATEI,VBSTATEI,TMPJA)='HW'
             VBATO(VBSTATEI,VBSTATEI,TMPKA)='HW'
             QVB(VBSTATEI,VBSTATEI,TMPJA)=0.417d0
             QVB(VBSTATEI,VBSTATEI,TMPKA)=0.417d0
             WATER=WATER+1
         ENDIF
      ENDDO

      !===========================================================================
      ! NOW NEED TO BUILD THE OTHER STATES
      !===========================================================================
      DO VBSTATEI=2,NVB
         DO IA=1,NAT
            VB(IA)=.FALSE.
         ENDDO
         MINDIST=10000.
         DO IA=1,NAT
            IF (INPSYM(IA).EQ.'O'.AND..NOT.USED(IA)) THEN
              DO JA=1,NAT
                IF (HYDR(JA).AND..NOT.USED(JA)) THEN
                  IF (DIST(IA,JA).LT.MINDIST) THEN
                    MINDIST=DIST(IA,JA)
                    TMPIA=IA
                    TMPJA=JA
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
         ENDDO
         USED(TMPJA)=.TRUE.
         USED(TMPIA)=.TRUE.
         MINDIST=10000.
         DO KA=1,NAT
           IF (INPSYM(KA).EQ.'H'.AND.KA.NE.TMPJA)THEN
              DO LA=1,NAT
                IF(INPSYM(LA).EQ.'H'.AND.(LA.NE.TMPJA.AND.LA.NE.KA))THEN
                  TMPDIST=DIST(TMPIA,KA)+DIST(TMPIA,LA)
                  IF (TMPDIST.LT.MINDIST) THEN
                    MINDIST=TMPDIST
                    TMPKA=KA
                    TMPLA=LA
                  ENDIF
                ENDIF
              ENDDO
           ENDIF
         ENDDO
         HYDR(TMPKA)=.TRUE.
         HYDR(TMPLA)=.TRUE.
         VB(TMPIA)=.TRUE.
         VB(TMPJA)=.TRUE.
         VB(TMPKA)=.TRUE.
         VB(TMPLA)=.TRUE.
         VBMOL(VBSTATEI,VBSTATEI,TMPIA)=1
         VBMOL(VBSTATEI,VBSTATEI,TMPJA)=1
         VBMOL(VBSTATEI,VBSTATEI,TMPKA)=1
         VBMOL(VBSTATEI,VBSTATEI,TMPLA)=1
         VBATO(VBSTATEI,VBSTATEI,TMPIA)='OP'
         VBATO(VBSTATEI,VBSTATEI,TMPJA)='HP'
         VBATO(VBSTATEI,VBSTATEI,TMPKA)='HP'
         VBATO(VBSTATEI,VBSTATEI,TMPLA)='HP'
         QVB(VBSTATEI,VBSTATEI,TMPIA)=-0.5d0
         QVB(VBSTATEI,VBSTATEI,TMPJA)=0.5d0
         QVB(VBSTATEI,VBSTATEI,TMPKA)=0.5d0
         QVB(VBSTATEI,VBSTATEI,TMPLA)=0.5d0
         WATER=2
         DO IA=1,NAT
            IF (INPSYM(IA).EQ.'O'.AND..NOT.VB(IA)) THEN
              VBMOL(VBSTATEI,VBSTATEI,IA)=WATER
              VBATO(VBSTATEI,VBSTATEI,IA)='OW'
              QVB(VBSTATEI,VBSTATEI,IA)=-0.834d0
              VB(IA)=.TRUE.
              MINDIST=10000.d0
              DO JA=1,NAT
                IF (INPSYM(JA).EQ.'H'.AND..NOT.VB(JA)) THEN
                  DO KA=1,NAT
                    IF(INPSYM(KA).EQ.'H'.AND.(.NOT.VB(KA).AND.KA.NE.JA)) THEN
                      TMPDIST=DIST(IA,JA)+DIST(IA,KA)
                      IF (TMPDIST.LT.MINDIST) THEN
                        MINDIST=TMPDIST
                        TMPJA=JA
                        TMPKA=KA
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
              VB(TMPJA)=.TRUE.
              VB(TMPKA)=.TRUE.
              VBMOL(VBSTATEI,VBSTATEI,TMPJA)=WATER
              VBMOL(VBSTATEI,VBSTATEI,TMPKA)=WATER
              VBATO(VBSTATEI,VBSTATEI,TMPJA)='HW'
              VBATO(VBSTATEI,VBSTATEI,TMPKA)='HW'
              QVB(VBSTATEI,VBSTATEI,TMPJA)=0.417
              QVB(VBSTATEI,VBSTATEI,TMPKA)=0.417
              WATER=WATER+1
            ENDIF
         ENDDO
      ENDDO

      !==================================================================================
      ! NOW NEED TO BUILD THE OFF DIAGONAL TERMS
      !==================================================================================
       DO VBSTATEI=1,NVB
         DO VBSTATEJ=1,NVB
           DIAG(VBSTATEI,VBSTATEJ)=.FALSE.
           DO IA=1,NAT
             VB(IA)=.FALSE.
             IF ((VBATO(VBSTATEI,VBSTATEI,IA).EQ.&
                  VBATO(VBSTATEJ,VBSTATEJ,IA))&
                  .AND.VBATO(VBSTATEI,VBSTATEI,IA).EQ.'HP'&
                  .AND.(VBSTATEI.NE.VBSTATEJ)) THEN
               VBMOL(VBSTATEI,VBSTATEJ,IA)=0
               VBATO(VBSTATEI,VBSTATEJ,IA)='HZ'
               QVB(VBSTATEI,VBSTATEJ,IA)=0.04002d0
               VB(IA)=.TRUE.
               DIAG(VBSTATEI,VBSTATEJ)=.TRUE.
             ENDIF
           ENDDO
           IF (DIAG(VBSTATEI,VBSTATEJ)) THEN
             DO IA=1,NAT
               IF ((VBATO(VBSTATEI,VBSTATEI,IA).EQ.'HP'.OR.  &
                   VBATO(VBSTATEJ,VBSTATEJ,IA).EQ.'HP').AND. &
                   .NOT.VBATO(VBSTATEI,VBSTATEJ,IA).EQ.'HZ') THEN
                VBMOL(VBSTATEI,VBSTATEJ,IA)=0
                 VBATO(VBSTATEI,VBSTATEJ,IA)='HM'
                 QVB(VBSTATEI,VBSTATEJ,IA)=0.04002d0
                 VB(IA)=.TRUE.
               ELSEIF(VBATO(VBSTATEI,VBSTATEI,IA).EQ.'OP'.OR. &
                   VBATO(VBSTATEJ,VBSTATEJ,IA).EQ.'OP') THEN
                 VBMOL(VBSTATEI,VBSTATEJ,IA)=0
                 VBATO(VBSTATEI,VBSTATEJ,IA)='OM'
                 QVB(VBSTATEI,VBSTATEJ,IA)=-0.10005d0
                 VB(IA)=.TRUE.
               ENDIF
             ENDDO
             WATER=3
             DO IA=1,NAT
               IF (VBATO(VBSTATEI,VBSTATEI,IA).EQ.'OW'   &
                  .AND..NOT.VB(IA)) THEN
                 VBMOL(VBSTATEI,VBSTATEJ,IA)=WATER
                 VBATO(VBSTATEI,VBSTATEJ,IA)='OW'
                 QVB(VBSTATEI,VBSTATEJ,IA)=-0.834d0
                 VB(IA)=.TRUE.
                 DO JA=1,NAT
                   IF (VBATO(VBSTATEI,VBSTATEI,JA).EQ.'HW'   &
                     .AND.VBMOL(VBSTATEI,VBSTATEI,JA).EQ.    &
                     VBMOL(VBSTATEI,VBSTATEI,IA)) THEN
                     VBMOL(VBSTATEI,VBSTATEJ,JA)=WATER
                     VBATO(VBSTATEI,VBSTATEJ,JA)='HW'
                     QVB(VBSTATEI,VBSTATEJ,JA)=0.417d0
                     VB(JA)=.TRUE.
                   ENDIF
                 ENDDO
                 WATER=WATER+1
               ENDIF
             ENDDO
           ENDIF
         ENDDO
       ENDDO

      !debug
      ! do state=1,nvb
      !    do stt=1,nvb
      !       write(*,*) 'state',state,'stt',stt,'diag',diag(state,stt)
      !       if (diag(state,stt)) then
      !          do ia=1,nat
      !             write(*,*) state,stt,ia,vbato(state,stt,ia), &
      !                   qvb(state,stt,ia),vbmol(state,stt,ia)
      !          enddo
      !       endif
      !    enddo
      ! enddo

      RETURN

   end subroutine INITEVB

   !***********************************************************************

   subroutine msevb(hevb)

      ! This is the main subroutine to calculate the potential
      ! using the Schmitt-Voth MSEVB potential for an extra
      ! proton in water

      implicit none
      
      real*8, intent(out) :: hevb(4,4)

      ! local variables
      integer, parameter :: two=4
      integer :: statei,statej,i,j
      integer :: ierr,ia,ato
      real*8  :: etmp
      REAL*8  :: r,RIJ,RIJSQ,VIJ,FIJ,DXIJ,DYIJ,DZIJ
      REAL*8  :: XI,YI,ZI,FXI,FYI,FZI
      REAL*8  :: DYI,DZI,RI,RISQ,FI,VI
      REAL*8  :: DSQKL,DXKL,DYKL,DZKL,FCXKL,FCYKL,FCZKL
      real*8  :: work1(4),work2(4)
      real*8  :: EEVB(4),CEVB(4,4)
      real*8  :: DCdx(natmax),DCdy(natmax),DCdz(natmax)
      real*8  :: VEPT,VET
      logical, save :: startms=.true.

      ! For debugging
      integer :: nitera
      integer :: DBGIA,DBGCOOR
      real*8  :: distAC,distAB
      real*8  :: distCD,distCE
      real*8  :: interAB,interCD,FracK1,FracK2
      real*8  :: colQ,tmpAB,tmpCD
      real*8  :: colQ1,colQ2
      real*8  :: QEA,QED
      real*8, save :: XED,YED,ZED,XEA,YEA,ZEA

      !PARAMETER(SOMEA=0.3*627.5)
      !PARAMETER(SOMEBETA=1.4)
      !if (startms) then
      !   write(*,*) 'enter delta G value to be used'
      !   read(*,*) deltag
      !   write(*,*) 'enter elctron transfer distance'
      !   read (*,*)distDA
      !   XED=-dist_ET/2.
      !   XEA=dist_ET/2.
      !   write(*,*) 'XED',XED,'XEA',XEA
      !   write (*,*) 'enter method for calculating VEPT'
      !   write(*,*) '(0) sets VEPT to the value for VET'
      !   write(*,*) '(1) sets VEPT to the value for VET*0.1'
      !   write(*,*) '(2) sets VEPT to the value for VET*10.0'
      !   read(*,*) VEPTMeth
      !endif

      XED = -dist_ET/2.d0
      XEA =  dist_ET/2.d0
      YED =  0.d0
      ZED =  0.d0
      YEA =  0.d0
      ZEA =  0.d0
      QEA = -1.d0
      QED = -1.d0

      ! Now build hamiltonian and first derivative of the hamiltonian

      do statei=1,nvb
        do statej=1,nvb
          call hijevb(statei,statej)
        enddo
        call hiievb(statei)
      enddo

      ! VET according to Cave paper

      VET=SOMEA*dexp(-SOMEBETA*dist_ET/2.d0)

      ! Now build states with electron interactions

      CALL COULOMBIC(ETMP,X0evb,Y0evb,Z0evb,DCdx,DCdy,DCdz,XED,YED,ZED,QED,1)
      HEVB(1,1)=VBHAM(1,1)+ETMP
      CALL COULOMBIC(ETMP,X0evb,Y0evb,Z0evb,DCdx,DCdy,DCdz,XED,YED,ZED,QED,2)
      HEVB(2,2)=VBHAM(2,2)+ETMP
      HEVB(1,2)=VBHAM(1,2)
      HEVB(2,1)=VBHAM(2,1)
      HEVB(1,3)=VET
      HEVB(2,4)=VET

      if (VEPTMeth.eq.0)then
        HEVB(1,4)=VET
        HEVB(2,3)=VET
      elseif(VEPTMeth.eq.1) then
        HEVB(1,4)=VET*0.1
        HEVB(2,3)=VET*0.1
      elseif(VEPTMeth.eq.2) then
        HEVB(1,4)=VET*10.
        HEVB(2,3)=VET*10.
      endif

      CALL COULOMBIC(ETMP,X0evb,Y0evb,Z0evb,DCdx,DCdy,DCdz,XEA,YEA,ZEA,QEA,1)
      HEVB(3,3)=VBHAM(1,1)+ETMP+deltag
      CALL COULOMBIC(ETMP,X0evb,Y0evb,Z0evb,DCdx,DCdy,DCdz,XEA,YEA,ZEA,QEA,2)
      HEVB(4,4)=VBHAM(2,2)+ETMP+deltag
      HEVB(3,4)=VBHAM(1,2)
      HEVB(4,3)=VBHAM(2,1)
      HEVB(3,1)=VET
      HEVB(4,2)=VET
      if (VEPTMeth.eq.0)then
        HEVB(4,1)=VET
        HEVB(3,2)=VET
      elseif(VEPTMeth.eq.1) then
        HEVB(4,1)=VET*0.1
        HEVB(3,2)=VET*0.1
      elseif(VEPTMeth.eq.2) then
        HEVB(4,1)=VET*10.
        HEVB(3,2)=VET*10.
      endif
      if (startms) then
         !write (*,*) 'sample gas phase evb matrix for first rp point'
         !do i=1,2*nvb
         !   write(*,*) HEVB(i,1),HEVB(i,2),HEVB(i,3),HEVB(i,4)
         !enddo
         startms=.false.
      endif

      RETURN
      
   end subroutine msevb

   !***********************************************************************

   !C------------------------------------------------------------------------------
   !C     Modification History of file: hii_evb.f
   !C------------------------------------------------------------------------------
   !C     When | Who |      What
   !C----------|-----|-------------------------------------------------------------
   !C  27-05-99| hd  |- final creation
   !C------------------------------------------------------------------------------
   !C This is a new subroutine to go with the new build evb steps.
   !
   !C This subroutine calculates the diagonal terms of the evb hamiltonian
   !C given the atom types of the evb states and the positions of the
   !C atoms
   !C The parameters are those developped by Voth except for water where
   !C the parameters coem from the TIP3P potential by Jorgensen

   subroutine hiievb(statei)

      implicit none
      integer, intent(in) :: statei

      ! local variables
      integer ia,ja,ka
      real*8 hiitot,hii1,hii2
      real*8 hii1a,hii1b,hii2a,hii2b
      real*8 hii3,hii4
      real*8 hii3a,hii3b,hii3c,hii4a,hii4b
      real*8 dhii1adx(natmax),dhii1bdx(natmax)
      real*8 dhii1ady(natmax),dhii1bdy(natmax)
      real*8 dhii1adz(natmax),dhii1bdz(natmax)
      real*8 dhii2adx(natmax),dhii2bdx(natmax)
      real*8 dhii2ady(natmax),dhii2bdy(natmax)
      real*8 dhii2adz(natmax),dhii2bdz(natmax)
      real*8 dhii3adx(natmax),dhii3bdx(natmax)
      real*8 dhii3ady(natmax),dhii3bdy(natmax)
      real*8 dhii3adz(natmax),dhii3bdz(natmax)
      real*8 dhii3cdx(natmax),dhii3cdy(natmax)
      real*8 dhii3cdz(natmax)
      real*8 dhii4adx(natmax),dhii4bdx(natmax)
      real*8 dhii4ady(natmax),dhii4bdy(natmax)
      real*8 dhii4adz(natmax),dhii4bdz(natmax)
      real*8 dhii3dx(natmax),dhii4dx(natmax)
      real*8 dhii1dx(natmax),dhii2dx(natmax)
      real*8 dhiidy(natmax),dhii1dy(natmax),dhii2dy(natmax)
      real*8 dhii3dy(natmax),dhii4dy(natmax)
      real*8 dhiidz(natmax),dhii1dz(natmax),dhii2dz(natmax)
      real*8 dhii3dz(natmax),dhii4dz(natmax)
      real*8 tmpfac,tmpfac1,tmpfac2,tmpfac3,tmpfac4
      logical debug,usedang(natmax,natmax)

      !external dst,angle,queue
      !external drdx1,drdx2,drdy1,drdy2,drdz1,drdz2
      !external dwdx1,dwdx2,dwdy1,dwdy2,dwdz1,dwdz2
      !external dqdx1,dqdx2,dqdy1,dqdy2,dqdz1,dqdz2
      !external dwdx3,dwdy3,dwdz3
      !external dqdx3,dqdy3,dqdz3

      debug=.FALSE.

      ! The best way to do this is to break up the terms contributing
      ! to the potential
      ! hii1 is the intramolecular H3O+ term
      ! hii2 is the intramolecular term over all the water molecules
      ! hii3 is the intermolecular term between the one H3O+ and the
      !      remaining waters in the system
      ! hii4 is the intermolecular term between all the waters
      ! all hiitot is the total diagonal hamiltonian term
      ! all hii1-hii4 and hiitot depend on the evb state occupied

      ! Part A: calculate hii1 for all the evb states

      ! calculate the three intramolecular angles for H3O+: HPOPHP
      hii1   = 0.d0
      hii1a  = 0.d0
      hii1b  = 0.d0
      hii2   = 0.d0
      hii2a  = 0.d0
      hii2b  = 0.d0
      hii3   = 0.d0
      hii3a  = 0.d0
      hii3b  = 0.d0
      hii3c  = 0.d0
      hii4   = 0.d0
      hii4a  = 0.d0
      hii4b  = 0.d0
      hiitot = 0.d0
      do ia=1,nat
        dhii1dx(ia)  = 0.d0
        dhii2dx(ia)  = 0.d0
        dhii3dx(ia)  = 0.d0
        dhii4dx(ia)  = 0.d0
        dhii1dy(ia)  = 0.d0
        dhii2dy(ia)  = 0.d0
        dhii3dy(ia)  = 0.d0
        dhii4dy(ia)  = 0.d0
        dhii1dz(ia)  = 0.d0
        dhii2dz(ia)  = 0.d0
        dhii3dz(ia)  = 0.d0
        dhii4dz(ia)  = 0.d0
        dhii1adx(ia) = 0.d0
        dhii2adx(ia) = 0.d0
        dhii3adx(ia) = 0.d0
        dhii4adx(ia) = 0.d0
        dhii1ady(ia) = 0.d0
        dhii2ady(ia) = 0.d0
        dhii3ady(ia) = 0.d0
        dhii4ady(ia) = 0.d0
        dhii1adz(ia) = 0.d0
        dhii2adz(ia) = 0.d0
        dhii3adz(ia) = 0.d0
        dhii4adz(ia) = 0.d0
        dhii1bdx(ia) = 0.d0
        dhii2bdx(ia) = 0.d0
        dhii3bdx(ia) = 0.d0
        dhii4bdx(ia) = 0.d0
        dhii1bdy(ia) = 0.d0
        dhii2bdy(ia) = 0.d0
        dhii3bdy(ia) = 0.d0
        dhii4bdy(ia) = 0.d0
        dhii1bdz(ia) = 0.d0
        dhii2bdz(ia) = 0.d0
        dhii3bdz(ia) = 0.d0
        dhii4bdz(ia) = 0.d0
        dhii3cdx(ia) = 0.d0
        dhii3cdy(ia) = 0.d0
        dhii3cdz(ia) = 0.d0
      enddo
      do ia=1,nat
        do ja=1,nat
          usedang(ia,ja)=.FALSE.
        enddo
      enddo


      ! This next section calculates the hii1 contribution

      do ia=1,nat
        if (vbato(statei,statei,ia).eq.'OP') then
          do ja=1,nat
            if (vbato(statei,statei,ja).eq.'HP') then
              hii1a=hii1a+msDOH*(1.-exp(-AOH*(dst(ia,ja)- R0OH)))**2
              tmpfac=(2.*AOH*msDOH*(1. - exp(-(aoh*(dst(ia,ja)-        &
                 R0OH)))))/exp(AOH*(dst(ia,ja)- R0OH))
              dhii1adx(ia)=dhii1adx(ia)+tmpfac*drdx1(ia,ja)
              dhii1ady(ia)=dhii1ady(ia)+tmpfac*drdy1(ia,ja)
              dhii1adz(ia)=dhii1adz(ia)+tmpfac*drdz1(ia,ja)
              dhii1adx(ja)=dhii1adx(ja)+tmpfac*drdx2(ia,ja)
              dhii1ady(ja)=dhii1ady(ja)+tmpfac*drdy2(ia,ja)
              dhii1adz(ja)=dhii1adz(ja)+tmpfac*drdz2(ia,ja)
              do ka=1,nat
                if (vbato(statei,statei,ka).eq.'HP'.and.ka.ne.ja &
                  .and.(.not.usedang(ka,ja))) then
                  usedang(ka,ja)=.true.
                  usedang(ja,ka)=.true.
                  hii1b=hii1b+0.5*kalpha*((angle(ja,ka,ia)-alpha0)**2)
                  tmpfac=Kalpha*(angle(ja,ka,ia)-alpha0)
                  dhii1bdx(ia)=dhii1bdx(ia)+tmpfac*dwdx3(ja,ka,ia)
                  dhii1bdy(ia)=dhii1bdy(ia)+tmpfac*dwdy3(ja,ka,ia)
                  dhii1bdz(ia)=dhii1bdz(ia)+tmpfac*dwdz3(ja,ka,ia)
                  dhii1bdx(ja)=dhii1bdx(ja)+tmpfac*dwdx1(ja,ka,ia)
                  dhii1bdy(ja)=dhii1bdy(ja)+tmpfac*dwdy1(ja,ka,ia)
                  dhii1bdz(ja)=dhii1bdz(ja)+tmpfac*dwdz1(ja,ka,ia)
                  dhii1bdx(ka)=dhii1bdx(ka)+tmpfac*dwdx2(ja,ka,ia)
                  dhii1bdy(ka)=dhii1bdy(ka)+tmpfac*dwdy2(ja,ka,ia)
                  dhii1bdz(ka)=dhii1bdz(ka)+tmpfac*dwdz2(ja,ka,ia)
                endif
              enddo
            endif
          enddo
        endif
      enddo
      hii1=hii1a+hii1b
      do ia=1,nat
        dhii1dx(ia)=dhii1adx(ia)+dhii1bdx(ia)
        dhii1dy(ia)=dhii1ady(ia)+dhii1bdy(ia)
        dhii1dz(ia)=dhii1adz(ia)+dhii1bdz(ia)
      enddo

      ! This next section calculates the hii2 contribution
      ! This is the intramolecular contribution from the water molecules
      ! identified in the evb state

      do ia=1,nat
        if (vbato(statei,statei,ia).eq.'OW') then
          do ja=1,nat
            if (vbato(statei,statei,ja).eq.'HW'.and.   &
                vbmol(statei,statei,ia).eq.            &
                vbmol(statei,statei,ja)) then
              do ka=ja+1,nat
                if (vbato(statei,statei,ka).eq.'HW'.and. &
                  vbmol(statei,statei,ia)                &
                  .eq.vbmol(statei,statei,ka)) then
                  tmpfac1=knsOWHW*(dst(ia,ja)-distHWOW)
                  tmpfac2=knsOWHW*(dst(ia,ka)-distHWOW)
                  tmpfac3=knsOWHW*(dst(ja,ka)-distHWOW)
                  hii2a=hii2a+ 0.5*knsOWHW*(dst(ia,ja)-distHWOW)**2  &
                        +0.5*knsOWHW*(dst(ia,ka)-distHWOW)**2
                  dhii2adx(ia)=dhii2adx(ia)+tmpfac1*drdx1(ia,ja)     &
                               +tmpfac2*drdx1(ia,ka)
                  dhii2ady(ia)=dhii2ady(ia)+tmpfac1*drdy1(ia,ja)  &
                               +tmpfac2*drdy1(ia,ka)
                  dhii2adz(ia)=dhii2adz(ia)+tmpfac1*drdz1(ia,ja)  &
                               +tmpfac2*drdz1(ia,ka)
                  dhii2adx(ja)=dhii2adx(ja)+tmpfac1*drdx2(ia,ja)
                  dhii2ady(ja)=dhii2ady(ja)+tmpfac1*drdy2(ia,ja)
                  dhii2adz(ja)=dhii2adz(ja)+tmpfac1*drdz2(ia,ja)
                  dhii2adx(ka)=dhii2adx(ka)+tmpfac2*drdx2(ia,ka)
                  dhii2ady(ka)=dhii2ady(ka)+tmpfac2*drdy2(ia,ka)
                  dhii2adz(ka)=dhii2adz(ka)+tmpfac2*drdz2(ia,ka)
                  hii2b=hii2b+0.5*knsHWOWHW*(angle(ja,ka,ia)-     &
                        anglHWOWHW)**2
                  tmpfac=knsHWOWHW*(angle(ja,ka,ia)-anglHWOWHW)
                  dhii2bdx(ia)=dhii2bdx(ia)+tmpfac*dwdx3(ja,ka,ia)
                  dhii2bdy(ia)=dhii2bdy(ia)+tmpfac*dwdy3(ja,ka,ia)
                  dhii2bdz(ia)=dhii2bdz(ia)+tmpfac*dwdz3(ja,ka,ia)
                  dhii2bdx(ja)=dhii2bdx(ja)+tmpfac*dwdx1(ja,ka,ia)
                  dhii2bdy(ja)=dhii2bdy(ja)+tmpfac*dwdy1(ja,ka,ia)
                  dhii2bdz(ja)=dhii2bdz(ja)+tmpfac*dwdz1(ja,ka,ia)
                  dhii2bdx(ka)=dhii2bdx(ka)+tmpfac*dwdx2(ja,ka,ia)
                  dhii2bdy(ka)=dhii2bdy(ka)+tmpfac*dwdy2(ja,ka,ia)
                  dhii2bdz(ka)=dhii2bdz(ka)+tmpfac*dwdz2(ja,ka,ia)
                endif
              enddo
            endif
          enddo
        endif
      enddo
      hii2=hii2a+hii2b
      do ia=1,nat
        dhii2dx(ia)=dhii2adx(ia)+dhii2bdx(ia)
        dhii2dy(ia)=dhii2ady(ia)+dhii2bdy(ia)
        dhii2dz(ia)=dhii2adz(ia)+dhii2bdz(ia)
      enddo

      ! This next section calculates the hii3 contribution
      ! This contribution is the intermolecular potential between
      ! the H3O+ with each water in the system

      do ia=1,nat
        if (vbato(statei,statei,ia).eq.'OP') then
          do ja=1,nat
            if (vbato(statei,statei,ja).eq.'OW') then
              hii3a=hii3a+4.*epsOP*((sgmOP**12./dst(ia,ja)**12.) &
              - (sgmOP**6./dst(ia,ja)**6.))
              tmpfac1=4.*epsOP*(-12.*((sgmOP/dst(ia,ja))**12 &
                  +6.*(sgmOP/dst(ia,ja))**6))/dst(ia,ja)
              dhii3adx(ia)=dhii3adx(ia)+tmpfac1*drdx1(ia,ja)
              dhii3ady(ia)=dhii3ady(ia)+tmpfac1*drdy1(ia,ja)
              dhii3adz(ia)=dhii3adz(ia)+tmpfac1*drdz1(ia,ja)
              dhii3adx(ja)=dhii3adx(ja)+tmpfac1*drdx2(ia,ja)
              dhii3ady(ja)=dhii3ady(ja)+tmpfac1*drdy2(ia,ja)
              dhii3adz(ja)=dhii3adz(ja)+tmpfac1*drdz2(ia,ja)

              hii3c=hii3c+capB*(1.-tanh(smlb*(dst(ia,ja)-d0OO)))
              tmpfac3=-(capB*smlb*(1./(cosh((-D0OO + &
                      dst(ia,ja))*smlb)**2)))
              dhii3cdx(ia)=dhii3cdx(ia)+tmpfac3*drdx1(ia,ja)
              dhii3cdy(ia)=dhii3cdy(ia)+tmpfac3*drdy1(ia,ja)
              dhii3cdz(ia)=dhii3cdz(ia)+tmpfac3*drdz1(ia,ja)
              dhii3cdx(ja)=dhii3cdx(ja)+tmpfac3*drdx2(ia,ja)
              dhii3cdy(ja)=dhii3cdy(ja)+tmpfac3*drdy2(ia,ja)
              dhii3cdz(ja)=dhii3cdz(ja)+tmpfac3*drdz2(ia,ja)
            endif
          enddo
        endif
      enddo
      do ia=1,nat
        if(vbmol(statei,statei,ia).eq.1) then
          do ja=1,nat
            if(vbmol(statei,statei,ja).ne.1) then
              hii3b=hii3b+ &
               QVB(statei,statei,ia)* &
               QVB(statei,statei,ja)/dst(ia,ja)
              tmpfac4=-(((QVB(statei,statei,ia))*0.76* &
               (QVB(statei,statei,ja))*Esqr)/(dst(ia,ja)**2))
              dhii3bdx(ia)=dhii3bdx(ia)+tmpfac4*drdx1(ia,ja)
              dhii3bdy(ia)=dhii3bdy(ia)+tmpfac4*drdy1(ia,ja)
              dhii3bdz(ia)=dhii3bdz(ia)+tmpfac4*drdz1(ia,ja)
              dhii3bdx(ja)=dhii3bdx(ja)+tmpfac4*drdx2(ia,ja)
              dhii3bdy(ja)=dhii3bdy(ja)+tmpfac4*drdy2(ia,ja)
              dhii3bdz(ja)=dhii3bdz(ja)+tmpfac4*drdz2(ia,ja)
            endif
          enddo
        endif
      enddo
      hii3b=hii3b*Esqr*0.76
      hii3=hii3a+hii3b+hii3c
      do ia=1,nat
        dhii3dx(ia)=dhii3adx(ia)+dhii3bdx(ia)+dhii3cdx(ia)
        dhii3dy(ia)=dhii3ady(ia)+dhii3bdy(ia)+dhii3cdy(ia)
        dhii3dz(ia)=dhii3adz(ia)+dhii3bdz(ia)+dhii3cdz(ia)
      enddo

      ! This next section calculates the contribution from the hii4
      ! This term accounts for intermolecular interactions between the water
      ! molecules only

      do ia=1,nat
        if (vbato(statei,statei,ia).eq.'OW') then
          do ja=ia+1,nat
            if (vbato(statei,statei,ja).eq.'OW') then
              hii4a=hii4a+Asq/(dst(ia,ja)**12)+Csq/(dst(ia,ja)**6)
              tmpfac1=-12.*Asq/(dst(ia,ja)**13)
              tmpfac2=-6.*Csq/(dst(ia,ja)**7)
              dhii4adx(ia)=dhii4adx(ia)+(tmpfac1+tmpfac2)*drdx1(ia,ja)
              dhii4ady(ia)=dhii4ady(ia)+(tmpfac1+tmpfac2)*drdy1(ia,ja)
              dhii4adz(ia)=dhii4adz(ia)+(tmpfac1+tmpfac2)*drdz1(ia,ja)
              dhii4adx(ja)=dhii4adx(ja)+(tmpfac1+tmpfac2)*drdx2(ia,ja)
              dhii4ady(ja)=dhii4ady(ja)+(tmpfac1+tmpfac2)*drdy2(ia,ja)
              dhii4adz(ja)=dhii4adz(ja)+(tmpfac1+tmpfac2)*drdz2(ia,ja)
            endif
          enddo
        endif
      enddo
      do ia=1,nat
        if (vbmol(statei,statei,ia).ne.1) then
          do ja=ia+1,nat
            if (vbmol(statei,statei,ja).ne.1.and. &
                vbmol(statei,statei,ja).ne.       &
                vbmol(statei,statei,ia)) then
              hii4b=hii4b+((QVB(statei,statei,ia)* &
              QVB(statei,statei,ja)*Esqr)/(dst(ia,ja)))
              tmpfac4=-(((QVB(statei,statei,ia))* &
               (QVB(statei,statei,ja))*Esqr)/(dst(ia,ja)**2))
              dhii4bdx(ia)=dhii4bdx(ia)+tmpfac4*drdx1(ia,ja)
              dhii4bdy(ia)=dhii4bdy(ia)+tmpfac4*drdy1(ia,ja)
              dhii4bdz(ia)=dhii4bdz(ia)+tmpfac4*drdz1(ia,ja)
              dhii4bdx(ja)=dhii4bdx(ja)+tmpfac4*drdx2(ia,ja)
              dhii4bdy(ja)=dhii4bdy(ja)+tmpfac4*drdy2(ia,ja)
              dhii4bdz(ja)=dhii4bdz(ja)+tmpfac4*drdz2(ia,ja)
           endif
          enddo
        endif
      enddo
      hii4=hii4a+hii4b
      do ia=1,nat
        dhii4dx(ia)=dhii4adx(ia)+dhii4bdx(ia)
        dhii4dy(ia)=dhii4ady(ia)+dhii4bdy(ia)
        dhii4dz(ia)=dhii4adz(ia)+dhii4bdz(ia)
      enddo

      ! Now add them up

      vbham(statei,statei)=hii1+hii2+hii3+hii4
      do ia=1,nat
         dhdx(statei,statei,ia)=dhii1dx(ia)+dhii2dx(ia)+dhii3dx(ia)+dhii4dx(ia)
         dhdy(statei,statei,ia)=dhii1dy(ia)+dhii2dy(ia)+dhii3dy(ia)+dhii4dy(ia)
         dhdz(statei,statei,ia)=dhii1dz(ia)+dhii2dz(ia)+dhii3dz(ia)+dhii4dz(ia)
      enddo

!===============some debugging========================================
!      if (statei.eq.1) then
!      write (*,*) 'state',statei
!      write(*,*) 'hii3a',hii3a,'hii3b',hii3b,'hii3c',hii3c
!      tmpfac1=4.*0.155*((3.164/dst(2,5))**12-(3.164/dst(2,5))**6)
!      write(*,*) 'hii3anum',tmpfac1
!      tmpfac2=332.177512*.76*((0.417*0.4/dst(3,7))+
!     ;       (0.417*0.5/dst(3,6))+(0.5*-0.834/dst(3,5))+
!     ;       (0.417*0.5/dst(1,7))+
!     ;       (0.417*0.5/dst(1,6))+(0.5*-0.834/dst(1,5))+
!     ;       (0.417*0.5/dst(4,7))+
!     ;       (0.417*0.5/dst(4,6))+(0.5*-0.834/dst(4,5))+
!     ;       (-0.5*0.417/dst(2,7))+
!     ;       (-0.5*0.417/dst(2,6))+(-0.5*-0.834/dst(2,5)))
!      write(*,*) 'hii3bnum',tmpfac2
!      tmpfac3=2.591*(1.-tanh(3.50*(dst(2,5)-2.50)))
!      write(*,*) 'hii3cnum',tmpfac3
!
!===============
!      elseif (statei.eq.2) then
!      write (*,*) 'state',statei
!      write(*,*) 'hii3a',hii3a,'hii3b',hii3b,'hii3c',hii3c
!      tmpfac1=4.*0.155*((3.164/dst(2,5))**12-(3.164/dst(2,5))**6)
!      write(*,*) 'hii3anum',tmpfac1
!      tmpfac2=332.17752*0.76*((0.4*0.5/dst(4,1))+
!     ;       (0.4*0.5/dst(4,3))+(0.5*(-0.8)/dst(4,2))+
!     ;       (0.4*0.5/dst(7,1))+
!     ;       (0.4*0.5/dst(7,3))+(0.5*(-0.8)/dst(7,2))+
!     ;       (0.4*0.5/dst(6,1))+
!     ;       (0.4*0.5/dst(6,3))+(0.5*(-0.8)/dst(6,2))+
!     ;       (-0.5*0.4/dst(5,1))+
!     ;       (-0.5*0.4/dst(5,3))+(-0.5*(-0.8)/dst(5,2)))
!      write(*,*) 'hii3bnum',tmpfac2
!      tmpfac3=2.591*(1.-tanh(3.50*(dst(2,5)-2.50)))
!      write(*,*) 'hii3cnum',tmpfac3
!
!      write(*,*) 'hii4',hii4
!
!
!      endif
!       if (statei.eq.1) then
!       write (*,*) 'diag statei',statei
!       write(*,*)'hii1',hii1
!       write(*,*)'hii1a',hii1a,'hii1b',hii1b
!       write(*,*)'hii2',hii2
!       write(*,*)'hii2a',hii2a,'hii2b',hii2b
!       write(*,*) 'hii1+hii2',hii1+hii2
!       write(*,*)'hii3',hii3
!       write(*,*)'hii3a',hii3a,'hii3b',hii3b,'hii3c',hii3c
!       write(*,*)'hii4',hii4
!       write(*,*)'total',hii1+hii2+hii3+hii4
!       endif

      RETURN

   end subroutine hiievb

   !***********************************************************************

   subroutine hijevb(statei,statej)

      !     Modification History of file: hij_evb.f
      !------------------------------------------------------------------------------
      !     When | Who |      What
      !----------|-----|-------------------------------------------------------------
      !  27-05-99| hd  |- final creation
      !------------------------------------------------------------------------------
      ! This is a new subroutine that will calculate the off diagonal terms
      ! of the evb Hamiltonian

      implicit none
      integer, intent(in) :: statei, statej

      ! local variables
      integer ia,ja,ka
      integer tmpia,tmpja,tmpka
      real*8 hijtot
      real*8 vijexch,fncF,fncG
      real*8 dvijexchdr,dFdr,dGdq
      real*8 f1,f2,df1dr,df2dr
      real*8 dhijdx(natmax),dhijdy(natmax),dhijdz(natmax)
      logical debug

      !external dst,angle,queue
      !external drdx1,drdx2,drdy1,drdy2,drdz1,drdz2
      !external dwdx1,dwdx2,dwdy1,dwdy2,dwdz1,dwdz2
      !external dqdx1,dqdx2,dqdy1,dqdy2,dqdz1,dqdz2
      !external dwdx3,dwdy3,dwdz3
      !external dqdx3,dqdy3,dqdz3


      ! First calculate the A(Roo,queue) function
      ! which is in turn a function of fncf(Roo),
      ! fncg(queue)

        fncF=0.
        dFdr=0.
        f1=0.
        f2=0.
        df1dr=0.
        df2dr=0.
        dGdq=0.
        fncG=0.
        hijtot=0.
        vijexch=0.
        do ia=1,nat
          dhijdx(ia)=0.
          dhijdy(ia)=0.
          dhijdz(ia)=0.
        enddo

        if (diag(statei,statej)) then
        
        do ia=1,nat
          if (vbato(statei,statej,ia).eq.'HZ') then
            do ja=1,nat
              if (vbato(statei,statej,ja).eq.'OM') then
                do ka=1,nat
                  if (vbato(statei,statej,ka).eq.'OM'.and.ka.ne.ja) then
                    fncG=exp(-gamma*((queue(ja,ka,ia)**2)))
                    dGdq=(-2.*gamma*queue(ja,ka,ia))*(exp(-gamma*(queue(ja,ka,ia)**2)))

                    f1=(1.+knsP*(exp(-knsk*(dst(ja,ka)-DOO)**2)))
                    f2=0.5*(1.-tanh(beta*(dst(ja,ka)-capR0OO)))+&
                           10.0*exp(-alpha*(dst(ja,ka)-smlr0OO))
                    df1dr=-2.*knsK*knsP*(dst(ja,ka)-DOO)* &
                            Exp(-knsK*((dst(ja,ka)-DOO)**2))
                    df2dr=-10.0*alpha*exp(-alpha*(dst(ja,ka)-smlr0OO)) &
                          -0.5*beta*(1./(cosh(beta*(-capR0OO + &
                           dst(ja,ka)))))**2
                    fncF=f1*f2
                    dfdr=f1*df2dr+f2*df1dr

                    tmpia=ia
                    tmpja=ja
                    tmpka=ka
                  endif
                enddo
              endif
            enddo
          endif
        enddo

        do ia=1,nat
          if (vbmol(statei,statej,ia).eq.0) then
            do ja=1,nat
              if (vbmol(statei,statej,ja).ge.3) then
                vijexch=vijexch+ &
                  (qvb(statei,statej,ia)*qvb(statei,statej,ja) &
                  *Esqr)/dst(ia,ja)
                dhijdx(ja)=dhijdx(ja)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdx2(ia,ja))
                dhijdy(ja)=dhijdy(ja)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdy2(ia,ja))
                dhijdz(ja)=dhijdz(ja)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdz2(ia,ja))
                dhijdx(ia)=dhijdx(ia)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdx1(ia,ja))
                dhijdy(ia)=dhijdy(ia)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdy1(ia,ja))
                dhijdz(ia)=dhijdz(ia)+(FncF*FncG* &
                           ((-QVB(statei,statej,ia)*Esqr* &
                           QVB(statei,statej,ja))/(dst(ia,ja)**2)) &
                           *drdz1(ia,ja))
              endif
            enddo
          endif
        enddo
        
        do ia=1,nat
          if(ia.eq.tmpja) then
              dhijdx(ia)=dhijdx(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdx1(ia,tmpka,tmpia)+FncG*dFdr* &
                         drdx1(ia,tmpka))
              dhijdy(ia)=dhijdy(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdy1(ia,tmpka,tmpia)+FncG*dFdr* &
                         drdy1(ia,tmpka))
              dhijdz(ia)=dhijdz(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdz1(ia,tmpka,tmpia)+FncG*dFdr* &
                         drdz1(ia,tmpka))
          elseif(ia.eq.tmpka) then
              dhijdx(ia)=dhijdx(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdx2(tmpja,ia,tmpia)+FncG*dFdr* &
                         drdx2(tmpja,ia))
              dhijdy(ia)=dhijdy(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdy2(tmpja,ia,tmpia)+FncG*dFdr* &
                         drdy2(tmpja,ia))
              dhijdz(ia)=dhijdz(ia)+(Vijkns+vijexch)*(FncF*dGdq* &
                         dqdz2(tmpja,ia,tmpia)+FncG*dFdr* &
                         drdz2(tmpja,ia))
          elseif(ia.eq.tmpia) then
              dhijdx(ia)=dhijdx(ia)+((Vijkns+vijexch)*(FncF*dGdq* &
                         dqdx3(tmpja,tmpka,ia)))
              dhijdy(ia)=dhijdy(ia)+((Vijkns+vijexch)*(FncF*dGdq* &
                         dqdy3(tmpja,tmpka,ia)))
              dhijdz(ia)=dhijdz(ia)+((Vijkns+vijexch)*(FncF*dGdq* &
                         dqdz3(tmpja,tmpka,ia)))
        endif
      enddo
      endif

      vbham(statei,statej)=(Vijkns+vijexch)*fncF*fncG

      !write (*,*) 'offd statei',statei,statej
      !write (*,*)'vijexch,fncf,fncg',vijexch,fncF,fncG
      !write(*,*) vbham(statei,statej),'vbham(statei,statej)'

      do ia=1,nat
        dhdx(statei,statej,ia)=dhijdx(ia)
        dhdy(statei,statej,ia)=dhijdy(ia)
        dhdz(statei,statej,ia)=dhijdz(ia)
      enddo
      
      !if (statei.ne.statej) then
      !   write (*,*) 'offd statei',statei,statej
      !   write (*,*)'vijexch,fncf,fncg',vijexch,fncF,fncG
      !   tmpfac1=0.
      !   tmpfac2=(1.+0.27*exp(-11.5*(dst(2,5)-2.875)**2))*
      !;       (0.5*(1.-tanh(4.50*(dst(2,5)-3.14)))+
      !;       exp(-15.0*(dst(2,5)-1.92)))
      !   tmpfac3=exp(-1.85*queue(2,5,4)**2)
      !   write (*,*)'vijnumh,fncf,fncg',tmpfac1,tmpfac2,tmpfac3
      !endif

      return

   end subroutine hijevb

   !***********************************************************************
   SUBROUTINE COULOMBIC(ETMP,XT,YT,ZT,dvijdx,dvijdy,dvijdz,XDA,YDA,ZDA,QDA,STATEI)

      implicit real*8 (a-h,o-z)

      real*8, intent(out) :: etmp
      real*8, intent(in),  dimension(NATMAX) :: XT,YT,ZT
      real*8, intent(out), dimension(NATMAX) :: dvijdx,dvijdy,dvijdz
      real*8, intent(in) :: XDA,YDA,ZDA,QDA
      integer, intent(in) :: STATEI

      ! local variables
      INTEGER :: IA,ATO
      REAL*8  :: RIJ,RIJSQ,VIJ,FIJ,DXIJ,DYIJ,DZIJ
      REAL*8  :: XI,YI,ZI
      REAL*8  :: DYI,DZI,RI,RISQ,FI,VI
      REAL*8  :: DVIJDRIJ
      ! SOME FOR THE TWO STATE EVB PART OF THE SUBROUTINE
      ! THAT PERTAINS TO THE ELECTRON
      REAL*8 :: DERF,DEXP

      PI = 4.0d0*DATAN(1.d0)

      etmp=0.d0

      ! COPY THE NECESSARY COORDINATES AND RESET FORCES
      dxij=0.d0
      dyij=0.d0
      dzij=0.d0
      vij=0.d0
      dvijdrij=0.d0
      rijsq=0.d0
      rij=0.d0
      fij=0.d0
      DO IA=1,NAT
        DXIJ=XT(IA)-XDA
        DYIJ=YT(IA)-YDA
        DZIJ=ZT(IA)-ZDA
        RIJSQ=DXIJ*DXIJ+DYIJ*DYIJ+DZIJ*DZIJ
        RIJ=DSQRT(RIJSQ)
        VIJ=QDA*(QVB(STATEI,STATEI,IA))*Esqr/RIJ
        dvijdx(ia)=-QDA*(QVB(STATEI,STATEI,IA))/(RIJ)**2 &
                   *(-((-X0evb(IA)+XDA)/DSQRT((-X0evb(IA)+XDA)**2 &
              +(-Y0evb(IA)+YDA)**2+(-Z0evb(IA)+ZDA)**2)))
        dvijdy(ia)=-QDA*(QVB(STATEI,STATEI,IA))/(RIJ)**2 &
                   *(-((-Y0evb(IA)+YDA)/DSQRT((-X0evb(IA)+XDA)**2 &
              +(-Y0evb(IA)+YDA)**2+(-Z0evb(IA)+ZDA)**2)))
        dvijdz(ia)=-QDA*(QVB(STATEI,STATEI,IA))/(RIJ)**2 &
                   *(-((-Z0evb(IA)+ZDA)/DSQRT((-X0evb(IA)+XDA)**2 &
              +(-Y0evb(IA)+YDA)**2+(-Z0evb(IA)+ZDA)**2)))
        ETMP=ETMP+VIJ
      ENDDO

      RETURN

   END subroutine COULOMBIC

   !***********************************************************************
   !=======================================================================
   ! SET OF USEFUL PARTIAL DERIVATIVES FUNCTIONS
   ! AND OTHER SUBROUTINES THAT ARE USED IN THE EVB STATES
   !=======================================================================

   !=========ANGLE BETWEEN THREE POINTS====================================
   !======WHERE KA IS BETWEEN IA AND JA====================================
   REAL*8 FUNCTION ANGLE(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      REAL*8 :: PI
      PI = 4.0*ATAN(1.0)
      ANGLE=DACOS((SQDST(IA,JA)-SQDST(IA,KA)-SQDST(JA,KA))/(-2.*DST(IA,KA)*DST(JA,KA)))
      !CORRECT IS ANGLE IS OVER 180
      !IF (ANGLE.GE.PI) THEN
      !   ANGLE=ANGLE-PI
      !ENDIF
   END function angle

   !=========TMP FUNCTION====================================
   REAL*8 FUNCTION TFNC(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      TFNC=-X0evb(IA)*X0evb(JA)+X0evb(IA)*X0evb(KA)+X0evb(JA)*X0evb(KA)-(X0evb(KA))**2   &
           -Y0evb(IA)*Y0evb(JA)+Y0evb(IA)*Y0evb(KA)+Y0evb(JA)*Y0evb(KA)-(Y0evb(KA))**2   &
           -Z0evb(IA)*Z0evb(JA)+Z0evb(IA)*Z0evb(KA)+Z0evb(JA)*Z0evb(KA)-(Z0evb(KA))**2
   END function TFNC

   !===========(D(ANGLE)/DX0(IA))==============================================
   REAL*8 FUNCTION DWDX1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
        DWDX1=-((-((X0evb(KA)-X0evb(IA))*TFNC(IA,JA,KA)/((DST(IA,KA)**3)*  &
              DST(JA,KA)))-(X0evb(KA)-X0evb(JA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                          &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdx1

   !===========(D(ANGLE)/DX0(JA))==============================================
   REAL*8 FUNCTION DWDX2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDX2=-((-((X0evb(KA)-X0evb(JA))*TFNC(IA,JA,KA)/((DST(IA,KA))*          &
              DST(JA,KA)**3))-(X0evb(KA)-X0evb(IA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                             &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdx2

   !===========(D(ANGLE)/DX0(KA))==============================================
   REAL*8 FUNCTION DWDX3(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDX3=-((((X0evb(KA)-X0evb(JA))*TFNC(IA,JA,KA))/((DST(IA,KA))*    &
              (DST(JA,KA)**3)))+((X0evb(KA)-X0evb(IA))*TFNC(IA,JA,KA))/ &
              ((DST(IA,KA)**3)*DST(JA,KA))-((X0evb(IA)+X0evb(JA)-       &
              2.*(X0evb(KA)))/(DST(IA,KA)*DST(JA,KA))))              &
              /(SQRT(1-((TFNC(IA,JA,KA)**2)/                      &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2)))))
   END function dwdx3

   !===========(D(ANGLE)/DY0(IA))==============================================
   REAL*8 FUNCTION DWDY1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDY1=-((-((Y0evb(KA)-Y0evb(IA))*TFNC(IA,JA,KA)/((DST(IA,KA)**3)*    &
              DST(JA,KA)))-(Y0evb(KA)-Y0evb(JA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                          &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdy1

   !===========(D(ANGLE)/DZ0(IA))==============================================
   REAL*8 FUNCTION DWDZ1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDZ1=-((-((Z0evb(KA)-Z0evb(IA))*TFNC(IA,JA,KA)/((DST(IA,KA)**3)*    &
              DST(JA,KA)))-(Z0evb(KA)-Z0evb(JA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                          &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdz1

   !===========(D(ANGLE)/DY0(JA))==============================================
   REAL*8 FUNCTION DWDY2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDY2=-((-((Y0evb(KA)-Y0evb(JA))*TFNC(IA,JA,KA)/((DST(IA,KA))*          &
              DST(JA,KA)**3))-(Y0evb(KA)-Y0evb(IA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                             &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdy2

   !===========(D(ANGLE)/DZ0(JA))==============================================
   REAL*8 FUNCTION DWDZ2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDZ2=-((-((Z0evb(KA)-Z0evb(JA))*TFNC(IA,JA,KA)/((DST(IA,KA))*          &
              DST(JA,KA)**3))-(Z0evb(KA)-Z0evb(IA))/(DST(IA,KA)*DST(JA,KA)))/ &
              (SQRT(1-((TFNC(IA,JA,KA)**2)/                             &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2))))))
   END function dwdz2

   !===========(D(ANGLE)/DY0(KA))==============================================
   REAL*8 FUNCTION DWDY3(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDY3=-((((Y0evb(KA)-Y0evb(JA))*TFNC(IA,JA,KA))/((DST(IA,KA))*    &
              (DST(JA,KA)**3)))+((Y0evb(KA)-Y0evb(IA))*TFNC(IA,JA,KA))/ &
              ((DST(IA,KA)**3)*DST(JA,KA))-((Y0evb(IA)+Y0evb(JA)-       &
              2.*(Y0evb(KA)))/(DST(IA,KA)*DST(JA,KA))))              &
              /(SQRT(1-((TFNC(IA,JA,KA)**2)/                      &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2)))))
   END  function dwdy3

   !===========(D(ANGLE)/DZ0(KA))==============================================
   REAL*8 FUNCTION DWDZ3(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DWDZ3=-((((Z0evb(KA)-Z0evb(JA))*TFNC(IA,JA,KA))/((DST(IA,KA))*    &
              (DST(JA,KA)**3)))+((Z0evb(KA)-Z0evb(IA))*TFNC(IA,JA,KA))/ &
              ((DST(IA,KA)**3)*DST(JA,KA))-((Z0evb(IA)+Z0evb(JA)-       &
              2.*(Z0evb(KA)))/(DST(IA,KA)*DST(JA,KA))))              &
              /(SQRT(1-((TFNC(IA,JA,KA)**2)/                      &
              ((DST(IA,KA)**2)*(DST(JA,KA)**2)))))
   END function dwdz3

   !=========DISTANCE BETWEEN TWO POINTS===================================
   REAL*8 FUNCTION DST(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DST=(X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2
      DST=DSQRT(DST)
   END function dst

   !=========SQUARE DISTANCE BETWEEN TWO POINTS===================================
   REAL*8 FUNCTION SQDST(IA,JA)
      INTEGER, intent(in) :: IA,JA
      SQDST=(X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2
   END function sqdst

   !===========(D(DST)/DX0(IA))================================================
   REAL*8 FUNCTION DRDX1(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDX1=-((-X0evb(IA)+X0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2   &
              +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2))
   END function drdx1

   !===========(D(DST)/DX0(JA))================================================
   REAL*8 FUNCTION DRDX2(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDX2=(-X0evb(IA)+X0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2    &
             +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2)
   END function drdx2

   !===========(D(DST)/DY0(IA))================================================
   REAL*8 FUNCTION DRDY1(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDY1=-(-Y0evb(IA)+Y0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2  &
              +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2)
   END function drdy1

   !===========(D(DST)/DY0(JA))================================================
   REAL*8 FUNCTION DRDY2(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDY2=(-Y0evb(IA)+Y0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2   &
             +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2)
   END function drdy2

   !=========(D(DST)/DZ0(IA))==================================================
   REAL*8 FUNCTION DRDZ1(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDZ1=-(-Z0evb(IA)+Z0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2  &
              +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2)
   END function drdz1

   !===========(D(DST)/DZ0(JA))================================================
   REAL*8 FUNCTION DRDZ2(IA,JA)
      INTEGER, intent(in) :: IA,JA
      DRDZ2=(-Z0evb(IA)+Z0evb(JA))/SQRT((-X0evb(IA)+X0evb(JA))**2   &
             +(-Y0evb(IA)+Y0evb(JA))**2+(-Z0evb(IA)+Z0evb(JA))**2)
   END function drdz2

   !===========ASSYMETRIC STRETCH (QUEUE) FOR THREE ATOMS =================
   !====================(KA BETWEEN JA AND IA)=============================
   !====================QUEUE=ABS((ROO/2)-ROH)=============================
   REAl*8 FUNCTION QUEUE(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      QUEUE=DABS(SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+   &
             (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+   &
             (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2))
   END function QUEUE

   !===========(D(QUEUE)/DX0(IA))==============================================
   REAL*8 FUNCTION DQDX1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDX1=((X0evb(IA)-X0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+      &
              (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))-         &
              (X0evb(IA)-X0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+         &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) THEN
         DQDX1=-DQDX1
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) THEN
          DQDX1=0.
      ENDIF
   END function dqdx1

   !===========(D(QUEUE)/DX0(JA))==============================================
   REAL*8 FUNCTION DQDX2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDX2=-(X0evb(IA)-X0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+  &
             (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) THEN
         DQDX2=-DQDX2
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) THEN
         DQDX2=0.
      ENDIF
   END function dqdx2

   !===========(D(QUEUE)/DX0(KA))==============================================
   REAL*8 FUNCTION DQDX3(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDX3=(X0evb(IA)-X0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+  &
            (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+           &
               (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) THEN
           DQDX3=-DQDX3
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) THEN
          DQDX3=0.
      ENDIF
   END function dqdx3

   !===========(D(QUEUE)/DY0(IA))==============================================
   REAL*8 FUNCTION DQDY1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDY1=((Y0evb(IA)-Y0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+   &
              (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))-      &
              (Y0evb(IA)-Y0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+      &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) THEN
          DQDY1=-DQDY1
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) THEN
          DQDY1=0.
      ENDIF
   END function dqdy1

   !===========(D(QUEUE)/DZ0(IA))==============================================
   REAL*8 FUNCTION DQDZ1(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDZ1=((Z0evb(IA)-Z0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+ &
              (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))-    &
              (Z0evb(IA)-Z0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+    &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) then
          DQDZ1=-DQDZ1
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) then
          DQDZ1=0.
      ENDIF
   END function dqdz1

   !===========(D(QUEUE)/DY0(JA))==============================================
   REAL*8 FUNCTION DQDY2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDY2=-(Y0evb(IA)-Y0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+   &
             (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) then
          DQDY2=-DQDY2
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) then
          DQDY2=0.
      ENDIF
   END function dqdy2

   !===========(D(QUEUE)/DZ0(JA))==============================================
   REAL*8 FUNCTION DQDZ2(IA,JA,KA)
      INTEGER, intent(in) :: IA,JA,KA
      DQDZ2=-(Z0evb(IA)-Z0evb(JA))/(2.*SQRT((X0evb(IA)-X0evb(JA))**2+ &
              (Y0evb(IA)-Y0evb(JA))**2+(Z0evb(IA)-Z0evb(JA))**2))
      IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+          &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) then
          DQDZ2=-DQDZ2
      ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+      &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) then
          DQDZ2=0.
      ENDIF
   END function dqdz2

   !===========(D(QUEUE)/DY0(KA))==============================================
    REAL*8 FUNCTION DQDY3(IA,JA,KA)
       INTEGER, intent(in) :: IA,JA,KA
       DQDY3=(Y0evb(IA)-Y0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+ &
             (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)
       IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+         &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) then
          DQDY3=-DQDY3
       ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+     &
             (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
             (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) then
          DQDY3=0.
       ENDIF
    END function dqdy3

   !===========(D(QUEUE)/DZ0(KA))==============================================
    REAL*8 FUNCTION DQDZ3(IA,JA,KA)
       INTEGER, intent(in) :: IA,JA,KA
       DQDZ3=(Z0evb(IA)-Z0evb(KA))/SQRT((X0evb(IA)-X0evb(KA))**2+ &
             (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)
       IF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+         &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).LT.(0.)) then
          DQDZ3=-DQDZ3
       ELSEIF((SQRT((X0evb(IA)-X0evb(JA))**2+(Y0evb(IA)-Y0evb(JA))**2+     &
              (Z0evb(IA)-Z0evb(JA))**2)/2.-SQRT((X0evb(IA)-X0evb(KA))**2+  &
              (Y0evb(IA)-Y0evb(KA))**2+(Z0evb(IA)-Z0evb(KA))**2)).EQ.(0.)) then
          DQDZ3=0.
       ENDIF
    END function dqdz3

end module msevb_water
