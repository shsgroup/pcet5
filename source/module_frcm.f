      module frcm

      !----------------------------------------------------------------|
      ! Author:      Ivan Rostov (module structure by A. V. Soudackov)
      ! Date:        11/05/2003
      ! Module name: frcm
      ! Description: contains all the FRCM related routines for
      !              the calculation of the outer sphere
      !              reorganization energy matrices
      ! Format:      fixed F77 format
      ! Notes:       All FRCM common blocks are replaced by module
      !              private variables; interfaces are added.
      !----------------------------------------------------------------|
      !
      !  $Author: souda $
      !  $Date: 2011-01-04 21:00:35 $
      !  $Revision: 5.4 $
      !  $Log: not supported by cvs2svn $
      !  Revision 5.3  2011/01/04 19:56:42  souda
      !  some unimportant changes, mostly cleaning (double precision constants etc.)
      !
      !  Revision 5.2  2010/10/28 21:29:36  souda
      !  First (working and hopefully bug-free) source of PCET 5.x
      !
      !
      !----------------------------------------------------------------|

      !----------------------------------------------------------------|
      ! Use statements
      !----------------------------------------------------------------|
      use pardim
      use parsol
      use elmnts
      use strings

      !----------------------------------------------------------------|
      ! Format, data, parameter, type, e.t.c. statements
      !----------------------------------------------------------------|
      implicit real(8) (a-h,o-z)
      private
      save

      !============================================
      ! Parameters for the FRCM portion of the code
      !--------------------------------------------
      real(8),  PARAMETER :: GATINC = 2.0D0
      real(8),  PARAMETER :: RTFQ   = 1.1D0
      real(8),  PARAMETER :: GSAV   = 0.9D0

      integer, PARAMETER :: MAXHEV=300, MAXLIT=300
      integer, PARAMETER :: NUMATM=MAXHEV+MAXLIT
      integer, PARAMETER :: MAXORB=4*MAXHEV+MAXLIT
      integer, PARAMETER :: MAXPAR=3*NUMATM
      integer, PARAMETER :: MPACK=MAXORB*(MAXORB+1)/2
      integer, PARAMETER :: NMECI=10
      integer, PARAMETER :: NSECMX=200
      integer, PARAMETER :: NSF=NUMATM*GATINC
      integer, PARAMETER :: NS=NUMATM*NSECMX
      integer, PARAMETER :: NMEMAT=NMECI*(NMECI+1)/2
      integer, PARAMETER :: NMESTR=(NMECI**4+3)/4
      integer, PARAMETER :: NMI=NMECI*NMECI/2+1
      integer, PARAMETER :: NINMS=NMI*(NMI+1)/2
      integer, PARAMETER :: NTQ=360
      integer, PARAMETER :: NFQ=720
      integer, PARAMETER :: NB=500
      integer, PARAMETER :: NT=50
      integer, PARAMETER :: NF=150
      integer, PARAMETER :: NBS=2000
      integer, PARAMETER :: NBGEN=100000
      integer, PARAMETER :: NBSQ=NBS*(NBS+1)/2
      integer, PARAMETER :: NTFQ=NTQ*NFQ
      integer, PARAMETER :: NTFQ1=NTFQ*RTFQ
      integer, PARAMETER :: NDNBMX=500
      integer, PARAMETER :: NDFMAX=5000
      integer, PARAMETER :: MXRAS=100
      integer, PARAMETER :: MSZSV=NTQ*NFQ*GSAV/2+28*NB+4*NBS


      !=============================================================
      ! COMMON blocks for the FRCM portion of the code
      !-------------------------------------------------------------

      character*241 :: KEYWRD                             ! /KEYWRD/
      character*81  :: KOMENT,TITLE                       ! /TITLES/
      character*8, dimension(MAXPAR) :: SIMBOL            ! /SIMBOL/
      character*2, dimension(108) :: ELEMNT               ! /ELEMTS/

      logical :: SF                                       ! /SF/
      logical :: ISOK                                     ! /OKMANY/

      integer, dimension(NTFQ1) :: NSET, NSEF             ! /NSTSF/
      integer :: LSECT(NTQ,NFQ,3), LNUSEC(NTQ,NFQ)        ! /LSECTM/

      integer :: NUMAT                                    ! /MOLKST/
      integer :: IREG                                     ! /REGSO/
      integer :: IJ,NW(NSF)                               ! /SFE2/
      integer :: IJEL,NWEL(NSF)                           ! /SFE2EL/

      integer :: NATOMS,LABELS(NUMATM)                    ! /GEOKST/

      integer :: NN,KFI,LTH                               ! /SFE/
      integer :: ITE1,NCFINR1,ITE2,NCFINR2                ! /MODF/
      integer :: NTETM,NFIM                               ! /TETFI/

      integer :: NDLATS(NSF),NDNEB(NDNBMX)                ! /CURLEN/

      real(8) :: CORE(107)                                 ! /CORE/
      real(8) :: QATOM(NUMATM)                             ! /QATOM/

      real(8) :: TSF1,TSF2,TSF3,TSFE,TCOLA,TCOLE,TCONNU    ! /TIMSF/
      real(8) :: TIME0                                     ! /TIMING/

      real(8), public :: XX(NSF),YY(NSF),ZZ(NSF)           ! /PR/
      real(8) :: RV(NSF),QSFE(NSF)                         ! /PR/
      real(8) :: COL2,COL1,COIN,COEL                       ! /PECAR/

      real(8) :: X0(NS),Y0(NS),Z0(NS),
     ,          AS(NS),QS(NS),QCOR(NS)                    ! /SOLMAT/

      real(8) :: X0EL(NS),Y0EL(NS),Z0EL(NS),
     ,          ASEL(NS),QSEL(NS),QCOREL(NS)              ! /SOLMAEL/

      real(8) :: RVEL(NSF),QSFEEL(NSF)                     ! /PREL/
      real(8) :: FACTOR,FACTOR2                            ! /FACTSF/

      real(8),  public :: COORD(3,NUMATM)                  !
      integer, public :: NHB(3)                           ! /GEOMXYZ/

      real(8)  :: EPS,EPSEL                                ! /EPSS

      real(8)  :: CHARGE,QA(NMECI,NUMATM)                  !
      integer :: N                                        ! /CHARGE/

      real(8)  :: RADD(NSF),NADD(NSF)                      !
      integer :: NUMADD                                   ! /ADDSF/

      real(8)  :: CON,VSOLV,RSOLV,SELFCR,CHDIFF            !
      integer :: ITSE                                     ! /VOL/

      real(8) :: RADDEL(NSF)                               ! /ADDSFEL/

      real(8)  :: DATTE0(MXRAS),DPRAM(MXRAS),DRASD,DSMIN   !
      integer :: ITMA                                     ! /RESINF/

      real(8) :: COTETM(NTQ),SITETM(NTQ)
      real(8) :: COFIM(NFQ),SIFIM(NFQ),FIM(NFQ)            ! /SICO/

      real(8) :: COSBT,RAT,COSBT1,RNE,RR,RDS,RDS1         ! /SCHCON/

      real(8) :: SITET0(NB),COTET0(NB),FI0(NB),SIFI0(NB),  !
     ,          COFI0(NB),COTECR(NB),SITECR(NB),          !
     ,          TETECR(NB),RRECR(NB),DATTET(MXRAS)        ! /RESCON/

      integer :: NSQRD,NSQDUM                             !
      real(8)  :: SQRMAS(1001)                             ! /SQRMS/

      integer :: NNBNEB(NDNBMX)
      real(8)  :: COTNB(NDNBMX), SITNB(NDNBMX),            !
     ,           COFNB(NDNBMX), SIFNB(NDNBMX),            !
     ,           COECNB(NDNBMX),SIECNB(NDNBMX)            ! /NEBCOR/

      real(8) :: XXOLD(NSF),YYOLD(NSF),ZZOLD(NSF),         !
     ,          RVOLD(NSF),QSFOLD(NSF)                    ! /PROLD/

      real(8)  :: X01(NS),Y01(NS),Z01(NS),                 !
     ,           AS1(NS),QCOR1(NS)                        !
      integer :: IJ0,NW0(NSF)                             ! /SOLMT1/

      real(8)  :: X01EL(NS),Y01EL(NS),Z01EL(NS),           !
     ,           AS1EL(NS),QCOR1EL(NS)                    !
      integer :: IJ0EL, NW0EL(NSF)                        ! /SOLMT1EL/

      real(8) :: CHARC, CHARE                              ! /CHH/
      real(8) :: EN, EPO, EVAC, ESCFOL                     ! /ENUCL/
      
      real(8) :: STIME, STIM0, STIM3

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Van der Waals radii
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(8), dimension(108) :: RVDW =                    ! /RVDWS/
     + (/
     +     1.20D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.60D0, 1.50D0,
     +     1.40D0, 1.30D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.90D0, 1.85D0, 1.81D0, 1.50D0, 1.33D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.96D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 2.20D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0, 1.50D0,
     +     1.50D0, 1.50D0, 1.50D0
     +  /)

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Standard atomic masses
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      real(8), dimension(108) :: AMS =                     ! /ISTOPE/
     + (/
     +    1.00790D0,   4.00260D0,   6.94000D0,   9.01218D0,
     +   10.81000D0,  12.01100D0,  14.00670D0,  15.99940D0,  18.99840D0,
     +   20.17900D0,  22.98977D0,  24.30500D0,  26.98154D0,  28.08550D0,
     +   30.97376D0,  32.06000D0,  35.45300D0,  39.94800D0,  39.09830D0,
     +   40.08000D0,  44.95590D0,  47.90000D0,  50.94150D0,  51.99600D0,
     +   54.93800D0,  55.84700D0,  58.93320D0,  58.71000D0,  63.54600D0,
     +   65.38000D0,  69.73500D0,  72.59000D0,  74.92160D0,  78.96000D0,
     +   79.90400D0,  83.80000D0,  85.46780D0,  87.62000D0,  88.90590D0,
     +   91.22000D0,  92.90640D0,  95.94000D0,  98.90620D0, 101.07000D0,
     +  102.90550D0, 106.40000D0, 107.86800D0, 112.41000D0, 114.82000D0,
     +  118.69000D0, 121.75000D0, 127.60000D0, 126.90450D0, 131.30000D0,
     +  132.90540D0, 137.33000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +  178.49000D0, 180.94790D0,
     +  183.85000D0, 186.20700D0, 190.20000D0, 192.22000D0, 195.09000D0,
     +  196.96650D0, 200.59000D0, 204.37000D0, 207.20000D0, 208.98040D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,
     +    0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0,   0.00000D0
     +  /)


      !----------------------------------------------------------------
      ! Member subprograms
      !----------------------------------------------------------------
      public :: frcminit, solint

      contains

      !================================================================!
      SUBROUTINE SOLINT(TK,TINFK,TRK,TRINFK)
      !================================================================!
      ! Calculates the reorganization energy matrices
      ! using FRCM model (Ivan Rostov)
      !-----------------------------------------------------------------
      !
      ! souda
      ! 2010/06/25 20:02:36
      ! 4.1
      ! Exp
      ! module_frcm.f,v 4.1 2010/06/25 20:02:36 souda Exp
      ! module_frcm.f,v
      ! Revision 4.1  2010/06/25 20:02:36  souda
      ! Release 4.1
      !
      ! Revision 1.1.1.1  2004/01/08 21:18:18  souda
      ! Initial PCET-4.0 Release
      !
      ! Revision 1.2  2003/05/19 14:01:16  souda
      ! CVS Tags added
      !
      !=================================================================

      IMPLICIT REAL(8) (A-H,O-Z)
      real(8), PARAMETER :: EV2KCAL = 23.061D0

      !include 'SIZES'
      !include 'parsol.h'

      real(8), intent(out), dimension(4,4) :: TK,TINFK,TRK,TRINFK
      logical, save :: FIRST=.true., PRNT, S12DR

      !COMMON /EPSS/ EPS,EPSEL
      !COMMON /FACTSF/ FACTOR,FACTOR2
      !COMMON /SOLMAT/ X0(NS),Y0(NS),Z0(NS),AS(NS),QS(NS),QCOR(NS)
      !COMMON /SOLMAEL/ X0EL(NS),Y0EL(NS),Z0EL(NS),
      !                 ASEL(NS),QSEL(NS),QCOREL(NS)
      !COMMON /SFE2  / IJ ,NW(NSF)
      !COMMON /SFE2EL/ IJEL,NWEL(NSF)
      !COMMON /QATOM/ QATOM(NUMATM)
      !COMMON /MOLKST/ NUMAT
      !COMMON /GEOMXYZ/COORD(3,NUMATM)
      !COMMON /CHARGE/CHARGE,QA(NMECI,NUMATM),N

      real(8), dimension(NMI,NMI) :: TIN,TEL,T
      real(8), dimension(4,4) :: TKBACK,TINFKBACK

      !SAVE PRNT,S12DR
      !DATA FIRST /.TRUE./

      IF (FIRST) THEN
          FIRST  = .FALSE.
          PRNT = (INDEX(KEYWRD,'SOLINT').NE.0)
          S12DR=(INDEX(KEYWRD,'S12DR').NE.0)
          TIME0 = SECOND()
      ENDIF

      ICHOUT=6

      N1=1

      DO J=1,N

C         N1=1

         DO K=1,NUMAT
            QATOM(K)=QA(J,K)
         ENDDO

         CALL FRCMDR(N1,COORD,T(J,J),TEL(J,J),TIN(J,J))

         T(J,J)=-T(J,J)
         TEL(J,J)=-TEL(J,J)
         TIN(J,J)=-TIN(J,J)

         IF (PRNT) THEN
            WRITE(6,'(81(''-'')/)')
            WRITE(6,'(''  Calculation of interaction energy for the '',
     *      ''given solute charge distribution''/''  with '',
     *      ''polarization field described by the following '',
     *      ''parameters:''/)')
            WRITE(6,'(7x,''Dielectric permittivity in Volume 3 ='',
     *      F7.3)')EPS
            WRITE(6,'(7x,''Dielectric permittivity in Volume 2 ='',
     *      F7.3)')EPSEL
            WRITE(6,'(36x,'' KAPPA ='',F7.3)')FACTOR
            WRITE(6,'(35x,'' DELTA ='',F7.3)')FACTOR2
            WRITE(6,'(/,24X,''The integral RO('',I2,'')*FI('',I2,'') =''
     *      ,F8.4,'' eV'')')J,J,T(J,J)
            WRITE(6,'(/81(''-''))')
         ENDIF

         DO K=1,N
            IF (K.EQ.J) GOTO 1
            DO L=1,NUMAT
               QATOM(L)=QA(K,L)
            ENDDO
            CALL CONN(COORD,X0,Y0,Z0,QS,IJ,TIN(J,K))
            TIN(J,K)=-TIN(J,K)
            CALL CONN(COORD,X0EL,Y0EL,Z0EL,QSEL,IJEL,TEL(J,K))
            TEL(J,K)=-TEL(J,K)
            IF (PRNT) THEN
               WRITE(6,*)'N1=',N1
               WRITE(6,*)'J,K,TIN(J,K):',J,K,TIN(J,K)
               WRITE(6,*)'J,K,TEL(J,K):',J,K,TEL(J,K)
            ENDIF
            IF(FACTOR2.EQ.0.0D0) THEN
               T(J,K)=TIN(J,K)
            ELSE
               T(J,K)=TIN(J,K)+TEL(J,K)
            ENDIF
    1       CONTINUE
         ENDDO

         !N1=2

         IF (S12DR) THEN

            IF (PRNT) THEN
               WRITE(6,'(/'' APPROXIMATION S_12=0 applied during '',
     *         ''solution of FRCM equations.'')')
               WRITE(6,'('' thus RO*FI_el =RO*FI_1, RO*FI_in = RO*FI_2''
     *         )')
            ENDIF

         ELSE

            COL2OLD=COL2
            COINOLD=COIN
            COL1OLD=COL1
            COL2=0.0D0
            COIN=0.0D0
            COL1=COEL
            EPS2=0.0D0
            DO K=1,NUMAT
               QATOM(K)=QA(J,K)
            ENDDO
            CALL FRCMDR(N1,COORD,DUMMY1,TEL(J,J),DUMMY2)
            TEL(J,J)=-TEL(J,J)

            IF (PRNT) THEN
               WRITE(6,'(81(''-'')/)')
               WRITE(6,'(''  Calculation of interaction energy for '',
     *         ''the given solute charge distribution''/''  with '',
     *         ''polarization field described by the following '',
     *         ''parameters:''/)')
               WRITE(6,'(7x,''Dielectric permittivity in Volume 3 ='',
     *         F7.3)')EPSEL
               WRITE(6,'(7x,''Dielectric permittivity in Volume 2 ='',
     *         F7.3)')EPS2
               WRITE(6,'(36x,'' KAPPA ='',F7.3)')FACTOR
               WRITE(6,'(35x,'' DELTA ='',F7.3)')FACTOR2
               WRITE(6,'(/,24X,''The integral RO('',I2,'')*FI_el('',I2,
     *         '') ='',F8.4,'' eV'')')J,J,TEL(J,J)
               WRITE(6,'(/81(''-''))')
            ENDIF

            DO K=1,N
               IF (K.EQ.J) GOTO 2
               DO L=1,NUMAT
                  QATOM(L)=QA(K,L)
               ENDDO
               CALL CONN(COORD,X0EL,Y0EL,Z0EL,QSEL,IJEL,TEL(J,K))
      	       TEL(J,K)=-TEL(J,K)
               IF (PRNT) THEN
                  WRITE(6,*)'N1=',N1
                  WRITE(6,*)'J,K,TEL(J,K):',J,K,TEL(J,K)
               ENDIF
    2          CONTINUE
            ENDDO

            COL2=COL2OLD
            COIN=COINOLD
            COL1=COL1OLD

         ENDIF

      ENDDO !J

      IF (PRNT) THEN
         WRITE(6,'(/81(''-''))')
         WRITE (6,'(/'' MATRIX T (Y-REPRESENTATION), eV'')')
         CALL PRMATR(ICHOUT,NMI,N,T)
         WRITE (6,'(/'' MATRIX T_EL (Y-REPRESENTATION), eV'')')
         CALL PRMATR(ICHOUT,NMI,N,TEL)
         WRITE(6,'(/81(''-''))')
      ENDIF
C
C     SYMMETRIZATION OF T-MATRICES AND CALCULATION OF T_in MATRIX
C
      DO I=1,N
         DO J=I+1,N

            T(I,J)=(T(I,J)+T(J,I))*0.5D0
            T(J,I)=T(I,J)

            TEL(I,J)=(TEL(I,J)+TEL(J,I))*0.5D0
            TEL(J,I)=TEL(I,J)

         ENDDO
         DO J=1,N
            TIN(I,J)=T(I,J)-TEL(I,J)
         ENDDO
      ENDDO

      IF (PRNT) THEN
         WRITE(6,'(/81(''-''))')
         WRITE(6,'(/('' AFTER SYMMETRIZATION:''/))')
         WRITE (6,'(/'' MATRIX T (Y-REPRESENTATION), eV'')')
         CALL PRMATR(ICHOUT,NMI,N,T)
         WRITE (6,'(/'' MATRIX T_EL (Y-REPRESENTATION), eV'')')
         CALL PRMATR(ICHOUT,NMI,N,TEL)
         WRITE (6,'(/'' MATRIX T_IN (Y-REPRESENTATION), eV'')')
         CALL PRMATR(ICHOUT,NMI,N,TIN)
         ERIN=(TIN(1,1)+TIN(2,2)-TIN(1,2)*2.D0)*0.5D0
         WRITE (6,'(/'' Reorganization energy ='',F10.4)')ERIN
         WRITE(6,'(/81(''-''))')
      ENDIF

      DO I=1,N
         DO J=1,N
	    T(I,J)=T(I,J)*EV2KCAL
	    TEL(I,J)=TEL(I,J)*EV2KCAL
	    TIN(I,J)=TIN(I,J)*EV2KCAL
	 ENDDO
      ENDDO

      WRITE(6,'(/('' IN KCAL/MOL:''/))')
      WRITE (6,'(/'' MATRIX T (Y-REPRESENTATION), KCAL/MOL'')')
      CALL PRMATR(ICHOUT,NMI,N,T)
      WRITE (6,'(/'' MATRIX T_EL (Y-REPRESENTATION), KCAL/MOL'')')
      CALL PRMATR(ICHOUT,NMI,N,TEL)
      WRITE (6,'(/'' MATRIX T_IN (Y-REPRESENTATION), KCAL/MOL'')')
      CALL PRMATR(ICHOUT,NMI,N,TIN)
      ERIN12=(TIN(1,1)+TIN(2,2)-TIN(1,2)*2.D0)*0.5D0
      ERIN34=(TIN(3,3)+TIN(4,4)-TIN(3,4)*2.D0)*0.5D0
      ERIN13=(TIN(1,1)+TIN(3,3)-TIN(1,3)*2.D0)*0.5D0
      ERIN24=(TIN(2,2)+TIN(4,4)-TIN(2,4)*2.D0)*0.5D0
      ERIN14=(TIN(1,1)+TIN(4,4)-TIN(1,4)*2.D0)*0.5D0
      ERIN23=(TIN(2,2)+TIN(3,3)-TIN(2,3)*2.D0)*0.5D0
      WRITE (6,'(/'' Table of reorganization energies for all possible''
     *,'' couples of diabatic states''//
     *''    E_s(1,2)  E_s(3,4)  E_s(1,3)  E_s(2,4)  E_s(1,4)  E_s(2,3)''
     *)')
      WRITE (6,'(2x,F9.4,5F10.4)')ERIN12,ERIN34,ERIN13,ERIN24,ERIN14,
     *ERIN23
      WRITE(6,'(/81(''-''))')
C
C     TRANSFER TO SASHAS ARRAYS
C
      DO I=1,4
         DO J=1,4
	    TK(I,J)=TIN(I,J)
	    TINFK(I,J)=TEL(I,J)
         ENDDO
      ENDDO

C==(HYD) 19 JULY 2001: NEW OPTION TO TAKE CARE OF SOLVATION MATRICES
C========(HYD) NOW TAKING CARE OF SOLVATION MATRIX OPTIONS==============
C======ALL MANIPULATIONS DONE USING ONE TRIANGLE OF MATRICES============
 
C=================================SYMT OPTION===========================
      IF (SYMT) THEN
C  first two line to ensure property amongst matrix elements respected
        TK(2,3)=TK(1,2)+TK(1,3)-TK(2,2)
        TK(1,4)=TK(1,2)+TK(1,3)-TK(1,1)
        TINFK(2,3)=TINFK(1,2)+TINFK(1,3)-TINFK(2,2)
        TINFK(1,4)=TINFK(1,2)+TINFK(1,3)-TINFK(1,1)
C  Now symmetrizing matrix for a symmetric system
        TK(4,4)=TK(1,1)
        TK(3,3)=TK(2,2)
        TK(3,4)=TK(1,2)
        TK(2,4)=TK(1,3)
        TINFK(4,4)=TINFK(1,1)
        TINFK(3,3)=TINFK(2,2)
        TINFK(3,4)=TINFK(1,2)
        TINFK(2,4)=TINFK(1,3)
       ENDIF
C================================SYMPT OPTION===========================
      IF (SYMPT) THEN
        TK(2,2)=TK(1,1)
        TINFK(2,2)=TINFK(1,1)
      ENDIF
C================================SYMET OPTION===========================
      IF (SYMET) THEN
        TK(3,3)=TK(1,1)
        TINFK(3,3)=TINFK(1,1)
      ENDIF
 
C========(HYD) NOW TAKING CARE OF REMOVING 2B DEPENDENSE================
C============USING A SET OF REDUCED DENSITIES===========================
      IF (REDDENS) THEN
        TK(1,4)=TK(1,2)+TK(1,3)-TK(1,1)
        TK(2,4)=TK(2,2)+TK(2,3)-TK(1,2)
        TK(3,4)=TK(2,3)+TK(3,3)-TK(1,3)
        TK(4,4)=TK(1,1)+TK(2,2)+TK(3,3) + 2.D0*(TK(2,3)-
     ;          TK(1,2)-TK(1,3))
        TINFK(1,4)=TINFK(1,2)+TINFK(1,3)-TINFK(1,1)
        TINFK(2,4)=TINFK(2,2)+TINFK(2,3)-TINFK(1,2)
        TINFK(3,4)=TINFK(2,3)+TINFK(3,3)-TINFK(1,3)
        TINFK(4,4)=TINFK(1,1)+TINFK(2,2)+TINFK(3,3)
     ;             + 2.D0*(TINFK(2,3)-TINFK(1,2)-TINFK(1,3))
      ENDIF
C==========(HYD) END SOLVATION MATRIX OPTIONS===========================


C============(HYD) NOW SYMMETRIZING RELATIVE TO DIAGONAL================
C============MUST ALWAYS BE DONE ELSE WHOLE APPROACH NOT VALIDATED======
C======NOSYMD SHOULD BE USED TO CHECK ONLY TO SEE HOW VALID APPROACH IS
      IF (NOSYMD) THEN
        GOTO 99
      ELSE
      DO I=1,4
        IF ((I+1).LT.4) THEN
          DO J=I+1,4
            TK(J,I)=TK(I,J)
            TINFK(J,I)=TINFK(I,J)
          ENDDO
        ENDIF
      ENDDO
      ENDIF
 
C (END HYD SECTION)=====================================================
 
C
C     The reduced matrices (TR,TRINF) are defined as
C
C     t'_{11} = t_{11}
C     t'_{1i} = t_{1i} - t_{11}                    (i .ne. 1)
C     t'_{ij} = t_{ij} + t_{11} - t_{i1} - t_{j1}  (i,j .ne. 1)
C
C     and
C
C     t_{ij} = T_{ii,jj} = -INT{\phi_{ii} \rho_{jj}}

C     Construction of reduced matrices [t']

   99 T00  = TK(1,1)
      T008 = TINFK(1,1)
      TRK(1,1) = T00
      TRINFK(1,1) = T008

      DO I=2,4
         TRK(I,1) = TK(I,1) - T00
         TRK(1,I) = TRK(I,1)
         TRINFK(I,1) = TINFK(I,1) - T008
         TRINFK(1,I) = TRINFK(I,1)
      ENDDO

      DO I=2,4
         DO J=2,4
            TRK(I,J) = TK(I,J) + T00 - TK(I,1) - TK(J,1)
            TRINFK(I,J) = TINFK(I,J) + T008 - TINFK(I,1) - TINFK(J,1)
         ENDDO
      ENDDO

      CLOSE (71)
      RETURN
      END SUBROUTINE SOLINT

!======================================================================!
      SUBROUTINE FRCMINIT(ISFILE,EPS0,EPS8,KAPPA,DELTA,CHRTOT,
     ,                    NAT,IPT,LAB,XYZ,CHR)
C======================================================================C
C     In case of FRCM solvation model:
C     - reads specific FRCM keywords
C     - initializes internal COMMON-blocks used in FRCM
C       calculations of the reorganization energy matrices
C======================================================================C

      IMPLICIT REAL(8) (A-H, O-Z)
      !include 'SIZES'
      !include 'pardim.h'
      !include 'elmnts.h'

      integer, intent(in)                     :: ISFILE, NAT
      integer, intent(in), dimension(:)       :: IPT, LAB
      real(8), intent(in)                     :: EPS0,EPS8,KAPPA,DELTA
      real(8), intent(inout), dimension(:,:)  :: XYZ
      real(8), intent(in), dimension(:,:)     :: CHR
      real(8), intent(in)                     :: CHRTOT

      !DIMENSION IPT(*),LAB(*),XYZ(3,*),CHR(NELST,*)
      !CHARACTER KEYWRD*241, KOMENT*81, TITLE*81
      !CHARACTER SPACE*1, SPACE2*2, CH*1, CH2*1, ECHO*80
      !CHARACTER ELEMNT*2
      !LOGICAL ISOK

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     FRCM-MOPAC COMMON-blocks
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !COMMON /KEYWRD/  KEYWRD
      !COMMON /TITLES/  KOMENT,TITLE
      !COMMON /ELEMTS/  ELEMNT(108)
      !COMMON /OKMANY/  ISOK
      !COMMON /GEOMXYZ/ COORD(3,NUMATM),NHB(3)
      !COMMON /GEOKST/  NATOMS,LABELS(NUMATM)
      !COMMON /FACTSF/  FACTOR,FACTOR2
      !COMMON /EPSS/    EPS,EPSEL
      !COMMON /MOLKST/  NUMAT
      !COMMON /CHARGE/  CHARGE,QA(NMECI,NUMATM),N

      !SAVE SPACE, SPACE2
      !DATA SPACE, SPACE2/' ','  '/

      !character*1, save :: space  = " "
      character*2, save :: space2 = "  "
      character*1 :: ch, ch2
      character*80 :: echo


      CALL GETTXT(ISFILE)

      IF(INDEX(KEYWRD,'ECHO').NE.0)THEN
         REWIND ISFILE
         ISOK=.FALSE.
         DO 50 I=1,1000
            READ(ISFILE,'(A)',END=60) ECHO
            DO 20 J=80,2,-1
   20       IF(ECHO(J:J).NE.' ') GOTO 30
            J=1
   30       DO 40 K=1,J
   40       IF(ICHAR(ECHO(K:K)).LT.32) ECHO(K:K)='*'
            WRITE(6,'(1X,A)') ECHO(1:J)
   50    CONTINUE
   60    CONTINUE
         REWIND ISFILE
         CALL GETTXT(ISFILE)
      ENDIF

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     FRCM specific keywords
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(INDEX(KEYWRD,'ECHO').NE.0)WRITE(6,'(''1'')')
      IF(KEYWRD(1:1) .NE. SPACE) THEN
         CH=KEYWRD(1:1)
         KEYWRD(1:1)=SPACE
         DO 70 I=2,239
            CH2=KEYWRD(I:I)
            KEYWRD(I:I)=CH
            CH=CH2
            IF(KEYWRD(I+1:I+2) .EQ. SPACE2) THEN
               KEYWRD(I+1:I+1)=CH
               GOTO 80
            ENDIF
   70    CONTINUE
         CH2=KEYWRD(240:240)
         KEYWRD(240:240)=CH
         KEYWRD(241:241)=CH2
   80    CONTINUE
      ENDIF

      CALL WRTKEY(KEYWRD)

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Comments
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (KOMENT(1:1) .NE. SPACE) THEN
         CH=KOMENT(1:1)
         KOMENT(1:1)=SPACE
         DO 90 I=2,79
            CH2=KOMENT(I:I)
            KOMENT(I:I)=CH
            CH=CH2
            IF(KOMENT(I+1:I+2) .EQ. SPACE2) THEN
               KOMENT(I+1:I+1)=CH
               GOTO 100
            ENDIF
   90    CONTINUE
         CH2=KOMENT(80:80)
         KOMENT(80:80)=CH
         KOMENT(81:81)=CH2
  100    CONTINUE
      ENDIF

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Title
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF(TITLE(1:1) .NE. SPACE) THEN
         CH=TITLE(1:1)
         TITLE(1:1)=SPACE
         DO 110 I=2,79
            CH2=TITLE(I:I)
            TITLE(I:I)=CH
            CH=CH2
            IF(TITLE(I+1:I+2) .EQ. SPACE2) THEN
               TITLE(I+1:I+1)=CH
               GOTO 120
            ENDIF
  110    CONTINUE
         CH2=TITLE(80:80)
         TITLE(80:80)=CH
         TITLE(81:81)=CH2
  120    CONTINUE
      ENDIF

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     WRITE HEADER
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      WRITE(6,'(81(''*'')/1X,''*'')')
      WRITE(6,'('' *  FRCM SOLVATION CODES WRITTEN BY I.V. ROSTOV'')')
      WRITE(6,'(1X,''*''/'' *  Charge distribution of the solute is'',
     *'' described as a sum of point charges'')')
      WRITE(6,'('' *  at the centers of atomic spheres'')')
      WRITE(6,'(1X,''*''/81(''*''))')

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Initialize cartesian coordinates and atomic labels
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NUMAT  = NAT
      NATOMS = NAT
      DO I=1,NUMAT
         DO J=1,3
            COORD(J,I) = XYZ(J,I)
         ENDDO
         LABELS(I) = LAB(I)
      ENDDO

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Initialize charges for EVB states
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CHARGE = CHRTOT
      N = NELST
      DO I=1,NUMAT
         DO J=1,N
            QA(J,I) = CHR(J,I)
         ENDDO
      ENDDO

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Initialize  PT interface triade
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO I=1,3
         NHB(I) = IPT(I)
      ENDDO
         
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Initialize dielectric constants and cavity parameters
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      EPS     = EPS0
      EPSEL   = EPS8
      FACTOR  = KAPPA
      FACTOR2 = DELTA

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Initialize atomic symbols
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO I=1,108
         ELEMNT(I) = ELSYM(I)(3:4)
      ENDDO

C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     build spheres around atoms
C     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL GETSFE(ISFILE,XYZ,NAT)

      RETURN
      END SUBROUTINE FRCMINIT

C======================================================================C
      SUBROUTINE FRCMDR (NITER_,COORD_,EPO_,EPOEL_,EPOIN_)
************************************************************************
*
*     THIS SUBPROGRAM IS THE DRIVER OF SOLVATION CALCULATIONS
C
C     SFERE STARTS MIERTUSH'S METHOD OF INTERACTION MOLECULES WITH
C     SOLVENT VIA DIELECTRIC CONTINUUM AND INITIALIZE PRINTING OF
C     POLARIZATION INFORMATION
C
C   ON INPUT
C             NUMAT  = NUMBER OF  ATOMS.
C    COOR(NUMATM,3)  = COORDINATE OF  ATOMS (IN ANGSTROMS).
C
************************************************************************
*
      IMPLICIT REAL(8) (A-H,O-Z)
      !INCLUDE 'SIZES'

      integer, intent(inout) :: NITER_
      real(8), intent(in), dimension(3,NUMATM) :: COORD_
      real(8), intent(inout) :: EPO_, EPOEL_, EPOIN_

      LOGICAL, save :: FIRST=.true., TIMES, PRNT
      LOGICAL :: BIGPRT

      !CHARACTER*241 KEYWRD
      !COMMON /PECAR/  COL2,COL1,COIN,COEL
      !COMMON /SOLMAT/ X0(NS),Y0(NS),Z0(NS),AS(NS),QS(NS),QCOR(NS)
      !COMMON /SOLMAEL/ X0EL(NS),Y0EL(NS),Z0EL(NS),ASEL(NS),
      !                 QSEL(NS),QCOREL(NS)
      !COMMON /PR   /XX(NSF),YY(NSF),ZZ(NSF),RV(NSF),QSFE(NSF)
      !COMMON /PREL/RVEL(NSF),QSFEEL(NSF)
      !COMMON /KEYWRD/ KEYWRD
      !COMMON /TIMING/ TIME0
      !COMMON /FACTSF/ FACTOR,FACTOR2
      !COMMON /TIMSF / TSF1,TSF2,TSF3,TSFE,TCOLA,TCOLE,TCONNU
      !COMMON /REGSO / IREG
      !COMMON /SFE   / NN,KFI,LTH
      !COMMON /SF    / SF
      !COMMON /SFE2  / IJ ,NW(NSF)
      !COMMON /SFE2EL/ IJEL,NWEL(NSF)
      !COMMON /MOLKST/ NUMAT
      !DIMENSION COORD_(3,NUMATM)


      IF ( FIRST ) THEN
          FIRST  = .FALSE.
          SF=(INDEX(KEYWRD,' SFERE') .NE. 0)
          PRNT=((INDEX(KEYWRD,'DEBUG').NE.0).AND.SF)
          TIMES  = (INDEX(KEYWRD,'TIMES') .NE. 0 )
          TSFE=0.0D0
      ENDIF
      BIGPRT=(INDEX(KEYWRD,'LARGE') .NE. 0 ) .AND. PRNT
      IF ( NITER_.EQ.1 ) THEN
          SF=.TRUE.
      ENDIF
      STIME=SECOND()
      TIMEA=SECOND()
C
C     SFERA1 STARTS ONLY AT FIRST CALL.
C
      IF (SF) THEN
         CALL SFERA1T
         WRITE(6,'(/'' Total number of tesserae produced:''/
     *   '' on inner cavity ='',I6/'' on outer cavity ='',I6/)')IJEL,IJ
         WRITE(6,'( //''  N.SF'',5X,''X'',9X,''Y'',9X,''Z'',8X,''RV1''
     *   ,7X,''RV2'',3X,''N.TES1'',1X,''N.TES2'')')
         WRITE(6,'(I4,5F10.4,2I6)')(I,XX(I),YY(I),ZZ(I),RVEL(I),RV(I),
     *   NWEL(I+1)-NWEL(I),NW(I+1)-NW(I),I=1,NN)
         CALL VOLSQU(X0EL,Y0EL,Z0EL,ASEL,XX,YY,ZZ,RVEL,NN,NWEL,
     *   SQC1,VMOLC1,PLENGTH1)
         WRITE (6,'(/7X,''SQUARE OF INTERNAL CAVITY  ='',F13.3,
     *    ''  ANGSTROM**2'')')SQC1
         WRITE (6,'(7X,''VOLUME OF INTERNAL CAVITY  ='',F13.3,
     *   ''  ANGSTROM**3'')')VMOLC1
         WRITE (6,'(5X,''MEAN RADIUS OF INNER CAVITY  ='',F13.3,
     *   ''  ANGSTROM'')')PLENGTH1
         IF (FACTOR2.NE.0.D0) THEN
            CALL VOLSQU(X0,Y0,Z0,AS,XX,YY,ZZ,RV,NN,NW,
     *      SQC2,VMOLC2,PLENGTH2)
            WRITE (6,'(/7X,''SQUARE OF EXTERNAL CAVITY  ='',F13.3,
     *      ''  ANGSTROM**2'')')SQC2
            WRITE (6,'(7X,''VOLUME OF EXTERNAL CAVITY  ='',F13.3,
     *      ''  ANGSTROM**3'')')VMOLC2
            WRITE (6,'(5X,''MEAN RADIUS OF OUTER CAVITY  ='',F13.3,
     *      ''  ANGSTROM'')')PLENGTH2
         ENDIF
         SF=.FALSE.
      END IF
      CALL CONNUDT(COORD_,X0,Y0,Z0,QCOR,AS,NW,COIN,XX,YY,ZZ,RV)
      CALL CONNUDT(COORD_,X0EL,Y0EL,Z0EL,QCOREL,ASEL,NWEL,COEL,
     *   XX,YY,ZZ,RVEL)
      IF(BIGPRT) THEN
         WRITE (6,'(//''  CAVITY FOR INERTIAL POLARIZATION'')')
         WRITE (6,'(''   COORDINATES AND CHARGES OF TESSERAE AFTER'',
     *   '' CONNUDT''/4X,''IJ'',6X,''X0'',9X,
     *   ''Y0'',9X,''Z0'',12X,''AS'',12X,''QCOR''/)')
         WRITE (6,'( 2X,I4,3F11.5,2F14.7)') (II,X0(II),Y0(II),Z0(II),
     *   AS(II),QCOR(II),II=1,IJ)
         WRITE (6,'(//''  CAVITY FOR INERTIALLESS POLARIZATION'')')
         WRITE (6,'(''   COORDINATES AND CHARGES OF TESSERAE2 AFTER'',
     *   '' CONNUDT''/2X,''IJEL'',5X,''X0EL'',7X,
     *   ''Y0EL'',7X,''Z0EL'',10X,''ASEL'',9X,''QCOREL''/)')
         WRITE (6,'( 2X,I4,3F11.5,2F14.7)') (II,X0EL(II),Y0EL(II),
     *   Z0EL(II),ASEL(II),QCOREL(II),II=1,IJEL)
      END IF

      DO I=1,NS
         QS(I)=QCOR(I)
      ENDDO

      IF (FACTOR2.EQ.0.D0.AND.COIN.NE.0.D0) THEN
         DO I=1,NS
            QSEL(I)=0.D0
         ENDDO
      ELSE
         DO I=1,NS
            QSEL(I)=QCOREL(I)
         ENDDO
      ENDIF

      CALL SFERA3T(X0,Y0,Z0,AS,QS,NW,IJ,X0EL,Y0EL,Z0EL,ASEL,
     *   QSEL,NWEL,IJEL)

      CALL CONN(COORD_,X0,Y0,Z0,QS,IJ,EPOIN_)
      CALL CONN(COORD_,X0EL,Y0EL,Z0EL,QSEL,IJEL,EPOEL_)
      IF (PRNT) THEN
         WRITE(6,*)' FRCMDR: EPOIN =',EPOIN_
         WRITE(6,*)' FRCMDR: EPOEL =',EPOEL_
      ENDIF
      IF(FACTOR2.EQ.0.0D0) THEN
         EPO_=EPOIN_
      ELSE
         EPO_=EPOIN_+EPOEL_
      ENDIF
      NITER_=NITER_+1
      TIMEB=SECOND()
      TSFE=TSFE+TIMEB-TIMEA
      IF(TIMES) THEN
          WRITE(6,'(''##### TIME TOTAL IS  '',F8.2,
     1    ''    SFERE   '',F8.2)')TIMEB-TIME0,TIMEB-TIMEA
      END IF
      RETURN
      END SUBROUTINE FRCMDR


!======================================================================!
      SUBROUTINE CONN (COORD_,X0_,Y0_,Z0_,QS_,IJ_,EPO_)
C*---------------------------------------------------------------------
C*     CONN CALCULATES ADDING ENEGY TO TOTAL ONES FROM CORES - TESERA
C*          INTERACTION (  IN EV  ).
C*----------------------------------------------------------------------
      IMPLICIT REAL(8) (A-H,O-Z)
      !INCLUDE 'SIZES'

      real(8),  intent(in), dimension(:,:) :: COORD_
      real(8),  intent(in), dimension(:) :: X0_, Y0_, Z0_, QS_
      integer, intent(in) :: IJ_
      real(8),  intent(inout) :: EPO_

      LOGICAL, save :: FIRST=.true.,PRINT
      real(8),  save :: second = 0.529167D00

      !COMMON /KEYWRD/ KEYWRD
      !COMMON /QATOM/ QATOM(NUMATM)
      !COMMON /MOLKST/ NUMAT
      !DIMENSION COORD_(3,*)
      !DIMENSION X0_(*),Y0_(*),Z0_(*),QS_(*)

      !SAVE PRINT
      !DATA SECOND/0.529167D00 /
      !DATA FIRST /.TRUE./

      IF (FIRST) THEN
          FIRST=.FALSE.
          PRINT=(INDEX(KEYWRD,'CONN') .NE. 0)
      END IF

      EPO_=0.D0

      DO M=1,NUMAT
C          NI=NAT(M)
          DO J=1,IJ_
              BX=COORD_(1,M)-X0_(J)
              BY=COORD_(2,M)-Y0_(J)
              BZ=COORD_(3,M)-Z0_(J)
              RA =DSQRT(BX*BX+BY*BY+BZ*BZ)
C
C             RA    IN  ANGSTREM !!!      RAU IN  ATOMIC  UNITS !!!
C
C              EPO_ = EPO_ + QS_(J)*CORE(NI)/RA
              EPO_ = EPO_ + QS_(J)*QATOM(M)/RA
         ENDDO
      ENDDO

      EPO_=EPO_*27.21D00*SECOND

      IF (PRINT) THEN
C          WRITE(6,'('' POLARIZATION ENERGY FROM CORES '',F14.5,
C     *    ''(EV)'')') EPO_
          WRITE(6,'('' CONN: POLARIZATION ENERGY ='',F14.5,
     *    ''(EV)'')') EPO_
      END IF

      RETURN
      END SUBROUTINE CONN


!======================================================================!
      SUBROUTINE CONNUDT(COORD_,X0_,Y0_,Z0_,QCOR_,AS_,NW_,CON_,
     ,                          XX_,YY_,ZZ_,RV_)
C*---------------------------------------------------------------------
C*     THE PROGRAM CALCULATES POTENTIAL FROM ATOMS AS POINT
C*     CHARGES POTENTIALS IS CALCULATED IN  A. U.
C*---------------------------------------------------------------------
      IMPLICIT REAL(8) (A-H,O-Z)
      !INCLUDE 'SIZES'

      real(8), intent(in), dimension(:,:) :: COORD_
      real(8), intent(in), dimension(:) :: X0_,Y0_,Z0_,AS_,
     ,                                    XX_,YY_,ZZ_,RV_
      real(8), intent(in) :: CON_
      integer, intent(in), dimension(:) :: NW_

      real(8), intent(inout), dimension(:) :: QCOR_

      !DIMENSION X0_(*),Y0_(*),Z0_(*),QCOR_(*),AS_(*),NW_(*)
      !DIMENSION XX_(*),YY_(*),ZZ_(*),RV_(*)

      !COMMON /MOLKST/ NUMAT
      !COMMON /QATOM/ QATOM(NUMATM)
      !COMMON /SFE   / NN
      !COMMON /KEYWRD/ KEYWRD
      !COMMON /TIMSF / TSF1,TSF2,TSF3,TSFE,TCOLA,TCOLE,TCONNU
      !COMMON /TIMING/ TIME0

      !DATA THRE/0.00001D0/,FIRST /.TRUE./
      real(8),  save :: THRE = 0.00001D0
      logical, save :: FIRST=.true.

      IF (FIRST ) THEN
         FIRST=.FALSE.
         TCONNU=0.0D0
      ENDIF

      TIMEA = SECOND()

      PI=DATAN(1.0D0)*4.0D0
      CONTH=CON_/(PI*4.0D0)
      CHARC_=0.0D0
      DO ISF=1,NN
          DO IL=NW_(ISF),NW_(ISF+1)-1
              QCOR_(IL)=0.0D00
              AX=X0_(IL)-XX_(ISF)
              AY=Y0_(IL)-YY_(ISF)
              AZ=Z0_(IL)-ZZ_(ISF)
              RVISF=RV_(ISF)
              DO IBT=1,NUMAT
                  BX=X0_(IL)-COORD_(1,IBT)
                  BY=Y0_(IL)-COORD_(2,IBT)
                  BZ=Z0_(IL)-COORD_(3,IBT)
C
                  RA2 = BX*BX+BY*BY+BZ*BZ
                  IF ( RA2.LT.THRE ) THEN
                      DIST=DSQRT(RA2)
                      WRITE (6,'(''  STOP IN CONNUD DISTANCE BETWEEN'',
     *                I6,'' SECTOR AND'',I6, ''ATOM  IS '',D12.3)')
     *                IL,IBT,DIST
                      STOP
                  END IF
                  COF=(AX*BX+AY*BY+AZ*BZ)/DSQRT(RA2)/RVISF
C                  QCOR_(IL)=QCOR_(IL)+CORE(NAT(IBT))/RA2*COF
                  QCOR_(IL)=QCOR_(IL)+QATOM(IBT)/RA2*COF
             ENDDO
             QCOR_(IL)=QCOR_(IL)*CONTH*AS_(IL)
             CHARC_=CHARC_+QCOR_(IL)
         ENDDO
      ENDDO
      TIMEB=SECOND()
      TCONNU=TCONNU+TIMEB-TIMEA

      RETURN
      END SUBROUTINE CONNUDT

!======================================================================!
      SUBROUTINE GETRAD(IREAD,NRAD,LABRAD,RADRAD)
      IMPLICIT REAL(8) (A-H,O-Z)
      !INCLUDE 'SIZES'
      
      integer, intent(in) :: IREAD
      integer, intent(out) :: NRAD
      integer, intent(inout), dimension(:) :: LABRAD
      real(8),  intent(inout), dimension(:) :: RADRAD

      !DIMENSION RADRAD(*),LABRAD(*)
************************************************************************
*
*   GETRAD READS IN THE VdW RADII. THE ELEMENT IS SPECIFIED BY IT'S
*          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
*
*  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
*
* ON OUTPUT LABRAD = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
*           RADRAD = WdV RADII OF ATOMS
*           NRAD   = THE NUMBER OF ATOMS READ
************************************************************************
      !CHARACTER KEYWRD*241, SIMBOL*8
      !COMMON /SIMBOL/ SIMBOL(MAXPAR)
      !COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
      !1                NA(NUMATM), NB(NUMATM), NC(NUMATM)
      !COMMON /ADDSF / RADD(NSF),NADD(NSF),NUMADD
      !COMMON /KEYWRD/ KEYWRD
      !CHARACTER ELEMNT*2
      !1COMMON /ELEMTS/  ELEMNT(108)

      integer, DIMENSION(40) :: ISTART
      LOGICAL :: LEADSP
      CHARACTER :: LINE*80, STRING*80, ELE*2
      character(1), save :: NINE="9", ZERO="0"

      !DATA COMMA,SPACE,NINE,ZERO/',',' ','9','0'/

      !TAB=CHAR(9)
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
      ICAPZ = ICHAR('Z')
      MAXTXT=0
      NRAD=0
      ISERR=0

      DO 10 I=1,MAXPAR
   10 SIMBOL(I)= '---'

   20 READ(IREAD,'(A)',END=130,ERR=240)LINE
      IF(LINE.EQ.' ') GO TO 130

*   CLEAN THE INPUT DATA
************************************************************************
      DO 30 I=1,80
         ILINE=ICHAR(LINE(I:I))
c         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
c            LINE(I:I)=CHAR(ILINE+ICAPA-ILOWA)
c         ENDIF
   30 CONTINUE
************************************************************************
      ICOMMA=ICHAR(COMMA)
      ITAB=ICHAR(TAB)
      DO 40 I=1,80
         KHAR=ICHAR(LINE(I:I))
         IF(KHAR.EQ.ICOMMA.OR.KHAR.EQ.ITAB)LINE(I:I)=SPACE
   40 CONTINUE
*
*   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      DO 50 I=1,10
   50 ISTART(I)=80
*
* FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
*     BY A CHARACTER AND STORE IN ISTART
      LEADSP=.TRUE.
      NVALUE=0
      DO 60 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   60 CONTINUE
*
* ESTABLISH THE ELEMENT'S NAME, CHECK FOR ERRORS OR E.O.DATA
*
      STRING=LINE(ISTART(1):ISTART(2)-1)
      IF( STRING(1:1) .GE. ZERO .AND. STRING(1:1) .LE. NINE) THEN
*  ATOMIC NUMBER USED
         LABEL=READA(STRING,1)
         IF (LABEL.EQ.0) GO TO 130
         IF (LABEL.LT.0.OR.LABEL.GT.107) THEN
            WRITE(6,'(''  ILLEGAL ATOMIC NUMBER'')')
            WRITE (6,'(A)') LINE
            GO TO 240
         END IF
         GO TO 80
      END IF
*  ATOMIC SYMBOL USED
      ELE=STRING(1:2)
*   CHECK FOR ERROR IN ATOMIC SYMBOL
      IF(ELE(1:1).EQ.'-'.AND.ELE(2:2).NE.'-')ELE(2:2)=' '
      DO 70 I=1,107
         IF(ELE.EQ.ELEMNT(I)) THEN
            LABEL=I
            GO TO 80
         END IF
   70 CONTINUE
      WRITE(6,'(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')')ELE
      GOTO 240
*
* ALL O.K.
*
   80 IF(LABEL.NE.99.AND.LABEL.NE.108) THEN
C
C  ENTERING INFORMATION ABOUT THE RADIUS OF THE ATOMS OF THE GIVEN TYPE
C
         NRAD=NRAD+1
         LABRAD(NRAD)=LABEL
         RADRAD(NRAD)=READA(LINE,ISTART(2))
      ELSE
C
C     ENTERING INFORMATION ABOUT THE RADIUS OF THE GIVEN ATOM
C
         IF(ISTART(3).EQ.80) THEN
            WRITE (6,'('' ERROR IN PROVIDING THE RADIUS OF THE '',
     *         ''SELECTED ATOM''/1X,A)') LINE
            STOP
         ENDIF
         NUMADD=NUMADD+1
         NADD(NUMADD)=READA(LINE,ISTART(2))
         RADD(NUMADD)=READA(LINE,ISTART(3))
         IF(NADD(NUMADD).GT.NATOMS) THEN
            WRITE (6,'('' ERROR IN PROVIDING THE NUMBER OF THE ATOM''/
     *         '' NAT='',I4,'' > NATOMS='',I4/1X,A)')
     *         NADD(NUMADD),NATOMS,LINE
         ENDIF
      ENDIF
      GOTO 20
  130 RETURN
  240 STOP
      END SUBROUTINE GETRAD

!======================================================================!
      SUBROUTINE GETSFE(ISFILE,COORD_,NUMAT_)

C     ***************************************************************
C     KEY WORDS FOR SOLVATION CALCULATIONS
C     ***************************************************************
C
C
C       KEY WORDS FOR CONTROLLING OUTPUT
C
C    'COLEND'  THE POLARIZATION PERTURBATION WILL BE PRINTED AT EACH
C              SCF CYCLE BY SUBPROGRAM COLEND
C
C    'CONN'    NUCLEAR ENERGY AND CONTRIBUTION OF CORES IN POLARIZATION
C              ENERGY WILL BE WRITTEN AFTER SFERE HAS BEEN IMPLEMENTED
C              AT EACH SCF CYCLE.
C
C    ' SF1'    WRITING INFORMATION ABOUT SPHERES IN SUBPROGRAM SFERA1
C              ( RV, XX, YY, ZZ ) ON FORMAT 4F10.6
C    ' SF2'    WRITING INFORMATION ABOUT SPHERES IN SUBPROGRAM SFERA2
C              ( N, XX, YY, ZZ, RV, QSFE, N.TESS ) WITH FORMAT
C              I4,5F12.4,I10  .
C    ' SF3'    INFO ABOUT SELF-POLARIZATION ITERATIONS WILL BE WRITTEN
C
C    'TESS1'   WRITING ALL INFORMATION ABOUT ALL TESSERAE IN
C              SUBPROGRAM SFERA1 ( I,XIN,X0,XOUT,YIN,Y0,YOUT,ZIN,
C              Z0,ZOUT,AS,PIMPOU ) ON FORMAT 2X,I4,2X,11F9.4  .
C    'TESS2'   WRITING ALL INFORMATION ABOUT ALL TESSERAE IN
C              SUBPROGRAM SFERA2 ( I,X0,Y0,Z0,AS,PIMPOU,QS)
C              ON FORMAT 4X,I4,9F8.3,2F9.5.
C    'TESS3'   WRITING SOME INFORMATION ABOUT ALL TESSERAE IN
C              SUBPROGRAM SFERA3 ( I,X0,Y0,Z0,QS,AS)
C              ON FORMAT 6X,I3,3F10.3,F11.5,F10.3   .
C     ***************************************************************
C     OPTIONS ...
C     ***************************************************************
C     _SMOOTH - SMOOTHING OF THE SURFACE BY THE INCORPORATION
C               OF ADDITIONAL SPHERES. (THE DEFAULT IN CALCULATIONS
C               WITHOUT GEOMETRY OPTIMIZATION).
C     NOSMOOTH - WITHOUT SMOOTHING OF THE SURFACE BY THE INCORPORATION
C                OF ADDITIONAL SPHERES. (THE DEFAULT IN CALCULATIONS
C                WITH GEOMETRY OPTIMIZATION).
C     SELFCR=N - DEFINITION OF MISFIT IN SELF-POLARIZATION ITERATIONS
C     MERT=N   - THE MIERTUSH METHOD OF TESSERAE PREPARATION TO BE USED
C     ITSE=N   - MAX. NUMBER OF POLARIZATION ITERATIONS IN SFERA3 (BY
C                DEFAULT EQUAL TO 15).
C     ***************************************************************
C     SOLRD=RR - RADIUS OF THE SOLVENT MOLECULE. BY DEFAULT  RR=1.0D0
C     EXVOL=V  - THE SUFFICIENT VOLUME TO AVOID INCORPORATION OF THE
C     ADDITIONAL SFERES BETWEEN NEIBOURING ATOMS. BY DEFAULT  V=1.0D0
C     MODFE=N  - DEGREE OF STRUCTURATION OF SURFACE ELEMENTS. BY
C     DEFAULT MODFE=6.
C     EPS=N      - THE DIELECTRIC PERMITTIVITY OF THE SOLVENT. BY
C     DEFAULT EPS=80.D0
C     EPSEL=N OPTIC PERMITTIVITI. BY DEFAULT EPSEL=2.D0.
C     INPUT DATA
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   LABEL(I) = THE ATOMIC NUMBER OF ATOM\I\.
C            = 99, THEN THE I-TH ATOM IS A DUMMY ATOM USED ONLY TO
C              SIMPLIFY THE DEFINITION OF THE MOLECULAR GEOMETRY.
      
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: ISFILE, NUMAT_
      real(8), intent(inout), dimension(:,:) :: COORD_

      !INCLUDE 'SIZES'
      !CHARACTER*241 KEYWRD, LINE *241,ELEMNT*2
      !LOGICAL SF,MERT
      !COMMON /GEOKST/  NATOMS,LABELS(NUMATM)  !Already initialized
      !COMMON /KEYWRD/  KEYWRD
      !COMMON /ELEMTS/  ELEMNT(108)
      !COMMON /PECAR/   COL2,COL1,COIN,COEL
      !COMMON /EPSS/    EPS,EPSEL       !Already initialized
      !COMMON /SF/      SF
      !COMMON /REGSO/   IREG
      !COMMON /PR/      XX(NSF),YY(NSF),ZZ(NSF),RV(NSF),QSFE(NSF)
      !COMMON /PREL/    RVEL(NSF),QSFEEL(NSF)
      !COMMON /SFE/     NN,KFI,LTH
      !COMMON /MODF  /  ITE1,NCFINR1,ITE2,NCFINR2
      !COMMON /VOL/     CON,VSOLV,RSOLV,SELFCR,CHDIFF,ITSE
      !COMMON /RVDWS/   RVDW(108)
      !COMMON /ADDSF/   RADD(NSF),NADD(NSF),NUMADD
      !COMMON /ADDSFEL/ RADDEL(NSF)
      !COMMON /FACTSF/  FACTOR,FACTOR2  !Already initialized
      !DIMENSION COORD_(3,*)
      !DIMENSION LABRAD(107),RADRAD(107)

      character(241) :: LINE
      logical :: MERT
      integer, dimension(107) :: LABRAD
      real(8),  dimension(107) :: RADRAD
C
C
C     INPUT VdW ATOMIC RADII
C
      I=INDEX(KEYWRD,'RADIUS')
      IF(I.NE.0)THEN
         WRITE(6,'(/''    VdW ATOMIC RADII BEING READ IN'')')
         REWIND ISFILE
         DO I=1,1000
            READ(ISFILE,'(A)')LINE
            IF (INDEX(LINE,'RADIUS').NE.0) GOTO 101
         ENDDO
  101    CONTINUE
         DO I=1,1000
            READ(ISFILE,'(A)',END=175) LINE
            IF (INDEX(LINE,'RADIUS').NE.0) GOTO 180
         ENDDO
  175    WRITE (6,'(/''   KEYWORD "RADIUS" AT THE BEGINNING OF '',
     *      ''THE LIST OF RADII IS ABSENT'')')
         STOP
  180    CALL GETRAD(ISFILE,NRAD,LABRAD,RADRAD)
         DO I=1,NRAD
            RVDW(LABRAD(I))=RADRAD(I)
         ENDDO
      ENDIF
C
      SF=.TRUE.
      MERT=.FALSE.
C
      IREG=2 ! Point charge representation and self-polarization
C
      NN=NUMAT_
      NNSE=INDEX(KEYWRD,'NSNN=')
      IF(NNSE.NE.0) NN=READA(KEYWRD,NNSE)
      IF(NN.GT.NSF) THEN
          WRITE (6,'( 3X,I3, ''  MAXIMUM NUMBER OF SPHERES EXCEEDED  '',
     *    3X,''STOP IN SUBPROGRAM GETSFE'')' ) NSF
          STOP
      END IF
      IF(NN.EQ.0) THEN
          WRITE(6,'(///9X,''NN EQUAL TO ZERO. STOP IN GETSFE'')')
          STOP
      END IF
C
C
      RSOLV=1.0D0
      ITS=INDEX(KEYWRD,'SOLRD=')
      IF(ITS.NE.0) RSOLV=READA(KEYWRD,ITS)
C
      VSOLV=1.0D0
      ITS=INDEX(KEYWRD,'EXVOL=')
      IF(ITS.NE.0) VSOLV=READA(KEYWRD,ITS)
C
      LTH=10
      KFI=20

      ITS=INDEX(KEYWRD,'MERT=')
      IF(ITS.NE.0) THEN
         MERT=.TRUE.
         IE=READA(KEYWRD,ITS)
         LTH=IE
         KFI=IE*2
      END IF

      ITE1=6
      ITE2=9
      NCFINR1=1
      NCFINR2=1
      SELFCR=2.D-4
      ITSE=15
      CHDIFF=1.0D-1
      IF(INDEX(KEYWRD,'PRECISE').NE. 0) THEN
         ITSE=30
         ITE1=9
         ITE2=12
         NCFINR1=2
         NCFINR2=2
         SELFCR=2.D-5
      ENDIF
      ITSELF=INDEX(KEYWRD,'SELFCR=')
      IF(ITSELF.NE.0) SELFCR=READA(KEYWRD,ITSELF)*1.0D-4
      ITS=INDEX(KEYWRD,'ITSE=')
      IF(ITS.NE.0) ITSE=READA(KEYWRD,ITS)
      ITSELF=INDEX(KEYWRD,'CHDIFF=')
      IF(ITSELF.NE.0) CHDIFF=READA(KEYWRD,ITSELF)

      IF(INDEX(KEYWRD,'MODFE=(').NE.0)THEN
         ITE2=READA(KEYWRD,INDEX(KEYWRD,'MODFE=(')+9)
         ITE1=READA(KEYWRD,INDEX(KEYWRD,'MODFE=(')+6)
      ELSEIF (INDEX(KEYWRD,'MODFE=').NE.0)THEN
         ITE1=READA(KEYWRD,INDEX(KEYWRD,'MODFE=')+6)
         ITE2=ITE1
      ENDIF

      IF (ITE1.GE.24) THEN
         NCFINR1=12
      ELSEIF (ITE1.GE.18) THEN
         NCFINR1=6
      ELSEIF (ITE1.GE.12) THEN
         NCFINR1=3
      ELSEIF (ITE1.GE.9) THEN
         NCFINR1=2
      ENDIF

      IF (ITE2.GE.24) THEN
         NCFINR2=12
      ELSEIF (ITE2.GE.18) THEN
         NCFINR2=6
      ELSEIF (ITE2.GE.12) THEN
         NCFINR2=3
C      ELSEIF (ITE1.GE.9) THEN
C         NCFINR2=2
      ENDIF

      IF(INDEX(KEYWRD,'NTETFI=(').NE.0)THEN
         NCFINR2=READA(KEYWRD,INDEX(KEYWRD,'NTETFI=(')+10)
         NCFINR1=READA(KEYWRD,INDEX(KEYWRD,'NTETFI=(')+7)
      ELSEIF (INDEX(KEYWRD,'NTETFI=').NE.0)THEN
         NCFINR1=READA(KEYWRD,INDEX(KEYWRD,'NTETFI=')+7)
         NCFINR2=NCFINR1
      ENDIF

      IF(FACTOR2.EQ.0.0D0) THEN
         IF (ITE2.NE.ITE1.OR.NCFINR2.NE.NCFINR1) THEN
            WRITE(6,'('' Because of DELTA = 0.0 parameters for outer'',
     *      '' sphere (MODFE, NTETFI) corrected to be equal to'',
     *      '' parameters''/''of inner sphere'')')
            ITE2=ITE1
            NCFINR2=ITE1
         ENDIF
      ENDIF
C

      IF ( INDEX(KEYWRD,'NSNN=').NE.0) THEN
          DO 23 I=1,NN
          READ (ISFILE, * ) RV(I),XX(I),YY(I),ZZ(I)
  23      CONTINUE
C***      ALL  IN ANGSTROMS .
      ELSE
         LW=0
         DO 70 I=1,NN
            XX(I)=COORD_(1,I)
            YY(I)=COORD_(2,I)
            ZZ(I)=COORD_(3,I)
   71       LW=LW+1
            IF (LABELS(LW).EQ.99.OR.LABELS(LW).EQ.108) GOTO 71
         RVEL(I)=RVDW(LABELS(LW))*FACTOR
   70    RV(I)=RVEL(I)+FACTOR2
CR         RV(I)=RVDW(LABELS(LW))*FACTOR
CR   70    RVEL(I)=RVDW(LABELS(LW))*FACTOR2
         DO 95 I=1,NN
          IF (RV(I).LT.0.08D0)  THEN
             WRITE (6,'( //6X,''THE MAGNITUDE OF THE RADIUS OF ATOM '',
     *       ''NUMBER '',I3/6X,''IS LESS THAN 0.08 A'',
     *       ''  STOP IN GETSFE (LABEL 70) '' )') I
          END IF
   95    CONTINUE

C-----------------------------------------------------------------------
CAS:     The following block (adding auxiliary spheres)
CAS:     is temporarily unavailable...
CAS:     (needs consulting with I.V.Rostov)
C-----------------------------------------------------------------------
C
C        ADD THE AUXILIARY SPHERES
C
C         IF(NUMADD.NE.0) THEN
C            DO I=1,NUMADD
C               RADDEL(I)=RADD(I)*FACTOR
C               RADD(I)=RADDEL(I)+FACTOR2
C            ENDDO
C            CALL ADDSFET(XX,YY,ZZ,RV,RVEL,NN)
C            write (6,'(/'' IN GETSFE'')')
C            WRITE (6,'('' COORDINATES AND RADII OF SPHERES'')')
C            WRITE (6,'(4X,''I '',6X,''XX'',8X,''YY'',8X,''ZZ'',8X,
C     *         ''RV'',7X,''RVEL'')')
C            DO I=1,NN
C                WRITE(6,8) I,XX(I),YY(I),ZZ(I),RV(I),RVEL(I)
C            ENDDO
C         ENDIF
C-----------------------------------------------------------------------

      END IF

C  120 FORMAT(    4F10.6)

      CON=-(EPS-1.0D00)/EPS
      IF (FACTOR2.EQ.0.0D0) THEN
         COL2=CON
         COL1=0.0D0
         IF (EPSEL.EQ.0.0D0) THEN
            COEL=0.0D0
         ELSE
            COEL=-(EPSEL-1.0D00)/EPSEL
         ENDIF
      ELSE
         IF (EPSEL.EQ.0.0D0) THEN
            COL2=CON
            COL1=0.0D0
            COEL=0.0D0
         ELSE
            COL2=-(EPS-EPSEL)/EPS
            COL1=-(EPSEL-1.0D00)/EPSEL
            COEL=-(EPSEL-1.0D00)/EPSEL
         ENDIF
      ENDIF
      IF (EPSEL.EQ.0.0D0.OR.FACTOR2.EQ.0.0D0) THEN
         COIN=CON
      ELSE
         COIN=-(EPS-EPSEL)/EPS/EPSEL
      ENDIF
      WRITE (6,'(//'' *********     PARAMETERS FOR THE SOLVATION'',
     *   '' CALCULATION     *********''/)')
      WRITE (6,'('' STATIC DIELECTRIC PERMITTIVITY,'',21X,''EPS = '',
     *F6.3)')EPS
      WRITE (6,'('' OPTICAL DIELECTRIC PERMITTIVITY,'',18X,''EPSEL = '',
     *F6.3)')EPSEL
      WRITE (6,'('' PARAMETER FOR INNER CAVITY SPHERES RADII,'',8X,
     *'' KAPPA = '',F6.3)')FACTOR
      WRITE (6,'('' PARAMETER FOR OUTER CAVITY SPHERES RADII,'',8X,
     *'' DELTA = '',F6.3,'' ANGSTROEMS'')')FACTOR2
      WRITE(6,'(/'' POLAR. PROCEDURE CRITERION,'',21X,'' SELFCR = '',
     *D9.2)')SELFCR
      WRITE (6,'('' MAXIMUM NUMBER OF ITERATION IN POLAR.'',
     *'' PROCEDURE,   ITSE = '',I3)')ITSE
      WRITE (6,'('' MAX. ALLOWED DEVIATION OF SURFACE CHARGE''/
     *'' BEFORE COMPENSATION,'',29X,''CHDIFF = '',D9.2)')CHDIFF

      IF (FACTOR2.EQ.0.0D0) THEN
         WRITE(6,'(/''!!!'')')
         WRITE(6,'('' THE CASE WITH DELTA=0.0 CAN NOT BE TREATED IN'',
     *  '' FRCM MODEL''/'' FOR THIS CASE ONE-CAVITY BKO MODEL WILL BE'',
     *  '' APPLIED '')')
         WRITE(6,'(''!!!''/)')
      ENDIF
      IF(.NOT.MERT) THEN
         WRITE (6,'(/'' INNER SPHERE PARAMETERS:'')')
         IF(INDEX(KEYWRD,'NOSMOOTH').EQ.0) THEN
            WRITE (6,'('' EFFECTIVE RADIUS OF THE SOLVENT MOLECULE,'',
     *      9X,''SOLRD = '',F6.3)') RSOLV
            WRITE (6,'('' MIN. UNREACHABLE VOLUME TO INCLUDE EXTRA'',
     *      '' SPHERES, EXVOL = '',F6.3)') VSOLV
         ENDIF
         IF(ITE1.EQ.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE1 = '',
     *     I2,''  (NORMAL)'')') ITE1
         ELSEIF(ITE1.EQ.9) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE1 = '',
     *     I2,''  (PRECISE USUAL)'')') ITE1
         ELSEIF(ITE1.GT.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE1 = '',
     *     I2,''  (PRECISE UNUSUAL)'')') ITE1
         ELSEIF(ITE1.EQ.3) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE1 = '',
     *     I2,''  (DRAFT USUAL)'')') ITE1
         ELSEIF(ITE1.LT.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE1 = '',
     *     I2,''  (DRAFT UNUSUAL)'')') ITE1
         ENDIF
         WRITE(6,'(39x,''PARAMETER NTETFI1 = '',I2)')NCFINR1

         WRITE (6,'(/'' OUTER SPHERE PARAMETERS:'')')
         IF(ITE2.EQ.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE2 = '',
     *      I2,''  (NORMAL)'')') ITE2
         ELSEIF(ITE2.EQ.9) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE2 = '',
     *      I2,''  (PRECISE USUAL)'')') ITE2
         ELSEIF(ITE2.GT.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE2 = '',
     *      I2,''  (PRECISE UNUSUAL)'')') ITE2
         ELSEIF(ITE2.EQ.3) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE2 = '',
     *      I2,''  (DRAFT USUAL)'')') ITE2
         ELSEIF(ITE2.LT.6) THEN
         WRITE (6,'('' PARAMETER OF TESSERA SIZES'',23X,''MODFE2 = '',
     *      I2,''  (DRAFT UNUSUAL)'')') ITE2
         ENDIF
         WRITE(6,'(39x,''PARAMETER NTETFI2 = '',I2)')NCFINR2
      ENDIF
    8 FORMAT(2X,I4,2X, 5F9.4)


      WRITE(6,'(//5X,''CARTESIAN COORDINATES, RADII, AND ATOMIC'',
     *'' CHARGES OF SOLUTE '',/)')
      WRITE(6,'(2X,''NO.'',3X,''ATOM'',6X,''X'',9X,''Y'',9X,''Z'',6X,
     *''RVDW'',5X,''RAD1'',5X,''RAD2''/)')
      L=0
      DO I=1,NUMAT_
         IF (LABELS(I).NE.99.AND.LABELS(I).NE.108) THEN
            L=L+1
            WRITE(6,'(I4,4X,A2,1X,3F10.4,3F9.4)')L,ELEMNT(LABELS(I)),
     *      (COORD_(J,L),J=1,3),RVDW(LABELS(I)),RVEL(I),RV(I)
         ENDIF
      ENDDO
      DO I=L+1,NN
         WRITE(6,'(I4,4X,''XX'',1X,3F10.4,3F9.4)')I,XX(I),YY(I),ZZ(I),
     2RV(I)/FACTOR,RVEL(I),RV(I)
      ENDDO

      RETURN
      END SUBROUTINE GETSFE

!======================================================================!
      SUBROUTINE GETTXT(ISFILE)

      implicit real(8) (a-h,o-z)
      integer, intent(in) :: ISFILE

      !COMMON /KEYWRD/ KEYWRD
      !COMMON /TITLES/ KOMENT,TITLE
      !DIMENSION IS(3)
      !CHARACTER KEYWRD*241, KOMENT*81, TITLE*81, CH*1, CH2*1, FILEN*50

      character (1) :: CH, CH2
      character(50) :: FILEN

      integer, dimension(3) :: IS

      IS(1)=161
      IS(2)=81
      IS(3)=1
      KEYWRD=' '
      KOMENT='    NULL  '
      TITLE ='    NULL  '
      READ(ISFILE,'(A)',END=100,ERR=90)KEYWRD(:80)
      CALL UPCASE(KEYWRD(1:80),80)
      IF(INDEX(KEYWRD,'SETUP').NE.0)THEN
         I=INDEX(KEYWRD,'SETUP=')
         IF(I.NE.0)THEN
            J=INDEX(KEYWRD(I:),' ')
            FILEN=KEYWRD(I+6:I+J-1)
         ELSE
            FILEN='SETUP'
         ENDIF
         OPEN(UNIT=4,FILE=FILEN,STATUS='UNKNOWN',FORM='FORMATTED')
         REWIND 4
         READ(4,'(A)',END=40,ERR=40)KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160),80)
         READ(4,'(A)',END=10,ERR=10)KEYWRD(161:240)
         CALL UPCASE(KEYWRD(161:240),80)
   10    CONTINUE
         READ(ISFILE,'(A)',END=100,ERR=90)KOMENT,TITLE
      ELSEIF(INDEX(KEYWRD(1:80),' +') .NE.0)THEN
C
C  READ SECOND KEYWORD LINE
C
         READ(ISFILE,'(A)',END=100,ERR=90)KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160),80)
         IF(INDEX(KEYWRD(81:160),'SETUP').NE.0)THEN
            I=INDEX(KEYWRD,'SETUP=')
            IF(I.NE.0)THEN
               J=INDEX(KEYWRD(I:),' ')
               FILEN=KEYWRD(I+6:I+J)
            ELSE
               FILEN='SETUP'
            ENDIF
            OPEN(UNIT=4,FILE=FILEN,STATUS='UNKNOWN',FORM='FORMATTED')
            REWIND 4
            READ(4,'(A)',END=20,ERR=20)KEYWRD(161:240)
            CALL UPCASE(KEYWRD(161:240),80)
   20       CONTINUE
         ELSEIF(INDEX(KEYWRD(81:160),' +') .NE.0)THEN
C
C  READ THIRD KEYWORD LINE
C
            READ(ISFILE,'(A)',END=100,ERR=90)KEYWRD(161:240)
            CALL UPCASE(KEYWRD(161:240),80)
         ENDIF
C
C  READ TITLE LINE
C
         READ(ISFILE,'(A)',END=100,ERR=90)KOMENT,TITLE
      ELSEIF(INDEX(KEYWRD(:80),'&').NE.0)THEN
         READ(ISFILE,'(A)',END=100,ERR=90)KEYWRD(81:160)
         CALL UPCASE(KEYWRD(81:160),80)
         IF(INDEX(KEYWRD(81:160),'SETUP').NE.0)THEN
            I=INDEX(KEYWRD,'SETUP=')
            IF(I.NE.0)THEN
               J=INDEX(KEYWRD(I:),' ')
               FILEN=KEYWRD(I+6:I+J)
            ELSE
               FILEN='SETUP'
            ENDIF
            OPEN(UNIT=4,FILE=FILEN,STATUS='UNKNOWN',FORM='FORMATTED')
            REWIND 4
            READ(4,'(A)',END=30,ERR=30)KEYWRD(161:240)
            CALL UPCASE(KEYWRD(161:240),80)
            READ(ISFILE,'(A)',END=100,ERR=90)TITLE
   30       CONTINUE
         ELSEIF(INDEX(KEYWRD(81:160),'&').NE.0)THEN
            READ(ISFILE,'(A)',END=100,ERR=90)KEYWRD(161:240)
         ELSE
            READ(ISFILE,'(A)',END=100,ERR=90)TITLE
         ENDIF
      ELSE
         READ(ISFILE,'(A)',END=100,ERR=90)KOMENT,TITLE
      ENDIF
      GOTO 50
   40 WRITE(6,'(A)')' SETUP FILE MISSING OR CORRUPT'
   50 DO 80 J=1,3
         IF(KEYWRD(IS(J):IS(J)) .NE. ' ') THEN
            CH=KEYWRD(IS(J):IS(J))
            KEYWRD(IS(J):IS(J))=' '
            DO 60 I=IS(J)+1,239
               CH2=KEYWRD(I:I)
               KEYWRD(I:I)=CH
               CH=CH2
               IF(KEYWRD(I+1:I+2) .EQ. '  ') THEN
                  KEYWRD(I+1:I+1)=CH
                  GOTO 70
               ENDIF
   60       CONTINUE
            WRITE(6,'(A,I2,A)')' LINE',J,' OF KEYWORDS DOES NOT HAVE ENO
     1UGH'
            WRITE(6,'(A)')' SPACES FOR PARSING.  PLEASE CORRECT LINE.'
            STOP
   70       CONTINUE
         ENDIF
   80 CONTINUE
      RETURN
   90 WRITE(6,'(A)')' ERROR IN READ OF FIRST THREE LINES'
  100 STOP
      END SUBROUTINE GETTXT

!======================================================================!
      SUBROUTINE MERTU(K,L,NN_,XX_,YY_,ZZ_,RV_,NW_,X0_,Y0_,Z0_,AS_,IJ_)
      IMPLICIT REAL(8) (A-H,O-Z)

      !INCLUDE 'SIZES'

      integer, intent( in) :: K, L, NN_
      integer, intent(out) :: IJ_

      real(8),  intent(in),    dimension(:) :: XX_,YY_,ZZ_,RV_
      real(8),  intent(inout), dimension(:) :: X0_,Y0_,Z0_,AS_
      integer, intent(inout), dimension(:) :: NW_

      !DIMENSION XX_(*),YY_(*),ZZ_(*),RV_(*),NW_(*)
      !DIMENSION X0_(NS),Y0_(NS),Z0_(NS),AS_( NS)

      PI=ATAN(1.0D0)*4.0D0
      FIR=PI/180.D0
      FI0_=360.D0/(2*K)*FIR
      THO=180.D0/(2*L)*FIR
      THD=2.D00*THO
      FID=2.D00*FI0_
      DSEN=SIN(THO)
      IJ_=0
      NW_(1)=1
      DO 1 I=1,NN_
          THSI=0.D0
          TH_=THO
          FIT=0.D0
          RVI=RV_(I)
          XXI=XX_(I)
          YYI=YY_(I)
          ZZI=ZZ_(I)
          CTHS=COS(THSI)
          STHS=SIN(THSI)
          DO 2 J=1,L
              TH_=TH_+THD
              IF(J.EQ.1) TH_=THO
              FI=FIT
              CTH=COS(TH_)
              STH=SIN(TH_)
              ZLS0=RVI*CTH
              RVSI=RVI*STH
              ZLC0=ZLS0*CTHS
              ZLS0=ZLS0*STHS
C ***
C  DEFINITION OF THE SQUARES OF THE TESSERAE .
C ***
              ASPPPP=4.D0*PI*RVI*RVI*DSEN*STH/K
              DO 3 M=1,K
                  FI=FI+FID
                  IJ_=IJ_+1
                  FIC=COS(FI)
                  FAS=SIN(FI)
                  XL00=RVSI*FIC
                  YL0=RVSI*FAS
                  XL0=XL00*CTHS-ZLS0
                  ZL0=XL00*STHS+ZLC0
                  X0_(IJ_)=XL0+XXI
                  Y0_(IJ_)=YL0+YYI
                  Z0_(IJ_)=ZL0+ZZI
                  DO 4 JN=1,NN_
                      IF(JN.EQ.I) GOTO 4
                      RDIF=(X0_(IJ_)-XX_(JN))**2+(Y0_(IJ_)-YY_(JN))
     *                **2+(Z0_(IJ_)-ZZ_(JN))**2
                      RV2=RV_(JN)*RV_(JN)
                      IF(RDIF.GT.RV2)GOTO 4
                      IJ_ = IJ_ - 1
                      GOTO 3
    4             CONTINUE
                  AS_(IJ_)=ASPPPP
    3         CONTINUE
    2     CONTINUE
      NW_(I+1)=IJ_+1
  1   CONTINUE
      IF ( IJ_.GT.NS ) THEN
           NNNNN=NS
           WRITE (6,'( //6X,''THE NUMBER OF SURFACE ELEMENTS IS'',I6,
     *     '' MORE THEN  NS  - '',I6//6X,''- ( THE DIMENSION'',
     *     ''  OF THE COMMON BLOCK )'')')  IJ_, NNNNN
           STOP
      END IF
      RETURN
      END SUBROUTINE MERTU

!======================================================================!
      SUBROUTINE MODFE(ITE,NCFINR)
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(inout) :: ITE
      integer, intent(in) :: NCFINR



      !INCLUDE 'SPHSIZES'
      integer, PARAMETER :: MXRAS1=MXRAS-2

      !COMMON /RESINF/DATTE0,DPRAM,DRASD,DSMIN,ITMA
      !COMMON /TETFI/NTETM,NFIM
      !COMMON /KEYWRD/KEYWRD
      !CHARACTER*241 KEYWRD
      !DIMENSION DATTE0(MXRAS),DPRAM(MXRAS)

      real(8), DIMENSION(MXRAS) :: DDTET
      real(8), DIMENSION(MXRAS) :: DIND
      LOGICAL, save :: PRNT, FIRST=.true.

      !SAVE PRNT,FIRST
      !DATA DIND/6.,9.,MXRAS1*10./,FIRST/.TRUE./

      DIND = 10.d0
      DIND(1:2) = (/6.d0, 9.d0/)

      IF (FIRST) THEN
         FIRST=.FALSE.
         PRNT=(INDEX(KEYWRD,'SECTS') .NE. 0) .AND.
     *      (INDEX(KEYWRD,'DEBUG') .NE. 0)
      ENDIF
      PI=ATAN(1.0D0)*4.0D0
      DO 2 I=1,MXRAS
    2 DPRAM(I)=DIND(I)
      DPRAM(1)=1.33D0*DIND(1)
      DPRAM(2)=1.28D0*DIND(2)
      DO 3 I=1,MXRAS
   3  DPRAM(I)=DPRAM(I)/DIND(I)
      DRASD=1.D0
      DSMIN=0.75D0
      NTETS=60*NCFINR
      NFIS=120*NCFINR
      NTETS=MIN(NTETS,NTQ)
      NFIS=MIN(NFIS,NFQ)
      NTETMS=1
      NFIMS=1
      NTETM=NTETMS*NTETS
      NFIM=NFIMS*NFIS
      IF(ITE.EQ.0) ITE=MXRAS
      COF=40.0D0/3.0D0/ITE*PI/180.D0
      ITMA=0
      SS=0.D0
      DO 1 I=1,MXRAS
      SS=SS+DIND(I)*COF
      IF(SS.GT.PI) GOTO 4
      DATTE0(I)=DIND(I)*COF
      ITMA=ITMA+1
      DDTET(I)=DATTE0(I)*180.D0/PI
   1  CONTINUE
      WRITE (6,'(/'' THE NUMBER OF LEVELS IS MORE THEN SIZE (MODFE)''/
     *   '' INCREASE PARAMETER MXRAS IN SPHSIZES''/)')
      STOP
   4  DATTE0(I)=DIND(I)*COF
      IF(PRNT) THEN
         PRINT 101,NTETM,NFIM,ITE,ITMA,
     *      (DDTET(I),I=1,ITMA)
         PRINT 102,DRASD,DSMIN,(DPRAM(I),I=1,ITMA)
C         PRNT=.FALSE.
      ENDIF
  101 FORMAT (
     *//' NO. OF SECTORS WRT TETA NTETM=',I4/,
     *' NO. OF SECTORS WRT FI NFIM=',I4/
     *' THE DESIRED NO. OF LARGE SECTORS WRT TETA ITE=',I4/
     *' NO. OF LARGE SECTORS WRT TETA ITMA=',I4/
     *' ANGLES WRT TETA IN DEGREES:'/1X,20F6.1)
  102 FORMAT (/
     *' MINIMUM DISTANCE BETWEEN SECTIONS SMALLER THAN THE STANDARD BY',
     *F6.2,' TIMES',
     *' COEFFICIENT OF THE SMALL SECTOR DSMIN=',F6.3/
     *' ARRAY OF THE EXCESS OF THE STEP WRT FI OVER ',
     *' THE STEP WRT TETA DPRAM:'/1X,20F6.2)
      RETURN
      END SUBROUTINE MODFE

!======================================================================!
      LOGICAL FUNCTION MYWORD(KEYWRD,TESTWD)
      CHARACTER(len=*), intent(in)    :: TESTWD
      CHARACTER(len=*), intent(inout) :: KEYWRD

      MYWORD=.FALSE.
   10 J=INDEX(KEYWRD,TESTWD)
      IF(J.NE.0)THEN
   20    IF(KEYWRD(J:J).NE.' ')GOTO 30
         J=J+1
         GOTO 20
   30    MYWORD=.TRUE.
         DO 60 K=J,241
            IF(KEYWRD(K:K).EQ.'='.OR.KEYWRD(K:K).EQ.' ') THEN
C
C     CHECK FOR ATTACHED '=' SIGN
C
               J=K
               IF(KEYWRD(J:J).EQ.'=')GOTO 50
C
C     CHECK FOR SEPARATED '=' SIGN
C
               DO 40 J=K+1,241
                  IF(KEYWRD(J:J).EQ.'=') GOTO 50
   40          IF(KEYWRD(J:J).NE.' ')GOTO 10
C
C    THERE IS NO '=' SIGN ASSOCIATED WITH THIS KEYWORD
C
               GOTO 10
   50          KEYWRD(J:J)=' '
C
C   THERE MUST BE A NUMBER AFTER THE '=' SIGN, SOMEWHERE
C
               GOTO 20
            ENDIF
   60    KEYWRD(K:K)=' '
      ENDIF
      RETURN
      END FUNCTION MYWORD

!======================================================================!
      SUBROUTINE PRMATR(IU,NDIM,IK,A)

C  WRITING OF MATRIX A(IK,IK) INTO UNIT IU

      IMPLICIT REAL(8) (A-H,O-Z)
      integer, intent(in) :: IU, NDIM, IK
      real(8), intent(in), dimension(NDIM,NDIM) :: A

      WRITE (IU,'(2X,8I10)') (I,I=1,IK)
      DO I=1,IK
         WRITE(IU,'(1X,I4,8F10.4)') I,(A(I,J),J=1,IK)
      ENDDO
      RETURN
      END SUBROUTINE PRMATR

!======================================================================!
      SUBROUTINE PRSHOR(LSECT_,NTETM_,NFIM_,NTQ_,NFQ_,N3_)
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: NTETM_, NFIM_, NTQ_, NFQ_, N3_
      integer, intent(in) :: LSECT_(NTQ_,NFQ_,3)

      character(1) :: IBLAN=' ',IVER='!',IGOR='-'
      character(1), dimension(120) :: ISTRK

      character(1), dimension(0:35) :: LIN
      character(36) :: LINN=' 123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      EQUIVALENCE (LINN,LIN(0))

      !CHARACTER *1 LIN,ISTRK(120),IBLAN,IVER,IGOR
      !CHARACTER *36 LINN
      !DIMENSION LIN(0:35)
      !EQUIVALENCE (LINN,LIN(0))
      !DATA LINN/' 123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      !DATA IBLAN/' '/,IVER/'!'/,IGOR/'-'/

      WRITE (6,597) ((I/10),I=1,60,10)
      WRITE (6,598) ((I-(I/10)*10),I=1,60)
  597 FORMAT (/6X,10(1X,I1,18X))
  598 FORMAT (6X,60I2/)
      DO 539 I=1,120
  539 ISTRK(I)=IBLAN
      DO 529 I=1,NTETM_
      DO 528 J=1,3
      IFF=120
      IFT=1
      IF(LSECT_(I,IFT,J).NE.LSECT_(I,IFF,J)) GOTO 527
  528 CONTINUE
      GOTO 529
  527 ISTRK(2*I)=IGOR
      ISTRK(2*I-1)=IGOR
      IF(I.GT.1) ISTRK(2*I-2)=IGOR
  529 CONTINUE
      WRITE (6,595) (ISTRK(I),I=1,120)
      DO 560 IFT=1,NFIM_
      DO 550,I=1,120
  550 ISTRK(I)=IBLAN
      LKT=0
      IS=1
      DO 570 I=1,NTETM_-1
      DO 571 J=1,3
      IF(LSECT_(I,IFT,J).NE.LSECT_(I+1,IFT,J)) GOTO 572
  571 CONTINUE
      LKT=LKT+1
      GOTO 570
  572 ISTRK(2*I)=IVER
      IF(LKT.EQ.0) GOTO 569
      ISTRK(IS)=LIN(LSECT_(I,IFT,1))
      IF(LKT.EQ.1) GOTO 569
      ISTRK(IS+1)=LIN(LSECT_(I,IFT,2))
      ISTRK(IS+2)=LIN(LSECT_(I,IFT,3))
  569 LKT=0
      IS=2*I+2
  570 CONTINUE
      IF(LKT.EQ.0) GOTO 568
      ISTRK(IS)=LIN(LSECT_(I,IFT,1))
      IF(LKT.EQ.1) GOTO 568
      ISTRK(IS+1)=LIN(LSECT_(I,IFT,2))
      ISTRK(IS+2)=LIN(LSECT_(I,IFT,3))
  568 CONTINUE
      WRITE (6,590) IFT,(ISTRK(I),I=1,120),IFT
  590 FORMAT (1X,I3,2X,120A1,2X,I3)
      DO 540 I=1,120
  540 ISTRK(I)=IBLAN
      DO 530 I=1,NTETM_
      DO 531 J=1,3
      IFF=IFT+1
      IF(IFT.EQ.NFIM_) IFF=1
      IF(LSECT_(I,IFT,J).NE.LSECT_(I,IFF,J)) GOTO 532
  531 CONTINUE
      GOTO 530
  532 ISTRK(2*I)=IGOR
      ISTRK(2*I-1)=IGOR
      IF(I.GT.1) ISTRK(2*I-2)=IGOR
  530 CONTINUE
      WRITE (6,595) (ISTRK(I),I=1,120)
  595 FORMAT (6X,120A1,2X)
  560 CONTINUE
      WRITE (6,597) ((I/10),I=1,60,10)
      WRITE (6,598) (I,I=1,60)
      RETURN
      END SUBROUTINE PRSHOR

!======================================================================!
      SUBROUTINE SECEXL(SSECT,KSECT,ISCFI,ITMX,IFMX,
     *LSECTK,LSECTI,LSECTJ,DFSTI,RAT_,NNEB,SPL)
***********************************************************************
*
*     THIS SUBPROGRAM ELIMINATES THE SMALLEST TESSERAE
*
***********************************************************************
*
      IMPLICIT REAL(8) (A-H,O-Z)
      !INCLUDE 'SPHSIZES'
C
C     NSRM - NUMBER OF ELEMENTS OF THE SMALL SECTOR
C            WHICH ARE TO BE DISTRIBUTED
C
      integer, PARAMETER :: NSRM=NTQ*NFQ/10
C
C     PARAMETER GRTN INCREASES MEMORY FOR THE ARRAYS NSET,NSEF
C     INCREASING THIS PARAMETER SPEEDS UP THE REMOVAL OF SMALL SECTORS
C
      real(8), PARAMETER :: GRTN=1.1D0
      integer, PARAMETER :: NTFQ1=NTFQ*GRTN

      real(8),    intent(inout) :: SSECT(NBGEN)
      integer,   intent(inout) :: KSECT(NBGEN)
      integer,   intent(inout) :: ISCFI(NB,NT)
      integer,   intent(inout) :: ITMX(NB), IFMX(NB,NT)
      integer,   intent(inout) :: LSECTK(NTQ,NFQ),LSECTI(NTQ,NFQ),
     ,                            LSECTJ(NTQ,NFQ)
      real(8),    intent(inout) :: DFSTI
      real(8),    intent(inout) :: RAT_
      integer,   intent(inout) :: NNEB
      real(8),    intent(inout) :: SPL

      ! ARRAYS FOR THE ELEMENTS OF A SMALL SECTOR,
      ! WHICH HAVE TO BE DISTRIBUTED
      integer, dimension(NSRM) :: MZK0,MZI0,MZJ0,MZTET,MZFI

      LOGICAL, save :: PRINT, FIRST=.true.
      logical :: flag
C
      !real(8) :: SITETM(NTQ), COTETM(NTQ)
      !real(8) :: FIM(NFQ),SIFIM(NFQ),COFIM(NFQ)
C
C     ARRAYS FOR THE DISTRIBUTION OF SMALL SECTORS WITHIN SECTORS
C
      integer :: LS(NBS)
      integer :: LSNM(NBGEN), LLS(NBS), LIS(NBS)
      real(8)  :: SECDL(NBS), SEDDS(MXRAS)
C
C     ARRAYS FOR THE INTERMEDIATE STORAGE OF INFORMATION ABOUT SECTORS
C     WHICH ARE TOO LARGE
C
      integer :: KBSIJT(NSRM), KBSDDM(NSRM)
C
C     ARRAYS FOR DEBUGGING PRINTING INFORMATION
C
      integer :: KOTL(NBGEN)
      real(8)  :: SOTL(NBGEN)
      integer, dimension(48) :: INNBIT=(/
     , 1,0,1,0,1,-1,-1,-1,2,0,2,1,2,0,2,1,2,-2,-2,-1,-2,-2,
     ,-1,-2,3,0,3,1,3,2,3,0,3,1,3,2,3,-3,-3,-1,-3,-2,-3,-3,-1,-3,-2,-3/)
      integer, dimension(48) :: INNBIF=(/
     , 0,1,1,-1,-1,0,1,-1,0,2,1,2,2,-2,-1,-2,-2,0,1,2,2,-1,
     ,-2,-2,0,3,1,3,2,3,3,-3,-1,-3,-2,-3,-3,0,1,3,2,3,3,-1,-3,-2,-3,-3/)
      integer, dimension( 4) :: INR=(/0,8,24,48/)


      !DIMENSION DATTE0(MXRAS),DPRAM(MXRAS)
      !COMMON /RESINF/DATTE0,DPRAM,DRASD,DSMIN,ITMA
      !COMMON /TETFI/NTETM,NFIM
      !COMMON /SICO/COTETM,SITETM,COFIM,SIFIM,FIM
      !COMMON /NSTSF/NSET,NSEF
      !COMMON /KEYWRD/ KEYWRD
      !SAVE PRINT

      !DATA INNBIT/1,0,1,0,1,-1,-1,-1,2,0,2,1,2,0,2,1,2,-2,-2,-1,-2,-2,
      !-1,-2,3,0,3,1,3,2,3,0,3,1,3,2,3,-3,-3,-1,-3,-2,-3,-3,-1,-3,-2,-3/
      !DATA INNBIF/0,1,1,-1,-1,0,1,-1,0,2,1,2,2,-2,-1,-2,-2,0,1,2,2,-1,
      !-2,-2,0,3,1,3,2,3,3,-3,-1,-3,-2,-3,-3,0,1,3,2,3,3,-1,-3,-2,-3,-3/
      !DATA INR/0,8,24,48/
      !DATA FIRST/.TRUE./
C
      IF (FIRST) THEN
         FIRST=.FALSE.
         PRINT=(INDEX(KEYWRD,'SECEXL').NE.0)
      ENDIF
      DSMAX=1.2D0
      DFMAXX=1.35D0
      DSAM=DFSTI*RAT_*RAT_*DSMIN
      DO 5 I=1,20
    5 SEDDS(I)=DATTE0(I)*DATTE0(I)*DPRAM(I)*DSAM
C
C     DEBUGGING RECORD OF THE STRUCTURE OF SECTORS
C     UP UNTIL THE EXCLUSION OF SMALL SECTORS
C
      IF(PRINT) THEN
         K1=0
         K2=0
         KK=0
         SDDM=1.D10
         DO 10 K=1,NNEB
         DO 10 I=1,ITMX(K)
         DDS=DATTE0(I)*DATTE0(I)*DPRAM(I)*DSAM
         ISFDL=ISCFI(K,I)
         DO 10 J=1,IFMX(K,I)
            IF(KSECT(ISFDL+J).LE.0) GOTO 10
            KOTL(ISFDL+J)=KSECT(ISFDL+J)
            SOTL(ISFDL+J)=SSECT(ISFDL+J)
            KK=KK+1
            IF(SSECT(ISFDL+J).GT.DDS) THEN
               WRITE (6,'('' KK='',I4,
     *            '' LARGE SECTOR  K='',I3,'' I='',I3,
     *            '' J='',I3,'' K1='',I3,'' N='',
     *            I4,'' S='',F7.4,'' S/S0='',F7.4)')
     *            KK,K,I,J,K1,KSECT(ISFDL+J),
     *            SSECT(ISFDL+J),SSECT(ISFDL+J)/DDS
            ELSE
               K2=K2+1
               WRITE (6,'('' KK='',I4,
     *            '' SMALL SECTOR  K='',I3,'' I='',I3,
     *            '' J='',I3,'' K2='',I3,'' N='',
     *            I4,'' S='',F7.4,'' S/S0='',F7.4)')
     *            KK,K,I,J,K2,KSECT(ISFDL+J),
     *            SSECT(ISFDL+J),SSECT(ISFDL+J)/DDS
            ENDIF
  10     CONTINUE
         K1=KK-K2
      ENDIF
      IPERS=1
C
C     INBR - LOOP COUNTER FOR DEBUGGING
C
      INBR=0
      NBGSE=0
   20 CONTINUE
      INBR=INBR+1
      IF(IPERS.EQ.1) THEN
         KDDM=0
         CALL SECSOR(KSECT,SSECT,ISCFI,ITMX,IFMX,
     *     LSECTK,LSECTI,LSECTJ,LS,LSNM,LLS,LIS,NSET,NSEF,
     *     SECDL,SEDDS,NNEB,KK)
         IF(INBR.EQ.1.AND.KK.EQ.0) RETURN
         IPERS=0
      ENDIF
C
C     IF THE TRANFER HAS BEEN COMPLETED OR THE SECTOR BEING CONSIDERED
C     HAS BEEN EXCLUDED
C

C(AVS)
      flag = .false.
      if (kddm.eq.0) then
         flag=.true.
      elseif (secdl(kddm).le.1.d-6) then
         flag = .true.
      else
         flag = .false.
      endif
C(AVS)


C     IF(KDDM.EQ.0.OR.SECDL(KDDM).LE.1.D-6) THEN
      if (flag) then
C
C        SEPARATING OUT THE SMALLEST SECTOR WHICH HAS TO BE REMOVED
C
C        EXCLUDING SECTORS WHICH ARE TOO LARGE
C
         DO 85 I=1,NBGSE
            KSECT(KBSIJT(I))=-ABS(KSECT(KBSIJT(I)))
            SECDL(KBSDDM(I))=-ABS(SECDL(KBSDDM(I)))
   85    CONTINUE
         NBGSE=0
         SDDM=1.D10
         DO 30 K=1,KK
            IF(SECDL(K).GT.1.D-6) THEN
               IF(SECDL(K).LT.SDDM) THEN
                  SDDM=SECDL(K)
                  KDDM=K
               ENDIF
            ENDIF
   30    CONTINUE
C
C        IF THERE ARE NO SMALL SECTORS
C
         IF(SDDM.GE.1) GOTO 80
      ENDIF
      INLL=0
C
C     DEFINITION OF COORDINATES K,I,J OF THE SECTOR BEING EXCLUDED
C
      ITETM=NSET(LLS(KDDM)+1)
C
C     IF ITETM=0 THEN THE GIVEN SECTOR IS ALREADY BEING EXCLUDED AND
C     ITS PREVIOUS COORDINATES CAN BE KEPT
C
      IF(ITETM.NE.0) THEN
         IFIM=NSEF(LLS(KDDM)+1)
         KNEB=LSECTK(ITETM,IFIM)
         INEB=LSECTI(ITETM,IFIM)
         JNEB=LSECTJ(ITETM,IFIM)
      ENDIF
      I=INEB
      J=JNEB
      KIJT=ISCFI(KNEB,I)+J
      DDS=SEDDS(I)
      DO 50 IPOL=LLS(KDDM)+1,LLS(KDDM)+LIS(KDDM)
         ITETM=NSET(IPOL)
C
C        IF THIS SECTOR HAS ALREADY BEEN EXCLUDED
C
         IF(ITETM.EQ.0) GOTO 50
         IFIM=NSEF(IPOL)
C
C        LOOKING FOR A SUITABLE SECTOR
C
C        LOOP OVER DIVERGING CIRCLES
C
         IUKRA=0
         SDKRA=1.D10
            DO 60 IROU=1,3
               INR1=INR(IROU)+1
               INR2=INR(IROU+1)
               DO 70 IRNB=INR1,INR2
                  ITPOP=INNBIT(IRNB)+ITETM
                  IFPOP=INNBIF(IRNB)+IFIM
C
C                 MAKING USE OF CYCLIC PERMUTATIONS
C
                  IF(ITPOP.GT.NTETM) THEN
                     ITPOP=2*NTETM+1-ITPOP
                     IFPOP=IFPOP+NFIM/2
                  ENDIF
                  IF(ITPOP.LE.0    ) THEN
                     ITPOP=1-ITPOP
                     IFPOP=IFPOP+NFIM/2
                   ENDIF
                   IF(IFPOP.GT.NFIM) IFPOP=IFPOP-NFIM
                   IF(IFPOP.LE.0   ) IFPOP=IFPOP+NFIM
                   KNE0=LSECTK(ITPOP,IFPOP)
                   IF(KNE0.EQ.0) GOTO 70
                   I0=LSECTI(ITPOP,IFPOP)
                   J0=LSECTJ(ITPOP,IFPOP)
                   ISFD0=ISCFI(KNE0,I0)
C
C                  STILL THE SAME SECTOR
C
                   IF(LSNM(ISFD0+J0).EQ.KDDM) GOTO 70
                   IF(KSECT(ISFD0+J0).LE.0) GOTO 70
                   SECDDS=ABS(SSECT(ISFD0+J0))/SEDDS(I0)
                   IF(SECDDS.GT.SDKRA) GOTO 70
                   IUKRA=1
                   I00=I0
                   J00=J0
                   K00=KNE0
                   SDKRA=SECDDS
   70           CONTINUE
                IF(IUKRA.EQ.1) THEN
C
C                  A SUITABLE SECTOR HAS BEEN FOUND
C
                   INLL=INLL+1
                   IF(INLL.GT.NSRM) THEN
                      WRITE (6,'(''    THE NUMBER OF SMALL SECTORS '',
     *                   ''IN SECEXL EXCEEDS ARRAY SIZE ''/
     *                   '' INLL='',I4,'' > NSRM='',I4/
     *                   '' INCREASE PARAMETER NSRM (IN SECEXL)'')') 
     *                   INLL,NSRM
                   STOP
                ENDIF
                MZK0(INLL)=K00
                MZI0(INLL)=I00
                MZJ0(INLL)=J00
                MZTET(INLL)=ITETM
                MZFI(INLL)=IFIM
C
C               REMOVING INFORMATION ABOUT THE EXCLUDED SECTOR
C
                NSET(IPOL)=0
                SECDL(KDDM)=SECDL(KDDM)-SPL*SITETM(ITETM)/DDS
                GOTO 50
             ENDIF
   60     CONTINUE
   50 CONTINUE
C
C     IF NOTHING HAS CHANGED, THEN THE GIVEN SECTOR
C     WILL NOT BE ANALYSED IN THE FUTURE
C
      IF(INLL.EQ.0) THEN
         KSECT(KIJT)=-ABS(KSECT(KIJT))
         SECDL(KDDM)=0.D0
         GOTO 20
      ENDIF
C
C     THE SMALL SECTOR CAN BE (PARTIALLY) DISTRIBUTED
C
      DO 90 INP=1,INLL
         ISFDL=ISCFI(KNEB,I)
         I0=MZI0(INP)
         ISFD0=ISCFI(MZK0(INP),I0)
         J0=MZJ0(INP)
         ITETM=MZTET(INP)
         IFIM=MZFI(INP)
         S0=SPL*SITETM(ITETM)
         SSECT(ISFD0+J0)=SSECT(ISFD0+J0)+SIGN(S0,SSECT(ISFD0+J0))
         SSECT(ISFDL+J)=SSECT(ISFDL+J)-SIGN(S0,SSECT(ISFDL+J))
         KSECT(ISFD0+J0)=KSECT(ISFD0+J0)+1
         KSECT(ISFDL+J)=KSECT(ISFDL+J)-1
         LSECTK(ITETM,IFIM)=MZK0(INP)
         LSECTI(ITETM,IFIM)=I0
         LSECTJ(ITETM,IFIM)=J0
C
C        DID THE SECTOR TURN OUT TO BE TOO LARGE?
C
         IF(ABS(SSECT(ISFD0+J0)/SEDDS(I0)).GT.DFMAXX/DSMIN) THEN
            SSECT(ISFD0+J0)=-ABS(SSECT(ISFD0+J0))
            KSECT(ISFD0+J0)=-ABS(KSECT(ISFD0+J0))
            SECDL(LSNM(ISFD0+J0))=-ABS(SECDL(LSNM(ISFD0+J0)))
         ENDIF
         IF(ABS(SSECT(ISFD0+J0))/SEDDS(I0).GT.DSMAX/DSMIN.
     *    AND.SSECT(ISFD0+J0).GT.0.D0) THEN
            NBGSE=NBGSE+1
            KBSIJT(NBGSE)=ISFD0+J0
            KBSDDM(NBGSE)=LSNM(ISFD0+J0)
            SSECT(ISFD0+J0)=-ABS(SSECT(ISFD0+J0))
         ENDIF
C
C        IS THERE OVERFLOW?
C
         IF(IPERS.EQ.0) THEN
            KN=LSNM(ISFD0+J0)
C
C           WILL THERE BE OVERFLOW
C
            IF(LLS(KN)+LIS(KN).GE.LLS(KN+1)) THEN
               IPERS=1
               GOTO 90
            ENDIF
            LIS(KN)=LIS(KN)+1
            LI=LLS(KN)+LIS(KN)
            NSET(LI)=ITETM
            NSEF(LI)=IFIM
            SECDL(KN)=SECDL(KN)+SIGN(S0/SEDDS(I0),SECDL(KN))
         ENDIF
   90 CONTINUE
C
C     COUNT OF LOOPS FOR DEBUGGING
C
      IF(INBR.GT.10000) GOTO 95
      GOTO 20
   95 WRITE (6,'('' THIS MASSAGE SHOULD NOT APPEAR (SECEXL)''/
     *   '' CONSULT PROGRAMMER''/)')
   80 CONTINUE
C
C     CHECK HOW MANY OF THE SMALL AND MAXIMALLY LARGE SECTORS REMAIN
C
      K1=0
      K2=0
      DO 25 K=1,NNEB
         DO 25 I=1,ITMX(K)
            ISFDL=ISCFI(K,I)
            DO 25 J=1,IFMX(K,I)
               SSECT(ISFDL+J)=ABS(SSECT(ISFDL+J))
               IF(KSECT(ISFDL+J).LT.0) THEN
                  IF(SSECT(ISFDL+J)/SEDDS(I).GT.1) THEN
                     K1=K1+1
                  ELSE
                     K2=K2+1
                  ENDIF
                  KSECT(ISFDL+J)=-KSECT(ISFDL+J)
               ENDIF
   25 CONTINUE
C
C     DEBUGGING RECORD OF THE STRUCTURE OF SECTORS AFTER THE EXCLUSION
C     OF SMALL SECTORS
C
      IF(PRINT) THEN
         K1=0
         K2=0
         KK=0
         SXA=0.D0
         SXB=0.D0
         DO 15 K=1,NNEB
         DO 15 I=1,ITMX(K)
         DDS=DATTE0(I)*DATTE0(I)*DPRAM(I)*DSAM
         ISFDL=ISCFI(K,I)
         DO 15 J=1,IFMX(K,I)
            IF(KSECT(ISFDL+J).LE.0) GOTO 15
            KK=KK+1
            IF(SSECT(ISFDL+J).GT.DDS) THEN
               WRITE (6,'('' KK='',I4,
     *            '' LARGE SECTOR K='',I3,'' I='',I3,
     *            '' J='',I3,'' K1='',I3,'' N='',
     *            I4,'' S='',F7.4,'' S/S0='',F7.4,
     *            '' INCREMENT N='',I3,'' DS='',F7.4,'' DS/S='',F7.4)')
     *            KK,K,I,J,K1,KSECT(ISFDL+J),
     *            SSECT(ISFDL+J),SSECT(ISFDL+J)/DDS,
     *            KSECT(ISFDL+J)-KOTL(ISFDL+J),
     *            SSECT(ISFDL+J)-SOTL(ISFDL+J),
     *            (SSECT(ISFDL+J)-SOTL(ISFDL+J))/SOTL(ISFDL+J)
               SASS=SSECT(ISFDL+J)/DDS
               SXA=SXA+SASS
               SXB=SXB+SASS*SASS
            ELSE
               K2=K2+1
               WRITE (6,'('' KK='',I4,
     *            '' SMALL SECTOR K='',I3,'' I='',I3,
     *            '' J='',I3,'' K2='',I3,'' N='',
     *            I4,'' S='',F7.4,'' S/S0='',F7.4,
     *            '' INCREMENT N='',I3,'' DS='',F7.4,'' DS/S='',F7.4)')
     *            KK,K,I,J,K1,KSECT(ISFDL+J),
     *            SSECT(ISFDL+J),SSECT(ISFDL+J)/DDS,
     *            KSECT(ISFDL+J)-KOTL(ISFDL+J),
     *            SSECT(ISFDL+J)-SOTL(ISFDL+J),
     *            (SSECT(ISFDL+J)-SOTL(ISFDL+J))/SOTL(ISFDL+J)
            ENDIF
  15     CONTINUE
         K1=KK-K2
         IF(K1.GT.0) THEN
            SSRED=SXA/K1
            SDISP1=SDISP1+SXA
            SDISP2=SDISP2+SXB
            NDISP=NDISP+K1
            IF(K1.GT.1) THEN
               SSDIS=SQRT((SXB-SXA*SXA/K1)/(K1-1))
            ENDIF
         ENDIF
      ENDIF
      STIM=SECOND()-STIME
      STIME=SECOND()
      STIM0=STIM0+STIM
      STIM3=STIM3+STIM
      RETURN
      END SUBROUTINE SECEXL

!======================================================================!
      SUBROUTINE SECINS(SECTTE,SECTFI,ISCFI,ITMX,IFMX,
     *     NNEB,LSECTK,LSECTI,LSECTJ,KSECT,SSECT,NMSBR,LNUSEC_,SPL)
***********************************************************************
*
*     THIS SUBPROGRAM CALCULATES TESSERAE AS A COMBINATION
*     OF SMALL SECTORS
*
***********************************************************************
*
      IMPLICIT REAL(8) (A-H,O-Z)

      ! input
      real(8),  intent(in) :: SECTTE(NB,NT),SECTFI(NBGEN)
      integer, intent(in) :: ISCFI(NB,NT)
      integer, intent(in) :: ITMX(NB),IFMX(NB,NT)
      integer, intent(in) :: NNEB
      integer, intent(in) :: NMSBR(NB)
      real(8),  intent(in) :: SPL

      ! output
      integer, intent(inout) :: LSECTK(NTQ,NFQ),LSECTI(NTQ,NFQ),
     ,                          LSECTJ(NTQ,NFQ)
      integer, intent(inout) :: KSECT(NBGEN)
      real(8),  intent(inout) :: SSECT(NBGEN)
      integer, intent(inout) :: LNUSEC_(NTQ,NFQ)

      ! local variables
      real(8) :: COTALL(NB)

      !INCLUDE 'SPHSIZES'
      !integer LSECTK,LSECTI,LSECTJ
      !integer LNUSEC_(NTQ,NFQ)
      !COMMON /TETFI/NTETM,NFIM
      !COMMON /SICO/COTETM,SITETM,COFIM,SIFIM,FIM
      !COMMON /RESCON/SITET0,COTET0,FI0,SIFI0,COFI0,
      !               COTECR,SITECR,TETECR,RRECR
      !COMMON /SQRMS/NSQRD,NSQDUM,SQRMAS(1)

      !DIMENSION SITETM(NTQ),COTETM(NTQ)
      !DIMENSION FIM(NFQ),SIFIM(NFQ),COFIM(NFQ)
      !DIMENSION SITET0(NB),COTET0(NB)
      !DIMENSION FI0(NB),SIFI0(NB),COFI0(NB)
      !DIMENSION COTECR(NB),SITECR(NB),TETECR(NB),RRECR(NB)

      !DIMENSION SECTTE(NB,NT),ISCFI(NB,NT),SECTFI(NBGEN)
      !DIMENSION ITMX(NB),IFMX(NB,NT)
      !DIMENSION LSECTK(NTQ,NFQ),LSECTI(NTQ,NFQ),LSECTJ(NTQ,NFQ)
      !DIMENSION NMSBR(NB)
      !DIMENSION KSECT(NBGEN),SSECT(NBGEN)
      !DIMENSION COTALL(NB)

      PI = DATAN(1.0D0)*4.0D0
      STIME=SECOND()
      STIM0=0.D0
C
C          LOOP OVER TETA
C
           DO 50 ITETM=1,NTETM
                COT=COTETM(ITETM)
                SIT=SITETM(ITETM)
C
C               LOOP OVER FI
C
                DO 51 IFIM=1,NFIM
                     COF=COFIM(IFIM)
                     SIF=SIFIM(IFIM)
C
C                    LOOP OVER NEIGHBOURING ATOMS TO CHECK
C                    IF THIS SECTOR IS ON THE SURFACE
C
                     DO 56 INEB=1,NNEB
                          COT0=COTET0(INEB)
                          SIT0=SITET0(INEB)
                          COF0=COFI0(INEB)
                          SIF0=SIFI0(INEB)
                          RR_=RRECR(INEB)
C
C                         CALCULATION OF COS(ALPHA)
C
                          COTAL=SIT*SIT0*(COF*COF0+SIF*SIF0)+COT*COT0
                          IF(COTAL.GT.COTECR(INEB)) THEN
                               LNUSEC_(ITETM,IFIM)=-NMSBR(INEB)
                               GOTO 51
                          ENDIF
                          COTALL(INEB)=COTAL
   56                CONTINUE
C
C                    LOOP OVER NEIGHBOURING ATOMS
C
                     GMX0=-1.D10
                     DO 55 INEB=1,NNEB
                          COTAL=COTALL(INEB)
C
C                         DEFINITION OF THE NEAREST NEIGHBOUR
C
                          ITAL=ABS(COTAL)*NSQRD+1
                          SITAL=SQRMAS(ITAL)
                          GMXT=COTAL*COTECR(INEB)+SITAL*SITECR(INEB)
                          IF(GMXT.LT.GMX0) GOTO 55
                          GMX0=GMXT
                          SIMX=SITAL
                          COMX=COTAL
                          KNEB=INEB
   55                 CONTINUE
                      COT0=COTET0(KNEB)
                      SIT0=SITET0(KNEB)
                      COF0=COFI0(KNEB)
                      SIF0=SIFI0(KNEB)
                      ITMN=ITMX(KNEB)-1
                      DO 60 I=1,ITMN
                           IF(COMX.GT.SECTTE(KNEB,I)) GOTO 61
  60                  CONTINUE
                      I=ITMN+1
                      J=1
                      ISFDL=ISCFI(KNEB,I)
                      GOTO 63
  61                  CONTINUE
                      IFMN=IFMX(KNEB,I)
                      IF(DABS(COT0).LT.1.D0-1.D-5) GOTO 64
                      CFF=COF
                      IF(FIM(IFIM).GT.PI) CFF=-CFF-2.D0
                      GOTO 66
   64                 STT=SIMX
                      CFF=(COT*SIT0-SIT*COT0*(COF*COF0+SIF*SIF0))/STT
                      FIF=FIM(IFIM)-FI0(KNEB)
                      IF((FIF.GE.0.D0.AND.FIF.LE.PI.AND.FI0(KNEB).LE.PI)
     *                .OR.
     *                ((FIF.GE.0.D0.OR.FIF.LE.-PI).AND.FI0(KNEB).GT.PI))
     *                CFF=-CFF-2.D0
   66                 CONTINUE
                      ISFDL=ISCFI(KNEB,I)
                      DO 62 J=1,IFMN
                         ISFDLJ=ISFDL+J
                         IF (ISFDLJ.GE.NBGEN) THEN
                            WRITE(6,101) ISFDLJ,NBGEN
                            WRITE(6,102)
                            STOP
                         ENDIF
                         IF(CFF.GT.SECTFI(ISFDL+J)) GOTO 63
   62                 CONTINUE
                      J=1
   63                 CONTINUE
                      LSECTK(ITETM,IFIM)=KNEB
                      LSECTI(ITETM,IFIM)=I
                      LSECTJ(ITETM,IFIM)=J
                      KSECT(ISFDL+J)=KSECT(ISFDL+J)+1
                      S=SPL*SIT
                      SSECT(ISFDL+J)=SSECT(ISFDL+J)+S
   51      CONTINUE
   50      CONTINUE
      STIM=SECOND()-STIME
      STIME=SECOND()
      STIM0=STIM0+STIM
      RETURN
  101 FORMAT(' SECINS: ISFDLJ = ',I8/
     *'  THIS IS LARGER THAN PARAMETER NBGEN =',I8)
  102 FORMAT(' SUCH SITUATION CAN NOT BE TREATED BY PROGRAM. IT IS',
     *' STOPPED'/
     *' SUGGESTION: INCREASE PARAMETER NBGEN, IT SUPPOSED TO BE LARGER'/
     *' THAT ANY POSSIBLE VALUE ISFDLJ')
      END SUBROUTINE SECINS

!======================================================================!
      SUBROUTINE SECNEB(XA,YA,ZA,RA,COFACT,NAT,NATNEW)
***********************************************************************
*
*     THIS SUBPROGRAM SMOOTHES THE SURFACE OF CAVITY
*     BY INCLUDING ADDITIONAL SHPERES
*
*     ON INPUT:
*
*     XA,YA,ZA - COORDINATES OF SPHERES
*     RA - RADII OF SPHERES
*     NAT - THE NUMBER OF SPHERES
*
*     ON OUTPUT
*
*     XA,YA,ZA - COORDINATES OF SPHERES
*     RA - RADII OF SPHERES
*     NATNEW - THE NEW NUMBER OF SPHERES
*
***********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: NAT
      real(8),  intent(inout), dimension(NSF) :: XA, YA, ZA, RA
      real(8),  intent(inout), dimension(NSF) :: COFACT
      integer, intent(out) :: NATNEW

      real(8),  dimension(NSF) :: TET0,SITET0_,COTET0_,FI0_,
     ,                           SIFI0_,COFI0_,RDIS,RUP
      integer, dimension(NSF) :: NMNI, INMB
      
      !INCLUDE 'SIZES'
      !INCLUDE 'SPHSIZES'
      !COMMON /VOL/   CON,VSOLV,RSOLV,SELFCR,CHDIFF,ITSE
      !COMMON /SCHCON/ COSBT,RAT,COSBT1,RNE,RR,RDS,RDS1
      !DIMENSION XA(NAT),YA(NAT),ZA(NAT),RA(NAT),COFACT(NAT)
      !DIMENSION TET0(NSF),SITET0_(NSF),COTET0_(NSF)
      !DIMENSION FI0_(NSF),SIFI0_(NSF),COFI0_(NSF)
      !DIMENSION RDIS(NSF),NMNI(NSF),RUP(NSF),INMB(NSF)
      !EXTERNAL VOLMIN

C
C     Cycle through atoms
C
      PI=ATAN(1.0D0)*4.0D0
      NATNEW=NAT
      NSFR=NAT

  50  NSFR0=NSFR

      DO IAT=1,NATNEW

         XAT=XA(IAT)
         YAT=YA(IAT)
         ZAT=ZA(IAT)
         RAT=RA(IAT)
C
C        Neighbour spheres definition
C
         NNEI=0
         NAT1=NSFR

         DO INEB=1,NAT1

            IF(INEB.EQ.IAT) GOTO 2
            XNE=XA(INEB)
            YNE=YA(INEB)
            ZNE=ZA(INEB)
            RNE=RA(INEB)
            XR=XNE-XAT
            YR=YNE-YAT
            ZR=ZNE-ZAT
            RR=DSQRT(XR*XR+YR*YR+ZR*ZR)
            NNEI=NNEI+1
            RDIS(NNEI)=RR
            RUP(NNEI)=RR
            INMB(NNEI)=INEB
C
C           Calculation of coordinates of neighbour sphere
C
            COTET0_(NNEI)=ZR/RR
            TET0(NNEI)=ACOS(ZR/RR)
            SITET0_(NNEI)=DSIN(TET0(NNEI))
            IF(DABS(XR)+DABS(YR).GE.1.D-5) THEN
               FI00=DATAN2(YR,XR)
               IF(FI00.LT.0.0D0) FI00=FI00+2.0D0*PI
            ELSE
               FI00=PI/2.0D0
            ENDIF
            FI0_(NNEI)=FI00
            COFI0_(NNEI)=DCOS(FI0_(NNEI))
            SIFI0_(NNEI)=DSIN(FI0_(NNEI))
    2       CONTINUE

         ENDDO
C
C        The number of neighbour spheres nneb
C
         NNEB=NNEI
C
C        Arrange of neighbour spheres by distances
C
         DO I=1,NNEB

            RM=1.0D10
            K=0

            DO J=1,NNEB
               IF(RUP(J).GT.RM) GOTO 4
               RM=RUP(J)
               K=J
    4          CONTINUE
            ENDDO

            NMNI(I)=K
            RUP(K)=1.0D10

         ENDDO

C         WRITE (6,'('' INMB'',40I3)') (INMB(I),I=1,NNEB)
C         WRITE (6,'('' NMNI'',40I3)') (NMNI(I),I=1,NNEB)
C         WRITE (6,'('' RDIS'',(1X,8F10.5))') (RDIS(I),I=1,NNEB)

         DO INEV=1,NNEB

            INEB=ABS(NMNI(INEV))
            RR=RDIS(INEB)
            RNE=RA(INMB(INEB))
            CRITRN=0.7D0
            RRED=RNE*CRITRN
            COCRIT=RR/SQRT(RR*RR+RRED*RRED)
            COT0=COTET0_(INEB)
            SIT0=SITET0_(INEB)
            COF0=COFI0_(INEB)
            SIF0=SIFI0_(INEB)

            DO INEN=INEV+1,NNEB

               INEI=NMNI(INEN)
               IF(INEI.LE.0) GOTO 6
               COT=COTET0_(INEI)
               SIT=SITET0_(INEI)
               COF=COFI0_(INEI)
               SIF=SIFI0_(INEI)
C
C              Cos(angle) for neighbour sphere
C
               CT000=SIT*(COF*COF0+SIF*SIF0)
               CTT=SIT0*CT000+COT*COT0
C
C              If angle < critical
C
               IF(CTT.GT.COCRIT) THEN
                  NMNI(INEN)=-INEI
               ENDIF
    6          CONTINUE

            ENDDO

            IF(NMNI(INEV).LE.0.OR.INMB(INEB).GT.NATNEW.
     *      OR.INMB(INEB).LE.IAT) GOTO 5
C
C           Calculation of radius of additional sphere
C
C           WRITE (6,'(/'' RR,RAT,RNE,RSOLV'',5F10.5)') RR,RAT,RNE,RSOLV
C
C           If the distance between spheres is more than solute size
C
            IF(RAT+RNE+2.0D0*RSOLV.LE.RR) GOTO 5

            IF(RR+RNE.LE.RAT.OR.RR+RAT.LE.RNE) THEN
               WRITE (6,'('' ONE SPHERE INCLUDES THE OTHER: ATOM'',
     *         I3,'' NEIGH.'',I3,''  R1='',F10.5,''  R2='',F10.5,
     *         ''  RR='',F10.5)') IAT,INMB(INEB),RAT,RNE,RR
               GOTO 5
            ENDIF

            COSBT=COSFN(RR,RAT+RSOLV,RNE+RSOLV)
            COSBT1=COSFN(RR,RNE+RSOLV,RAT+RSOLV)
            SINBT=SQRT(1.0D0-COSBT*COSBT)
            SINBT1=SQRT(1.0D0-COSBT1*COSBT1)
            XT=(RAT+RSOLV)*SINBT
            VIN=FINT(XT,RSOLV,RSOLV*COSBT)+FINT(XT,RSOLV,RSOLV*COSBT1)
            RA1=RAT*SINBT
            RB1=RNE*SINBT1

            IF(RAT+RNE.GT.RR) THEN
C
C              Spheres are crossed
C
               COSAL=COSFN(RR,RAT,RNE)
               COSAL1=COSFN(RR,RNE,RAT)
               RAB=RAT*SQRT(1.0D0-COSAL*COSAL)
               HA=RAT*(COSAL-COSBT)
               HB=RNE*(COSAL1-COSBT1)
               VSF=VSEC(RA1,RAB,HA)+VSEC(RB1,RAB,HB)
               VEXC=VIN-VSF
C
C              Check, if additional sphere really need
C
               IF(VEXC.LT.VSOLV) GOTO 5

               COSGM=COSBT*COSAL+SINBT*SQRT(1.0D0-COSAL*COSAL)
               XX_=DLENF(RAT+RSOLV,RAT,COSGM)
               XXX=SQRT(XX_)
               COSDL=COSFN(XXX,RAT+RSOLV,RAT)
               SINDL=SQRT(1.0D0-COSDL*COSDL)
               SINEP=SINDL*COSBT+COSDL*SINBT
               XF=(RSOLV+RAT)*SINBT/SINEP
               COF=(RSOLV+RAT)*SINDL/(SINEP*RR)

            ELSE
C
C              Spheres are not crossed
C
               HA=RAT*(1.0D0-COSBT)
               HB=RNE*(1.0D0-COSBT1)
               VSF=VSEC(RA1,0.D0,HA)+VSEC(RB1,0.D0,HB)
               IF(XT.LT.RSOLV) THEN
                  XV=SQRT((RSOLV*RSOLV)-XT*XT)
                  VDP=2.0D0*FINT(XT,RSOLV,XV)
                  VIN=VIN-VDP
               ENDIF
               VEXC=VIN-VSF
C
C              Check, if additional sphere is really needed
C
               IF(VEXC.LT.VSOLV) GOTO 5

               COSAL=COSFN(RR,RNE+RSOLV,RAT+RSOLV)
               DL1=DLENF(RAT,RAT+RSOLV,COSBT)
               DL2=DLENF(RNE,RNE+RSOLV,COSAL)
               XF=0.50D0*SQRT(2.0D0*DL1+2.0D0*DL2-(RR-RAT-RNE)**2)
               IF(2.0D0*(XF-RSOLV)+RAT+RNE.GT.RR) THEN
C
C                 One additional sphere
C
                  COF=(RAT+(RR-RAT-RNE)/2.0D0)/RR

               ELSE
C
C                 Two additional spheres
C
                  COST=COSBT
                  RADT=RAT
                  CALL SEALINV(0.D0,0.D0,XX_,FF,1.D-3,RAT/2.0D0)
                  COF=XX_/RR
                  RDOP=RDS
                  XF=RDS+RSOLV
                  NSFR=NSFR+1
                  COFACT(NSFR)=COF
                  RA(NSFR)=RDOP
                  XA(NSFR)=XAT+(XA(INMB(INEB))-XAT)*COF
                  YA(NSFR)=YAT+(YA(INMB(INEB))-YAT)*COF
                  ZA(NSFR)=ZAT+(ZA(INMB(INEB))-ZAT)*COF
                  COF=(RR-XX_*RNE/RAT)/RR
                  XF=RDS1+RSOLV

               ENDIF
            ENDIF

            RDOP=XF-RSOLV
            NSFR=NSFR+1
            COFACT(NSFR)=COF
            RA(NSFR)=RDOP
            XA(NSFR)=XAT+(XA(INMB(INEB))-XAT)*COF
            YA(NSFR)=YAT+(YA(INMB(INEB))-YAT)*COF
            ZA(NSFR)=ZAT+(ZA(INMB(INEB))-ZAT)*COF
            IF(NSF.LT.NSFR) THEN
               WRITE (6,'(/'' SECNEB: the number of spheres NSFR = '',
     *         I4,'' > value of paramter NSF = '',I4/
     *         '' Program stopped.''/
     *         '' Suggestion: increase parameter NSF in SIZES''/)')
     *         NSFR,NSF
               STOP
            ENDIF

    5       CONTINUE

         ENDDO

      ENDDO

      NATNEW=NSFR
C
C     Check whether new spheres were added
C
      IF(NSFR.NE.NSFR0) GOTO 50
C
C     Check spheres which don't give contribution to cavity surface
C
      DO IAT=NAT+1,NSFR

         XAT=XA(IAT)
         YAT=YA(IAT)
         ZAT=ZA(IAT)
         RAT=RA(IAT)

         DO INEB=1,NSFR

            IF(INEB.EQ.IAT) GOTO 7
            XNE=XA(INEB)
            YNE=YA(INEB)
            ZNE=ZA(INEB)
            RNE=RA(INEB)
            XR=XNE-XAT
            YR=YNE-YAT
            ZR=ZNE-ZAT
            RR=DSQRT(XR*XR+YR*YR+ZR*ZR)
            IF(RR+RAT.LE.RNE) RA(IAT)=0.D0
    7       CONTINUE

         ENDDO

      ENDDO
C
C     Deleting inactive spheres
C
      I=NAT
      DO IAT=NAT+1,NSFR
         I=I+1
         IF(RA(I).EQ.0.D0) THEN
            DO J=I,NATNEW
               XA(I)=XA(I+1)
               YA(I)=YA(I+1)
               ZA(I)=ZA(I+1)
               RA(I)=RA(I+1)
            ENDDO
            NATNEW=NATNEW-1
            I=I-1
         ENDIF
      ENDDO

      NSFR=NATNEW

      IF(NAT.NE.NATNEW)
     *    WRITE (6,'(/'' THE NUMBER OF SHPERES WAS INCREASED FROM'',I4,
     *   ''  TO'',I4)') NAT,NATNEW
C      WRITE (6,'((1X,I5,4F10.5))') (I,XA(I),YA(I),ZA(I),RA(I),
C     *I=1,NSFR)

      RETURN
      END SUBROUTINE SECNEB

!======================================================================!
      SUBROUTINE SECRES(XA,YA,ZA,RA,IAT,SECTTE,SECTFI,ISCFI,
     *                  ITMX,IFMX,DFSTI,NAT,NNEB,NMSBR,NMSABS)
C
C     DATTE0 - ARRAY OF ANGLES WRT TETA
C     DPRAM  - ARRAY OF INCREASES OF THE STEP IN FI OVER THE STEP
C              IN TETA FOR THE SECTOR IN TETA # I
C     DRASD  - COEFFICIENT OF THE SMALLEST DISTANCE
C              BETWEEN SECTIONS IN FI
C     DSMIN  - COEFFICIENT OF THE MINIMUM SIZE OF THE BIGGEST SECTOR
C
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: NAT, IAT
      real(8),  intent(in) :: XA(NAT),YA(NAT),ZA(NAT),RA(NAT)

      real(8),  intent(out) :: SECTTE(NB,NT),SECTFI(NBGEN)
      integer, intent(out) :: ISCFI(NB,NT)

      integer, intent(inout) :: ITMX(NB),IFMX(NB,NT)

      real(8),  intent(out) :: DFSTI
      integer, intent(out) :: NNEB
      integer, intent(out) :: NMSBR(NB),NMSABS(NB)

      ! local variables
      real(8) :: FRES(NB), SECT00(NF)
      real(8), save :: PI
      LOGICAL, save :: PRINT, FIRST=.true.
c
      !INCLUDE 'SPHSIZES'
      !LOGICAL PRINT, FIRST
      !CHARACTER    KEYWRD*241
      !COMMON /RESINF/DATTE0,DPRAM,DRASD,DSMIN,ITMA
      !COMMON /RESCON/SITET0,COTET0,FI0,SIFI0,COFI0,
      !               COTECR,SITECR,TETECR,
      !RRECR,DATTET
      !COMMON /KEYWRD/ KEYWRD
      !DIMENSION XA(NAT),YA(NAT),ZA(NAT),RA(NAT)
      !DIMENSION SITET0(NB),COTET0(NB)
      !DIMENSION FI0(NB),SIFI0(NB),COFI0(NB)
      !DIMENSION NMSBR(NB),NMSABS(NB)
      !DIMENSION COTECR(NB),SITECR(NB),TETECR(NB),RRECR(NB)
      !DIMENSION SECTTE(NB,NT),ISCFI(NB,NT),SECTFI(NBGEN),FRES(NB)
      !DIMENSION DATTE0(MXRAS),DPRAM(MXRAS),DATTET(MXRAS)
      !DIMENSION ITMX(NB),IFMX(NB,NT)
      !DIMENSION SECT00(NF)
      !SAVE PI,PRINT
      !DATA FIRST/.TRUE./
C
      IF (FIRST) THEN
         FIRST=.FALSE.
         PRINT=(INDEX(KEYWRD,'SECRES').NE.0)
         PI=ATAN(1.0D0)*4.0D0
      ENDIF

      XAT=XA(IAT)
      YAT=YA(IAT)
      ZAT=ZA(IAT)
      RAT_=RA(IAT)
      DO I=1,NAT
         NMSABS(I)=0
      ENDDO
C
C     DEFINITION OF NEIGHBOURS
C     NNEI - VARYING NUMBER OF THE NEIGHBOUR
C
      NNEI=0
C
C     FIRST PASS - DEFINING THE ANGLE OF INTERSECTION OF SPHERES
C     FOR ARRANGING THEM IN ORDER
C
      DO 15 INEB=1,NAT

         IF(INEB.EQ.IAT) GOTO 15
         XNE=XA(INEB)
         YNE=YA(INEB)
         ZNE=ZA(INEB)
         RNE_=RA(INEB)
         XR=XNE-XAT
         YR=YNE-YAT
         ZR=ZNE-ZAT
         RR_=SQRT(XR*XR+YR*YR+ZR*ZR)
C
C        IF THE ATOMS DO NOT INTERSECT
C
         IF(RR_.GE.RAT_+RNE_) GOTO 15
         IF(RR_+RNE_.LE.RAT_) GOTO 15
C
C        IF THE GIVEN ATOM LIES INSIDE ANOTHER
C
         IF(RR_+RAT_.LE.RNE_) THEN
            NNEB=-1
            RETURN
         ENDIF
         NNEI=NNEI+1
         NMSBR(NNEI)=INEB
C
C        CALCULATING THE ANGLE OF INTERSECTION OF SPHERES
C
         RRECR(NNEI)=RR_
         COTECR(NNEI)=(RAT_*RAT_-RNE_*RNE_+RR_*RR_)/(2.D0*RAT_*RR_)

   15 CONTINUE
      JMAXF=1

      DO 50 INEB=1,NNEI-1
         COTR=2.D0
         DO 55 JNEB=INEB,NNEI
            IF(COTR.LE.COTECR(JNEB)) GOTO 55
            COTR=COTECR(JNEB)
            JTR=JNEB
   55    CONTINUE
         COTECR(JTR)=COTECR(INEB)
         COTECR(INEB)=COTR
         RR_=RRECR(INEB)
         RRECR(INEB)=RRECR(JTR)
         RRECR(JTR)=RR_
         NN_=NMSBR(INEB)
         NMSBR(INEB)=NMSBR(JTR)
         NMSBR(JTR)=NN_
         IF(INEB.EQ.1) JMAXF=JTR
   50 CONTINUE

      DO 60 INEB=1,NNEI

         INEV=NMSBR(INEB)
         XNE=XA(INEV)
         YNE=YA(INEV)
         ZNE=ZA(INEV)
         RNE_=RA(INEV)
         XR=XNE-XAT
         YR=YNE-YAT
         ZR=ZNE-ZAT
         RR_=RRECR(INEB)
C
C        DEFINING THE COORDINATES OF THE NEIGHBOUR
C
         COTET0(INEB)=ZR/RR_
         XSQRT=1.D0-COTET0(INEB)**2
         IF(XSQRT.LT.0.D0) XSQRT=0.D0
         SITET0(INEB)=SQRT(XSQRT)
         IF(DABS(XR)+DABS(YR).GE.1.D-5) THEN
            FI00=DATAN2(YR,XR)
            IF(FI00.LT.0.D0) FI00=FI00+2.D0*PI
         ELSE
            FI00=PI/2.D0
         ENDIF
         FI0(INEB)=FI00
         DRD=XR*XR+YR*YR
         IF(DRD.GE.1.D-10) THEN
            DRS=SQRT(DRD)
            COF=XR/DRS
            SIF=YR/DRS
         ELSE
            COF=0.D0
            SIF=1.D0
         ENDIF
         COFI0(INEB)=COF
         SIFI0(INEB)=SIF
C
C        CALCULATING THE ANGLE OF INTERSECTION OF SPHERES
C
         TETECR(INEB)=ACOS(COTECR(INEB))
         XSQRT=1.D0-COTECR(INEB)*COTECR(INEB)
         IF(XSQRT.LT.0.D0) XSQRT=0.D0
         SITECR(INEB)=SQRT(XSQRT)
         IF(PRINT) WRITE (6,'('' NEIGBBOUR #'',I3,''  ATOM #'',I3,
     *                ''  TETA='',F6.1)') INEB,INEV,R2D(TETECR(INEB))
   60 CONTINUE
C
C     THE NUMBER OF NEIGHBOURS NNEB
C
      NNEB=NNEI
      IF(NNEI.GT.NB) THEN
         WRITE (6,'(/'' THE NUMBER OF NEIGHBOURS TO ATOM #'',I3,
     *           '' IS TOO LARGE''/'' NNEB='',I3,'' > NB='',I3,
     *           ''  INCREASE PARAMETER NB'')') IAT,NNEB,NB
         STOP
      ENDIF
C
C     EXCLUDING SPHERES, THE INTERSECTION OF WHICH LIES INSIDE
C     ANOTHER INTERSECTION
C
      DO INEB=1,NNEI

         IF(INEB.GT.NNEB) cycle
         COT=COTET0(INEB)
         SIT=SITET0(INEB)
         COF=COFI0(INEB)
         SIF=SIFI0(INEB)
         COTEF=COTECR(INEB)
         INEI=INEB

         DO INEV=INEB+1,NNEI

            INEI=INEI+1
            IF (INEI.GT.NNEB) cycle
            COT0=COTET0(INEI)
            SIT0=SITET0(INEI)
            COF0=COFI0(INEI)
            SIF0=SIFI0(INEI)
            COTEF0=COTECR(INEI)
            SITEF0=SITECR(INEI)
            COTAL=SIT*SIT0*(COF*COF0+SIF*SIF0)+COT*COT0
            XSQRT=1.D0-COTAL*COTAL
            IF (XSQRT.LT.0.D0) XSQRT=0.D0
            COTS=COTAL*COTEF0-SQRT(XSQRT)*SITEF0
            IF (COTS.LT.COTEF) cycle
C
C           SPHERE IS EXCLUDED
C
            NNEB=NNEB-1

            DO I=INEI,NNEB

               COTET0(I)=COTET0(I+1)
               SITET0(I)=SITET0(I+1)
               COFI0(I)=COFI0(I+1)
               SIFI0(I)=SIFI0(I+1)
               FI0(I)=FI0(I+1)
               RRECR(I)=RRECR(I+1)
               COTECR(I)=COTECR(I+1)
               TETECR(I)=TETECR(I+1)
               SITECR(I)=SITECR(I+1)
               NMSBR(I)=NMSBR(I+1)

            ENDDO

            INEI=INEI-1

         ENDDO

      ENDDO

C
C     CONSTRUCTING THE RAYS OF SECTORS, CORRESPONDING
C     TO EACH OF THE NEIGHBOURS
C     SECTORS WRT TETA - I - SECTTE(INEB,I)
C     SECTORS WRT FI - J - SECFI(INEB,I,J)
C
      IONEAT=0

      IF(NNEB.EQ.0) THEN
C
C        THE ATOM HAS NO NEIGHBOURS - EMULATE ONE NEIGHBOUR
C
         IONEAT=1
         NNEB=1
         RRECR(1)=2.D0*RAT_
         COTECR(1)=1.D0
         SITECR(1)=0.D0
         TETECR(1)=0.D0
         COTET0(1)=0.D0
         SITET0(1)=1.D0
         FI0(1)=0.D0
         COFI0(1)=1.D0
         SIFI0(1)=0.D0
      ENDIF

      STFIMS=0.D0
      NTFIMS=0
      ISFDL=0

      DO 20 INEB=1,NNEB
C
C        CALCULATING THE SECTORS OF ATOM INEB
C
C        MAKING THE SECTORS A WHOLE NUMBER, BY FACTOR AA
C
         TETERR=TETECR(INEB)
         TETOS= PI-TETERR
         PD2=1.0D0/SQRT(PI)
         TESS=0.D0

         IF(IONEAT.NE.1) THEN

            DO 31 I=1,ITMA
               TT=DATTE0(I)*PD2
               IF(TESS+TT.GT.TETOS) GOTO 32
   31       TESS=TESS+DATTE0(I)
   32       A2=TETOS/(TESS+TT)

            TT2=TT
            II=I

            IF(I.GT.1) THEN
               TESS=TESS-DATTE0(I-1)
               TT=DATTE0(I-1)*PD2
               TT1=TT
               A1=TETOS/(TESS+TT)
               IF(DABS(A2-1.D0).LT.DABS(1.D0-A1)) THEN
                  AA=A2
                  TT=TT2
                  II=I
               ELSE
                  AA=A1
                  TT=TT1
                  II=I-1
               ENDIF
            ELSE
               AA=A2
            ENDIF

            AA=AA*1.01D0

            DO 33 J=1,ITMA
   33       DATTET(J)=DATTE0(J)*AA
            DATTET(II)=TT*AA

         ELSE

C
C           NO NEIGHBOURS - EMULATE
C
            PD=PI/DATTE0(3)
            I=PD-2.D0*PD2
            A1=PD/(2.D0*PD2+I)
            A2=PD/(2.D0*PD2+I+1)
            TT=PD2*DATTE0(3)

            IF(ABS(1.D0-A2).LT.ABS(A1-1.D0)) THEN
               AA=A2
               II=I+3
            ELSE
               AA=A1
               II=I+2
            ENDIF

            AA=AA*1.01D0
            DO 44 J=2,ITMA
   44          DATTET(J)=DATTE0(3)*AA

            DATTET(1)=TT
            DATTET(II)=TT

         ENDIF

         IF(PRINT) THEN
            WRITE (6,'('' NEIGHBOUR #'',I4,''  AA='',F8.3,
     *                 ''  NO. OF SECTORS WRT TETA'',I4)') INEB,AA,II
            WRITE (6,'('' TETOS,TESS,TT,A1,A2,AA'',6F8.2/
     *          '' DATTET(J)'',9F8.2)') R2D(TETOS),R2D(TESS),R2D(TT),
     *          A1,A2,AA,(R2D(DATTET(J)),J=1,II)
         ENDIF

         ITMXX=0
         DTTE=0.D0
C
C        IF THERE IS ONLY ONE NEIGHBOUR
C
         IF(NNEB.EQ.1) THEN
            ISTSEC=1
            GOTO 22
         ENDIF

         COT0=COTET0(INEB)
         SIT0=SITET0(INEB)
         COF0=COFI0(INEB)
         SIF0=SIFI0(INEB)
         STIM=SECOND()-STIME
         STIME=SECOND()
         STIM0=STIM0+STIM
C
C        NUMBER OF SECTIONS WRT FI
C
         MNE=0
         DO 21 INEI=1,NNEB

            IF(INEI.EQ.INEB) GOTO 21
            COT=COTET0(INEI)
            SIT=SITET0(INEI)
            COF=COFI0(INEI)
            SIF=SIFI0(INEI)
C
C           COS OF THE ANGLE OF DIRECTION TO THE NEIGHBOUR
C
            CT000=SIT*(COF*COF0+SIF*SIF0)
            CTT=SIT0*CT000+COT*COT0
C
C           IF THE ANGLE IS > 120'
C
            IF(CTT.LT.-.5D0) GOTO 21
C
C           IF THE NEIGHBOURS ARE IN ONE LINE
C
            IF(CTT.GT.1.D0-1.D-3) GOTO 21
            MNE=MNE+1
            IF(DABS(COT0).LT.1.D0-1.D-5) GOTO 19
            FFI=FI0(INEI)
            GOTO 18
C
C           SIN OF THE ANGLE OF DIRECTION TO THE NEIGHBOUR
C
  19        XSQRT=1.D0-CTT*CTT
            IF(XSQRT.LE.0.D0) XSQRT=1.D-16
            STT=SQRT(XSQRT)
C
C           THE ANGLE FFI TO THE NEIGHBOUR
C
            CFF=(COT*SIT0-COT0*CT000)/STT
            IF(ABS(CFF).GT.1.D0) CFF=DSIGN(1.D0,CFF)
            FFI=ACOS(CFF)
C
C                    DEFINING THE SIGN OF FFI
C
            FIF=FI0(INEI)-FI0(INEB)
            IF(FIF.GE.0.D0.AND.FIF.LE.PI.AND.FI0(INEB).LE.PI.OR
     *         .(FIF.GE.0.D0.OR.FIF.LE.-PI).AND.FI0(INEB).GT.PI)
     *         FFI=2.D0*PI-FFI

   18       FRES(MNE)=FFI

   21    CONTINUE

         IF(MNE.EQ.0) THEN
            ISTSEC=1
            GOTO 22
         ENDIF
C
C        ARRANGING THE ANGLES FI TO THE NEIGHBOURS IN ORDER
C
         IF(MNE.EQ.1) GOTO 27
         MNEN=MNE-1

         DO 25 I=1,MNEN
            JJ=1
            S=FRES(1)
            MNN=MNE+1-I
            DO 26 J=2,MNN
               IF(FRES(J).LE.S) GOTO 26
               JJ=J
               S=FRES(J)
   26       CONTINUE
            S1=FRES(MNN)
            FRES(MNN)=FRES(JJ)
            FRES(JJ)=S1
   25    CONTINUE
C
C        EXCLUDING CLOSE SUBSECTIONS
C
         DATT=DATTET(1)*DPRAM(1)*DRASD
         DFFR=FRES(MNE)-2.D0*PI
         I=1

   28    IF(I.GT.MNE) GOTO 27
         DFRES=-DFFR+FRES(I)

         IF(DFRES.GT.DATT) THEN
            DFFR=FRES(I)
            I=I+1
            GOTO 28
         ELSE
            MNE=MNE-1
            DO 29 J=I,MNE
   29          FRES(J)=FRES(J+1)
         ENDIF
         GOTO 28
   27    CONTINUE
C
C        IF THERE ARE NO SUBSECTIONS
C
         IF(MNE.EQ.0) THEN
            ISTSEC=1
            GOTO 22
         ENDIF
C
C        CONSTRUCTION OF SECTORS, RELATING TO THE GIVEN NEIGHBOUR
C        CONSTRUCTION OF SECTORS WRT TETA SECTTE(INEB,I)
C
C        WE NEED THE ARRAY DATTET(1;ITMA) -
C        THE ANGLES TETA WRT THE NUMBERS
C
C        SUBSECTIONS ARE ACCOUNTED FOR ONLY
C        IN THE FIRST TWO LAYERS
C
         ISTSMX=2
         DO 30 I=1,ISTSMX
 
            TTER=TETERR+DTTE+DATTET(I)
            IF(TTER.GT.PI) GOTO 34
            ITMXX=ITMXX+1
            SECTTE(INEB,I)=DCOS(TTER)
C
C           CONSTRUCTING THE SECTORS WRT FI SECTFI(INEB,I,J)
C
            TETLN=TETERR+DTTE+DATTET(I)/2.D0
            DFI=DATTET(I)/DSIN(TETLN)*DPRAM(I)
            FINO1=FRES(MNE)-2.D0*PI
            JJ=0
            IJNN=0
 
            DO 35 J=1,MNE

               FINO2=FRES(J)
               FIST=FINO2-FINO1
               NFIST=FIST/DFI+0.5D0
               IF(NFIST.LE.0) NFIST=1
               STFIST=FIST/NFIST
C
C              CORRECTION OF THE AVERAGE SIZE OF THE SECTOR
C
               STFIMS=STFIMS+STFIST/DFI*NFIST
               NTFIMS=NTFIMS+NFIST

               DO 36 JK=1,NFIST
                  JJ=JJ+1
                  SE0000=FINO1+STFIST*(JK-1)
                  IF(SE0000.GT.0.D0) GOTO 37
                  IJNN=IJNN+1
                  SE0000=SE0000+2.D0*PI
   37             CONTINUE
                  IF(SE0000.LE.PI) THEN
                     SECT00(JJ)=DCOS(SE0000)
                  ELSE
                     SECT00(JJ)=-DCOS(SE0000)-2.D0
                  ENDIF
   36          CONTINUE

               FINO1=FINO2

   35       CONTINUE

            ISCFI(INEB,I)=ISFDL
C
C           SORTING THE ARRAY SECTFI IN ORDER
C
            IJNN0=JJ-IJNN
            IF(IJNN0.EQ.0) GOTO 43
            IMA1=ISFDL+IJNN0+IJNN
            IF (IMA1.GT.NBGEN) THEN
               WRITE(6,101)IMA1,NBGEN
               WRITE(6,102)
               STOP
            ENDIF

            DO 40 J=1,IJNN0
  40           SECTFI(ISFDL+J)=SECT00(J+IJNN)

  43        IF(IJNN.EQ.0) GOTO 41
            DO 42 J=1,IJNN
  42           SECTFI(ISFDL+J+IJNN0)=SECT00(J)
  41        CONTINUE
C
C           THE NUMBER OF SECTORS WRT FI
C
            IFMX(INEB,I)=JJ
            ISFDL=ISFDL+JJ

            IF(IFMX(INEB,I).GT.NF) THEN
               WRITE (6,'(/'' NO. OF PARTS OF THE LARGE SECT'',
     *                 ''ORS BEING PREPARED IS TOO LARGE''/'' ISF'',
     *                 ''DL='',I4,'' > NBGEN='',I4)') ISFDL,NBGEN
               STOP
            ENDIF

            IF(IFMX(INEB,I).GT.NF) THEN
               WRITE (6,'(/'' NO. OF LARGE SECTORS WRT FI IS'',
     *               '' TOO LARGE''/'' IFMX='',I3,'' > NF='',I3)')
     *                  IFMX(INEB,I),NF
               STOP
            ENDIF

            DTTE=DTTE+DATTET(I)

   30    CONTINUE
C
C        REMAINING SECTORS WITHOUT SUBSECTIONS
C
         ISTSEC=ISTSMX+1
         GOTO 22
C
C        THE LAST SECTOR WRT TETA
C
   34    CONTINUE
C
C        THE NUMBER OF SECTORS WRT TETA
C
         ITMXX=ITMXX+1
         IFMX(INEB,ITMXX)=1
         ITMX(INEB)=ITMXX
         ISCFI(INEB,ITMXX)=ISFDL
         ISFDL=ISFDL+1

         IF(ITMXX.GT.NT) THEN
            WRITE (6,'(/'' NO. OF LARGE SECTORS WRT TETA'',
     *                  '' IS TOO LARGE''/'' ITMXX='',I3,
     *                  '' > NT='',I3)') ITMXX,NT
            STOP
         ENDIF

         GOTO 20

C
C        IF NO ADDITIONAL SUBSECTIONS
C

   22    CONTINUE

         DO 23 I=ISTSEC,ITMA

            TTER=TETERR+DTTE+DATTET(I)
            IF(TTER.GT.PI) GOTO 23
            ITMXX=ITMXX+1
            SECTTE(INEB,I)=DCOS(TTER)
            TETLN=TETERR+DTTE+DATTET(I)/2.D0
 
            IF(IONEAT.NE.1) THEN
               DFI=DATTET(I)/DSIN(TETLN)*DPRAM(I)
            ELSE
               IF(I.NE.1) THEN
                  DFI=DATTET(I)/DSIN(TETLN)
               ELSE
                  DFI=2.D0*PI
               ENDIF
            ENDIF

            NFIST=2.D0*PI/DFI+0.5D0
            IF(NFIST.LE.0) NFIST=1
            STFIST=2.D0*PI/NFIST
C
C           CORRECTION OF THE AVERAGE SIZE OF THE SECTOR
C
            STFIMS=STFIMS+STFIST/DFI*NFIST
            NTFIMS=NTFIMS+NFIST
            ISCFI(INEB,I)=ISFDL
            IMA1=ISFDL+NFIST

            IF (IMA1.GT.NBGEN) THEN
               WRITE(6,103)IMA1,NBGEN
               WRITE(6,102)
               STOP
            ENDIF

            DO 24 J=1,NFIST
               SE0000=STFIST*(J-1)
               IF(SE0000.LE.PI) THEN
                  SECTFI(ISFDL+J)=DCOS(SE0000)
               ELSE
                  SECTFI(ISFDL+J)=-DCOS(SE0000)-2.D0
               ENDIF
   24       CONTINUE

            IFMX(INEB,I)=NFIST
            ISFDL=ISFDL+NFIST

            IF(IFMX(INEB,I).GT.NF) THEN
               WRITE (6,'(/'' NO. OF PARTS OF THE LARGE SECT'',
     *                     ''ORS BEING PREPARED IS TOO LARGE''/'' ISF'',
     *                     ''DL='',I4,'' > NBGEN='',I4)') ISFDL,NBGEN
               STOP
            ENDIF

            IF(IFMX(INEB,I).GT.NF) THEN
               WRITE (6,'(/'' NO. OF LARGE SECTORS WRT FI IS'',
     *            '' TOO LARGE''/'' IFMX='',I3,'' > NF='',I3)')
     *             IFMX(INEB,I),NF
               STOP
            ENDIF

            JJ=NFIST

   23       DTTE=DTTE+DATTET(I)

         ITMXX=ITMXX+1
         IFMX(INEB,ITMXX)=1
         ITMX(INEB)=ITMXX
         ISCFI(INEB,ITMXX)=ISFDL
         ISFDL=ISFDL+1

         IF(ITMXX.GT.NT) THEN
            WRITE (6,'(/'' NO. OF LARGE SECTORS WRT TETA'',
     *                  '' IS TOO LARGE''/'' ITMXX='',I3,
     *                  '' > NT='',I3)') ITMXX,NT
            STOP
         ENDIF

   20 CONTINUE

C
C     CORRECTION OF THE AVERAGE SIZE OF THE SECTOR
C
      IF(NTFIMS.NE.0) THEN
         DFSTI=STFIMS/NTFIMS
      ELSE
         DFSTI=1.D0
      ENDIF
C
C     COMPLETION OF THE ARRAY OF ADDRESSES
C
      DO 80 I=1,NNEB
   80    NMSABS(NMSBR(I))=I

      STIM=SECOND()-STIME
      STIME=SECOND()
      STIM0=STIM0+STIM
      RETURN

  101 FORMAT(' SECRES: ISFDL+IJNN0+IJNN',12X,'=',I8/
     *' THIS IS LARGER THAN PARAMETER NBGEN =',I8)
  102 FORMAT(' SUCH SITUATION CAN NOT BE TREATED BY PROGRAM. IT IS',
     *' STOPPED'/
     *' SUGGESTION: INCREASE PARAMETER NBGEN, IT SUPPOSED TO BE LARGER'/
     *' THAT ANY POSSIBLE VALUE OF ABOVE EXPRESSION')
  103 FORMAT(' SECRES: ISFDL+NFIST',12X,'=',I8/
     *' THIS IS LARGER THAN PARAMETER NBGEN =',I8)

      END SUBROUTINE SECRES

!======================================================================!
      SUBROUTINE SECSOR(KSECT,SSECT,ISCFI,ITMX,IFMX,
     *LSECTK,LSECTI,LSECTJ,LS,LSNM,LLS,LIS,NSET_,NSEF_,
     *SECDL,SEDDS,NNEB,KK)
***********************************************************************
*
*     THIS SUBPROGRAM ARRANGES SMALL SECTORS FOR THEM ELIMINATION
*
***********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      real(8),  PARAMETER :: GRTN=1.1d0
      integer, PARAMETER :: NTFQ1=NTFQ*GRTN

      integer, intent(in) :: KSECT(NBGEN)
      real(8),  intent(in) :: SSECT(NBGEN)
      integer, intent(in) :: ISCFI(NB,NT)
      integer, intent(in) :: ITMX(NB),IFMX(NB,NT)
      real(8),  intent(in) :: SEDDS(MXRAS)
      integer, intent(in), dimension(NTQ,NFQ) :: LSECTK,LSECTI,LSECTJ

      integer, intent(inout), dimension(:) :: NSET_,NSEF_

      integer, intent(out) :: LS(NBS),LSNM(NBGEN),LLS(NBS),LIS(NBS)
      real(8),  intent(out) :: SECDL(NBS)

      !INCLUDE 'SPHSIZES'
      !INTEGER LSECTK,LSECTI,LSECTJ
C
C     PARAMETER GRTN INCREASES REQUIRED MEMORY FOR ARRAYS NSET,NSEF
C     BUT DECREASES TIME OF ELIMINATING
C
      !PARAMETER (GRTN=1.1)
      !PARAMETER (NTFQ1=NTFQ*GRTN)
C
      !COMMON /TETFI/NTETM,NFIM
C
      !DIMENSION ISCFI(NB,NT)
      !DIMENSION ITMX(NB),IFMX(NB,NT)
      !DIMENSION LSECTK(NTQ,NFQ),LSECTI(NTQ,NFQ),LSECTJ(NTQ,NFQ)
      !DIMENSION KSECT(NBGEN),SSECT(NBGEN)
C
C     ARRAYS FOR DISTRIBUTION OF SMALL SECTORS IN TESSERAE
C
      !DIMENSION LS(NBS)
      !DIMENSION LSNM(NBGEN),LLS(NBS),LIS(NBS)
      !DIMENSION SECDL(NBS),SEDDS(MXRAS)

C
      KK=0
      KFUL=0
      DO 10 K=1,NNEB
      DO 10 I=1,ITMX(K)
      ISFDL=ISCFI(K,I)
      DDS=SEDDS(I)
      DO 10  J=1,IFMX(K,I)
         IF(KSECT(ISFDL+J).LE.0) GOTO 10
         KK=KK+1
         LS(KK)=KSECT(ISFDL+J)
         SECDL(KK)=SSECT(ISFDL+J)/DDS
         LSNM(ISFDL+J)=KK
         KFUL=KFUL+KSECT(ISFDL+J)
   10 CONTINUE
      IF(KFUL.EQ.0) THEN
C
C        IF WE HAVE NO SECTORS - RETURN
C
         KK=0
         RETURN
      ENDIF
      IF(KK.GE.NBS) THEN
           WRITE (6,'(/'' THE NUMBER OF LARGE SECTORS '',
     *        ''EXCEEDES ARRAY SIZE (SECSOR)''/
     *        '' K='',I4,'' > NBS='',I4/
     *        '' INCREASE PARAMETER NBS'')') KK,NBS-1
         STOP
      ENDIF
      IF(KFUL.GT.NTFQ1) THEN
           WRITE (6,'(/'' THE NUMBER OF SMALL SECTORS '',
     *        ''EXCEEDES ARRAY SIZE (SECSOR)''/
     *        '' KFUL='',I4,'' > NTFQ1='',I4/
     *        '' INCREASE PARAMETER GRTN'')') KFUL,NTFQ1
         STOP
      ENDIF
C
C     APR DEFINES THE ADDITIONAL SPACE IN ARRAYS NSET, NSEF
C
      APR=NTFQ1*1.D0/KFUL
C
C     ARRANGE OF SECTORS
C
      LLS(1)=0
      DO 20 KS=1,KK
C
C        THE FIRST POINTS OF SECTORS
C
         LLS(KS+1)=LLS(KS)+LS(KS)*APR
         LIS(KS)=0
   20 CONTINUE
      LIS(KK)=0
      DO 30 ITETM=1,NTETM
      DO 30 IFIM=1,NFIM
         KNEB=LSECTK(ITETM,IFIM)
         IF(KNEB.EQ.0) GOTO 30
         I=LSECTI(ITETM,IFIM)
         J=LSECTJ(ITETM,IFIM)
         ISFDL=ISCFI(KNEB,I)
         IF(KSECT(ISFDL+J).LE.0) GOTO 30
         KS=LSNM(ISFDL+J)
         LIS(KS)=LIS(KS)+1
         NSET_(LLS(KS)+LIS(KS))=ITETM
         NSEF_(LLS(KS)+LIS(KS))=IFIM
   30 CONTINUE
C
C     REARRANGE ARRAY LS USING FREE SPACE
C
      DO 40 I=1,KK
   40 LS(I)=LS(I)*APR
      RETURN
      END SUBROUTINE SECSOR

!======================================================================!
      SUBROUTINE SECTEG(XA,YA,ZA,RA,SNS,XNS,YNS,ZNS,NW_,NAT,KAT,KSST)
***********************************************************************
*
*     THIS SUBPROGRAM IS DRIVER FOR CALCULATING OF CARTESIAN COORDINATES
*     AND SQUARES OF TESSERAE ON THE SURFACE OF CROSSED SPHERES
*
*     ON INPUT:
*
*     XA,YA,ZA - ARRAYS OF CARTESIAN COORDINATES OF CENTERS OF SPHERES
*     RA - ARRAY OF RADII OF SPHERES
*     NAT - THE NUMBER OF SPHERES
*
*     ON OUTPUT
*
*     XNS, YNS, ZNS - ARRAYS OF CARTESIAN COORDINATES OF TESSERAE
*     SNS - ARRAY OF SQUARES OF TESSERAE
*     NW_ - ARRAY OF ABSOLUTE NUMBERS OF FIRST TESSERAE ON EACH SPHERE
*
***********************************************************************
C
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in)  :: NAT, KAT
      real(8),  intent(in)  :: XA(NAT),YA(NAT),ZA(NAT),RA(NAT)
      real(8),  intent(out) :: SNS(KAT),XNS(KAT),YNS(KAT),ZNS(KAT)
      integer, intent(out) :: NW_(NAT+1)
      integer, intent(out) :: KSST

      ! local variables and arrays
      logical :: PRINT, PRTMAP, BIGPRT

      integer, dimension(NTQ,NFQ)  :: LSECTK,LSECTI,LSECTJ
      integer, dimension(NB,NT,NF) :: LSNM

      integer, dimension(NB)    :: NMSBR,NMSABS,ITMX
      integer, dimension(NB,NT) :: IFMX,ISCFI
      integer, dimension(NBS)   :: LS
      integer, dimension(NBGEN) :: KSECT

      real(8), DIMENSION(NTQ)   :: TETM
      real(8), DIMENSION(NB,NT) :: SECTTE
      real(8), DIMENSION(NBS)   :: CRELRD
      real(8), DIMENSION(NBGEN) :: SSECT,XSECT,YSECT,ZSECT,SECTFI
      real(8), DIMENSION(NBS)   :: SS,XS,YS,ZS

      !INCLUDE 'SIZES'
      !INCLUDE 'SPHSIZES'

      !LOGICAL PRINT,PRTMAP,BIGPRT
      !CHARACTER    KEYWRD*241

      !COMMON /KEYWRD/ KEYWRD
      !COMMON /RESINF/DATTE0,DPRAM,DRASD,DSMIN,ITMA
      !COMMON /TETFI/NTETM,NFIM
      !COMMON /SICO/COTETM,SITETM,COFIM,SIFIM,FIM
      !COMMON /RESCON/SITET0,COTET0,FI0,SIFI0,COFI0,
      !               COTECR,SITECR,TETECR,RRECR,DATTET
      !COMMON /SQRMS/NSQRD,NSQDUM,SQRMAS(1001)
      !COMMON /NEBCOR/NNBNEB,COTNB,SITNB,COFNB,SIFNB,COECNB,SIECNB
      !COMMON /CURLEN/NDLATS,NDNEB,NDFSEC,NDFFIR,NDFLAS
      !COMMON /LSECTM/LSECT,LNUSEC
C
      !DIMENSION XA(NAT),YA(NAT),ZA(NAT),RA(NAT),SNS(KAT),
      !*XNS(KAT),YNS(KAT),ZNS(KAT),NW_(NAT)

      !DIMENSION NDNEB(NDNBMX)
      !DIMENSION NDFSEC(NDFMAX),NDFFIR(NDFMAX),NDFLAS(NDFMAX)
      !DIMENSION NDLATS(NSF),NNBNEB(NDNBMX),COTNB(NDNBMX),SITNB(NDNBMX)
      !DIMENSION COFNB(NDNBMX),SIFNB(NDNBMX)
      !DIMENSION COECNB(NDNBMX),SIECNB(NDNBMX)
C
      !DIMENSION TETM(NTQ)
      !DIMENSION SITETM(NTQ),COTETM(NTQ)
      !DIMENSION FIM(NFQ),SIFIM(NFQ),COFIM(NFQ)
      !DIMENSION SITET0(NB),COTET0(NB)
      !DIMENSION FI0(NB),SIFI0(NB),COFI0(NB)
      !DIMENSION COTECR(NB),SITECR(NB),TETECR(NB),RRECR(NB)
      !DIMENSION SECTTE(NB,NT),ISCFI(NB,NT),SECTFI(NBGEN)
      !DIMENSION DATTE0(MXRAS),DPRAM(MXRAS),DATTET(MXRAS)
      !DIMENSION ITMX(NB),IFMX(NB,NT)
      !DIMENSION LSECTK(NTQ,NFQ),LSECTI(NTQ,NFQ),LSECTJ(NTQ,NFQ)
      !DIMENSION LSECT(NTQ,NFQ,3)
      !DIMENSION LNUSEC(NTQ,NFQ)
      !DIMENSION CRELRD(NBS)
      !DIMENSION KSECT(NBGEN),SSECT(NBGEN),XSECT(NBGEN),
      !*YSECT(NBGEN),ZSECT(NBGEN)
      !DIMENSION SS(NBS),XS(NBS),YS(NBS),ZS(NBS),LS(NBS)
      !DIMENSION NMSBR(NB),NMSABS(NB)
      !DIMENSION LSNM(NB,NT,NF)
C
      !EQUIVALENCE (LSECT(1,1,1),LSECTK(1,1))
      !EQUIVALENCE (LSECT(1,1,2),LSECTI(1,1))
      !EQUIVALENCE (LSECT(1,1,3),LSECTJ(1,1))

      ! DON't FORGET TO ASSIGN IT BACK !!!!!!!!!!!!

      LSECTK = LSECT(:,:,1)
      LSECTI = LSECT(:,:,2)
      LSECTJ = LSECT(:,:,3)
C
      PI=DATAN(1.0D0)*4.0D0
      STIME=SECOND()
      STIM0=0.0D0
      STIM1=0.0D0
      STIM2=0.0D0
      STIM3=0.0D0
      STIM4=0.0D0
      STIM5=0.0D0
      NSQRD=1000

      PRINT=(INDEX(KEYWRD,'SECTEG').NE.0).OR.
     *   (INDEX(KEYWRD,'SECEXL').NE.0)
      PRTMAP=(INDEX(KEYWRD,'SECMAP').NE.0)
      BIGPRT=PRTMAP.AND.(INDEX(KEYWRD,'LARGE').NE.0)
C
C     Preparation of arrays TETM, SITETM (SIN), COTETM (COS)
C
      STTET=PI/NTETM

      IF(NTETM.GT.NTQ.OR.NFIM.GT.NFQ) THEN
      WRITE (6,'(/'' THE NUMBER OF STEPS MORE THEN SIZE (SECTEG)''/
     *   '' NTETM,NFIM='',2I4,'' > NTQ,NFQ='',2I4/
     *   '' INCREASE PARAMETERS NTQ,NFQ IN SPHSIZES''/)')
     *   NTETM,NFIM,NTQ,NFQ
         STOP
      ENDIF
      DO I=1,NTETM
         TETMM=STTET*(I*1.D0-0.5D0)
         TETM(I)=TETMM
         SITETM(I)=DSIN(TETMM)
         COTETM(I)=DCOS(TETMM)
      ENDDO
C
C     Preparation of arrays FIM, SIFIM (SIN), COFIM (COS)
C
      STFI=2.0D0*PI/NFIM
      DO I=1,NFIM
         FIMM=STFI*(I*1.0D0-0.5D0)
         FIM(I)=FIMM
         SIFIM(I)=DSIN(FIMM)
         COFIM(I)=DCOS(FIMM)
      ENDDO
      STSQR=1.0D0/NSQRD
      DO I=1,NSQRD+1
         X=STSQR*(I-1)
         IF(X.GT.1.D0) THEN
            X=1.D0
C           WRITE (6,*) 'X=',1.D0-X
         ENDIF
         SQRMAS(I)=SQRT(1.D0-X*X)
      ENDDO
      KSST=0
      NTKAT=0
      NTKDL=0
      NNEBZ=0
      NDUGS=0
C
C     Main cycle per spheres
C

      DO IAT=1,NAT

         XAT=XA(IAT)
         YAT=YA(IAT)
         ZAT=ZA(IAT)
         RAT_=RA(IAT)
         IF(PRINT)
     *   WRITE (6,'(/'' ATOM #'',I3,''  RAT='',F6.3/)') IAT,RAT_
C
C        SPL is a factor for calculating of tesserae square
C
         DSI=4.0D0*PI*DSIN(STTET/2.0D0)/NFIM
         SPL=DSI*RAT_*RAT_
C
C        Definition the structure of tesserae
C
         CALL SECRES(XA,YA,ZA,RA,IAT,SECTTE,SECTFI,ISCFI,
     *   ITMX,IFMX,DFSTI,NAT,NNEB,NMSBR,NMSABS)
C
C        Check if this sphere does not give contribution in the surface
C
         IF(NNEB.LT.0) THEN
            KSS=0
            GOTO 15
         ENDIF
C
C        Clean arrays KSECT,SSECT,XSECK,YSECT,ZSECT
C
         DO J=1,NBGEN
            XSECT(J)=0.D0
            YSECT(J)=0.D0
            ZSECT(J)=0.D0
            KSECT(J)=0
            SSECT(J)=0.D0
         ENDDO
         DO ITETM=1,NTETM
            DO IFIM=1,NFIM
               LSECTK(ITETM,IFIM)=0
               LSECTI(ITETM,IFIM)=0
               LSECTJ(ITETM,IFIM)=0
            ENDDO
         ENDDO

         LSECT(:,:,1) = LSECTK
         LSECT(:,:,2) = LSECTI
         LSECT(:,:,3) = LSECTJ

         STIM=SECOND()-STIME
         STIME=SECOND()
         STIM0=STIM0+STIM
         STIM1=STIM1+STIM
C
C        Calculation ot original tesserae
C
         CALL SECINS(SECTTE,SECTFI,ISCFI,ITMX,IFMX,NNEB,
     *   LSECTK,LSECTI,LSECTJ,KSECT,SSECT,NMSBR,LNUSEC,SPL)
         LSECT(:,:,1) = LSECTK
         LSECT(:,:,2) = LSECTI
         LSECT(:,:,3) = LSECTJ
         STIM=SECOND()-STIME
         STIME=SECOND()
         STIM0=STIM0+STIM
         STIM2=STIM2+STIM
         IF (BIGPRT) CALL PRSHOR(LSECT,NTETM,NFIM,NTQ,NFQ,3)
C
C        Eliminates of the smallest tesserae
C
         CALL SECEXL(SSECT,KSECT,ISCFI,ITMX,IFMX,LSECTK,LSECTI,LSECTJ,
     *   DFSTI,RAT_,NNEB,SPL)
         LSECT(:,:,1) = LSECTK
         LSECT(:,:,2) = LSECTI
         LSECT(:,:,3) = LSECTJ
         IF(PRTMAP) CALL PRSHOR(LSECT,NTETM,NFIM,NTQ,NFQ,3)
         STIM=SECOND()-STIME
         STIME=SECOND()
         STIM0=STIM0+STIM
         STIM3=STIM3+STIM
C
C        Calculation of coordinates and squares of tesserae
C
         DO ITETM=1,NTETM

            COT=COTETM(ITETM)
            SIT=SITETM(ITETM)

            DO IFIM=1,NFIM
               COF=COFIM(IFIM)
               SIF=SIFIM(IFIM)
               KNEB=LSECTK(ITETM,IFIM)
               IF(KNEB.EQ.0) GOTO 1
               I=LSECTI(ITETM,IFIM)
               J=LSECTJ(ITETM,IFIM)
               S=SPL*SIT
               X=SIT*COF*S
               Y=SIT*SIF*S
               Z=COT*S
               ISFDL=ISCFI(KNEB,I)
               XSECT(ISFDL+J)=XSECT(ISFDL+J)+X
               YSECT(ISFDL+J)=YSECT(ISFDL+J)+Y
               ZSECT(ISFDL+J)=ZSECT(ISFDL+J)+Z
    1          CONTINUE
            ENDDO

         ENDDO
C
C        Reordering of tesserae
C
         KS=0
         DO INEB=1,NNEB
            ITMXX=ITMX(INEB)
            DO I=1,ITMXX
               IFMXX=IFMX(INEB,I)
               ISFDL=ISCFI(INEB,I)

               DO J=1,IFMXX

                  IF(KSECT(ISFDL+J).EQ.0) GOTO 2
                  KS=KS+1
                  SS(KS)=SSECT(ISFDL+J)
                  XT=XSECT(ISFDL+J)
                  YT=YSECT(ISFDL+J)
                  ZT=ZSECT(ISFDL+J)
                  RR_=DSQRT(XT*XT+YT*YT+ZT*ZT)
C
C                 RAT - sphere radius
C
                  XS(KS)=XT/RR_*RAT_+XAT
                  YS(KS)=YT/RR_*RAT_+YAT
                  ZS(KS)=ZT/RR_*RAT_+ZAT
C
C                 CRELRD - array for geometry optimization
C
                  CRELRD(KS)=RAT_/RR_
                  LS(KS)=KSECT(ISFDL+J)
                  LSNM(INEB,I,J)=KS
    2             CONTINUE

                ENDDO

            ENDDO
         ENDDO

C
C        KSS - the number of tesserae on this sphere
C
         KSS=KS
C
C        Preparation of final arrays
C
   15    DO KS=1,KSS
            KSST=KSST+1
            SNS(KSST)=SS(KS)
            XNS(KSST)=XS(KS)
            YNS(KSST)=YS(KS)
            ZNS(KSST)=ZS(KS)
         ENDDO
         NW_(IAT)=KSST-KSS+1
         STIM=SECOND()-STIME
         STIME=SECOND()
         STIM0=STIM0+STIM
         STIM4=STIM4+STIM
      ENDDO !IAT=1,NAT
      NDLATS(NAT+1)=NTKAT+1
      NDNEB(NTKAT+1)=NTKDL+1
      NW_(NAT+1)=KSST+1
      RETURN
      END SUBROUTINE SECTEG

!======================================================================!
      SUBROUTINE SFERA1T
      IMPLICIT REAL(8) (A-H,O-Z)

      logical :: TIMES,MERT,NSMO,PRNT,BIGPRT
      logical, save :: FIRST=.true.
      real(8), DIMENSION(NSF) :: CCC

      !INCLUDE 'SIZES'
      !CHARACTER*241 KEYWRD
      !LOGICAL TIMES,MERT,NSMO,PRNT,BIGPRT,FIRST
      !COMMON /FACTSF/ FACTOR,FACTOR2
      !COMMON /PR   /XX(NSF),YY(NSF),ZZ(NSF),RV(NSF),QSFE(NSF)
      !COMMON /PREL/RVEL(NSF),QSFEEL(NSF)
      !COMMON /PROLD/XXOLD(NSF),YYOLD(NSF),ZZOLD(NSF),RVOLD(NSF),
      !*QSFOLD(NSF)
      !COMMON /SOLMAT/ X0(  NS),Y0(  NS),Z0(  NS),AS(  NS),
      !1  QS(  NS),QCOR(  NS)
      !COMMON /SOLMAEL/
      !*X0EL(NS),Y0EL(NS),Z0EL(NS),ASEL(NS),QSEL(NS),QCOREL(NS)
      !COMMON /SOLMT1/ X01(NS),Y01(NS),Z01(NS),AS1(NS),QCOR1(NS),IJ0,
      !*   NW0(NSF)
      !COMMON /SOLMT1EL/X01EL(NS),Y01EL(NS),Z01EL(NS),AS1EL(NS),
      !*   QCOR1EL(NS),IJ0EL,NW0EL(NSF)
      !COMMON /CHH   / CHARC, CHARE
      !COMMON /SFE   / NN,KFI,LTH
      !COMMON /SFE2  / IJ,NW(NSF)
      !COMMON /SFE2EL/ IJEL,NWEL(NSF)
      !COMMON /TIMSF / TSF1,TSF2,TSF3,TSFE,TCOLA,TCOLE,TCONNU
      !COMMON /TIMING/ TIME0
      !COMMON /ENUCL / EN, EPO,EVAC,ESCFOL
      !COMMON /MOLKST/ NUMAT
      !COMMON /MODF  /  ITE1,NCFINR1,ITE2,NCFINR2
      !COMMON /KEYWRD/ KEYWRD
      !DIMENSION CCC(NSF)
      !DATA FIRST/.TRUE./

      TIMEA=SECOND()
      IF (FIRST ) THEN
         FIRST=.FALSE.
         TSF1=0.0D0
      ENDIF
      TIMES  =(INDEX(KEYWRD,'TIMES') .NE. 0 )
      PRNT=(INDEX(KEYWRD,' SF1') .NE. 0) .AND.
     *   (INDEX(KEYWRD,'DEBUG') .NE. 0)
      BIGPRT=(INDEX(KEYWRD,'LARGE') .NE. 0 ) .AND. PRNT
      MERT   = (INDEX(KEYWRD,'MERT=') .NE. 0 )
      NSMO   = (INDEX(KEYWRD,'NOSMOOTH') .NE. 0 )
C
C ***
C  CALCULATION OF THE CARTESIAN COORDINATES AND SQUARES OF THE TESSERA.
C ***
C
C        SMOOTHING OF THE SURFACE BY INCORPORATING ADDITIONAL SPHERES
C
      IF (.NOT.NSMO) THEN
         IF (PRNT) WRITE (6,'('' SMOOTHING'')')
         CALL SECNEB (XX,YY,ZZ,RVEL,CCC,NN,NNEL)
         DO I=NN+1,NNEL
            RV(I)=RVEL(I)+FACTOR2
         ENDDO
         NN=NNEL
	 IF(PRNT) THEN
         WRITE (6,'(/'' TOTAL NUMBER OF SPHERES FOR BOTH CAVITIES''/
     *   '' AFTER ADDING EXTRASPHERES FOR SMOOTHING ='',I5)') NN
         ENDIF
      ENDIF
      IF (PRNT) THEN
         WRITE (6,'(/'' IN SFERA1T'')')
         WRITE (6,'('' COORDINATES AND RADII OF SPHERES'')')
         WRITE (6,'(4X,''I '',6X,''XX'',8X,''YY'',8X,''ZZ'',8X,
     *      ''RV'',7X,''RVEL'')')
         DO I=1,NN
             WRITE(6,8) I,XX(I),YY(I),ZZ(I),RV(I),RVEL(I)
         ENDDO
      ENDIF
      IF ( .NOT.MERT ) THEN
         IF (PRNT) WRITE (6,'('' SECTEG'')')
         CALL MODFE(ITE1,NCFINR1)
         CALL SECTEG(XX,YY,ZZ,RVEL,ASEL,X0EL,Y0EL,Z0EL,NWEL,NN,NS,IJEL)
         CALL MODFE(ITE2,NCFINR2)
         CALL SECTEG(XX,YY,ZZ,RV,AS,X0,Y0,Z0,NW,NN,NS,IJ)
         IF(PRNT) THEN
            WRITE(6,*)'SFERA1T: NW(NN)=',(NW(I),I=1,NN)
            WRITE(6,*)'SFERA1T: NW(NN)=',(NWEL(I),I=1,NN)
         ENDIF
      ELSE
         IF (PRNT) WRITE (6,'('' MERTU'')')
         CALL MERTU(KFI,LTH,NN,XX,YY,ZZ,RVEL,NWEL,
     ,              X0EL,Y0EL,Z0EL,ASEL,IJEL)
         CALL MERTU(KFI,LTH,NN,XX,YY,ZZ,RV,NW,X0,Y0,Z0,AS,IJ)
      END IF
      IF ( IJ.GT.NS.OR.IJEL.GT.NS ) THEN
           WRITE (6,'( //6X,''NUMBER OF SURFACE ELEMENTS IS'',I6,
     *     '' AND'',I6,
     *     ''  IS MORE THEN  NS  - '',I6//6X,''- ( THE DIMENSION'',
     *     ''  OF THE COMMON BLOCK )''/'' INCREASE PARAMETER NSECMX'')')
     *         IJ,IJEL,NS
           STOP
      END IF
      TIMEB=SECOND()
      TSF1=TSF1+TIMEB-TIMEA
      IF(TIMES) THEN
           WRITE(6,'(''##### TIME TOTAL IS  '',F8.2,
     1     ''    SFERA1  '',F8.2)')TIMEB-TIME0,TIMEB-TIMEA
      END IF
      RETURN
    8 FORMAT(2X,I4,2X, 5F9.4)
      END SUBROUTINE SFERA1T

!======================================================================!
      SUBROUTINE SFERA3T(X0_,Y0_,Z0_,AS_,QS_,NW_,IJ_,X0EL_,Y0EL_,Z0EL_,
     *                   ASEL_,QSEL_,NWEL_,IJEL_)
***********************************************************************
*
*     THIS SUBPROGRAM CALCULATES THE SELF-POLARIZATION OF
*     SURFACE ELEMENTS
*
***********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)

      logical, save :: FIRST=.true., SF, TIMES, UNIFORM1

      real(8), intent(in), dimension(:) :: X0_,Y0_,Z0_,AS_
      real(8), intent(in), dimension(:) :: X0EL_,Y0EL_,Z0EL_,ASEL_
      real(8), intent(inout), dimension(:) :: QS_,QSEL_
      integer, intent(in), dimension(:) :: NW_, NWEL_
      integer, intent(in) :: IJ_,IJEL_

      ! local arrays
      real(8), dimension(NS) :: QSO,QS1,QSOEL,QS1EL,QS12,QS21

      !INCLUDE 'SIZES'
      !LOGICAL FIRST,SF,TIMES,UNIFORM1
      !CHARACTER*241 KEYWRD
      !COMMON /VOL/     CON,VSOLV,RSOLV,SELFCR,CHDIFF,ITSE
      !COMMON /PECAR/  COL2,COL1,COIN,COEL
      !COMMON /PR   /XX(NSF),YY(NSF),ZZ(NSF),RV(NSF),QSFE(NSF)
      !COMMON /PREL/RVEL(NSF),QSFEEL(NSF)
      !COMMON /REGSO / IREG
      !COMMON /TIMSF / TSF1,TSF2,TSF3,TSFE,TCOLA,TCOLE,TCONNU
      !COMMON /TIMING/ TIME0
      !COMMON /KEYWRD/ KEYWRD
      !COMMON /FACTSF/ FACTOR,FACTOR2
      !DIMENSION X0_(*),Y0_(*),Z0_(*),AS_(*),NW_(*),QS_(*)
      !DIMENSION X0EL_(*),Y0EL_(*),Z0EL_(*),ASEL_(*),NWEL_(*),QSEL_(*)
      !DIMENSION QSO(NS),QS1(NS)
      !DIMENSION QSOEL(NS),QS1EL(NS)
      !DIMENSION QS12(NS),QS21(NS)
      !SAVE FIRST,SF,TIMES,UNIFORM1
      !DATA FIRST /.TRUE./

      TIMEA=SECOND()
      IF ( FIRST ) THEN
          FIRST  = .FALSE.
          SF=(INDEX(KEYWRD,' SFERE') .NE. 0)
          TIMES  =(INDEX(KEYWRD,'TIMES') .NE. 0 )
          TSF3=0.0D0
          UNIFORM1=(INDEX(KEYWRD,'UNI1') .NE. 0 )
      END IF
      IF ((INDEX(KEYWRD,'UNIFORM').NE.0)) IREG=3
      IF (IREG.EQ.3) GOTO 55
      PI=ATAN(1.0D0)*4.0D0
      DMP=0.65D00
      TLM=1.0D00-DMP
      DCC2=4.D0*PI
      DCC1=PI*2.D00
      DCC=COL2/DCC2
      DCCT=COL1/DCC2
      DCC3=1.D0+DCC1*DCC
      DCC3T=1.D0+DCC1*DCCT
      DQSS=0.D0
      DQSST=0.D0
      DO I=1,IJ_
         DQSS=DQSS+AS_(I)
      ENDDO
      DO I=1,IJEL_
         DQSST=DQSST+ASEL_(I)
      ENDDO
C
CR      QSO=QS
C
      DO I=1,IJ_
         QSO(I)=QS_(I)
      ENDDO
      DO I=1,IJEL_
         QSOEL(I)=QSEL_(I)
      ENDDO
C
CR  11  QS1=QS_
C
      DO M=1,ITSE

  11     DO I=1,IJ_
           QS1(I)=QS_(I)
           QS21(I)=0.D0
           QS_(I)=0.D0
         ENDDO
         DO I=1,IJEL_
           QS1EL(I)=QSEL_(I)
           QS12(I)=0.D0
           QSEL_(I)=0.D0
         ENDDO
         CALL SFERT0(X0_,Y0_,Z0_,AS_,QS1,QS_,NW_,XX,YY,ZZ,RV)
         CALL SFERT0(X0EL_,Y0EL_,Z0EL_,ASEL_,QS1EL,QSEL_,NWEL_,
     *      XX,YY,ZZ,RVEL)
         CALL SFERT2(X0_,Y0_,Z0_,QS1EL,QS21,NW_,XX,YY,ZZ,RV,
     *      X0EL_,Y0EL_,Z0EL_,IJEL_)
         CALL SFERT2(X0EL_,Y0EL_,Z0EL_,QS1,QS12,NWEL_,XX,YY,ZZ,RVEL,
     *      X0_,Y0_,Z0_,IJ_)
         QQ1=0.0D0   !Total charge on inner cavity initialization
         QQ2=0.0D0   !Total charge on outer cavity intialization
         IF (COIN.NE.0.0D0) THEN
           IF (UNIFORM1.OR.FACTOR2.EQ.0.0D0) THEN
C             WRITE (6,'('' UNIFORM APPROXIMATION NO.1 WILL BE USED'')')
              DO I=1,IJ_
                  QS_(I)=(QSO(I)+QS_(I)*AS_(I)*DCC)/DCC3
                  QS_(I)=QS_(I)*DMP+QS1(I)*TLM
                  QQ2=QQ2+QS_(I)
              ENDDO
           ELSE
              DO I=1,IJ_
                 QS_(I) = (QSO(I)*COL2/COIN + (QS_(I) +
     +                     QS21(I))*AS_(I)*DCC)/DCC3
                 QS_(I) = QS_(I)*DMP+QS1(I)*TLM
                 QQ2 = QQ2+QS_(I)
              ENDDO
           ENDIF
         ELSE
           DO I=1,IJ_
              QS_(I)=0.0D0
           ENDDO
         ENDIF
         IF (COL1.NE.0.0D0) THEN
           IF (UNIFORM1.OR.FACTOR2.EQ.0.0D0) THEN
C             WRITE (6,'('' UNIFORM APPROXIMATION NO.1 WILL BE USED'')')
              DO I=1,IJEL_
                  QSEL_(I)=(QSOEL(I)+QSEL_(I)*ASEL_(I)*DCCT)/DCC3T
                  QSEL_(I)=QSEL_(I)*DMP+QS1EL(I)*TLM
                  QQ1=QQ1+QSEL_(I)
              ENDDO
           ELSE
              DO I=1,IJEL_
                 QSEL_(I) = (QSOEL(I) + (QSEL_(I) +
     +                       QS12(I))*ASEL_(I)*DCCT)/DCC3T
                 QSEL_(I) = QSEL_(I)*DMP + QS1EL(I)*TLM
                 QQ1=QQ1+QSEL_(I)
              ENDDO
           ENDIF
         ELSE
           DO I=1,IJEL_
              QSEL_(I)=0.0D0
           ENDDO
         ENDIF

CR       DEL=QS1-QS_
         DQSQ=0.D0
         DQSQT=0.D0
         DO I=1,IJ_
            IF(ABS(AS_(I)).GT.1.D-4) THEN
               DQSQ=DQSQ+(QS_(I)-QS1(I))*(QS_(I)-QS1(I))/AS_(I)
            ENDIF
         ENDDO
         DO I=1,IJEL_
            IF(ABS(ASEL_(I)).GT.1.D-4) THEN
               DQSQT=DQSQT+(QSEL_(I)-QS1EL(I))*
     *                     (QSEL_(I)-QS1EL(I))/ASEL_(I)
            ENDIF
         ENDDO
         DQNN1=SQRT(DQSQT/DQSST)
         DQNN2=SQRT(DQSQ/DQSS)

         IF (SF) THEN
            WRITE(6,'(/1X,''. ITER'',I2,'' Sum(sigma_1)='',F10.6,
     *      '' MISFIT='',D10.3)')M,QQ1,DQNN1
            WRITE(6,'( 1X,''. ITER'',I2,'' Sum(sigma_2)='',F10.6,
     *      '' MISFIT='',D10.3)')M,QQ2,DQNN2
         ENDIF

CR      IF DEL>CRITER GOTO 11
        IF(DABS(DQNN1).LE.SELFCR*0.1D0 .AND. DABS(DQNN2).LE.SELFCR)
C        IF(DABS(DQNN1).LE.SELFCR .AND. DABS(DQNN2).LE.SELFCR)
     *  GOTO 55

      ENDDO ! OVER M - NUMBER OF ITERATIONS

      WRITE(6,100)ITSE
      STOP

  55  IF(COIN.NE.0.0D0)CALL SFERT1(QS_,NW_,IJ_,COIN,CHDIFF,QSFE)
      IF(COL1.NE.0.0D0)CALL SFERT1(QSEL_,NWEL_,IJEL_,COL1,CHDIFF,QSFEEL)
      TIMEB=SECOND()
      TSF3=TSF3+TIMEB-TIMEA
      IF (TIMES ) THEN
         WRITE(6,'(''##### TIME TOTAL IS  '',F8.2,
     1   ''    SFERA3  '',F8.2)')TIMEB-TIME0,TIMEB-TIMEA
      END IF
      RETURN
  100 FORMAT(/'SFERA3T: Cavity surface charges are not converged',
     *' within the maximum allowed '/
     *' number of iterations in polarization procedure ITSE = ',I3/
     *' Program stopped at this point'//
     *' Suggestions:'//
     *' It seems that something wrong with parameters specifying the',
     *' number'/
     *' and sizes of tesserae on cavity surfaces. Try to provide the',
     *' same number of '/
     *' tesserae on inner and outer cavity'//
     *' Or, if your specified very tough convergence criteria, allowed',
     *' number of '/
     *' iterations is not enough to reach the convergence. Increase',
     *' the number of '/
     *'iterations by using keyword ITSE=')
      END SUBROUTINE SFERA3T

!======================================================================!
      SUBROUTINE SFERT0(X0_,Y0_,Z0_,AS_,QS1,QS_,NW_,XX_,YY_,ZZ_,RV_)
***********************************************************************
*
*     THIS SUBPROGRAM CALCULATES THE SELF-POLARIZATION OF
*     SURFACE ELEMENTS
*
***********************************************************************
*
      IMPLICIT REAL(8) (A-H,O-Z)
      
      logical, save :: FIRST=.true.,SF3
      
      real(8),  intent(in),    dimension(:) :: X0_,Y0_,Z0_,AS_,QS1
      real(8),  intent(in),    dimension(:) :: XX_,YY_,ZZ_,RV_
      integer, intent(in),    dimension(:) :: NW_
      real(8),  intent(inout), dimension(:) :: QS_

      !INCLUDE 'SIZES'
      !DIMENSION X0_(*),Y0_(*),Z0_(*),AS_(*),NW_(*),QS_(*),QS1(*)
      !DIMENSION XX_(*),YY_(*),ZZ_(*),RV_(*)
      !COMMON /KEYWRD/ KEYWRD
      !COMMON /SFE   / NN
      !CHARACTER*241 KEYWRD
      !LOGICAL FIRST,SF3
      !SAVE SF3
      !DATA FIRST /.TRUE./

      IF (FIRST) THEN
          FIRST=.FALSE.
          SF3=(INDEX(KEYWRD,' SF3') .NE. 0)
C           RSELF - IF THE DISTANCE BETWEEN SURFACE ELEMENTS
C           LESS THEN RSELF (IN SFERA3) THEIR INTERACTION WILL NOT
C           BE TAKEN INTO ACCOUNT.
      ENDIF
      PI=ATAN(1.0D0)*4.0D0
      DCC2=4.D0*PI
      DCC1=PI*2.D00
      RSELF=.026D0
      QQ=0.0D0
      DO 120 IAT=1,NN
          XIT=XX_(IAT)
          YIT=YY_(IAT)
          ZIT=ZZ_(IAT)
          RIT=RV_(IAT)
          RITT=RIT*RIT
          RIT2=2.D0*RIT
          DO 110 IL=NW_(IAT),NW_(IAT+1)-1
C
             XSG=X0_(IL)
             YSG=Y0_(IL)
             ZSG=Z0_(IL)
             QS1IL=QS1(IL)
C
C      THE TESSERAE ARE ON THE SAME SPHERES
C
  70         DO 80 JL=IL+1,NW_(IAT+1)-1
                  RRR=SQRT((XSG-X0_(JL))*(XSG-X0_(JL))+(YSG-Y0_(JL))*
     *               (YSG-Y0_(JL))+(ZSG-Z0_(JL))*(ZSG-Z0_(JL)))
                  COF=1.D0/(2.D0*RIT*RRR)
C
                  QS_(JL)=QS_(JL)+COF*QS1IL
                  QS_(IL)=QS_(IL)+COF*QS1(JL)
   80        CONTINUE
C
C         SELFPOLARIZATION OF THE TESSERAE
C
             IF(AS_(IL).GT.1.D-6) THEN
               COF=DCC1* SQRT(AS_(IL)/(DCC2*RIT*RIT)) /AS_(IL)
               QS_(IL)=QS_(IL)+COF*QS1IL
             ELSEIF(ABS(AS_(IL)).LT.1.D-6) THEN
C
C            IF THE SQUARE IS 0, DO NOTHING
C
             ELSE
C
C            ERROR - NEGATIVE SQUARE
C
                WRITE (6,'('' ERROR IN SFERA3. TESSER'',I4,'' OF ATOM'',
     *          I3,'' HAS THE SQUARE S='',1P,D12.4,'' < 0'')')
     *          IL,IAT,AS_(IL)
             ENDIF
             DO 100 JAT=IAT+1,NN
                 RAT_=RV_(JAT)
                 XAT=XX_(JAT)
                 YAT=YY_(JAT)
                 ZAT=ZZ_(JAT)
                 RATT=RAT_*RAT_
                 RAT2=2.D0*RAT_
                 RTO=(XSG-XAT)*(XSG-XAT)+(YSG-YAT)*(YSG-YAT)+
     *           (ZSG-ZAT)*(ZSG-ZAT)
                 DO 60 JL=NW_(JAT),NW_(JAT+1)-1
                     XTS=X0_(JL)
                     YTS=Y0_(JL)
                     ZTS=Z0_(JL)
                     RRK=(XSG-XTS)*(XSG-XTS)+(YSG-YTS)*
     *               (YSG-YTS)+(ZSG-ZTS)*(ZSG-ZTS)
                     IF (RRK.LT.RSELF) GO  TO  60
                     RRG=RRK*SQRT(RRK)
                     RTU=(XTS-XIT)*(XTS-XIT)+(YTS-YIT)*(YTS-YIT)+
     *               (ZTS-ZIT)*(ZTS-ZIT)
                     COF=(RRK+RATT-RTO)/(RAT2*RRG)
                     COG=(RRK+RITT-RTU)/(RIT2*RRG)
C
                     QS_(JL)=QS_(JL)+COF*QS1IL
                     QS_(IL)=QS_(IL)+COG*QS1(JL)
   60            CONTINUE
  100        CONTINUE
             QQ=QQ+QS_(IL)
  110     CONTINUE
  120 CONTINUE
      IF (SF3)WRITE (6,'('' IN SFERT0 QQ='',F12.7)')QQ
      RETURN
      END SUBROUTINE SFERT0

!======================================================================!
      SUBROUTINE SFERT1(QS_,NW_,IJ_,CON_,CHDIFF_,QSFE_)
C
C    TESSERAL  CHARGE'S   COMPENSATION
C
      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: IJ_
      real(8),  intent(in) :: CON_, CHDIFF_
      integer, intent(in),    dimension(:) :: NW_
      real(8),  intent(inout), dimension(:) :: QS_
      real(8),  intent(inout), dimension(:) :: QSFE_

      ! local variables
      logical, save :: FIRST=.true.,SFERT1P

      !INCLUDE 'SIZES'
      !LOGICAL  FIRST,SFERT1P
      !CHARACTER*241 KEYWRD
      !DIMENSION NW_(*),QS_(*)
      !DIMENSION QSFE_(*)
      !COMMON /KEYWRD/ KEYWRD
      !COMMON /CHARGE/ CHARGE
      !COMMON /SFE   / NN
      !SAVE SFERT1P
      !DATA FIRST /.TRUE./

      IF ( FIRST ) THEN
          FIRST  = .FALSE.
          SFERT1P = (INDEX(KEYWRD,'SFERT1' ) .NE. 0 )
      END IF
      QQ2=0.0D0
      DO I=1,IJ_
         QQ2=QQ2+QS_(I)
      ENDDO
      IF (SFERT1P) THEN
         WRITE (6,'(/'' TOTAL CHARGE OF SOLUTE'',13X,''='',F9.5)')CHARGE
         WRITE (6,'('' SURFACE CHARGE BEFORE COMPENSATION ='',F9.5)')QQ2
      ENDIF
      CHARGE_=CHARGE*CON_
      DC=ABS(QQ2-CHARGE_)
C      IF (ABS(CHARGE_).GT.1.D0) DC=DC/ABS(CHARGE_)
      IF (DC.GT.CHDIFF_) THEN
         WRITE (6,'(//'' CATION!!!''/
     *   '' SFERT1:  CHARGE ON CAVITY TOO FAR FROM EXPECTED'')')
         WRITE (6,'('' EXPECTED CHARGE ON SURFACE'',9X,''='',F9.5)')
     *   CHARGE_
         WRITE (6,'('' SURFACE CHARGE BEFORE COMPENSATION ='',F9.5)')QQ2
         WRITE (6,'('' MOST LIKELY SOMETHING WRONG IN INPUT FILES'')')
         WRITE (6,'('' FURTHER JOB DOES NOT HAVE A SENSE. IT STOPPED'')
     *      ')
         WRITE(6,100)
         STOP
      ENDIF
      IF (ABS(CHARGE).LT.1.0D-5) THEN

C        System is considered neutral, CHARGE=0.0

         QQN=0.0D0
         DO I=1,IJ_
            IF(QS_(I).LT.0.0D0) QQN=QQN+QS_(I)
         ENDDO
         QQP=QQ2-QQN
         IF(ABS(QQ2).LE.1.D-10) GOTO 5 ! Surface charges are good enough
C
C     MERTUSH'S CHARGE COMPENSATION.
C
         ZQN=(1.D00-.5D00*QQ2/QQN)
         ZQP=(1.D00-.5D00*QQ2/QQP)
         DO I=1,IJ_
            IF(QS_(I).LT.0.D0) THEN
                QS_(I)=QS_(I)*ZQN
            ELSE
                QS_(I)=QS_(I)*ZQP
            ENDIF
         ENDDO
      ELSE
C         WRITE (*,'('' CON='',F9.5)') CON_
         IF (SFERT1P)
     *   WRITE (6,'('' EXPECTED CHARGE ON SURFACE'',9X,''='',F9.5)')
     *   CHARGE_
         AA=CHARGE_/QQ2
         DO I=1,IJ_
            QS_(I)=QS_(I)*AA
         ENDDO
      ENDIF
    5 CHARGE1=0.0D0
      DO I=1,NN
          QSFE_(I)=0.0D00
          L=NW_(I)
          M=NW_(I+1)-1
          DO JJ=L,M
             QSFE_(I)=QSFE_(I)+QS_(JJ)
          ENDDO
          CHARGE1=CHARGE1+QSFE_(I)
      ENDDO
      IF (SFERT1P)
     *WRITE (6,'('' SURFACE CHARGE AFTER COMPENSATION  ='',F9.5/)')
     *CHARGE1
      RETURN
  100 FORMAT(' Suggestions:'/
     *' It seems that something wrong with parameters specifying the',
     *' number'/
     *' and sizes of tesserae on cavity surfaces. Try to provide the',
     *' same number of '/
     *' tesserae on inner and outer cavity'/)
      END SUBROUTINE SFERT1

!======================================================================!
      SUBROUTINE SFERT2(X0_,Y0_,Z0_,QS2,QS_,NW_,XX_,YY_,ZZ_,RV_,
     *   X02,Y02,Z02,IJ2)

      IMPLICIT REAL(8) (A-H,O-Z)

      integer, intent(in) :: IJ2
      integer, intent(in), dimension(:) :: NW_
      real(8), intent(in), dimension(:) :: X0_,Y0_,Z0_,XX_,YY_,ZZ_,RV_
      real(8), intent(in), dimension(:) :: X02,Y02,Z02,QS2
      real(8), intent(inout), dimension(:) :: QS_

      ! local variables
      logical, save :: FIRST=.true.,SF3
      real(8),  save :: RSELF

      !INCLUDE 'SIZES'
      !DIMENSION X0_(*),Y0_(*),Z0_(*),NW_(*),QS_(*)
      !DIMENSION XX_(*),YY_(*),ZZ_(*),RV_(*)
      !DIMENSION X02(*),Y02(*),Z02(*),QS2(*)
      !COMMON /SFE   / NN
      !COMMON /KEYWRD/ KEYWRD
      !CHARACTER*241 KEYWRD
      !LOGICAL FIRST,SF3
      !SAVE RSELF,SF3
      !DATA FIRST /.TRUE./

      IF (FIRST) THEN

          FIRST=.FALSE.
          SF3=(INDEX(KEYWRD,' SF3') .NE. 0)
C
C           RSELF - if the distance between surface elements is
C           less then RSELF (in SFERA3) their interaction will not
C           be taken into account.
C
          RSELF=.026D0

      ENDIF

      QQ=0.0D0

      DO 120 IAT=1,NN

         XIT=XX_(IAT)
         YIT=YY_(IAT)
         ZIT=ZZ_(IAT)
         RIT=RV_(IAT)
         RITT=RIT*RIT
         RIT2=2.D0*RIT

         DO 110 IL=NW_(IAT),NW_(IAT+1)-1

            XSG=X0_(IL)
            YSG=Y0_(IL)
            ZSG=Z0_(IL)

            DO 100 JL=1,IJ2

                  XTS=X02(JL)
                  YTS=Y02(JL)
                  ZTS=Z02(JL)
                  RRK=(XSG-XTS)*(XSG-XTS)+(YSG-YTS)*
     *            (YSG-YTS)+(ZSG-ZTS)*(ZSG-ZTS)
                  RRKSR=SQRT(RRK)
                  IF (RRK.LT.RSELF) GOTO 100

C                  IF (RRK.LT.RSELF) THEN
C                     WRITE(6,'('' DISTANSE BETWEEN TESSERA SMALLER'',
C     *               ''MINIMAL IN SFERT2 (S12 AND S21 OPERATORS)'')')
C                  WRITE(6,1010)IL,JL,RRKSR
C                     GOTO 100
C                  ENDIF

                  RRG=RRK*SQRT(RRK)
                  RTU=(XTS-XIT)*(XTS-XIT)+(YTS-YIT)*(YTS-YIT)+
     *            (ZTS-ZIT)*(ZTS-ZIT)
                  CSA=(RRK+RITT-RTU)/(RIT2*RRKSR)
                  IF (DABS(CSA-1.D0).LE.1.0D-5) GOTO 100
                  COG=(RRK+RITT-RTU)/(RIT2*RRG)
                  QS_(IL)=QS_(IL)+COG*QS2(JL)
C                 WRITE(6,1000)IL,JL,RRKSR,CSA,COG,COG*QS2(JL)

  100       CONTINUE

              QQ=QQ+QS_(IL)

  110    CONTINUE

  120 CONTINUE

      IF (SF3) WRITE (6,'('' IN SFERT2 QQ='',F12.7)') QQ

 1000 FORMAT(' IL,JL,RRKSR,COG,COG*QS2(JL):',2I5,5X,4G11.5)
 1010 FORMAT(' IL,JL,RRKSR:',2I5,5X,3G11.5)

      END SUBROUTINE SFERT2

!======================================================================!
      REAL(8) FUNCTION VOLMIN(X)
      IMPLICIT REAL(8) (A-H,O-Z)

      !COMMON /SCHCON/ COSBT,RAT,COSBT1,RNE,RR,RDS,RDS1
      !COMMON /VOL/     CON,VSOLV,RSOLV,SELFCR,CHDIFF,ITSE
C
      X1=X*RNE/RAT
      DL=SQRT(DLENF(RAT+RSOLV,X,COSBT))
      DL1=SQRT(DLENF(RNE+RSOLV,X1,COSBT1))
      RDS=DL-RSOLV
      RDS1=DL1-RSOLV
      IF(RDS.LE.0.D0.OR.X.GT.RAT+RDS.OR.RDS1.LE.0.D0.OR.X1.GT.RNE+RDS)
     1THEN
         VOLMIN=0.D0
         RETURN
      ENDIF
      COSS=COSFN(RAT,X,RDS)
      COSS1=COSFN(RNE,X1,RDS1)
      SINS=SQRT(1.D0-COSS*COSS)
      SINS1=SQRT(1.D0-COSS1*COSS1)
      RA=RAT*SINS
      RA1=RNE*SINS1
      IF(X+X1+RDS+RDS1.GT.RR) THEN
         COSG=COSFN(RDS,RR-X-X1,RDS1)
         RG=RDS*SQRT(1.D0-COSG*COSG)
         HR=RDS*COSG
         HR1=RR-X-X1-HR
      ELSE
         COSG=1.D0
         RG=0.D0
         HR=RDS
         HR1=RDS1
      ENDIF
      HA=RAT*(1.D0-COSS)
      HA1=RNE*(1.D0-COSS1)
      HB=HR+X-RAT*COSS
      HB1=HR1+X1-RNE*COSS1
      VOLMIN=-((VSEC(RA,RG,HB)-VSEC(RA,0.D0,HA))+
     *(VSEC(RA1,RG,HB1)-VSEC(RA1,0.D0,HA1)))
C     WRITE (6,'('' VOLMIN: X='',1PD12.5,''  V='',D12.5)') X,-VOLMIN
C     WRITE (6,'('' VOLMIN: V='',1PD12.5)') -VOLMIN
C     WRITE (6,'('' COSG,RG'',2F10.5)') COSG,RG
C     WRITE (6,'('' DL,RDS,COSS,SINS/RA,HA,HB,V1,V2'',4F10.5/
C    *1X,5F10.5)') DL,RDS,COSS,SINS,RA,HA,HB,VSEC(RA,0.D0,HA),
C    *VSEC(RA,RG,HB)
C     WRITE (6,'('' DL1,RDS1,COSS1,SINS1/RA1,HA1,HB1,V1,V2'',4F10.5/
C    *1X,5F10.5)') DL1,RDS1,COSS1,SINS1,RA1,HA1,HB1,VSEC(RA1,0.D0,HA1),
C    *VSEC(RA1,RG,HB1)
      RETURN
      END FUNCTION VOLMIN

!======================================================================!
      REAL(8) FUNCTION COSFN(X,Y,Z)
      REAL(8) X,Y,Z
      COSFN=(X*X+Y*Y-Z*Z)/(2.0D0*X*Y)
      END FUNCTION COSFN

!======================================================================!
      REAL(8) FUNCTION DLENF(X,Y,C)
      REAL(8) X,Y,C
      DLENF=X*X+Y*Y-2.0D0*X*Y*C
      END FUNCTION DLENF

!======================================================================!
      REAL(8) FUNCTION FINT(H,R,X)
      REAL(8) H,R,X,PI
      PI=ATAN(1.0D0)*4.0D0
      FINT=PI*(H*H*X+R*R*X-X*X*X/3.0D0-H*X*SQRT(R*R-X*X)
     *-H*R*R*ASIN(X/R))
      END FUNCTION FINT

!======================================================================!
      REAL(8) FUNCTION VSEC(R1,R2,H)
      REAL(8) R1,R2,H,PI
      PI=DATAN(1.0D0)*4.0D0
      VSEC=PI/6.0D0*H*(3.0D0*R1*R1+3.0D0*R2*R2+H*H)
      END FUNCTION VSEC

!======================================================================!
      REAL(8) FUNCTION R2D(X) !transfer from radians to degrees
      REAL(8) X,PI
      PI=DATAN(1.0D0)*4.0D0
      R2D=X*180.0D0/PI
      END FUNCTION R2D

!======================================================================!
      SUBROUTINE VOLSQU
     ,   (X0_,Y0_,Z0_,AS_,XX_,YY_,ZZ_,RV_,NN_,NW_,SQC,VMOLC,PLENGTH)
C
C   Subroutine calculates square, volume and mean radius of molecule
C
C   Output:   SQC     - square
C             VMOLC   - volume
C             PLENGTH - mean radius
C
      IMPLICIT REAL(8) (A-H,O-Z)

      real(8), intent(in),  dimension(:) :: X0_,Y0_,Z0_,AS_
      real(8), intent(in),  dimension(:) :: XX_,YY_,ZZ_,RV_
      integer, intent(in), dimension(:) :: NW_
      integer, intent(in) :: NN_
      real(8), intent(out) :: SQC,VMOLC,PLENGTH

      !INCLUDE 'SIZES'
      !DIMENSION X0_(*),Y0_(*),Z0_(*),AS_(*)
      !DIMENSION XX_(*),YY_(*),ZZ_(*),RV_(*),NW_(*)
      !COMMON /SOLMAEL/ X0EL(NS),Y0EL(NS),Z0EL(NS),ASEL(NS),
      !1QSEL(NS),QCOREL(NS)

      PI=DATAN(1.0D0)*4.0D0
      SQC=0.d0
C   Calculation of  cavity surface square
C   Loop over all atoms of solute molecule
      DO 3 I=1,NN_
         S=0.d0
         DO 4 J=NW_(I),NW_(I+1)-1
    4       S=S+AS_(J)
         SQC=SQC+S
    3 CONTINUE
C
C   Calculation of  cavity volume
C
      VMOLC=0.d0
      I=1
C  Coordinates of atom I of soluted molecule
C
      XI=XX_(I)
      YI=YY_(I)
      ZI=ZZ_(I)
C
C   Loop over all atoms of solute molecule
C
      DO 1 K=1,NN_
           XK=XX_(K)
           YK=YY_(K)
           ZK=ZZ_(K)
           RK=RV_(K)
C     WRITE (6,'('' RK='',D10.5)') RK
C     P2 - square of the distance between the nuclei
           P2=(XI-XK)*(XI-XK)+(YI-YK)*(YI-YK)+(ZI-ZK)*(ZI-ZK)
C
C     Loop over all tessera of surface of atom K
C
           DO 2 J=NW_(K),NW_(K+1)-1
C     Coordinates of tessera
                XJ=X0_(J)
                YJ=Y0_(J)
                ZJ=Z0_(J)
C     WRITE (6,'(1X,''XJ,YJ,ZJ  '',3D15.10)')XJ,YJ,ZJ
                SJ=AS_(J)
                RO=SQRT((XI-XJ)*(XI-XJ)+(YI-YJ)*(YI-YJ)+(ZI-ZJ)*(ZI-ZJ))
C     WRITE (6,'('' RO='',D10.5)') RO
C     Cosinus gamma
                COM=(RO*RO+RK*RK-P2)/(2.D0*RK*RO)
                VMOLC=VMOLC+COM*RO*SJ/3.D0
    2      CONTINUE
    1 CONTINUE
      PLENGTH=(VMOLC*3.D0/4.D0/PI)**(1.D0/3.D0)
      RETURN
      END SUBROUTINE VOLSQU

!======================================================================!
      SUBROUTINE WRTKEY(KEYWRD_)
***********************************************************************
*
*  WRTKEY CHECKS ALL KEY-WORDS AND PRINTS THOSE IT RECOGNIZES.  IF IT
*  FINDS A WORD IT DOES NOT RECOGNIZE THE PROGRAM WILL BE STOPPED.
*
***********************************************************************
      IMPLICIT REAL(8) (A-H,O-Z)
      
      CHARACTER(len=*), intent(in) :: KEYWRD_
      
      character(241) :: ALLKEY
      character(  1) :: CH
      logical :: TP

      !INCLUDE 'SIZES'
      !CHARACTER*241 KEYWRD_, ALLKEY
      !CHARACTER CH*1
      !LOGICAL MYWORD,TP

      ALLKEY=KEYWRD_

C     DUMMY IF STATEMENT TO REMOVE AMPERSAND AND PLUS SIGNS, IF PRESENT

      IF(MYWORD(ALLKEY(160:),' SETUP'))I=1
      IF(MYWORD(ALLKEY,'&'))I=2
      IF(MYWORD(ALLKEY,' +'))I=3

      IF (MYWORD(ALLKEY,'TIMES') )WRITE(6,280)
  280 FORMAT(' *  TIMES    - TIMES OF VARIOUS STAGES TO BE PRINTED')
      IF (MYWORD(ALLKEY,'LARGE') ) WRITE(6,340)
  340 FORMAT(' *  LARGE    - EXPANDED OUTPUT TO BE PRINTED')
      IF (MYWORD(ALLKEY,' XYZ') ) WRITE(6,580)
  580 FORMAT(' *   XYZ     - CARTESIAN COORDINATE SYSTEM TO BE USED')
      IF (MYWORD(ALLKEY,'ECHO') ) WRITE(6,600)
  600 FORMAT(' *  ECHO     - ALL INPUT DATA TO BE ECHOED BEFORE RUN')

      IF (MYWORD(ALLKEY,'DEBUG ') ) WRITE(6,660)
  660 FORMAT(' *  DEBUG    - DEBUG OPTION TURNED ON')
      IF (MYWORD(ALLKEY,'CHARGE') )
     1 WRITE(6,850)NINT(READA(KEYWRD_,INDEX(KEYWRD_,'CHARGE')))
  850 FORMAT(1(' *',/),' *',15X,'  CHARGE ON SYSTEM =',I3,1(/,' *'))

      IF (MYWORD(ALLKEY,'PREC') ) WRITE(6,1110)
 1110 FORMAT(' *  PRECISE  - CRITERIA TO BE INCREASED BY 100 TIMES')
      IF (MYWORD(ALLKEY,'NOINTER') ) WRITE(6,1130)
 1130 FORMAT(' *  NOINTER  - INTERATOMIC DISTANCES NOT TO BE PRINTED')
      IF (MYWORD(ALLKEY,'SCFCRT') ) WRITE(6,1180)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'SCFCRT'))
 1180 FORMAT(' *  SCFCRT   - DEFAULT SCF CRITERION REPLACED BY',G12.3)
      IF (MYWORD(ALLKEY,'NOXYZ') ) WRITE(6,1200)
 1200 FORMAT(' *  NOXYZ    - CARTESIAN COORDINATES NOT TO BE PRINTED')
      IF (MYWORD(ALLKEY,'BYPASS') ) WRITE(6,1201)
 1201 FORMAT(' *  BYPASS   - no cartesian to internal conversion')
C
C     SOLVATATION KEY WORDS
C
      IF (MYWORD(ALLKEY,'RADIUS') ) WRITE(6,2190)
 2190 FORMAT(' *  RADIUS   - USER''S VAN-DER-WAALS RADII TO BE USED')
      IF (MYWORD(ALLKEY,'KAPPA')  )WRITE(6,2110)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'KAPPA'))
 2110 FORMAT(' *  KAPPA=   - FACTOR FOR FIRST RADII OF ATOMS IS ',F5.2)
      IF (MYWORD(ALLKEY,'DELTA')  )WRITE(6,2111)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'DELTA'))
 2111 FORMAT(' *  DELTA=  - ADDITIONAL FACTOR FOR THE 2D RADII IS ',
     *   F5.2)
      IF (MYWORD(ALLKEY,'EPS')  )WRITE(6,2010)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'EPS'))
 2010 FORMAT(' *  EPS=     - DIELECTRIC PERMITTIVITY OF SOLVENT',
     1' EPS=',F6.2)
      IF (MYWORD(ALLKEY,'EPSEL')  )WRITE(6,2020)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'EPSEL'))
 2020 FORMAT(' *  EPSEL=   - OPTICAL PERMITTIVITY OF SOLVENT ',
     1 'EPSEL=',F6.3)
      IF (MYWORD(ALLKEY,'EXVOL')  )WRITE(6,2030)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'EXVOL'))
 2030 FORMAT(' *  EXVOL=   - MINIMAL UNREACHABLE VOLUME TO ',
     1'INCLUDE ADDITIONAL SPHERES IS',G10.3)
      IF (MYWORD(ALLKEY,'SOLRD')  )WRITE(6,2140)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'SOLRD'))
 2140 FORMAT(' *  SOLRD=   - EFFECTIVE RADIUS OF THE SOLVENT MOLECULE ',
     1'IS ',F5.2)
      TP=MYWORD(ALLKEY,'MODFE')
      IF (TP) THEN
         J=INDEX(KEYWRD_,'MODFE=(')
         IF(J.NE.0)THEN
            WRITE(6,2059)INT(READA(KEYWRD_,INDEX(KEYWRD_,
     1         'MODFE=(')+6)),
     2         INT(READA(KEYWRD_,INDEX(KEYWRD_,'MODFE=(')+9))
         ELSE
            WRITE(6,2060)INT(READA(KEYWRD_,INDEX(KEYWRD_,'MODFE=')+6))
         ENDIF
      ENDIF
 2059 FORMAT(' *  MODFE=   - VALUE OF TESSERAE SIZE FOR 1ST CAVITY',
     *' IS',I3/15X,'VALUE OF TESSERAE SIZE FOR  2D CAVITY IS',I3)
 2060 FORMAT(' *  MODFE=   - VALUE OF TESSERAE SIZE FOR BOTH CAVITIES',
     *' IS ',I3)
      IF (MYWORD(ALLKEY,'MERT')  )WRITE(6,2050)
     1 NINT(READA(KEYWRD_,INDEX(KEYWRD_,'MERT')))
 2050 FORMAT(' *  MERT=    - SIMPLIFIED METHOD FOR TESSERAE',
     1 ' PREPARING TO BE USED, N=',I3)
      IF (MYWORD(ALLKEY,'NOSMOOTH') ) WRITE(6,2090)
 2090 FORMAT(' *  NOSMOOTH - SURFACE SMOOTHING NOT TO BE USED')
      IF (MYWORD(ALLKEY,' SMOOTH') ) WRITE(6,2130)
 2130 FORMAT(' *  _SMOOTH  - SURFACE SMOOTHING TO BE USED')
      IF (MYWORD(ALLKEY,'NSNN')  )WRITE(6,2100)
     1 NINT(READA(KEYWRD_,INDEX(KEYWRD_,'NSNN')))
 2100 FORMAT(' *  NSNN=    - THE NUMBER OF SPHERES IS ',I3)
      IF (MYWORD(ALLKEY,'ITSE')  )WRITE(6,2040)
     1 NINT(READA(KEYWRD_,INDEX(KEYWRD_,'ITSE')))
 2040 FORMAT(' *  ITSE=    - DO A MAXIMUM OF',I3,' ITERATIONS IN ',
     1'SELFPOLARIZATION')
      IF (MYWORD(ALLKEY,'SELFCR')  )WRITE(6,2120)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'SELFCR'))*1.D-4
 2120 FORMAT(' *  SELFCR=  - DEFAULT SELFPOLARIZATION CRITERION',
     1' REPLACED BY ',G10.3)
      IF (MYWORD(ALLKEY,'CHDIFF')  )WRITE(6,2121)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'CHDIFF'))*1.D0
 2121 FORMAT(' *  CHDIFF=  - DEFAULT MAX. COMPENSATION CRITERION',
     1' REPLACED BY ',G10.3)
      IF (MYWORD(ALLKEY,'SOLCRT')  )WRITE(6,2230)
     1 READA(KEYWRD_,INDEX(KEYWRD_,'SOLCRT'))
 2230 FORMAT(' *  SOLCRT   - DEFAULT CRITERION OF MEDIUM AND SOLUTE',
     1 'POLARIZATION REPLACED BY',G10.3)
      IF (MYWORD(ALLKEY,'NOCOMP') ) WRITE(6,2200)
 2200 FORMAT(' *  NOCOMP   - SURFACE CHARGES NORMALIZATION TO BE',
     1 ' SUPPRESSED')

      TP=MYWORD(ALLKEY,'NTETFI')
      IF (TP) THEN
         J=INDEX(KEYWRD_,'NTETFI=(')
         IF(J.NE.0)THEN
            WRITE(6,2219)INT(READA(KEYWRD_,INDEX(KEYWRD_,
     1         'NTETFI=(')+7)),
     2         INT(READA(KEYWRD_,INDEX(KEYWRD_,'NTETFI=(')+10))
         ELSE
            WRITE(6,2220)INT(READA(KEYWRD_,INDEX(KEYWRD_,'NTETFI=')+7))
         ENDIF
      ENDIF
 2219 FORMAT(' *  NTETFI   - THE NUMBER OF SECTORS WRT THE POLAR',
     1'ANGLES FOR 1ST CAVITY IS',I3/15X,
     2'THE NUMBER OF SECTORS WRT THE POLAR ANGLES FOR 2D CAVITY IS',I3)
 2220 FORMAT(' *  NTETFI   - THE NUMBER OF SECTORS WRT THE POLAR',
     1 'ANGLES FOR BOTH CAVITIES IS',I3)

      IF (MYWORD(ALLKEY,'S12DR') ) WRITE(6,2300)
 2300 FORMAT(' *  S12DR   - THE INFLUENCE OF CHARGES ON SECOND SURFACE',
     *' ON CHARGES'/14X,'ON THE FIRST  SURFACE IS NEGLECTED')
C
C         END OF SOLVATATION KEY WORDS
C
      IF (ALLKEY.NE.' ') THEN

         J=0
         DO 50 I=1,240
            IF(ALLKEY(I:I).NE.' '.OR.ALLKEY(I:I+1).NE.'  ')THEN
               J=J+1
               CH=ALLKEY(I:I)
               ALLKEY(J:J)=CH
            ENDIF
   50    CONTINUE

         IF(ALLKEY(241:241).NE.' ')THEN
            J=J+1
            CH=ALLKEY(241:241)
            ALLKEY(J:J)=CH
         ENDIF

         J=MAX(1,J)
         L=INDEX(KEYWRD_,'DEBUG')

         IF(L.NE.0)THEN
            WRITE(6,'('' *  DEBUG KEYWORDS USED:  '',A)')ALLKEY(:J)
         ELSE
            WRITE(6,'(///10X,''WARNING: UNKNOWN KEYWORDS: ('',A,'')'')')
     1      ALLKEY(:J)
C-----------WRITE(6,'(///10X,''CALCULATION STOPPED TO AVOID WASTING TIME
C-----1     .'')')
            WRITE(6,'(///10X,''IF THESE ARE DEBUG KEYWORDS, ADD THE KEYW
     1      ORD "DEBUG"'')')
C-----------STOP
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE WRTKEY

!======================================================================!
      SUBROUTINE SEALINV(X,FX,Y,FY,FTIN,STEP_IN)
      IMPLICIT REAL(8) (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      STEP = STEP_IN
      FTOL=FTIN
      IP=1
      NP=6
      IEXIT=0
      NTOL=0
      FTOL2=FTOL/100.D0
      FA=FX
      FB=FX
      FC=FX
      DA=0.D0
      DB=0.D0
      DC=0.D0
      K=-2
      M=0
      D=STEP
    1 Y=X+D
      F=VOLMIN(Y)
      K=K+1
      IF(F-FA)5,3,6
    3 Y=X+DA
      IF(IP.EQ.1) PRINT 2100
      FY=FA
 2100 FORMAT(' SEARCH FAILED. FUNCTION ',
     *'VALUE INDEPENDENT OF SEARCH DIRECTION')
      GOTO 326
    5 FC=FB
      FB=FA
      FA=F
      DC=DB
      DB=DA
      DA=D
      D=2.D0*D+STEP
      GOTO 1
    6 IF(K) 7,8,9
    7 FB=F
      DB=D
      D=-D
      STEP=-STEP
      GOTO 1
    8 FC=FB
      FB=FA
      DC=DB
      FA=F
      DB=DA
      DA=D
      GOTO 21
    9 DC=DB
      DB=DA
      DA=D
      FC=FB
      FB=FA
      FA=F
   10 D=0.5D0*(DA+DB)
      Y=X+D
      F=VOLMIN(Y)
   12 IF((DC-D)*(D-DB)) 15,13,18
   13 Y=X+DB
      FY=FB
      IF(IEXIT.EQ.1) GOTO 32
      IF(IP.EQ.1) WRITE(NP,2200)
 2200 FORMAT(' SEARCH FAILED. LOCATION OF',
     *' MINIMUM IS LIMITED BY ROUNDING')
      GOTO 325
   15 IF(F-FB) 16,13,17
   16 FC=FB
      FB=F
      DC=DB
      DB=D
      GOTO 21
   17 FA=F
      DA=D
      GOTO 21
   18 IF(F-FB) 19,13,20
   19 FA=FB
      FB=F
      DA=DB
      DB=D
      GOTO 21
   20 FC=F
      DC=D
   21 A=FA*(DB-DC)+FB*(DC-DA)+FC*(DA-DB)
      IF(A) 22,30,22
   22 D=0.5D0*((DB*DB-DC*DC)*FA+(DC*DC-
     *DA*DA)*FB+(DA*DA-DB*DB)*FC)/A
      IF((DA-D)*(D-DC)) 13,13,23
   23 Y=X+D
      F=VOLMIN(Y)
      IF(DABS(FB)-FTOL2)25,25,26
   25 A=1.D0
      GOTO 27
   26 A=1.D0/FB
   27 IF(DABS((FB-F)*A)-FTOL) 28,28,12
   28 IEXIT=1
      IF(F-FB) 29,13,13
   29 FY=F
      GOTO 32
   30 IF(M) 31,31,13
   31 M=M+1
      GOTO 10
   32 IF(Y.NE.X) GOTO 325
      GOTO 33
  325 IF(NTOL.NE.0.AND.IP.EQ.1)
     *WRITE(NP,3000) NTOL
 3000 FORMAT(' TOLERANCE REDUCED ',I1,
     *' TIME(S)')
  326 IF(FY.LT.FX) RETURN
      WRITE(NP,5000)
 5000 FORMAT(' SEARCH FAILED.  JOB ',
     *' TERMINATED')
      RETURN
   33 IF(NTOL.EQ.3) GOTO 34
      IEXIT=0
      NTOL=NTOL+1
      FTOL=FTOL/10.D0
      GOTO 12
   34 IF(IP.EQ.1) WRITE(NP,2000)
 2000 FORMAT(' A POINT BETTER THAN ',
     *' ENTERING POINT CANNOT BE FOUND')
      RETURN
      END SUBROUTINE SEALINV

      end module frcm
