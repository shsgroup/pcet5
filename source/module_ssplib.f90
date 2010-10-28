module ssplib

!---------------------------------------------------------------------
! Some routines from SSPLIB
!---------------------------------------------------------------------
!
!  $Author: souda $
!  $Date: 2010-10-28 21:29:36 $
!  $Revision: 5.2 $
!  $Log: not supported by cvs2svn $
!
!---------------------------------------------------------------------

   implicit none

contains
                                                                      
      SUBROUTINE DQH64(FCT,Y)                                           

      !..................................................................
      !                                                                       
      !   PURPOSE                                                        
      !      TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM    
      !                          -INFINITY TO +INFINITY).                
      !                                                                       
      !   USAGE                                                          
      !      CALL DQH64 (FCT,Y)                                          
      !      PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                
      !                                                                       
      !   DESCRIPTION OF PARAMETERS                                      
      !      FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  
      !               SUBPROGRAM USED.                                   
      !      Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     
      !                                                                       
      !   REMARKS                                                        
      !      NONE                                                        
      !                                                                       
      !   SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
      !      THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    
      !      MUST BE FURNISHED BY THE USER.                              
      !                                                                       
      !   METHOD                                                         
      !      EVALUATION IS DONE BY MEANS OF 64-POINT GAUSSIAN-HERMITE    
      !      QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       
      !      FCT(X) IS A POLYNOMIAL UP TO DEGREE 127.                    
      !      FOR REFERENCE, SEE                                          
      !      SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    
      !      CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     
      !      GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       
      !      TR00.1100 (MARCH 1964), PP.213-214.                         
      !                                                                       
      !..................................................................

      real(kind=8) :: X, Y, Z
      real(kind=8), external :: FCT                                        

      X=.10526123167960546D2                                            
      Z=-X                                                              
      Y=.55357065358569428D-48*(FCT(X)+FCT(Z))                          
      X=.9895287586829539D1                                             
      Z=-X                                                              
      Y=Y+.16797479901081592D-42*(FCT(X)+FCT(Z))                        
      X=.9373159549646721D1                                             
      Z=-X                                                              
      Y=Y+.34211380112557405D-38*(FCT(X)+FCT(Z))                        
      X=.8907249099964770D1                                             
      Z=-X                                                              
      Y=Y+.15573906246297638D-34*(FCT(X)+FCT(Z))                        
      X=.8477529083379863D1                                             
      Z=-X                                                              
      Y=Y+.25496608991129993D-31*(FCT(X)+FCT(Z))                        
      X=.8073687285010225D1                                             
      Z=-X                                                              
      Y=Y+.19291035954649669D-28*(FCT(X)+FCT(Z))                        
      X=.7689540164040497D1                                             
      Z=-X                                                              
      Y=Y+.7861797788925910D-26*(FCT(X)+FCT(Z))                         
      X=.7321013032780949D1                                             
      Z=-X                                                              
      Y=Y+.19117068833006428D-23*(FCT(X)+FCT(Z))                        
      X=.69652411205511075D1                                            
      Z=-X                                                              
      Y=Y+.29828627842798512D-21*(FCT(X)+FCT(Z))                        
      X=.66201122626360274D1                                            
      Z=-X                                                              
      Y=Y+.31522545665037814D-19*(FCT(X)+FCT(Z))                        
      X=.62840112287748282D1                                            
      Z=-X                                                              
      Y=Y+.23518847106758191D-17*(FCT(X)+FCT(Z))                        
      X=.59556663267994860D1                                            
      Z=-X                                                              
      Y=Y+.12800933913224380D-15*(FCT(X)+FCT(Z))                        
      X=.56340521643499721D1                                            
      Z=-X                                                              
      Y=Y+.52186237265908475D-14*(FCT(X)+FCT(Z))                        
      X=.53183252246332709D1                                            
      Z=-X                                                              
      Y=Y+.16283407307097204D-12*(FCT(X)+FCT(Z))                        
      X=.50077796021987682D1                                            
      Z=-X                                                              
      Y=Y+.39591777669477239D-11*(FCT(X)+FCT(Z))                        
      X=.47018156474074998D1                                            
      Z=-X                                                              
      Y=Y+.7615217250145451D-10*(FCT(X)+FCT(Z))                         
      X=.43999171682281376D1                                            
      Z=-X                                                              
      Y=Y+.11736167423215493D-8*(FCT(X)+FCT(Z))                         
      X=.41016344745666567D1                                            
      Z=-X                                                              
      Y=Y+.14651253164761094D-7*(FCT(X)+FCT(Z))                         
      X=.38065715139453605D1                                            
      Z=-X                                                              
      Y=Y+.14955329367272471D-6*(FCT(X)+FCT(Z))                         
      X=.35143759357409062D1                                            
      Z=-X                                                              
      Y=Y+.12583402510311846D-5*(FCT(X)+FCT(Z))                         
      X=.32247312919920357D1                                            
      Z=-X                                                              
      Y=Y+.8788499230850359D-5*(FCT(X)+FCT(Z))                          
      X=.29373508230046218D1                                            
      Z=-X                                                              
      Y=Y+.51259291357862747D-4*(FCT(X)+FCT(Z))                         
      X=.26519724354306350D1                                            
      Z=-X                                                              
      Y=Y+.25098369851306249D-3*(FCT(X)+FCT(Z))                         
      X=.23683545886324014D1                                            
      Z=-X                                                              
      Y=Y+.10363290995075777D-2*(FCT(X)+FCT(Z))                         
      X=.20862728798817620D1                                            
      Z=-X                                                              
      Y=Y+.36225869785344588D-2*(FCT(X)+FCT(Z))                         
      X=.18055171714655449D1                                            
      Z=-X                                                              
      Y=Y+.10756040509879137D-1*(FCT(X)+FCT(Z))                         
      X=.15258891402098637D1                                            
      Z=-X                                                              
      Y=Y+.27203128953688918D-1*(FCT(X)+FCT(Z))                         
      X=.12472001569431179D1                                            
      Z=-X                                                              
      Y=Y+.58739981964099435D-1*(FCT(X)+FCT(Z))                         
      X=.9692694230711780D0                                             
      Z=-X                                                              
      Y=Y+.10849834930618684D0*(FCT(X)+FCT(Z))                          
      X=.69192230581004458D0                                            
      Z=-X                                                              
      Y=Y+.17168584234908370D0*(FCT(X)+FCT(Z))                          
      X=.41498882412107868D0                                            
      Z=-X                                                              
      Y=Y+.23299478606267805D0*(FCT(X)+FCT(Z))                          
      X=.13830224498700972D0                                            
      Z=-X                                                              
      Y=Y+.27137742494130398D0*(FCT(X)+FCT(Z))                          
      RETURN                                                            
      END subroutine dqh64

                                                                       
      SUBROUTINE DQSF(H,Y,Z,NDIM)                                       
                                                                       
      !..................................................................
      !                                                                       
      !   PURPOSE                                                        
      !      TO COMPUTE THE VECTOR OF INTEGRAL VALUES FOR A GIVEN        
      !      EQUIDISTANT TABLE OF FUNCTION VALUES.                       
      !                                                                  
      !   USAGE                                                          
      !      CALL DQSF (H,Y,Z,NDIM)                                      
      !                                                                  
      !   DESCRIPTION OF PARAMETERS                                      
      !      H      - DOUBLE PRECISION INCREMENT OF ARGUMENT VALUES.     
      !      Y      - DOUBLE PRECISION INPUT VECTOR OF FUNCTION VALUES.  
      !      Z      - RESULTING DOUBLE PRECISION VECTOR OF INTEGRAL      
      !               VALUES. Z MAY BE IDENTICAL WITH Y.                 
      !      NDIM   - THE DIMENSION OF VECTORS Y AND Z.                  
      !                                                                  
      !   REMARKS                                                        
      !      NO ACTION IN CASE NDIM LESS THAN 3.                         
      !                                                                  
      !   SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
      !      NONE                                                        
      !                                                                  
      !   METHOD                                                         
      !      BEGINNING WITH Z(1)=0, EVALUATION OF VECTOR Z IS DONE BY    
      !      MEANS OF SIMPSONS RULE TOGETHER WITH NEWTONS 3/8 RULE OR A  
      !      COMBINATION OF THESE TWO RULES. TRUNCATION ERROR IS OF      
      !      ORDER H**5 (I.E. FOURTH ORDER METHOD). ONLY IN CASE NDIM=3  
      !      TRUNCATION ERROR OF Z(2) IS OF ORDER H**4.                  
      !      FOR REFERENCE, SEE                                          
      !      (1) F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,     
      !          MCGRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP.71-76.   
      !      (2) R.ZURMUEHL, PRAKTISCHE MATHEMATIK FUER INGENIEURE UND   
      !          PHYSIKER, SPRINGER, BERLIN/GOETTINGEN/HEIDELBERG, 1963, 
      !          PP.214-221.                                             
      !                                                                    
      !..................................................................

      integer, intent(in) :: NDIM
      real(8), intent(in) :: H
      real(8), dimension(NDIM), intent(in)  :: Y
      real(8), dimension(NDIM), intent(out) :: Z

      ! local variables
      integer :: I
      real(8) :: HT, SUM1, SUM2, AUX, AUX1, AUX2                 

      HT = 0.33333333333333333d0*H                                         
      IF (NDIM-5) 7, 8, 1

      ! NDIM IS GREATER THAN 5. PREPARATIONS OF INTEGRATION LOOP          

    1 SUM1=Y(2)+Y(2)                                                    
      SUM1=SUM1+SUM1                                                    
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          
      AUX1=Y(4)+Y(4)                                                    
      AUX1=AUX1+AUX1                                                    
      AUX1=SUM1+HT*(Y(3)+AUX1+Y(5))                                     
      AUX2=HT*(Y(1)+3.875d0*(Y(2)+Y(5))+2.625d0*(Y(3)+Y(4))+Y(6))       
      SUM2=Y(5)+Y(5)                                                    
      SUM2=SUM2+SUM2                                                    
      SUM2=AUX2-HT*(Y(4)+SUM2+Y(6))                                     
      Z(1)=0.d0                                                         
      AUX=Y(3)+Y(3)                                                     
      AUX=AUX+AUX                                                       
      Z(2)=SUM2-HT*(Y(2)+AUX+Y(4))                                      
      Z(3)=SUM1                                                         
      Z(4)=SUM2                                                         

      IF(NDIM-6) 5, 5, 2                                                   

      ! INTEGRATION LOOP                                                  

    2 DO 4 I=7,NDIM,2
         SUM1=AUX1                                                         
         SUM2=AUX2                                                         
         AUX1=Y(I-1)+Y(I-1)                                                
         AUX1=AUX1+AUX1                                                    
         AUX1=SUM1+HT*(Y(I-2)+AUX1+Y(I))                                   
         Z(I-2)=SUM1                                                       
         IF(I-NDIM) 3,6,6                                                   
    3    AUX2=Y(I)+Y(I)                                                    
         AUX2=AUX2+AUX2                                                    
         AUX2=SUM2+HT*(Y(I-1)+AUX2+Y(I+1))                                 
    4    Z(I-1)=SUM2                                                       

    5 Z(NDIM-1)=AUX1                                                    
      Z(NDIM)=AUX2                                                      
      RETURN                                                            

    6 Z(NDIM-1)=SUM2                                                    
      Z(NDIM)=AUX1                                                      
      RETURN                                                            

      ! END OF INTEGRATION LOOP

    7 IF(NDIM-3) 12, 11, 8

      ! NDIM IS EQUAL TO 4 OR 5                                           

    8 SUM2=1.125d0*HT*(Y(1)+Y(2)+Y(2)+Y(2)+Y(3)+Y(3)+Y(3)+Y(4))         
      SUM1=Y(2)+Y(2)                                                    
      SUM1=SUM1+SUM1                                                    
      SUM1=HT*(Y(1)+SUM1+Y(3))                                          
      Z(1)=0.d0                                                         
      AUX1=Y(3)+Y(3)                                                    
      AUX1=AUX1+AUX1                                                    
      Z(2)=SUM2-HT*(Y(2)+AUX1+Y(4))                                     
      IF(NDIM-5)10,9,9                                                  
    9 AUX1=Y(4)+Y(4)                                                    
      AUX1=AUX1+AUX1                                                    
      Z(5)=SUM1+HT*(Y(3)+AUX1+Y(5))                                     
   10 Z(3)=SUM1                                                         
      Z(4)=SUM2                                                         
      RETURN                                                            

      ! NDIM IS EQUAL TO 3                                                

   11 SUM1=HT*(1.25d0*Y(1)+Y(2)+Y(2)-.25d0*Y(3))                        
      SUM2=Y(2)+Y(2)                                                    
      SUM2=SUM2+SUM2                                                    
      Z(3)=HT*(Y(1)+SUM2+Y(3))                                          
      Z(1)=0.d0                                                         
      Z(2)=SUM1                                                         
   12 RETURN                                                            

      END subroutine DQSF

end module ssplib
