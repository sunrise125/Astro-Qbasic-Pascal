  program sol_equinox
!     program  equiprim2.pas ----- 04 Gen 2021 ------
!
!     gfortran equprim2.f90 2020.for jpl_430.f sky9.for -o equprim2.exe
!      
      IMPLICIT NONE

      INTEGER year,IY,MO,ID,IH,IM,FINAL_,EOP_final,M1,M2,M3,M4
      DOUBLE PRECISION L1,L2,L3,L4,SEC,TDB,GELOAPD
      CHARACTER*1 NEW
   
5     FINAL_ = 59555     ! 2021/12/07   reported in MJD IERS FINAL_ data (for update) 

      EOP_final = 59152  ! 2020/10/30 jump from IERS EOP14 C04 to IERS FINAL_ data (IAU2000). 


      WRITE (*,*) " ********************************************************* "       
      WRITE (*,*) " *         ---------   EQUINOX   ----------              * "
      WRITE (*,*) " * QUESTO PROGRAMMA CALCOLA I TEMPI DEI SOLSTIZI E DEGLI * "
      WRITE (*,*) " *   EQUINOZI RELATIVI AGLI ANNI POSTERIORI AL 1962      * "
      WRITE (*,*) " *    Allo scopo viene usato il database JPL DE430t      * "
      WRITE (*,*) " *               e software SOFA.                        * "
      WRITE (*,*) " *  gfortran equinox.f90 2020.for jpl_430.f sky9.for     * "
      WRITE (*,*) " ********************************************************* " 
      WRITE (*,*) "       Data e Tempo  UTC"    
      WRITE (*,*) "   =============================="
      WRITE (*,*)
      WRITE (*,*) " Input : Anno  " 
      READ *, year    

      L1 = 360D0
      L2 = 90D0
      L3 = 180D0
      L4 = 270D0
      M1 = 3
      M2 = 6
      M3 = 9
      M4 = 12 

!      CALL equinox( year,FINAL_,EOP_final,L1,M1,IY,MO,ID,IH,IM,SEC,TDB,GELOAPD) 
!     WRITE(*,*)"  ===================================================================================================="
!     WRITE(*,10)IY,MO,ID,IH,IM,SEC,TDB,GELOAPD
!     WRITE(*,*)"  ===================================================================================================="
!10   FORMAT(2X,"-- Equinozio di primavera ",4x,I4,"/",I2,"/",I2,2x,I2,":",I2,":",F4.1,x,"UTC /",2x,"JD",x,F16.7," TDB",2x, &
!                  "Long ",F11.7," gradi") 
  
      CALL equinox( year,FINAL_,EOP_final,L2,M2,IY,MO,ID,IH,IM,SEC,TDB,GELOAPD)
     WRITE(*,20)IY,MO,ID,IH,IM,SEC,TDB,GELOAPD
     WRITE(*,*)"  ===================================================================================================="
20   FORMAT(2X,"-- Solstizio d'estate     ",4x,I4,"/",I2,"/",I2,2x,I2,":",I2,":",F4.1,x,"UTC /",2x,"JD",x,F16.7," TDB",2x, &
                  "Long ",F11.7," gradi")

      CALL equinox( year,FINAL_,EOP_final,L3,M3,IY,MO,ID,IH,IM,SEC,TDB,GELOAPD)
     WRITE(*,30)IY,MO,ID,IH,IM,SEC,TDB,GELOAPD
     WRITE(*,*)"  ===================================================================================================="
30   FORMAT(2X,"-- Equinozio d'autunno    ",4x,I4,"/",I2,"/",I2,2x,I2,":",I2,":",F4.1,x,"UTC /",2x,"JD",x,F16.7," TDB",2x, &
                  "Long ",F11.7," gradi") 
 
      CALL equinox( year,FINAL_,EOP_final,L4,M4,IY,MO,ID,IH,IM,SEC,TDB,GELOAPD)
     WRITE(*,40)IY,MO,ID,IH,IM,SEC,TDB,GELOAPD
     WRITE(*,*)"  ===================================================================================================="
40   FORMAT(2X,"-- Solstizio d'inverno    ",4x,I4,"/",I2,"/",I2,2x,I2,":",I2,":",F4.1,x,"UTC /",2x,"JD",x,F16.7," TDB",2x, &
                  "Long ",F11.7," gradi") 
     WRITE (*,*)
  
      WRITE(*,*)
      WRITE(*,*)" ***************************************************" 
      WRITE(*,*)" Per continuare  digita 'Y' oppure  'N' per uscire. "
      WRITE(*,*)" ***************************************************" 
      READ *, NEW
      IF(NEW == "y")NEW ="Y"
      IF(NEW == "n")NEW ="N" 
             
      IF ( NEW == "Y" .or. NEW == "y") THEN
        PRINT*, "  SCELTA ... ",NEW," ===> CONTINUA " 
        GOTO 5
       ELSE IF (NEW == "N" .or. NEW == "n" ) THEN
        PRINT*, "  SCELTA ... ",NEW," ===> TERMINA " 
        GOTO 50
      END IF

50    CONTINUE
  end program sol_equinox

      subroutine equinox ( year,FINAL_,EOP_final,LG,MO,IY,M0,ID,IH,IM,SEC,ET,GELOAPD) 

      IMPLICIT NONE
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )      
! Arcseconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER ( AS2R = 4.848136811095359935899141D-6 )
! Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 ) 

      INTEGER IY, MO, ID, J, I,NTARG,NCENT,ISTEP,NT,FINAL_
      INTEGER IH, IM, K, JJ,IHMSF(4),NSTEP,N,year,M0,EOP_final
      DOUBLE PRECISION OBL,AU,GELO,GELA,OBL6,DATE_,UTC1, UTC2,UTC,UT,SEC,MJD,DTR 
      DOUBLE PRECISION DATE, TIME,angledif,MU,RB(3,3),GMB(3,2),PVE(3,2),FD,XYZ(3)
      DOUBLE PRECISION RMASSE,DJMJD0,RH,TC,WL,BODY_G(3),BODY_H(3),TAI1, TAI2,LG,AE        
      DOUBLE PRECISION GAST,LONG ,LAT,HM,julian,DSQRT,T0,ET,XET,OB(3),OV(3),GP(3)
      DOUBLE PRECISION DS, AUDAY, JD,C, DT, MM,GV(3),EB(3),EPSA,GN,de_S,VALS(400)       
      DOUBLE PRECISION MUC,BP(3),HE0(3,2),SUNB(3,2),EPV(3,2),GPB(3),u0(3),V1(3),U,V
      DOUBLE PRECISION THD,DET,HD,GD,HMOD,DJM0,DJM,UT11,UT12,HPE(3),BVE(3),BEP(3),EO
      DOUBLE PRECISION PV(3,2),PVH(3,2),PVB(3,2),BEV(3,2),PHPA,VC,iau_DTDB,TT1,TT2,MO1              
      DOUBLE PRECISION XPO,YPO,DUT,LOD,DX0,DY0,ut1_notid, TDB1, TDB2,TDB,DP00,DE00
      DOUBLE PRECISION RPB(3,3),RN(3,3),RNPB(3,3), RP(3,3),RBP(3,3),DDP00,DDE00,TT
      DOUBLE PRECISION GELOAP,GELAAP,GELOAPD,GELAAPD,iau_OBL06,DPSI,DEPS,EH,ETF,ET0
      DOUBLE PRECISION DELTA,DX00,DY00,DAT,longit,latit,height,V2(3),u1(3),u3(3),u4(3)
      DOUBLE PRECISION tol, passo
      
      common HPE
      CHARACTER*6  NAMS(400)

      tol= 1D-08  ! tolleranza di calcolo   
      passo= 1D0  ! passo iniziale (1 giorno)
     
      longit= 0d0; latit=0d0; height=0d0 
      DJMJD0 = 2400000.5D0
      
      IY = year
      MO = MO         !  default values       
!      ID = ID
      IH = 0
      IM = 0
      SEC = 0d0

      IF( MO == 3) THEN
         ID = 19
      ELSE IF ( MO == 6 ) THEN
         ID = 20
      ELSE IF ( MO == 9 ) THEN
         ID = 21
      ELSE IF ( MO == 12 ) THEN
         ID = 21
      END IF  
  
      TIME = ( 60D0*(60D0*DBLE(IH) + DBLE(IM)) + SEC ) / 86400D0
    
      CALL iau_DAT ( IY, MO, ID, TIME, DAT, J ) 
       IF(J == -2 .OR. J == -3 .OR. J == -4) THEN
        WRITE(*,*) " BAD DATE ! (out of range) "
       END IF             
    
!--------- CONSTANTS
      AE = 6378136.6                       ! Equatorial radius for Earth (meters)  
      C  = 299792458D0                     ! Speed of light  m/s 
      AU = 149597870700D0                  ! Astronomical Unit (meter)(TDB)
      DS = 86400D0                         ! n° seconds  Day
      K = 0.01720209895D0                  ! K  heliocentric gravitational costant (UA/Day)
      AUDAY = (C * DS / AU)                ! Speed of light (AU per day  "TDB")
      VC = AU / DS / C                     ! Speed of body in unit of C.
      MU = 1.32712440017987D20             ! Heliocentric gravitational costant (1976 system)   
      MUC = 2D0 * (MU / (C * C))/AU        ! Gravitational mass Sun delecting light
      RMASSE = 1D0/332946.050895D0         ! Gravitational mass Earth deflecting light
      OBL = 84381.406D0 * AS2R             ! Ecliptic obliquity IAU 2010
!----------------------------------------------------------------------------------------------

!------Geodetic coordinate observer

      LONG = longit * DD2R                         ! Long. neg. to Weast,(radians)
      LAT  = latit * DD2R                          ! Latit. (radians)
      HM  = height                                 ! ASL. (meter)

!-- Transform geodetic to geocentric coordinates of the site.(reference ellipsoid  WGS84)
    
       CALL iau_GD2GC ( 1, LONG, LAT, HM, XYZ, J )   ! Geodetic to Geocentric  coordinate
       
       IF ( J.NE.0 )STOP
       U = SQRT ( XYZ(1)*XYZ(1) + XYZ(2)*XYZ(2) )    ! U = distance from Earth spin axis (km)
       V = XYZ(3)                                    ! V = distance north of equatorial plane (km)
!-----------------------------------------------------------------------------------------------

      call iau_CAL2JD( IY, MO, ID, DJM0, DJM, J)
      MJD = DJM + TIME                          ! UTC or TT in MJD (Julian Day modified) at date     

!---------Transform into internal format.
      CALL iau_DTF2D ( "UTC", IY, MO, ID, IH, IM, SEC, UTC1, UTC2, J )   !UTC1 + UTC2 = Julian Day
      IF ( J.NE.0 ) STOP                                                 !UTC1 is the Julian Day number and
      
      UTC = UTC1 + UTC2-DJMJD0                                            ! UTC   is the MJD number
      DATE_ = UTC1 - DJMJD0                                               ! DATE_ is the MJD number
   
      IF( MJD > FINAL_ ) THEN
       de_S = (MJD - 51544.5D0) * DS                ! Time in sec. from J2000 (JD 2451545.0 TT)
       MM = 6.239996D0 + de_S * 1.99096871D-7
       DELTA = 32.184 + DAT + 0.001657D0 *sin(MM + 0.01671D0 * sin(MM))
      END IF                                                 

!----------Call the interpolation routine for per XP,YP,DUT,dX,dY
   
      CALL iers_calc(UTC,DATE_,FINAL_,EOP_final,XPO,YPO,DUT,LOD, &
                     DX0,DY0,ut1_notid)

      DX00 = DX0/1000D0 * AS2R           ! mas---> radians 
      DY00 = DY0/1000D0 * AS2R 
           
      IF( MJD > FINAL_ )  then 
       DUT = 0D0 
       DX00 = 0D0
       DY00 = 0D0    
      ELSE
       DUT = DUT
      END IF
!-------------------------------------------------------------------------------------------------------
!--------- UTC -> UT1.
      CALL iau_UTCUT1 ( UTC1, UTC2, DUT, UT11, UT12, J )
      IF ( J.NE.0 ) STOP   
!--------- Extract fraction for TDB-TT calculation, later.
      UT = MODULO ( MODULO(UT11,1D0)+MODULO(UT12,1D0), 1D0 )
!---------- UTC -> TAI -> TT -> TDB.
      CALL iau_UTCTAI ( UTC1, UTC2, TAI1, TAI2, J )
      IF ( J.NE.0 ) STOP
      CALL iau_TAITT ( TAI1, TAI2, TT1, TT2, J )
       IF ( J.NE.0 ) STOP
!---------- TDB-TT (using TT as a substitute for TDB).
      DTR = iau_DTDB ( TT1, TT2, UT, LONG, U, V )
      CALL iau_TTTDB ( TT1, TT2, DTR, TDB1, TDB2, J )
      IF ( J.NE.0 ) STOP
     
      OBL6 = iau_OBL06(TT1,TT2)    
      TT  = TT1 + TT2 - DJMJD0   !Time show in MJD  

! =========================================================
!            IAU 2000A, equinox based
! =========================================================
     
!-------- Nutation, IAU 2000A.
      CALL iau_NUT00A ( DJMJD0, TT, DP00, DE00 )

!-------- Precession-nutation quantities, IAU 2000.
      CALL iau_PN00 ( DJMJD0, TT, DP00, DE00, &
           EPSA, RB, RP, RPB, RN, RNPB )
!-------- Transform dX,dY corrections from GCRS to mean of date.
      V1(1) = DX00
      V1(2) = DY00
      V1(3) = 0D0
      CALL iau_RXP ( RNPB, V1, V2 )
      DDP00 = V2(1) / SIN ( EPSA )
      DDE00 = V2(2)
   
!------- Corrected nutation.
      DPSI = DP00 + DDP00
      DEPS = DE00 + DDE00
!------- Build the rotation matrix.
      CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )
!------- Combine the matrices: PN = N x P.
      CALL iau_RXR ( RN, RPB, RNPB )

!============================================================
            
      TDB = TDB1+TDB2            !TDB in JD  
!      print*,"TDB di partenza = ",TDB
!========================================================================================================

!------- Use subroutine  PLEPH ( DE430)  
      EH = 0 
      ET = TDB 
    DO          
      ET = ET      
      ISTEP=1
      NSTEP=1
      NTARG = 11
      CALL PLEPH ( ET, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector
      CALL PLEPH ( ET, 3, 11, PVH)              ! PVH = Earth heliocentric state vector 
      DT = 0D0
      CALL iau_PVUP (DT,PVH,HPE)                 ! HPE Heliocentric position Earth vector 
      CALL PLEPH ( ET, 3, 12, PVE)               ! PVE = Earth barycentric state vector
      CALL PLEPH ( ET, NTARG, 11, BEV)           ! BEV = Body's heliocentric state vector 
      DO  I=1,3
        BVE(I) = PVE(I,1) * DT + PVE(I,2)
      END DO                                     ! BVE Barycentric velocity Earth vector
      CALL iau_PVUP (DT,PVE,EB)                  ! EB Barycentric  position Earth vector
      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP = Body's Heliocentric position
      CALL iau_PVMPV(PVB,PVE,GMB)                ! GMB  Body's geocentric mean state vector
      CALL iau_PVUP (DT,GMB,GPB)                 ! GPB  Body's geocentric mean J2000 position    

!----------------------------------------------------------------------------------------------------------
      
      CALL L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,HE0,THD,EPV)
      
!------RELATIVISTIC DEFLECTION OF LIGHT 
   
      CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
    
!-------  RELATIVISTIC ANNUAL ABERRATION   ( for apparent GCRS coordinate)

      CALL ABERR (AUDAY,u1,BVE,u0)
      
      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 

!------------------------------------------------------------------------------------------------------------

     CALL   iau_RXP   ( RNPB, u3, u4 )
!      equinox = "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2000 EQUINOX BASED)" 
     CALL GEPV(u4,OBL6,1,GELO,GELA)    ! GELO = Geoc.Ecliptic longit; GELA = Geoc.Eclipic latit. (radians)
    
      IF ( GELO < 0D0) THEN
       GELO = GELO + 2D0 * PI
      END IF
     
     GELOAP = GELO 
     GELAAP = GELA                          
     GELOAPD = GELOAP * DR2D
     GELAAPD = GELAAP * DR2D 
    
     angledif = abs( LG - GELOAPD)

! --------------- Schema dicotomia usato su Python (oltre che su Pascal, in passato) -----
!    print ('Tempi progressivi=',t2,' (Lon-LG)=',Lon-LG, 'passo=',passo)
!    if Lon > LG:
!        t1= t1-passo
!        passo=oneday*abs(LG-Lon)/2.0 # nuovo passo per dimezzamento dello scarto
!    else:
!        t1 = t1 + passo
!        if abs(LG-Lon) < eps: # or kount >=100:
!            break # chiude il loop infinito di (while true) appena raggiunta la convergenza desiderata
! ------------------------------------------------------------------------------------------

! ============== Metodo di convergenza per DICOTOMIA ============== 
    IF (GELOAPD > LG) then
         ET = ET - passo
         passo= abs(LG-GELOAPD)/2.0
    else
         ET = ET + passo
         IF (abs(LG-GELOAPD) < tol) then
            goto 1000 
         end if    
    end if
! ==================================================================     

     END DO
1000 continue
!      print*," DEVIAZIONE = arcsec ", angledif *3600d0
       
      ET0 = INT(ET)
      ETF = ET - ET0

      CALL iau_TDBTT ( ET0, ETF, DTR, TT1, TT2, J )
       IF ( J.NE.0 ) STOP
      CALL iau_TTTAI ( TT1, TT2, TAI1, TAI2, J )
       IF ( J.NE.0 ) STOP
      CALL iau_TAIUTC ( TAI1, TAI2, UTC1, UTC2, J )
       IF ( J.NE.0 ) STOP
      CALL iau_JD2CAL( UTC1, UTC2, IY, M0, ID, FD, J )
      CALL time_hms ( FD, IH, IM, SEC)      

      end subroutine equinox
!-------------------------------------------------------------------------------------------------------      


 

