 program planets

!     THIS PROGRAN COMPUTES THE APPARENT DIRECTION OF SOLAR SYSTEM
!     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM
!
!    Compilation: gfortran -w rts.f90 sky77.for 1919.for jpl_430A.f -o rts.exe
!
!
!  - - - - - - - - - - -
      IMPLICIT NONE

! Pi
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )
      
! Arcseconds to radians
      DOUBLE PRECISION AS2R
      PARAMETER ( AS2R = 4.848136811095359935899141D-6 )
! 2Pi
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI = 6.283185307179586476925287D0 )
! Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )
!  Radians to arcseconds
      DOUBLE PRECISION DR2AS
      PARAMETER ( DR2AS = 206264.8062470963551564734D0 )
!  Ratio between solar and sidereal time
      DOUBLE PRECISION SOLSID
      PARAMETER (SOLSID=1.002737909350795D0)
      DOUBLE PRECISION DJMJD0
      PARAMETER (DJMJD0=2400000.50000000001d0)  


      
      INTEGER IY, MO, ID, IH, IM, J, I,mode,IREFR,SGN,ios,NVS,NTARG,NCENT,ISTEP,NSTEP
      INTEGER IHMSF1(4),IHMSF2(4),IHMSF3(4),IHMSF0(4),IDMSF0(4), Rtoday,YN
      INTEGER FINAL_,JW,EQTM,IY1, IM1, ID1,IDMSF6(4),conta,EOP_final
      INTEGER IYF, MOF, IDF,N0,MS,PT,NUM,JJ,IHMSF6(4),MES,YA
      integer,dimension(8) :: VALUES
      DOUBLE PRECISION SEC,XP,YP,DUT,TT,TCB1,TCB2,T,AU,PS,DEL,HA,TDB,Ltime,TDB1,TDB2,TDB0,OBL,s_s,JTRS
      DOUBLE PRECISION  UT11, UT12, UT,TAI1, TAI2, TT1, TT2, TCG1, TCG2,RVETT,UTC1, UTC2,UTC 
      DOUBLE PRECISION DDP80, DDE80, DX00, DY00, DX06, DY06,DELTA_T,DE_T,u01(3),Rdot,GN,Mtime,ETT,TRSJ
      DOUBLE PRECISION DATE_, TIME, DAT,TCB,TCG,POS2(3),DOT,TAI, TUT, UT1, DP80, DE80, MJD0,DJM,DTA,JTRS1
      DOUBLE PRECISION RA3,DE3,RA4,DE4,RA5,DE5,RA6,DE6,RA7,DE7,RA8,DE8,RA9,DE9,DJM0,RA10,DE10,ut1_notid
      DOUBLE PRECISION RMASSE,SUNLONG,COSELON,EDEC,EVET,GLON,GLAT,XEL,DIST,DIST4,ECC,dra,TDBA,TDBB         
      DOUBLE PRECISION DPSI, DEPS, EPSA, ERA, DP00, DE00, DDP00, DDE00 ,ELON, DELTAT, EARTHVEC,JTRS2 
      DOUBLE PRECISION GAST,LONG ,LAT,HM,DTR,MO1,FD,julian,DSQRT,HLON,HLAT,radius,T0,TDBJD,ET,XET
      DOUBLE PRECISION  X, Y, S , DT, ELONG, u02, MO2, longit, DELTA, DS, AUDAY, XPO, YPO, DIAM,AE,JD1,JD2      
      DOUBLE PRECISION latit,UJD,MJD,ANGLE,PANGLE,PDIAM,AZ,height,ZD,MAGN,AUM,Dtime,PVV(2),DJMF,Day
      DOUBLE PRECISION RA7H, DE7D,JD,RVET,RVET0,EPOS,EPSA0,DP000,DE000,MM,iau_OBL80,B(8,3),OBL6
      DOUBLE PRECISION DX0,DY0,UT_0,TT_0,GAST0,PL,SD,RIS,SET,TRS,AZM1,AZM2,HAL,MUC,BP,de_S,TIME1
      DOUBLE PRECISION C,MU,u2,RAR,DER,K,TAU,THD,TGD,ILL,COSI,DET,HD,GD,HMOD,DIST1,DIST2,DIST3
      DOUBLE PRECISION RC2IT(3,3),RC2I(3,3),RPB(3,3),RNPB0(3,3),RN(3,3),RNPB(3,3), RP(3,3),RBP(3,3)
      DOUBLE PRECISION PV(3,2),PVH(3,2),PVB(3,2),BEV1(3,2),BEV(3,2),RB(3,3),GMB(3,2),PVE(3,2),HSV(3,2)
      DOUBLE PRECISION HE0(3,2),SUNB(3,2),EPV(3,2),GPB(3),HV(3),POS(3),VEL(3),XYZ(3),V1(3),RA66,RA,RA0
      DOUBLE PRECISION EB(3),u0(3),u03(3),u1(3),u3(3),u4(3),u5(3),u6(3),U(3),V(3),V2(3),p(3),RAsun,DE0 
      DOUBLE PRECISION VEL2(3),VEL3(3),POS3(3),POS_H(3),GE(3),HPE(3),BVE(3),BEP(3),TSL,EO ,BGV(3,2)
      DOUBLE PRECISION BODY_G(3),BODY_H(3),HE(3),RS(3),OB(3),OV(3),GP(3),GV(3),POS_T(3),EO0,eta,PI_
      DOUBLE PRECISION iau_ANP,iau_GMST00,iau_S06,iau_DTDB,iau_GMST06,iau_GST06A,iau_GST00A,CHI,p_,A_,B_    
      DOUBLE PRECISION iau_ERA00, iau_EE00,iau_EE06A,iau_EORS,iau_EE00A,deldot,C_,OBL0,uu(3),SAB(3,2)
      DOUBLE PRECISION AOB,ZOB,HOB,DOB,ROB,PHPA,TC,RH,WL,ZOT,AOT, HOT, DOOT,ROT,EQEQ,GMST0,ZD0,GEL,atlon  
      DOUBLE PRECISION Sun_lon,Earth_lon,TDBJ,iau_FAE03,EQT,EQTS,EQT0,LOD,EQTD,LatLu,AGE,DJ2,Rsun 
      DOUBLE PRECISION DE_GMST,UTC0,UTC01,UTC02,RA06,EQTX,BODY_sun(3),R_sun,JHE,MJ,AX,JTT 
      DOUBLE PRECISION BODY_H1,GD1,HD1,XET1,PV1,HE01,THD1,E_P(3,2),D_sun,DISTKm,LoLu,LoLu_ap,OMEGA,OME0
      DOUBLE PRECISION JAP,DLR,DJ3,PH,TLong,TLat,LTL,LTB,PA,iau_OBL06,LON0,LAT0,T_0,a0,b0,c0,THL,THB 
      DOUBLE PRECISION BGE,Distmoon,R_sun1,VC,CVE(3),BM1,BCV,er,Tcy,PHE,KA,Lsun,Bsun,AZ1,SGE(3,2)
      DOUBLE PRECISION RVET1,Rdot1,Ltime1,deldot1,Mtime1,HD0,BEV0(3,2),LEA(6),JUB(3,2),JH0,K0,KU
      DOUBLE PRECISION  R(6)
      DOUBLE PRECISION  SS(3)
      DOUBLE PRECISION  VALS(400)
      DIMENSION body(11),radius(10),dayweek(7)
      
      CHARACTER*8 DATE
      CHARACTER*10 TIMEX
      CHARACTER*6 ZONE,data_b
      CHARACTER(len=1):: sign,sign1,sign2,aS   
      CHARACTER(len=6):: NAMS(400),TX,choose,yes,no,distance,EL
      CHARACTER(len=16):: XDIAM,dayweek,Dday,IDEM 
      CHARACTER(len=26):: body,planet,LDUT,LDET,EQ_label,Slod
      CHARACTER(len=2):: lS,bS,update,NEW,str,AS1,AS2,AS3
      CHARACTER(len=95):: equinox,equinox2,equinox3,commentR,commentS,commentT,commentU,labelL
      CHARACTER(len=95):: labelR,labelS,labelT,labelM,labelN,labelO,labelX,labelY,commentM,commentN
      CHARACTER(len=7):: APE 
      CHARACTER(len=14):: comment 
      common SGE,RVET0,HSV,NTARG                   

      DATA body / "MERCURY","VENUS ","EARTH ","MARS ","JUPITER","SATURN ","URANUS ","NEPTUNE","PLUTO ", &
                  "MOON ","SUN "/
                  
      DATA radius / 2440D0, 6051.8D0, 6378.137D0, 3397D0, 71492D0, &
                    60268D0, 25559D0, 24766D0, 1151D0, 1737.53D0 /

      DATA dayweek / "SUNDAY","MONDAY","TUESDAY","WEDNESDAY","THURSDAY","FRIDAY", &
                     "SATURDAY" /

      !  **************************************************************************
      !  *                                                                        *
      !  *      --------------------rts_planet.f90 -------- (May.20,2020) --------*
      !  *                                                                        *
      !  *     THIS PROGRAN COMPUTES RISE-TRANSIT-SET OF SOLAR SYSTEM             *
      !  *     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM      *
      !  *            for this purpose we use JPL DE430t database                 *
      !  *                                                                        *
      !  *     Compiling command:                                                 *
      !  *     gfortran -w rts.f90 sky7.for 1919.for jpl_430A.f                      *
      !  *     which generates the executable file: a.out                         *
      !  *                                                                        *
      !  *     Range of calendar dates:                                           *
      !  *      a) Not less 1962-01-01                                            *
      !  *      b) Not over MJD of line 128 (or nearby), i.e. FINAL_ = 58492 of   *
      !  *         main file .                                                    *
      !  *      c) Ephemerides computations may be forced up to 2100/01/01 upper  *                 
      !  *         limit of JPL database,but with less precision.                 *
      !  *         ------------------------------------------------------------   * 
      !  *         For FINAL_-iers.txt updating, please refer to : www.iers.org/  *
      !  *       - Data/Products - Standard EOP data files - Finals data IAU2000  *  
      !  *         and change the FINAL_ value accordingly, say MJD shown 3 lines *
      !  *         prior the last full row, idem for database EOP14 C04           * 
      !  *                                                                        *
      !  **************************************************************************
      
                   
!       ****  We recommend that you update at least every 15 days, checking ****  
!       ****      the FINAL_ dates of the respective IERS database           ****  
!        ========================================================================= 

!       Enter the last updated dates (MJD) of the IERS databases 
!       Finals data (IAU2000) and EOP14 C04 (IAU2000A) or uses   
!       default dates (1 or 2).  

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
      TSL = 20.0                           ! Is the approximate sea-level air temperature in K"
      RH = 0.5                             ! relative humidity at the observer (range 0-1)
      WL = 0.55                            ! wavelength (micrometers)
      KA = 20.49552d0                      ! Costant of aberration 
!---------------------------------------------------------------------------------------  


      FINAL_ = 59342       ! 2021/05/04   reported in MJD IERS FINAL_ data (for update) 

      EOP_final = 58940  ! 2020/04/01 jump from IERS EOP14 C04 to IERS FINAL_ data (IAU2000).
    

      !"========================================================================"              
      !"                                                                        "
      !"        THE NUMBERING CONVENTION FOR 'NTARG' IS:                        "
      !"                                                                        "
      !"     1 = MERCURY          7 = URANUS                                    "
      !"     2 = VENUS            8 = NEPTUNE                                   "
      !"     3 = EARTH            9 = PLUTO                                     "
      !"     4 = MARS            10 = MOON                                      "
      !"     5 = JUPITER         11 = SUN                                       "
      !"     6 = SATURN                                                         "
      !"                                                                        "
      !"========================================================================"
      
!      WRITE(*,*) 'Enter Target Body ' 
!      READ *, NTARG
      NTARG = 10

      DO I = 1,11
       IF ( NTARG == I) THEN
        planet = body(I)
       END IF
      END DO          

      WRITE(*,*)""
      WRITE(*,*) "  Observer geodetic coordinate (degrees) " 
         longit = 15.0d0
         latit =  42.0d0 
         height = 0.0 ! (m)
         TC = 10d0
      WRITE(*,5) longit,latit,height
5     FORMAT(3x,"Long",F9.5,2x"Latit",F9.5,2x,"height",F6.1,"m")  
    
!--------- INPUT TIME (UTC)
       IY = 1993; MO = 2; ID = 1
       IH =0d0 ; IM = 0d0; SEC = 0d0        
      WRITE (*,10) "Initial UTC time  ",IY,MO,ID,IH,IM,SEC
       IYF = 1993; MOF =3; IDF = 31
       IH =0d0 ; IM = 0d0; SEC = 0d0     
      WRITE (*,10) "Final UTC time    ", IYF,MOF,IDF,IH,IM,SEC
10    FORMAT(3x,A16,x,I4,"/",I2,"/",I2,4x,I2,":",I2,":",F6.3)

      CALL iau_CAL2JD(  IY, MO, ID,DJM0, DJM, J )
        IF ( J.NE.0 ) STOP
      CALL iau_CAL2JD(  IYF, MOF, IDF,DJM0, DJMF, J )
        IF ( J.NE.0 ) STOP
      TIME1 = ( 60D0*(60D0*DBLE(IH) + DBLE(IM)) + SEC ) / 86400D0
      JD1 = DJM0 + DJM + TIME1                      ! JD start
      JD2 = DJM0 + DJMF + TIME1                     ! JD stop
      NUM = JD2 - JD1
      JD = JD1
     

      WRITE (*,*)
      WRITE (*,*)" Celestial Body: ",planet
      WRITE (*,*)" Geocentric apparent RA - DEC 0.0h UT true equinox and ecliptic of the date." 
      WRITE (*,*)" TRANSIT (RA-DEC)Apparent topocentric coordinate equinox of date."
      WRITE (*,*) 
       
     IF (NTARG == 10)THEN 
      WRITE (*,*)"    TIME     MOONS'  ILL.%   RISE (UT)     TRANSIT (UT)    SET (UT)     Azimuth   Elev.   Azimuth      &   
                   RA 0.0h     DEC 0.0h      RA(transit)    DEC(transit)     "
      WRITE (*,*)"   0h UT      AGE              --              --             --          RISE   TRANSIT    SET      &            
                 (app.geoc)    (app. geoc)        UT             UT  " 
      WRITE (*,*)"===========================================================================================================&
                 ==================================================="
     ELSE IF(NTARG /= 10) THEN
      WRITE (*,*)"    TIME       RISE (UT)     TRANSIT (UT)   SET (UT)     Azimuth  Elev.   Azimuth     RA 0.0h      &
                  DEC 0.0h    RA(transit)   DEC(transit)   "
      WRITE (*,*)"     UT                                                   RISE   TRANSIT   SET       (app.geoc)   &
                 (app.geoc)       UT             UT  " 
      WRITE (*,*)"============================================================================================================&
                 ===================================="
    END IF  

     DO JJ = 1, NUM   
      
      JD = JD
      DJ2 = JD - DJMJD0
      CALL iau_JD2CAL ( DJMJD0, DJ2, IY, MO, ID, FD, J )
      CALL time_hms ( FD, IH, IM, SEC)      
      CALL iau_DTF2D ( "UTC", IY, MO, ID, IH, IM, SEC, UTC1, UTC2, J )   !UTC1 + UTC2 = Julian Day
      
      JD = UTC1 + UTC2                         ! UTC or TT in Julian Day at date
      MJD = JD - DJMJD0                        ! UTC or TT in MJD (Julian Day modified) at date  

!-------------------------------------------------------------------------------------------------

      CALL iau_DAT ( IY, MO, ID, TIME, DAT, J ) 
       IF(J == -2 .OR. J == -3 .OR. J == -4) THEN
        WRITE(*,*) " BAD DATE ! (out of range) "
       END IF           

!------Geodetic coordinate observer

      LONG = longit * DD2R                         ! Long. neg. to Weast,(radians)
      LAT  = latit * DD2R                          ! Latit. (radians)
      HM  = height                                 ! ASL. (meter)
      PHPA = 1013.25 * exp ( -HM / ( 29.3 * TSL ) )! pressure at the observer (hPa = mB

      IF (longit > 0D0) THEN
       lS = "E"
       ELSE IF (longit < 0D0) THEN
       lS = "W"
      END IF
      
      IF ( latit > 0D0 ) THEN 
       bS = "N"
       ELSE IF ( latit < 0D0 ) THEN
       bS = "S"  
      END IF 
     
!-- Transform geodetic to geocentric coordinates of the site.(reference ellipsoid  WGS84)
      
       CALL iau_GD2GC ( 1, LONG, LAT, HM, XYZ, J )   ! Geodetic to Geocentric  coordinate
       
       IF ( J.NE.0 )STOP
       U = SQRT ( XYZ(1)*XYZ(1) + XYZ(2)*XYZ(2) )    ! U = distance from Earth spin axis (km)
       V = XYZ(3)                                    ! V = distance north of equatorial plane (km)         
           
      IF( MJD > FINAL_ -3 ) THEN
       de_S = (MJD - 51544.5D0) * DS                ! Time in sec. from J2000 (JD 2451545.0 TT)
       MM = 6.239996D0 + de_S * 1.99096871D-7
       DELTA = 32.184 + DAT + 0.001657D0 *sin(MM + 0.01671D0 * sin(MM))
      END IF                                                 
       
!---------Transform into internal format.
175   CALL iau_DTF2D ( "UTC", IY, MO, ID, IH, IM, SEC, UTC1, UTC2, J )   !UTC1 + UTC2 = Julian Day
      IF ( J.NE.0 ) STOP                                                 !UTC1 is the Julian Day number and
      CALL iau_DTF2D ( "UTC", IY, 1, 1, 0, 0, 0D0, UTC01, UTC02, J )     !UTC2 is the fraction of a day.   
    
      UTC = UTC1 + UTC2-DJMJD0
      DATE_ = UTC1 - DJMJD0                                               ! DATE_ is the MJD number
      UTC0 = UTC01 + UTC02                                               ! UTC JD at 0.0h YYYY/01/01  
!----------Call the interpolation routine for per XP,YP,DUT,dX,dY
       
      CALL iers_calc(UTC,DATE_,FINAL_,EOP_final,XPO,YPO,DUT,LOD,DX0,DY0,ut1_notid)
      
      IF(MJD > (EOP_final-4)) then
       data_b = "IERS_A" 
       IF (MJD > FINAL_ .and. MJD > Rtoday) then
       data_b = "******" 
       END IF
      ELSE IF(MJD <= (EOP_final-4)) then 
       data_b = "EOP14 C04"
       IF (MJD < 37665.0) then
       data_b = "EOP C01" 
       END IF
      END IF      
     
      DUT = DUT                          !  sec.
      TUT = UTC2 + DUT/86400D0           !  in fraction of day          
      LOD = LOD
      XP = XPO * AS2R                    ! arcsec---> radians
      YP = YPO * AS2R   
      
      DX00 = DX0/1000D0 * AS2R           ! mas---> radians 
      DY00 = DY0/1000D0 * AS2R 

      IF( MJD > FINAL_ -4) THEN             !
       DUT = 0D0
       LOD = 0D0                         !
       TUT = UTC2                        !
       XP = 0D0                          !
       YP = 0D0                          !
       DX0 = 0D0                         !
       DY0 = 0D0                         ! IF MJD > FINAL_
       DX00 = 0D0                        !
       DY00 = 0D0                        !
      END IF                             !

!--------- UTC -> UT1.
      CALL iau_UTCUT1 ( UTC1, UTC2, DUT, UT11, UT12, J )
      IF ( J.NE.0 ) STOP

!--------- Extract fraction for TDB-TT calculation, later.
      UT = MODULO ( MODULO(UT11,1D0)+MODULO(UT12,1D0), 1D0 )
    
!---------- UTC -> TAI -> TT -> TCG.
      CALL iau_UTCTAI ( UTC1, UTC2, TAI1, TAI2, J )
      IF ( J.NE.0 ) STOP
      CALL iau_TAITT ( TAI1, TAI2, TT1, TT2, J )
      IF ( J.NE.0 ) STOP
      CALL iau_TTTCG ( TT1, TT2, TCG1, TCG2, J )
      IF ( J.NE.0 ) STOP
    
!---------- TDB-TT (using TT as a substitute for TDB).
      DTR = iau_DTDB ( TT1, TT2, UT, LONG, U, V )

!----------- TT -> TDB -> TCB.
      CALL iau_TTTDB ( TT1, TT2, DTR, TDB1, TDB2, J )
      IF ( J.NE.0 ) STOP
        
      CALL iau_TDBTCB ( TDB1, TDB2, TCB1, TCB2, J )
      IF ( J.NE.0 ) STOP

      TT  = DATE_+TT2          ! 
      TDB = DATE_+TDB2         !
      TAI = DATE_+TAI2         ! Time show in MJD
      TCG = DATE_+TCG2         ! 
      TCB = DATE_+TCB2         !
      UT1 = DATE_+UT12         !
      UTC = DATE_+UTC2
      IF( MJD > FINAL_ ) THEN
       DET = DELTA
       DELTAT = DELTA/86400D0
       GOTO 195
      ELSE IF ( MJD <= FINAL_ ) THEN
       DET = 32.184D0 + DAT - DUT           ! in sec.
       DELTAT = DET / 86400D0               ! in JD.        
      END IF       
      
195      UT_0 = INT (MJD) + DJMJD0            ! JD Time to 0.0 h UT1. 
     
         TT_0 = UT_0 + DELTAT                 ! JD Time to 0.0 h TT the day.

! ======================================================
!   IAU  2006/2000A,  CIO    based,   using  X,Y series
! ======================================================
 
!   CIP  and  CIO,  IAU 2006A. 
      CALL iau_XY06 ( DJMJD0, TT, X, Y )
     S = iau_S06 ( DJMJD0, TT, X, Y )
    
!   Add CIP corrections.
       X = X + DX00
       Y = Y + DY00

!   GCRS  to  CIRS  matrix. 
       CALL   iau_C2IXYS  (   X,  Y,  S, RC2I   ) 
      
!   Earth  rotation  angle. 
       ERA  =  iau_ERA00  (   UT11,  UT12    ) 
                          
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
     
!========================================================================================================

!---------- Nutation, IAU 2000A. to 0.0h TT
      CALL iau_NUT00A ( TT_0, 0D0, DP00, DE00 )
      
!---------- Precession-nutation quantities, IAU 2000. to 0.0h TT
      CALL iau_PN00 ( TT_0, 0D0, DP00, DE00, &
           EPSA, RB, RP, RPB, RN, RNPB0 )
!========================================================================================================
      
!--------- Greenwich apparent sidereal time (IAU 2000) for 0.0 h UT1.
      
      GAST0 = iau_ANP (iau_GMST00 ( UT_0, 0D0, TT_0, 0D0 ) &
          + iau_EE00 ( TT_0, 0D0, EPSA, DPSI ) )
                     
!--------- Greenwich apparent sidereal time (IAU 2000).
 
      GAST = iau_ANP ( iau_GMST00 ( UT11, UT12, TT1, TT2 ) &
          + iau_EE00 ( DJMJD0, TT, EPSA, DPSI ) )

      EO = ERA - GAST                                ! EO0  Equation of the origin.
      EO0= iau_EORS(RNPB,S)

!--------- Greenwich mean sidereal time (IAU 2000). 
     
      GMST0 = iau_GMST00 ( UT11, UT12, TT1, TT2 )
      
!--------- Equation of equinox (EQEQ  radians )  

      EQEQ = iau_EE00A ( TT1, TT2 )
      EQEQ = EQEQ * DR2D                            ! EQEQ in degs.

!=====================================================================================================
      CALL GEOCPV ( LONG,LAT,HM,GAST,POS,VEL,J)     ! Geocenter position and velocity observer. 
      IF ( J.NE.0 ) STOP

!==========================================================================================================
      
      ET = TDB1+TDB2
      T0 = 2451545.0D0
      Tcy = (ET - T0) / 36525d0     ! Time in century from Epoch J2000

!---------- Eccentricity of Earth orbit (ECC)

      ECC = 0.016708634d0 - 0.000042037d0 * Tcy - 0.0000001267d0 * Tcy**2d0     

!---------- Longitude of the Perihelion of Earth orbit (PHE radians)

      PHE = (102.93735d0 + 1.71946d0 * Tcy + 0.00046d0 * Tcy**2d0) * DD2R
      
!-------- Nutation - Precession POS. - VEL. for Topocentric coordinate 
      CALL NUTATION (-ET, DJMJD0, TT, POS, POS2)
      CALL NUTATION (-ET, DJMJD0, TT, VEL, VEL2)
      
!     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   GP )

!     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   GV )
      
!---------------------------------------------------------------------------------------------------
!------- Use subroutine  PLEPH ( DE430)    
  
      ISTEP=1
      NSTEP=1
    
!      CALL PLEPH ( ET, 14, 0, NUT)              ! 1980 IAU nutation angles. (Long,Obliq,Long. rate, Obliq rate ) 
        
      CALL PLEPH ( ET, 15, 0, LEA)              ! Lunar Euler angles: phi, theta, psi 
    
      CALL PLEPH ( ET, 17, 0, PVV)              ! PVV(1) = TT-TDB in sec./ PVV(2)=rate of change of TT-TDB in sec/day
   
      CALL PLEPH ( ET, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
      CALL PLEPH ( ET, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( ET, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector

      CALL PLEPH ( ET, NTARG, 11, BEV)          ! BEV = Body's heliocentric state vector 

      CALL PLEPH ( ET, 11, 12, SUNB )           ! SUNB = Sun barycentric state vector    

      CALL PLEPH ( ET, 11, 3, SGE )             ! SGE = Sun geocentric state vector  
      
      CALL PLEPH ( ET, NTARG, 3, BGV)           ! BGV = Body geocentric state vector
         
!------Update a PVH-vector, discarding the velocity component.

      DT = 0.0D0   
   
      CALL iau_PVUP (DT,PVH,HPE)                 ! HPE Heliocentric position Earth vector 

      CALL iau_PVUP (DT,PVE,EB)                  ! EB Barycentric  position Earth vector

      DO  I=1,3
        BVE(I) = PVE(I,1) * DT + PVE(I,2)
      END DO                                     ! BVE Barycentric velocity Earth vector
                
      CALL iau_PPP (EB,GP,OB)                    ! OB = EB + GP (barycentric position of observer)
 
      CALL iau_PPP (BVE,GV,OV)                   ! OV = BVE + GV (barycentric velocity of observer)

      CALL iau_PVUP (DT,PVB,BP)                  ! BP Body's Barycentric position 

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP Body's Heliocentric position

      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus
  
      CALL iau_PVMPV(PVB,SUNB,HSV)               ! HSV  Heliocentric state vector

      CALL iau_PVMPV(PVB,PVE,GMB)                ! GMB  Body's geocentric mean state vector
      
      CALL iau_PVUP (DT,GMB,GPB)                 ! GPB  Body's geocentric mean J2000 position    

      CALL  iau_C2S   ( GPB, RA4, DE4 )       
      RA4 = iau_ANP(RA4)                         ! range angle (0,2!pi)
      
      DIST2 = DSQRT(GPB(1)*GPB(1)+GPB(2)*GPB(2)+GPB(3)*GPB(3))     ! Geocentric Distance 

      EARTHVEC = DSQRT(HPE(1)*HPE(1)+HPE(2)*HPE(2)+HPE(3)*HPE(3))  ! Distance Sun - Earth
      
!---------- Astrometric coordinate      
    
      CALL L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,HE0,THD,EPV) ! vedi EPV
   
      BEV1 = PV - PVE                            ! Geocentric state vector
           
!-----------Form  astrometric position RA3,DE3
! ......... Call for direction cosines to spherical coordinates

      CALL  iau_C2S ( BODY_G, RA3, DE3 )  
     
      RA3 = iau_ANP(RA3)                         ! range angle (0,2!pi)
!----------- Transform RA3, DE3 in astrometric-topocentric 
!--------------Form RVET0 (vector ray mean J2000 )and  Rdot.
      
      RVET0 = HD  !DSQRT(BODY_H(1)**2D0+BODY_H(2)**2D0+BODY_H(3)**2D0)
      
      Rdot = (HE0(1,2)*HE0(1,1)+HE0(2,2)*HE0(2,1)+HE0(3,2)*HE0(3,1))/ HD  

      Rdot = Rdot*AU/86400D0/1000D0                 ! Rdot Km/sec.              
   
!-------------Form " deldot "            
       
      DOT = (BEV1(1,2)*BEV1(1,1)+BEV1(2,2)*BEV1(2,1)+BEV1(3,2)*BEV1(3,1))/GD 

      deldot = DOT*AU/86400D0/1000D0                ! deldot Km/sec. 

      BODY_G = BODY_G
      DIST1 = GD                                 ! true distance (astrometric)
      HD = HD                                    ! True vector ray 
      Ltime = DIST1 * 8.316746395045D0           ! Light time (min.)
      Mtime = DIST1 * 8.316746395045D0 * 60D0    ! Moon Light time (sec)
      Distmoon = GD*AU/1000d0/6378.137d0         ! Moon's distance in Earth-radius

!------------- Compute the phase angle of body

      RVET = RVET0
      COSI = (RVET**2D0+DIST1**2D0-EARTHVEC**2D0) / (2D0*RVET*DIST1)
      ANGLE = acos(COSI) * DR2D                     ! ANGLE = Phase angle

      ILL = ((RVET + DIST1)**2D0 - EARTHVEC**2D0) / (4D0 * RVET * DIST1) *100D0 ! % Illuminated fraction disk   

!-----------LON0 = ELIOCENTRIC ECLIPTIC LONGITUDE ; LAT0 = ELIOCENTRIC ECLIPTIC LATITUDE  (true apparent J2000)
      OBL0 = OBL 
      OBL6 = iau_OBL06(TT1,TT2)
      CALL GEPV(HSV,OBL0,1,LON0,LAT0)
     
!----------- Reduce Heliocentric Coordinate mean J2000 to mean of the date.(Meeus Astronomical Algorithms; rigorous method)
     
      T_0 = (JD - 2451545.0d0)/36525d0                                               ! time in centuries          

      eta = (47.0029d0 * T_0 - 0.03302d0 * T_0**2 + 0.00006d0 * T_0**3) * AS2R       ! ( arcsec to rad.)
      PI_ = (174.876384d0 - (869.8089d0 * T_0 +0.03536d0 * T_0**2) /3600d0) * DD2R   ! (arcsec to degs. to rad.)
      p_  = (5029.0966d0 * T_0 + 1.11113d0 * T_0**2 - 0.000006d0 * T_0**3) * AS2R    ! ( arcsec to rad.)

      A_ = cos(eta) * cos(LAT0) * sin(PI_ - LON0) - sin(eta)*sin(LAT0)
      B_ = cos(LAT0) * cos(PI_ - LON0)
      C_ = cos(eta) * sin(LAT0) + sin(eta) * cos(LAT0) * sin(PI_ - LON0)

      atlon = atan2(A_,B_)                                                             ! tan A_/B_ = (p_ + PI_ - HLON)   
      HLON = p_ + PI_ - atlon                                ! HLON = ELIOCENTRIC ECLIPTIC LONGITUDE (mean of the date)
      HLAT = asin(C_)                                        ! HLAT = ELIOCENTRIC ECLIPTIC LATITUDE  (mean of the date)

      IF ( LON0 < 0D0) THEN
       LON0 = LON0 + D2PI
      END IF
!-------------------------------------------------------------------------------------------------------------------           
      IF (NTARG == 11) THEN
        HLON = 0D0; HLAT = 0D0
      END IF
!----------- Form Magnitude and angular diameter 

      PANGLE = ANGLE
      CALL MAGNIT (ET,NTARG,ANGLE,RVET,DIST1,PVH,EARTHVEC,OBL,HLON,HLAT,ILL,MAGN) !!

      DO I = 1,10
       IF ( NTARG == I) THEN
        DIAM = 2D0 * radius(I)
       END IF
      END DO  

      IF ( NTARG == 5 ) PDIAM = DIAM *(1D0 - 0.06487D0)  
      IF ( NTARG == 6 ) PDIAM = DIAM *(1D0 - 0.09796D0)
      IF ( NTARG == 7 ) PDIAM = DIAM *(1D0 - 0.02293D0)
      IF ( NTARG == 8 ) PDIAM = DIAM *(1D0 - 0.0171D0)

      DIAM = DIAM * DR2AS /( DIST1 * AU/1000D0 )
      PDIAM = PDIAM * DR2AS /( DIST1 * AU/1000D0 )

      IF ( NTARG == 11 ) THEN
       DIAM = 1919.29D0 / DIST1
       ILL = 1.0D0 *100D0
      END IF
!--------------------------------------------------------------------------------------------------
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun,Saturn and Jupiter ( for apparent geocentric coordinate)
 
      B(1,1) = 0.00028574D0      ! Saturn parameter
      B(2,1) = 3D-10
      B(3,1) = SAB(1,1)   
      B(4,1) = SAB(2,1)
      B(5,1) = SAB(3,1)
      B(6,1) = SAB(1,2)
      B(7,1) = SAB(2,2)
      B(8,1) = SAB(3,2)

      B(1,2) = 0.00095435D0      ! Jupiter parameter 
      B(2,2) = 3D-9
      B(3,2) = JUB(1,1)
      B(4,2) = JUB(2,1)
      B(5,2) = JUB(3,1)
      B(6,2) = JUB(1,2)
      B(7,2) = JUB(2,2)
      B(8,2) = JUB(3,2)

      B(1,3) = 1D0              ! Sun parameter
      B(2,3) = 6D-6
      B(3,3) = SUNB(1,1)
      B(4,3) = SUNB(2,1)
      B(5,3) = SUNB(3,1)
      B(6,3) = SUNB(1,2)
      B(7,3) = SUNB(2,2)
      B(8,3) = SUNB(3,2)

      CALL iau_PN(BODY_G,MO1,uu)

      CALL iau_LDN( 3, B, EB, uu , u1) 

      CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

      CALL  iau_C2S   ( p, RA5, DE5 )       
      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
      
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent GCRS coordinate)

      CALL iau_SXP ( VC,BVE,CVE)    ! CVE body velocity vector in unit of C 
      
      BCV = SQRT(CVE(1)**2d0+CVE(2)**2d0+CVE(3)**2d0) !form module of vector CVE

      BM1 = SQRT(1d0-abs(BCV)**2d0)            ! for use with routine SOFA iau_AB()  

      CALL ABERR (AUDAY,u1,BVE,u0)
     
      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 
   
! ......... FINAL_ GCRS (RA6,DE6)  body apparent geocentric coordinate (RA6,DE6)

      CALL   iau_RXP   ( RNPB, u3, u4 )

! ......... Call for direction cosines to spherical coordinates
      CALL   iau_C2S   ( u4, RA6, DE6 )
 
         RA6= iau_ANP(RA6)                       ! range angle (0,2!pi)
         RA06 = RA6
         RA66= RA6-EO
 
      CALL iau_A2TF ( 1 ,RA6 ,sign , IHMSF6 )
      CALL iau_A2AF ( 1 ,DE6 ,sign2 , IDMSF6 )

!===================================================================================================
!--------- Compute  -TDB0- JD Time to 0.0 h TDB the day (for routine TRANSIT)
       
       TDB0 = INT(MJD) + DJMJD0   
       AUM = AU/1000D0
       PL = ASIN(6378.137/(DIST1 * AUM)) * DR2D  ! Lunar horizontal parallax 
       SD = (DIAM / 3600D0 ) / 2D0               ! Lunar semidiameter in arcmin
               
       CALL TRANSIT(NTARG,TDB0,GAST0,AUDAY,DELTAT,PL,SD,LAT,LONG, &
            RNPB0,RIS,SET,TRS,azm1,azm2,HAL,aS,N0,MS,PT)
            
      TRSJ = TRS/24d0 + JD                             ! Transit in UTC
      JTRS = (TRSJ + (32.184D0 + DAT)/86400D0)         ! Transit in JD (TT)
      JTRS1 = INT(JTRS) 
      JTRS2 = JTRS - JTRS1  
!--------- Form UT1-TAI.
      DTA = DUT - DAT                                      ! DTA in sec.
!----------- TT -> TAI 
      CALL iau_TTTAI ( JTRS1, JTRS2, TAI1, TAI2, J )       ! Transit in TAI  
!----------- TAI -> UT1 
      CALL iau_TAIUT1 ( TAI1, TAI2, DTA, UT11, UT12, J )   ! Transit UT1

!--------- Extract fraction for TDB-TT calculation, later.
      UT = MODULO ( MODULO(UT11,1D0)+MODULO(UT12,1D0), 1D0 )

       DTR = iau_DTDB ( JTRS1, JTRS2, UT, LONG, U, V )

!----------- TT -> TDB 
      CALL iau_TTTDB ( JTRS1, JTRS2, DTR, TDBA, TDBB, J )  ! Transit in TDB

      ETT = TDBA + TDBB
      T0 = 2451545.0D0
      JTT = JTRS - DJMJD0
!====================================================================================================
      GAST = iau_ANP ( iau_GMST00 ( UT11, UT12, JTRS1, JTRS2 ) &
          + iau_EE00 ( DJMJD0, TT, EPSA, DPSI ) )

!=====================================================================================================
      CALL GEOCPV ( LONG,LAT,HM,GAST,POS,VEL,J)     ! Geocenter position and velocity observer. 
      IF ( J.NE.0 ) STOP

!====================================================================================================

!   print*,"ETT,JTT----",ETT,JTT 
!-------- Nutation - Precession POS. - VEL. for Topocentric coordinate 
      CALL NUTATION (-ETT, DJMJD0, JTT, POS, POS2)
      CALL NUTATION (-ETT, DJMJD0, JTT, VEL, VEL2)
      
!     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ETT, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   GP )

!     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ETT, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   GV )
      ETT = ETT
!--------------------------------------------------------------------------------------------
 
      CALL PLEPH ( ETT, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
      CALL PLEPH ( ETT, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( ETT, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector

      CALL PLEPH ( ETT, NTARG, 11, BEV)          ! BEV = Body's heliocentric state vector

!------Update a PVH-vector, discarding the velocity component.

      DT = 0.0D0   
   
      CALL iau_PVUP (DT,PVH,HPE)                 ! HPE Heliocentric position Earth vector 

      CALL iau_PVUP (DT,PVE,EB)                  ! EB Barycentric  position Earth vector

      DO  I=1,3
        BVE(I) = PVE(I,1) * DT + PVE(I,2)
      END DO                                     ! BVE Barycentric velocity Earth vector
                
      CALL iau_PPP (EB,GP,OB)                    ! OB = EB + GP (barycentric position of observer)
 
      CALL iau_PPP (BVE,GV,OV)                   ! OV = BVE + GV (barycentric velocity of observer)

      CALL iau_PVUP (DT,PVB,BP)                  ! BP Body's Barycentric position 

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP Body's Heliocentric position

      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus

!----------------------------------------------------------------------------------------------
  
      CALL L_TIME(ETT,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST,HD0,XET,PV,HE0,THD,EPV) ! RECALL for Topocentric data
   
      BEV0(1,1) = PV(1,1) - OB(1)      ! Body baric.state vector at time T + dT minus Baric. posiction observer
      BEV0(2,1) = PV(2,1) - OB(2)
      BEV0(3,1) = PV(3,1) - OB(3)
      BEV0(1,2) = PV(1,2) - OV(1)
      BEV0(2,2) = PV(2,2) - OV(2)
      BEV0(3,2) = PV(3,2) - OV(3)      

      POS_T = POS_T        
      DIST4 = DIST
      CALL  iau_C2S ( POS_T, RA0, DE0 )       

      RA0 = iau_ANP(RA0)                  ! range angle (0,2!pi) 
     

      Ltime1 = DIST * 8.316746395045D0           ! Light time (min.)
      Mtime1 = DIST * 8.316746395045D0 * 60D0    ! Moon Light time (sec)
      Distmoon = DIST*AU/1000d0/6378.137d0         ! Moon's distance in Earth-radius
      
!====================================================================================
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun (for Topocentric coordinate)

      CALL DEFLIGHT_0(POS_T,BEP,HPE,u01)
      
!--------- Add in deflection due to Earth

      CALL DEFLIGHT_T(u01,BEP,HPE,u01)

      CALL iau_PN(u01,u02,p)          ! p = unit vector in order MU/C^2, u02 = modulus

!---------------------------------------------------------------------------       
!         RELATIVISTIC ANNUAL ABERRATION  (for Topocentric coordinate)

      CALL ABERR (AUDAY,u01,OV,u03)

      CALL iau_PN(u03,MO2,u5)         !Form unit vector of u03 = u5 ) 

      CALL  iau_C2S   ( u5, RA10, DE10 )       
      RA10 = iau_ANP(RA10)                 ! range angle (0,2!pi)

!---------- Form apparent topocentric coordinate RA7,DE7)

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

 
      UJD = UT11 + UT12               ! UJD = UT1  Julian day       
                                      ! IREFR = 1 whit atmoswferic refraction
                                      ! IREFR = 0 no atmosferic refraction 
                        !  INPUT -----! XPO ,YPO polar motion (arcsec)
                                      ! longit,latit (degrees)  
                                      ! HM = heigh (meter) / 0.0 no refr.
                                      ! RA7H, DE7D (hour,degrees)

                                      ! ZD Zenith distance (degrees) 
                        !  OUTPUT-----! AZ Azimuth (degrees)
                                      ! RA8  RA corrected for refraction (degrees)
                                      ! DE8  DEC corrected for refraction (degrees)
!------------------------------------------------------------------------------------------------------------------
  
!       equinox2 = "Apparent Topocentric position no refraction (IAU 2000 EQUINOX BASED)"
       CALL   iau_RXP  ( RNPB, u5, u6 )
! ......... Call for direction cosines to spherical coordinates
       CALL   iau_C2S   ( u6, RA7, DE7 )
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)

         RA7H = RA7 * DR2D / 15D0                ! RA7 in hour
         DE7D = DE7 * DR2D

       CALL SETDT ( DET )
       
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 0, &
                        ZD, AZ, RA8, DE8 )

       RA8 = RA8 * 15D0* DD2R
       DE8 = DE8 * DD2R 
     
      CALL iau_A2TF ( 1 ,RA8 ,sign , IHMSF0 )
      CALL iau_A2AF ( 1 ,DE8 ,sign1 , IDMSF0 )
       
      IF ( HOT < 0d0)then       !carry the hour angle HOT and HOB in the positive direction (West).
       HOT = HOT + 360d0 
      ELSE 
       HOT = HOT
      END IF
      IF ( HOB < 0d0)then 
       HOB = HOB + 360d0 
      ELSE 
       HOB = HOB
      END IF
      
      ZD0 = 90d0 - ZD0
      ZOB = 90d0 - ZOB

!================================   Search  RISE, SET, TRANSIT
       
       RIS = RIS / 24D0
       SET = SET / 24D0
       TRS = TRS / 24D0 
       
       CALL iau_D2TF ( 1, RIS, SIGN, IHMSF1)     
       CALL iau_D2TF ( 1, SET, SIGN, IHMSF2)
       CALL iau_D2TF ( 1, TRS, SIGN, IHMSF3)

         IF ( N0 == 1 )THEN
           AS1 = "* "
           IHMSF1 = 0.0
          ELSEIF ( N0 == 2) THEN
           AS1 = "**"
          ELSEIF ( N0 == 0) THEN
           AS1 = "  " 
         END IF
         IF ( MS == 1 )THEN
           AS2 = "* "
           IHMSF2 = 0.0
          ELSEIF ( MS == 2 )THEN
           AS2 = "**"
          ELSEIF ( MS == 0) THEN
           AS2 = "  " 
         END IF
         IF ( PT == 1 )THEN
           AS3 = "* "
           IHMSF3 = 0.0
           IHMSF0 = 0.0
           IDMSF0 = 0.0
          ELSEIF ( PT == 2) THEN
           AS3 = "**" 
          ELSEIF ( PT == 0) THEN
           AS3 = " " 
         END IF             
       
!-----------------------------------------------------------------------------------       
      
      IF (NTARG == 10) THEN
       DISTKm = DIST1 * AUM
      END IF
!       CALL moon ( ET,IY,MO,ID,DPSI,RA6,DE6,DISTKm,R_sun,D_sun,RA7,DE7,OBL,GD1,PL,HLON,HLAT, &
!                OMEGA,OME0,AGE,JHE,CHI,JAP,APE,DLR,TLong,TLat,LTL,LTB,PA,THL,THB)
!    print*,"AGE---",AGE 

!------ calculating of the lunar phase of New Moon
  
      If (MO ==  1) Day = 15.5
      If (MO ==  3) Day = 15.5
      If (MO ==  4) Day = 15
      If (MO ==  5) Day = 15.5
      If (MO ==  6) Day = 15
      If (MO ==  7) Day = 15.5
      If (MO ==  8) Day = 15.5
      If (MO ==  9) Day = 15
      If (MO ==  10) Day = 15.5
      If (MO ==  11) Day = 15
      If (MO ==  12) Day = 15.5
      If (MO ==  2 .And. IY == 1700 .Or. IY == 1800 .Or. IY == 1900 &
        .Or. IY == 2100) Day = 14
      If (MO ==  2 .And. Mod (IY, 4 ) .NE. 0 )Day = 14
      If (MO ==  2 .And. Mod (IY, 4 ) == 0 ) Day = 14.5

      KU = ((IY + MO / 12D0 - Day / 365.25d0) - 2000d0) * 12.3685d0    ! K0 = aproximate parameter 
      K0 = FLOOR(KU)
!      K1 = K0 + 0.25d0 
!      K2 = K0 + 0.5d0
!      K3 = K0 + 0.75d0
!      print*,"---K0-KU-",K0,KU
      CALL TIMEPHASE (K0,0,JH0)      ! JH0  time for New Moon  
    
      AGE = ET - JH0               !AGE ,age of last New Moon (day)
         
      If (AGE < 0) Then
       AGE = 29.4 + AGE
      ElseIf (AGE > 29.530589) Then
       AGE = AGE - 29.530589
      End if
      IF (NTARG == 10) THEN
       AGE = AGE
      ELSE IF (NTARG /= 10) then
       AGE = 0.00 
      END IF 
!!-----------------------------   OUTPUT  ----------------------------------------
      IF (NTARG == 10) THEN
       GOTO 40
      ELSE IF (NTARG /= 10) THEN
       GOTO 50
      END IF

40    WRITE (*,70)IY, MO, ID,AGE,ILL,IHMSF1,AS1,IHMSF3,AS3,IHMSF2,AS2,azm1,HAL,azm2,IHMSF6,sign2,IDMSF6,IHMSF0,AS3,sign1,IDMSF0,AS3 
      
      GOTO 55

50    WRITE (*,80)IY, MO, ID,IHMSF1,AS1,IHMSF3,AS3,IHMSF2,AS2,azm1,HAL,azm2,IHMSF6,sign2,IDMSF6,IHMSF0,AS3,sign1,IDMSF0,AS3                   
      
55       IF ( JJ == NUM ) THEN
          GOTO 60
         END IF
      JD = JD + 1
    END DO
60        WRITE(*,*)  ' '
          WRITE(*,*)  '-- (0:00:00.0*) - Event occurs the next day'
          WRITE(*,*)  ' '
          WRITE(*,*) "--- END PROGRAM ---" 


70    FORMAT ( x,I4,"/",I2,"/",I2,2x,F6.2,2x,F5.2,2x,I3,2(":",I2.2),".",I1.1,a2,2x,I3,2(":",I2.2),".",I1.1,a2,2X,I3,2(":",I2.2),  &
          ".",I1.1,a2,2X,F6.2,"°",3X,F5.2,"°",3X,F6.2,"°",2x,I3,2(":",I2.2),".",I1.1,3x,a,I3,2(":",I2.2),".",I1.1,3x,I3,2 &
          (":",I2.2),"." ,I1.1,a,2x,a,I3,2(":",I2.2),".",I1.1,a )

80    FORMAT ( x,I4,"/",I2,"/",I2,2x,I3,2(":",I2.2),".",I1.1,a2,2x,I3,2(":",I2.2),".",I1.1,a2,2X,I3,2(":",I2.2),  &
          ".",I1.1,a2,2X,F6.2,"°",2X,F5.2,"°",2X,F6.2,"°",2x,I3,2(":",I2.2),".",I1.1,2x,a,I3,2(":",I2.2),".",I1.1,2x,I3,2 &
          (":",I2.2),"." ,I1.1,a,2x,a,I3,2(":",I2.2),".",I1.1,a )
                     

 END program planets
      
