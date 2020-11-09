 program planets

!     THIS PROGRAN COMPUTES THE APPARENT DIRECTION OF SOLAR SYSTEM
!     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM
!
!     gfortran -w planet430F.f90 1919.for jpl_430.f sky9.for
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
      PARAMETER (DJMJD0=2400000.5d0)  
      DOUBLE PRECISION DJ00
      PARAMETER (DJ00=2451545d0)  


      
      INTEGER IY, MO, ID, IH, IM, J, I,mode,IREFR,SGN,ios,modetime,NVS,NTARG,NCENT,ISTEP,NSTEP,EOP_final
      INTEGER IHMSF1(4),IHMSF2(4),IHMSF3(4),IHMSF(4),IDMSF(4), IDMSF0(4),IHMSF4(4),IHMSF5(4),IHMSF6(4),Rtoday,YN
      INTEGER FINAL_,JW,EQTM,IY4, IM4, ID4,IDMSF1(4),d_d,h_h,m_m,conta,Dtoday,Mtoday,Ytoday,MJDtoday
      INTEGER IY1, IM1, ID1,IY2, IM2, ID2, IY3, IM3, ID3,IHMSF7(4),IHMSF8(4),NR,MS,PT,LX,LM,LN
      INTEGER IHMSF10(4),IHMSF11(4),IHMSF12(4),IHMSF13(4),IHMSF14(4),IHMSF15(4),N1,N2,N3,numday
      integer,dimension(8) :: VALUES
      DOUBLE PRECISION SEC,XP,YP,DUT,TT,TCB1,TCB2,T,AU,PS,DEL,HA,TDB,Ltime,TDB1,TDB2,TDB0,OBL,s_s
      DOUBLE PRECISION  UT11, UT12, UT,TAI1, TAI2, TT1, TT2, TCG1, TCG2,RVETT,xs,ys,pas,UTC1, UTC2,UTC 
      DOUBLE PRECISION DDP80, DDE80, DX00, DY00, DX06, DY06,DELTA_T,DE_T,u01(3),Rdot,GN,Mtime,newtime
      DOUBLE PRECISION DATE_, TIME, DAT,TCB,TCG,POS2(3),DOT,TAI, TUT, UT1, DP80, DE80, MJD0,DJM,JD0,day0
      DOUBLE PRECISION RA1,DE1,RA2,DE2,RA3,DE3,RA4,DE4,RA5,DE5,RA6,DE6,RA7,DE7,RA8,DE8,RA9,DE9,DJM0
      DOUBLE PRECISION RA10,DE10,RMASSE,SUNLONG,COSELON,EDEC,EVET,GLON,GLAT,DIST,DIST4,ECC,JDtoday         
      DOUBLE PRECISION DPSI, DEPS, EPSA, ERA, DP00, DE00, DDP00, DDE00 ,ELON, DELTAT, EARTHVEC,R_sun,JH0 
      DOUBLE PRECISION GAST,LONG ,LAT,HM,DTR,MO1,FD,julian,DSQRT,HLON,HLAT,radius,T0,ET,XET,BIS,MJ,AX
      DOUBLE PRECISION  X, Y, S , DT, ELONG, u02, MO2, longit, DELTA, DS, AUDAY, XPO, YPO, DIAM,AE,IYdate      
      DOUBLE PRECISION latit,UJD,MJD,ANGLE,PANGLE,PDIAM,AZ,height,ZD,MAGN,AUM,PVV(2),LAST,fract_year
      DOUBLE PRECISION RA7H, DE7D,JD,RVET,RVET0,EPOS,EPSA0,MM,LEarth,LSole,PEarth,iau_FAPA03,GELO_,GELA_
      DOUBLE PRECISION DX0,DY0,UT_0,TT_0,GAST0,PL,SD,RIS,SET,TRS,AZM1,AZM2,HAL,MUC,BP,de_S,UTC0,UTC01,UTC02
      DOUBLE PRECISION C,MU,u2,RAR,DER,K,THD,TGD,ILL,COSI,DET,HD,GD,HMOD,DIST1,DIST2,DIST3,EQTX,BODY_sun(3)
      DOUBLE PRECISION RC2IT(3,3),RC2I(3,3),RPB(3,3),RNPB0(3,3),RN(3,3),RNPB(3,3), RP(3,3),RBP(3,3),iau_EQEQ94
      DOUBLE PRECISION PV(3,2),PVH(3,2),PVB(3,2),BEV1(3,2),BEV(3,2),RB(3,3),GMB(3,2),PVE(3,2),HSV(3,2)
      DOUBLE PRECISION HE0(3,2),SUNB(3,2),EPV(3,2),GPB(3),HV(3),POS(3),VEL(3),XYZ(3),V1(3),RA66,RA,RA0
      DOUBLE PRECISION EB(3),u0(3),u03(3),u1(3),u3(3),u4(3),u5(3),u6(3),U(3),V(3),V2(3),p(3),RAsun,DE0 
      DOUBLE PRECISION VEL2(3),VEL3(3),POS3(3),POS_H(3),GE(3),HPE(3),BVE(3),BEP(3),TSL,EO ,BGV(3,2),NUT(4)
      DOUBLE PRECISION BODY_G(3),BODY_H(3),HE(3),RS(3),OB(3),OV(3),GP(3),GV(3),POS_T(3),EO0,eta,PI_
      DOUBLE PRECISION iau_ANP,iau_GMST00,iau_S06,iau_DTDB,iau_GMST06,iau_GST06A,iau_GST00A,CHI,p_,A_,B_    
      DOUBLE PRECISION iau_ERA00, iau_EE00,iau_EE06A,iau_EORS,iau_EE00A,dpsi_,deps_,deldot,C_,OBL0,JDT,TIME62
      DOUBLE PRECISION AOB,ZOB,HOB,DOB,ROB,PHPA,TC,RH,WL,ZOT,AOT, HOT, DOOT,ROT,EQEQ,GMST0,ZD0,GEL,atlon  
      DOUBLE PRECISION Sun_lon,Earth_lon,TDBJ,iau_FAE03,CORR,EQT,EQTS,EQT0,LOD,EQTD,AGE,DJ2,Rsun,GST  
      DOUBLE PRECISION BODY_H1,GD1,HD1,XET1,PV1,HE01,THD1,E_P(3,2),D_sun,DISTKm,OMEGA,OME0,OBL6
      DOUBLE PRECISION JAP,DLR,DJ3,PH,TLong,TLat,LTL,LTB,PA,iau_OBL06,LON0,LAT0,T_0,a0,b0,c0,THL,THB,LAV(6) 
      DOUBLE PRECISION LOGE,BGE1,Distmoon,R_sun1,GELO,GELA,GELOAP,GELOD,GELAD,GELOAPD,VC,CVE(3),BM1,BCV,er
      DOUBLE PRECISION Tcy,PHE,KA,Lsun,Bsun,SLO,SLA,SGP(3),DLO,DLA,GELAAPD,GELAAP,AZ1,SGE(3,2),iau_GMST82
      DOUBLE PRECISION RVET1,Rdot1,Ltime1,deldot1,Mtime1,HD0,BEV0(3,2),HSVT(3,2),PVBT(3,2),LEA(6),iau_OBL80,EE
      DOUBLE PRECISION sine,cose,nsine,ec(3,3),HEPV(3,2),ut1_notid,JH1,DJ4,JHP,JHU,DJ5,DJ6,SSH(3,2),SGV(3,2),tempo
      DOUBLE PRECISION SAB(3,2),JUB(3,2),B(8,3),ET0(2),uu(3),BCTW,ECTW,BNTW,ENTW,BASTW,EASTW,step,XX
      DOUBLE PRECISION  R(6)
      DOUBLE PRECISION  SS(3)
      DOUBLE PRECISION  VALS(400)
      DIMENSION body(17),radius(10),dayweek(7)
      
      CHARACTER*8 DATE
      CHARACTER*10 TIMEX
      CHARACTER*6 ZONE,data_b
      CHARACTER(len=8):: NAMS(400),TX,choose,yes,no,distance,EL
      CHARACTER(len=16):: XDIAM,dayweek,Dday,IDEM 
      CHARACTER(len=26):: body,planet,LDUT,LDET,EQ_label,Slod
      CHARACTER(len=1):: sign,aS,lS,bS,update,NEW,str,AS1,AS2,AS3
      CHARACTER(len=95):: equinox,equinox2,equinox3,commentR,commentS,commentT,commentU,labelL
      CHARACTER(len=95):: commentM,commentN,commentA,commentP
      CHARACTER(len=7):: APE 
      CHARACTER(len=14):: comment 
      common SGE,RVET0,HSV,BEP                   

      DATA body / "MERCURY","VENUS ","EARTH ","MARS ","JUPITER","SATURN ","URANUS ","NEPTUNE","PLUTO ", &
                  "MOON ","SUN ","SOLAR-SYSTEM BARYCENTER","EARTH-MOON BARYCENTER","NUTATIONS", &
                  "Lunar Euler angles"," Lunar libration","TT - TDB" / 
                  
      DATA radius / 2440D0, 6051.8D0, 6378.137D0, 3397D0, 71492D0, &
                    60268D0, 25559D0, 24766D0, 1151D0, 1737.53D0 /

      DATA dayweek / "SUNDAY","MONDAY","TUESDAY","WEDNESDAY","THURSDAY","FRIDAY", &
                     "SATURDAY" /

      WRITE (*,*)"  **************************************************************************"
      WRITE (*,*)"  *     Last release of Planet, renamed planet430F, to avoid mismatches    *"
      WRITE (*,*)"  *     New SOFA routines added into 1919.for , JPL DE430t database,       *" 
      WRITE (*,*)"  *     IERS FINAL_ data (IAU2000), and EOP C01 1900-now modified          *"           
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *      --------------------planet430F.f90 -------- (Oct.25,2020) --------*"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     THIS PROGRAN COMPUTES THE APPARENT DIRECTION OF SOLAR SYSTEM       *"
      WRITE (*,*)"  *     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM      *"
      WRITE (*,*)"  *            for this purpose we use JPL DE430t database                 *"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     Compiling command:                                                 *"
      WRITE (*,*)"  *     gfortran planet430F.f90 sky9.for 1919.for jpl_430.f                *"
      WRITE (*,*)"  *     which generates the executable file: a.out                         *"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     Range of calendar dates:                                           *"
      WRITE (*,*)"  *      a) Not less 1900-03-01                                            *"
      WRITE (*,*)"  *      b) Not over MJD of line 150 (or nearby), i.e. FINAL_ = 58492 of   *"
      WRITE (*,*)"  *         main file planet430F.f90.                                      *"
      WRITE (*,*)"  *      c) Ephemerides computations may be forced up to 2100/01/01 upper  *"                 
      WRITE (*,*)"  *         limit of JPL database,but with less precision.                 *"
      WRITE (*,*)"  *         ------------------------------------------------------------   *" 
      WRITE (*,*)"  *         For FINAL_-iers.txt updating, please refer to : www.iers.org/  *"
      WRITE (*,*)"  *       - Data/Products - Standard EOP data files - Finals data IAU2000 -*"  
      WRITE (*,*)"  *         and change the FINAL_ value accordingly, say MJD shown 3 lines *"
      WRITE (*,*)"  *         prior the last full row.                                       *" 
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  **************************************************************************"
      WRITE (*,*)
                   
      WRITE (*,*)" ****  We recommend that you update at least every 15 days, checking ****  "
      WRITE (*,*)" ****      the FINAL_ dates of the respective IERS database           ****  "
      WRITE (*,*)" ========================================================================= "

      WRITE (*,*)" Enter the last updated dates (MJD) of the IERS databases "
      WRITE (*,*)" Finals data (IAU2000) and EOP14 C04 (IAU2000A) or uses   "
      WRITE (*,*)" default dates (1 or 2)"  

      READ*,YN
                 
      IF ( YN == 1) THEN
        PRINT*, "  ENTER  FINAL_ ... " 
        READ*,FINAL_
        PRINT*, "  ENTER EOP_final ..."
        READ*,EOP_final
        GOTO 3                     
      ELSEIF ( YN == 2) THEN
        GOTO 2
      END IF                          

2     continue    
      FINAL_ = 59338       ! 2021/05/04   reported in MJD IERS FINAL_ data (for update) 

      EOP_final = 59059  ! 2020/07/29 jump from IERS EOP14 C04 to IERS FINAL_ data (IAU2000).       

3     continue            
10    WRITE(*,*)
      WRITE(*,*) " Input : Observer geodetic coordinate (degrees) " 
      WRITE(*,*) " Longitude  +/-yyy.xxx..(negat. at WEST)  " 
        READ *, longit
       IF (longit > 180D0 .OR. longit < -180D0) THEN
        PRINT*," *** ERROR *** Longitude out of range "
        GOTO 10
       END IF  
12      WRITE(*,*) " Latitude  +/-yy.xxx...(max. +/- 70 deg. - posit. at NORD)  " 
        READ *, latit
       IF (latit > 75D0 .OR. latit < -75D0) THEN
        PRINT*," *** ERROR *** Latitude out of range "
        GOTO 12
       END IF  
      WRITE(*,*) " Height (in meter above sea level)     "
      READ *, height 
      WRITE(*,*)" TC  ambient temperature at the observer (deg C)"  
      READ*, TC
      WRITE(*,*) "==========================================="

      WRITE (*,*)
      100   WRITE (*,*)" Choose Time input!  'UTC' = 1  or 'TT' = 2 ( Prior the 1962.0 use 'UT1'= 1 or 'TDB' = 2 )" 
      READ *, modetime

      IF(modetime > 2)THEN
       GOTO 100
      ELSE IF(modetime == 1) THEN
       TX = "UTC/UT1"
       GOTO 200
      ELSE IF (modetime == 2) THEN
       TX = "TT/TDB"
       GOTO 200 
      END IF 
 
200   WRITE (*,*)
      WRITE (*,*) "       Date and Time " ,TX   
      WRITE (*,*) "   =============================="
210   WRITE (*,*)
      WRITE (*,*) " Input : Year , Month, Day (yyyy,mm,dd) "
      READ *, IY, MO, ID
       IF(IY < 1900  .OR. IY > 2100 ) THEN              
        WRITE(*,*)"*****  YEAR out of range ! *****"
        GOTO 210
       ELSE IF(MO < 1  .OR. MO > 12 ) THEN              
        WRITE(*,*)"*****  ERROR ! WRONG MONTH ! *****"
        GOTO 210
       ELSE IF(ID < 1  .OR. ID > 31 ) THEN              
        WRITE(*,*)"*****  ERROR ! WRONG DAY ! *****"
        GOTO 210
       ELSE  
        GOTO 215
       END IF 
                
215   WRITE (*,*)
      WRITE (*,*) " Input : Hour, Minutes, Seconds (hh,mm,ss.xxx..) "
      READ *, IH, IM, SEC
       IF(IH < 0  .OR. IH > 23 ) THEN              
        WRITE(*,*)"*****  ERROR ! WRONG HOUR ! *****"
         GOTO 215
       ELSE IF(IM < 0  .OR. IM > 59 ) THEN              
        WRITE(*,*)"*****  ERROR ! WRONG MINUTS ! *****"
        GOTO 215
       ELSE IF(SEC < 0  .OR. SEC > 60 ) THEN              
        WRITE(*,*)"*****  ERROR ! WRONG SECONDS ! *****"
        GOTO 215
       END IF 

      conta = 0 
5100  TIME = ( 60D0*(60D0*DBLE(IH) + DBLE(IM)) + SEC ) / 86400D0
      conta = conta
      CALL iau_DAT ( IY, MO, ID, TIME, DAT, J ) 
       IF(J == -2 .OR. J == -3 .OR. J == -4) THEN
        WRITE(*,*) " BAD DATE ! (out of range) "
        GOTO 210
        ELSE 
        GOTO 220
       END IF
       
220   continue      
           
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
      TSL = 25.0                           ! Is the approximate sea-level air temperature in K"
      RH = 0.5                             ! relative humidity at the observer (range 0-1)
      WL = 0.55                            ! wavelength (micrometers)
      KA = 20.49552d0                      ! Costant of aberration 
!---------------------------------------------------------------------------------------  
      
!------- TODAY date    
   CALL DATE_AND_TIME(DATE, TIMEX, ZONE, VALUES) 
      Ytoday = VALUES(1)
      Mtoday = VALUES(2)
      Dtoday = VALUES(3)

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

!---------- Compute JD (Julian day) today .--------------------------

      GN = (100D0*Ytoday+Mtoday-190002.5)                    !I find the sign of the numerical expression
      IF (GN < 0D0) THEN
       SGN = -1
      ELSE IF (GN > 0D0) THEN
       SGN = +1
      END IF

      JDtoday = 367*Ytoday - (7*(Ytoday+((Mtoday+9)/12)) / 4)+(275*Mtoday/9) + Dtoday + 1721013.5d0  & 
             - 0.5d0 * SGN +0.5d0                ! JD TODAY  

      MJDtoday = INT(JDtoday - DJMJD0)           ! MJD today 

!------------------------------------------------------------------------------------------------------
      call iau_CAL2JD( IY, MO, ID, DJM0, DJM, J)
      MJD = DJM + TIME                          ! UTC or TT in MJD (Julian Day modified) at date     

      JD = DJMJD0 + MJD                         ! UTC or TT in Julian Day at date
             
!----------Compute the day of week---------------------------------

      JW = INT((JD + 1.5D0) - 7D0 * Int((JD + 1.5D0) / 7D0)) 
      
      DO I = 1,7
       IF ( JW == I-1) THEN
        Dday = dayweek(I)
       END IF
      END DO  
!-------------------------------------------------------------------    
 
      IF (MJD > 37666) then           ! MJD 37666 = final database EOP C01
       GOTO 178                       ! MJD 37667 = 1962/01/03 ... start database EPO14
      ELSE IF (MJD <= 37667)THEN      ! MJD 37665 -- 37671 initial range interpolation  
       GOTO 1640      
      END IF

178   WRITE(*,*)
      WRITE(*,*) "   Choose mode Equinox Based (1)or CIO Based (2) "    
      WRITE(*,*) "   =============================================="

      WRITE(*,*) "  Input number : 1 or 2 "
      READ *, mode
      WRITE(*,*) 


      WRITE(*,*)"========================================================================"              
      WRITE(*,*)"                                                                        "
      WRITE(*,*)"        THE NUMBERING CONVENTION FOR 'NTARG' IS:                        "
      WRITE(*,*)"                                                                        "
      WRITE(*,*)"     1 = MERCURY         10 = MOON                                      "
      WRITE(*,*)"     2 = VENUS           11 = SUN                                       "
      WRITE(*,*)"     3 = EARTH           12 = SOLAR-SYSTEM BARYCENTER                   "
      WRITE(*,*)"     4 = MARS            13 = EARTH-MOON BARYCENTER                     "
      WRITE(*,*)"     5 = JUPITER         14 = Nutations (Longitude and Obliquity) (*)   "
      WRITE(*,*)"     6 = SATURN          15 = Lunar Euler angles: phi, theta, psi (**)  "
      WRITE(*,*)"     7 = URANUS          16 = Lunar angular velocity: omegax,           "
      WRITE(*,*)"     8 = NEPTUNE              omegay, omegaz  (*)           "               
      WRITE(*,*)"     9 = PLUTO           17 = TT - TDB (respect geocenter) (**)         "
      WRITE(*,*)"                                                                        "
      WRITE(*,*)"     (*)  No parameter in this ephemeris database.                                "
      WRITE(*,*)"     (**) Parameters calculated internally.                             "
      WRITE(*,*)"                                                                        "
      WRITE(*,*)"========================================================================"
      WRITE(*,*)

      WRITE(*,*) 'Enter Target Body ' 
      READ *, NTARG

      DO I = 1,18
       IF ( NTARG == I) THEN
        planet = body(I)
        GOTO 166
       END IF
      END DO  
  
166 continue        
      
      IF( MJD > FINAL_ ) THEN
       de_S = (MJD - 51544.5D0) * DS                ! Time in sec. from J2000 (JD 2451545.0 TT)
       MM = 6.239996D0 + de_S * 1.99096871D-7
       DELTA = 32.184 + DAT + 0.001657D0 *sin(MM + 0.01671D0 * sin(MM))
      END IF                                                 
!      print*,"DELTA",DELTA 
!----------------------------------------------------------------------iau_GST00A
! Compute TT and trasform TT to UTC
      
      IF (modetime == 1) GOTO 175
      IF (modetime == 2) GOTO 170 
!----------TT (MJD).

170   TT = MJD
!     CALL iau_CAL2JD ( IY, MO, ID, DJMJD0, DATE, J )              ! input TT (IY,MO,ID)
!      TT = DATE + TIME                                             ! TT value (MJD)
                                                                                   
!----------UTC,TAI (MJD)      
       
      UTC = TT - (32.184D0 + DAT)/DS

      CALL iau_JD2CAL ( DJMJD0, UTC, IY, MO, ID, FD, J )
      CALL time_hms ( FD, IH, IM, SEC)

!---------Transform into internal format.
175   CALL iau_DTF2D ( "UTC", IY, MO, ID, IH, IM, SEC, UTC1, UTC2, J )   !UTC1 + UTC2 = Julian Day
      IF ( J.NE.0 ) STOP                                                 !UTC1 is the Julian Day number and
      CALL iau_DTF2D ( "UTC", IY, 1, 1, 0, 0, 0D0, UTC01, UTC02, J )     !UTC2 is the fraction of a day.   
    
      UTC = UTC1 + UTC2-DJMJD0
      DATE_ = UTC1 - DJMJD0                                               ! DATE_ is the MJD number
      UTC0 = UTC01 + UTC02 
                                                   ! UTC JD at 0.0h YYYY/01/01  
!----------Call the interpolation routine for per XP,YP,DUT,dX,dY
      
      CALL iers_calc(UTC,DATE_,FINAL_,EOP_final,XPO,YPO,DUT,LOD, &
                     DX0,DY0,ut1_notid)
      
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

      IF( MJD > FINAL_ ) THEN             !
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
!      UTC = DATE_+UTC2
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
           EPSA0, RB, RP, RPB, RN, RNPB0 )
!-------------------------------------------------------------------------------------------------------

!--------- Greenwich apparent sidereal time (IAU 2000) for 0.0 h UT1.
      
      GAST0 = iau_GST06A ( UT_0, 0D0, TT_0, 0D0 )
                     
!--------- Greenwich apparent sidereal time (IAU 2000-2006).

      GAST = iau_GST06A ( UT11, UT12, TT1, TT2 )

!--------- Local apparent sideral time (IAU 2000)

      LAST = GAST + LONG

!--------- Greenwich mean sidereal time (IAU 2000). 
     
      GMST0 = iau_GMST00 ( UT11, UT12, TT1, TT2 )

!-------------------------------------------------------------------------------------------------------      
      EO = ERA - GAST                                ! EO0  Equation of the origin.
      EO0= iau_EORS(RNPB,S)
      
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
  
!       ET0(1) = TDB1
!       ET0(2) = TDB2 
    
!      CALL PLEPH ( ET, 14, 0, NUT)              ! 1980 IAU nutation angles. (Long,Obliq,Long. rate, Obliq rate ) 
      
      CALL PLEPH ( ET, 15, 0, LEA)              ! Lunar Euler angles: phi, thETa, psi 
    
      CALL PLEPH ( ET, 17, 0, PVV)              ! PVV(1) = TT-TDB in sec./ PVV(2)=rate of change of TT-TDB in sec/day
   
      CALL PLEPH ( ET, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
      CALL PLEPH ( ET, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( ET, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector

      CALL PLEPH ( ET, NTARG, 11, BEV)          ! BEV = Body's heliocentric state vector 
   
      CALL PLEPH ( ET, 11, 12, SUNB )           ! SUNB = Sun barycentric state vector    

      CALL PLEPH ( ET, 11, 3, SGE )             ! SGE = Sun geocentric state vector  
      
      CALL PLEPH ( ET, NTARG, 3, BGV)           ! BGV = Body geocentric state vector

      CALL PLEPH ( ET, 6, 12, SAB)              ! SAB = Saturn barycentric state vector
      
      CALL PLEPH ( ET, 5, 12, JUB)              ! jUB = Jupiter barycentric state vector      
         
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

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP = Body's Heliocentric position
   
      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus
  
!      CALL iau_PVMPV(SUNB,PVE,SGV)               ! SGV  Sun Geocentric mean state vector

      CALL iau_PVMPV(PVB,PVE,GMB)                ! GMB  Body's geocentric mean state vector
      
      CALL iau_PVUP (DT,GMB,GPB)                 ! GPB  Body's geocentric mean J2000 position    

      CALL  iau_C2S   ( GPB, RA4, DE4 )       
      RA4 = iau_ANP(RA4)                         ! range angle (0,2!pi)
      
      DIST2 = DSQRT(GPB(1)*GPB(1)+GPB(2)*GPB(2)+GPB(3)*GPB(3))     ! Geocentric Distance 

      EARTHVEC = DSQRT(HPE(1)*HPE(1)+HPE(2)*HPE(2)+HPE(3)*HPE(3))  ! Distance Sun - Earth
      
!---------- Astrometric coordinate      
    
      CALL L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,HE0,THD,EPV) ! vedi EPV

      BEV1 = PV - PVE                            ! Geocentric state vector
           
!-----------Form  Geocentric astrometric position RA3,DE3
! ......... Call for direction cosines to spherical coordinates

      CALL  iau_C2S ( BODY_G, RA3, DE3 )  
     
      RA3 = iau_ANP(RA3)                         ! range angle (0,2!pi)

!----------- Transform RA3, DE3 in astrometric-topocentric 

      
!--------------Form RVET0 (vector ray mean J2000 )and  Rdot.
      
!      RVET0 = DSQRT(BODY_H(1)**2D0+BODY_H(2)**2D0+BODY_H(3)**2D0)

      RVET0 = HD 
      
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

!----------- HEPV = HELIOCENTRIC ECLIPTIC POS-VEL BODY

      OBL0 = OBL 
      OBL6 = iau_OBL06(TT1,TT2)
      cose = COS(OBL0)
      sine = SIN(OBL0)
      nsine = -1d0 * sine

      ec(1,1) = 1d0;  ec(1,2) = 0d0;   ec(1,3) = 0d0
      ec(2,1) = 0d0;  ec(2,2) = cose;  ec(2,3) = sine
      ec(3,1) = 0d0;  ec(3,2) = nsine; ec(3,3) = cose

      CALL iau_RXPV( ec, PVB, HEPV)!( ec, HSV, HEPV)

!-----------LON0 = ELIOCENTRIC ECLIPTIC LONGITUDE ; LAT0 = ELIOCENTRIC ECLIPTIC LATITUDE  (J2000)
           
      CALL GEPV(BEV,OBL0,1,LON0,LAT0)
     
!----------- Reduce Heliocentric Coordinate mean J2000 to mean of the date.(Meeus Astronomical Algorithms; rigorous method pag.136)
     
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
        ILL = 1.0D0 *100D0
      END IF
!===================================================================================================
      CALL L_TIME(ET,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST,HD0,XET,PV,HE0,THD,EPV) ! RECALL for Topocentric data
   
      BEV0(1,1) = PV(1,1) - OB(1)      ! Body baric.state vector at time T + dT minus Baric. posiction observer
      BEV0(2,1) = PV(2,1) - OB(2)
      BEV0(3,1) = PV(3,1) - OB(3)
      BEV0(1,2) = PV(1,2) - OV(1)
      BEV0(2,2) = PV(2,2) - OV(2)
      BEV0(3,2) = PV(3,2) - OV(3) 
     
!---------- ( POS_T = Topocentric Astrometric position vector )

      POS_T = POS_T        
      DIST4 = DIST
      CALL  iau_C2S ( POS_T, RA0, DE0 )          ! RA0, DE0: Topocentric Astrometric coordinate      

      RA0 = iau_ANP(RA0)                         ! range angle (0,2!pi) 
 
      Ltime1 = DIST * 8.316746395045D0           ! Light time (min.)
      Mtime1 = DIST * 8.316746395045D0 * 60D0    ! Moon Light time (sec)
      Distmoon = DIST*AU/1000d0/6378.137d0       ! Moon's distance in Earth-radius

      IF(NTARG == 10) THEN
      DIAM = (1737.53D0 * 2D0 * DR2AS /( DIST * AU )* 1000d0/60d0) ! arcmin-Topoc.Moon Diam. 
      END IF
      IF ( NTARG == 11 ) THEN
       DIAM = 1919.29D0 / (DIST * 60d0)                            ! arcmin-Topoc. Sun Diam.
      END IF
      
!--------------Form RVET1 (vector ray mean J2000 )and  Rdot.
      
      RVET1 = DSQRT(POS_H(1)**2D0+POS_H(2)**2D0+POS_H(3)**2D0)
      
      Rdot1 = (HE0(1,2)*HE0(1,1)+HE0(2,2)*HE0(2,1)+HE0(3,2)*HE0(3,1))/ HD0  

      Rdot1 = Rdot1*AU/86400D0/1000D0                 ! Rdot Km/sec.              
   
!-------------Form " deldot "            
       
      DOT = (BEV0(1,2)*BEV0(1,1)+BEV0(2,2)*BEV0(2,1)+BEV0(3,2)*BEV0(3,1))/DIST   

      deldot1 = DOT*AU/86400D0/1000D0                ! deldot Km/sec. 

!----------------------------------------------------------------------------

      IF(NTARG > 7 .AND. NTARG <= 9 ) THEN
        GOTO 555
       ELSE IF (NTARG <= 7 .OR. NTARG >= 10) THEN
        GOTO 666
      END IF 

!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun,Saturn and Jupiter ( for apparent geocentric coordinate)
 
555   B(1,1) = 0.00028574D0      ! Saturn parameter
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
 
      GOTO 777
666   CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
777   CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

      CALL  iau_C2S   ( p, RA5, DE5 )       
      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
       
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent GCRS coordinate)

      CALL iau_SXP ( VC,BVE,CVE)    ! CVE body velocity vector in unit of C 
      
      BCV = SQRT(CVE(1)**2d0+CVE(2)**2d0+CVE(3)**2d0) !form module of vector CVE

      BM1 = SQRT(1d0-abs(BCV)**2d0)            ! for use with routine SOFA iau_AB()  

      CALL ABERR (AUDAY,u1,BVE,u0)
     
      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 

   
! ......... FINAL_ GCRS (RA10,DE10) of body 
! ......... Call for direction cosines to spherical coordinates
 
      CALL  iau_C2S   ( u3, RA10, DE10 )       
      RA10 = iau_ANP(RA10)                 ! range angle (0,2!pi)
   
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
     
!========================================================================================================
!--------- Compute  -TDB0- JD Time to 0.0 h TDB the day (for routine TRANSIT)


       TDB0 = INT(MJD) + DJMJD0  
      IF(NTARG == 10) THEN 
       AUM = AU/1000D0
       PL = ASIN(6378.137/(DIST1 * AUM)) * DR2D  ! Lunar horizontal parallax 
       SD = (DIAM / 60D0 ) / 2D0               ! Lunar semidiameter in degs
      END IF 
               
       CALL TRANSIT(NTARG,TDB0,GAST0,AUDAY,DELTAT,PL,SD,LAT,LONG,DX00,DY00,RIS,SET,TRS,  &
            AZM1,AZM2,HAL,As,NR,MS,PT,LX,LM,LN) 

     IF(NTARG == 11) THEN
      CALL TWILIGHT(NTARG,TDB0,GAST0,AUDAY,DELTAT,LAT,LONG,DX00,DY00,Ha, &
           BCTW,ECTW,BNTW,ENTW,BASTW,EASTW,N1,N2,N3)
     END IF  
      
!----------------------------------------------------------------------------------------------

!---------- Form apparent geocentric coordinate (RA6,DE6)       

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

      IF (mode == 1) then 
       CALL   iau_RXP   ( RNPB, u3, u4 )
       equinox = "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2000 EQUINOX BASED)" 
      ELSE IF (mode == 2) then 
       CALL   iau_RXP   ( RC2I, u3, u4 )
       equinox = "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2006 CIO BASED)"  
      END IF 
! ......... Call for direction cosines to spherical coordinates CIRS (RA66 , DE6)

      CALL   iau_C2S   ( u4, RA66, DE6 )
            
      IF (mode == 1) then
         RA6 = RA66 - (0.0524D0 * AS2R)    ! 52mas = Diff. RA origin
         RA6= iau_ANP(RA6)                   ! range angle (0,2!pi)
         DE6 = DE6 
       ELSE IF (mode == 2) then                     
         RA6= iau_ANP(RA66-EO)         
      END IF

         DIST3 = DIST1                              ! Distance 
      
!-------------------------------------------------------------------------------------------------

       IF( MJD > FINAL_ ) THEN
        DET = DELTA
        LDUT =" value not determinated"
        LDET =" value aproximate "
        XPO = 0D0
        YPO = 0D0
       ELSE IF (MJD <= FINAL_) THEN
        LDUT =" " 
        LDET =" "
        XPO = XPO
        YPO = YPO
       END IF 

      if (LOD == 0 )then
       Slod = " value not available "
      else
       Slod = " "
      end if        

!---------- Form apparent topocentric coordinate RA7,DE7)

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

 
      UJD = 2400000.5D0 + UT1         ! UJD = UT1  Julian day       
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
  
     IF (mode == 1) then
       equinox2 = "Apparent Topocentric position no refraction (IAU 2000 EQUINOX BASED)"
       CALL   iau_RXP  ( RNPB, u5, u6 )
! ......... Call for direction cosines to spherical coordinates
       CALL   iau_C2S   ( u6, RA7, DE7 )
         RA7 = RA7 !- 2.54527183D-7       ! 2.54527183D-7 Rad. = Diff. RA origin
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)
      
         RA7H = RA7 * DR2D / 15D0                ! RA7 in hour
         DE7D = DE7 * DR2D
       CALL SETDT ( DET )
  
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 0, &
                        ZD, AZ, RA8, DE8 )  
        AZ = AZ
       ZD0 = ZD
       HOT = (LONG + GMST0 - RA7)* DR2D          ! HOT = Hour angle 
       RA8 = RA8 * 15D0* DD2R -(0.052D0 * AS2R)    ! 52mas = Diff. RA origin
       DE8 = DE8 * DD2R 
       
      ELSE IF (mode == 2) then

       equinox2 = "Apparent Topocentric position no refraction (IAU 2006 CIO BASED)"
!       CALL   iau_RXP  ( RC2I, u5, u6 )
!       CALL   iau_C2S   ( u6, RA7, DE7 )
!         RA7= iau_ANP(RA7 - 2.54527183D-7)                       ! range angle (0,2!pi)

!      INPUT CIRS (RA66,DE6)     

       CALL iau_ATIO13 ( RA66, DE6, UTC1, UTC2, DUT, LONG, LAT, HM, &      ! AOT = Azimut (CIO Based)
                XP, YP, 0D0, 0D0, 0D0, 0D0,AOT, ZOT, HOT, DOOT, ROT, J )  ! ZOT = Dist.Zenith(CIO Based)
       IF (J .NE. 0 ) STOP                                                ! HOT = Hour angle (CIO Based)

       AZ = AOT * DR2D;  ZD0 = ZOT * DR2D; HOT = HOT *DR2D 
       RA8 = ROT 
       DE8 = DOOT
      END IF        
!
!--------Transform topocentric right ascension and declination to zenith distance and azimuth.
!        and compute atmosferic refraction.
      
     IF (mode == 1) then
              equinox3 = "Topocentric observed position with refraction.(IAU 2000 EQUINOX BASED)"

!         RA7H = RA6 * DR2D / 15D0                ! RA7 in hour
!         DE7D = DE6 * DR2D

       CALL REFDAT (PHPA,TC,0.0,WL)
       CALL SETDT ( DET )
       
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 1, &
                        ZD, AZ, RA9, DE9 )
      
       RA9 = RA9 * 15D0* DD2R - (0.052D0 * AS2R)    ! 52mas = Diff. RA origin
       DE9 = DE9 * DD2R 
       HOB = (LONG + GMST0 - RA9)* DR2D
       ZOB = ZD 
       AZ1 = AZ

     ELSE IF (mode == 2) then
              equinox3 = "Topocentric observed position with refraction.(IAU 2006 CIO BASED)"
      
!      INPUT CIRS (RA66,DE6)
       CALL iau_ATIO13 ( RA66, DE6, UTC1, UTC2, DUT,LONG, LAT, HM, XP, &   !RA7,DE7
                        YP,PHPA,TC,RH,WL, AOB, ZOB, HOB, DOB, ROB, J ) 
       IF (J .NE. 0 ) STOP                                             
      AZ1 = AOB *DR2D; ZOB = ZOB * DR2D; HOB = HOB *DR2D        ! AOB = Azimut (CIO Based)                     
      RA9 = ROB                                                 ! ZOB = Dist.Zenith(CIO Based)
      DE9 = DOB                                                 ! HOB = Hour angle (CIO Based)
                                                                !PHPA = pressure at the observe(hPa=mB)
     END IF                                                     !  TC = ambient temperature(deg C)
                                                                !  RH = relative humidity (range 0-1)
                                                                !  WL = wavelength (micrometers
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

!====================================================================================================================
      CALL L_TIME(ET,11,EB,BVE,HMOD,BODY_sun,BODY_H1,GD1,HD1,XET1,PV1,HE01,THD1,E_P) ! search RA and DEC Sun (R_sun,D_sun)  

      CALL RA_SUN (mode,AUDAY,BODY_sun,RNPB,RC2I,BEP,HPE,BVE,MO1,R_sun,D_sun)          
       
!--------- Compute the Equation of time.(EQT)------------------------------------------------

      TDBJ = ((TDB1+TDB2) - 2451545D0) / 365250D0     ! TDB by 2000.0 JD in millenium  
      
      Sun_lon = 280.4664567D0+360007.6982779D0*TDBJ &   ! Sun mean Longit. at time TDBJ in degrees
           +0.3032028D0*TDBJ**2D0+TDBJ**3D0/49931D0 &   ! Meuss equation use TDB in millenium
           -TDBJ**4D0/15300D0-TDBJ**5D0/2000000D0

      Sun_lon = mod(Sun_lon,360D0)
      IF (Sun_lon < 0D0) THEN
        Sun_lon = Sun_lon + 360D0
      END IF

      R_sun1 = R_sun  !RA6 * DR2D

      CORR = 0.0057183D0                                         ! correction for long.aberration (degrees)

      EQT0 = Sun_lon - CORR - R_sun1 + (DPSI * COS(OBL))*DR2D     ! value of EQT in degrees
   
      EQT0 = mod(EQT0,360D0)
      EQT = EQT0 * 4D0                                           ! value of EQT in min. of time
      
      EQT = MOD(EQT,60D0)                                        ! value of EQT in the range 0-60 min. 
      EQTX = ABS(EQT)

      EQTM = INT(EQTX)                                           ! Part Equation of time in minutes.
      EQTS = (EQTX - EQTM) * 60D0                                ! Part Equation of time in seconds.
      
      IF (EQT < 0D0) THEN                                        ! choose the sign of EQT
        str="-" 
        EQTM= ABS(EQTM) 
       ELSE IF (EQT > 0D0) THEN 
        str="+"                 
      END IF
    
      IF(mode == 1)THEN
       EQ_label = "(Use Equinox Based)"
      ELSE IF(mode == 2)THEN
       EQ_label = "(Use CIO Based)"
      END IF
   
!---------  Choose unit of  distance ----------------------------------------------------------------

      IF ( NTARG == 10 ) THEN
        distance = " Km "
      ELSE
        distance = " AU "
      END IF  
           
!------------ Compute body's elongation to Sun

!-------- Compute Astrometric Sun Right ascension J2000
  
      Rsun = sqrt( SGE(1,1)**2d0 + SGE(2,1)**2d0 + SGE(3,1)**2d0)   ! Distance Earth - Sun 
     
      xs = SGE(1,1);  ys = SGE(2,1)

      RAsun = datan2( ys,xs)                  !Astrometric Sun Right ascension J2000(rad)
      RAsun =iau_ANP(RAsun)
      RAsun =RAsun * DR2D                     !Astrometric Sun Right ascension J2000(degrees)      

      ELONG = acos((Rsun**2d0 + DIST1**2d0 - RVET**2d0) / (2d0 * Rsun * DIST1))
      ELONG = ELONG * DR2D  
    
      RA = RA3 * DR2D          !Trasform "rad " to degrees astrometric position body                          
      
      EVET = DSQRT( EPV(1,1)**2D0 + EPV(2,1)**2D0 + EPV(3,1)**2D0)   ! Dist. Sun-Earth at XET time
      EPOS = ATAN2(EPV(2,1),EPV(1,1)) 
      EDEC = ASIN (EPV(3,1) / EVET)
      
      SUNLONG = ATAN2((sin(EPOS)*cos(OBL)+tan(EDEC)*sin(OBL)),(cos(EPOS))) + PI  !Geocentric Longit. Sun  
            
      GLON = ATAN2((sin(RA6)*cos(OBL)+tan(DE6)*sin(OBL)),(cos(RA6)))
      GLAT = ASIN(sin(DE6)*cos(OBL) - cos(DE6)*sin(OBL)*sin(RA6))
      COSELON = cos(GLAT) * cos(GLON - SUNLONG)
      
      GLON = GLON * DR2D
      SUNLONG = SUNLONG * DR2D

!------- Search relative position body to Sun ( E = East  or W = West) 

      if(RA >= RAsun) then
       if((RA-RAsun)>=180)then
        EL="West"
        else if((RA-RAsun)<180)then
        EL="East"
       end if
      else if(RA<RAsun) then
        if((RAsun-RA)>=180)then
        EL="East"
        else if((RAsun-RA)<180)then
        EL="West"
       end if
      end if
      
      IF (NTARG == 11) ELONG = 0.0D0    

!------Search  RISE, SET, TRANSIT
       
       RIS = RIS / 24D0
       SET = SET / 24D0
       TRS = TRS / 24D0 
       
       CALL iau_D2TF ( 1, RIS, SIGN, IHMSF1)     
       CALL iau_D2TF ( 1, SET, SIGN, IHMSF2)
       CALL iau_D2TF ( 2, TRS, SIGN, IHMSF3)

!------Search  civil,nautical,astronomical TWILIGHT   
    
      IF(NTARG == 11)THEN 
       CALL iau_D2TF ( 1, BCTW, SIGN, IHMSF10)     
       CALL iau_D2TF ( 1, ECTW, SIGN, IHMSF11)
       CALL iau_D2TF ( 1, BNTW, SIGN, IHMSF12)
       CALL iau_D2TF ( 1, ENTW, SIGN, IHMSF13)     
       CALL iau_D2TF ( 1, BASTW, SIGN, IHMSF14)
       CALL iau_D2TF ( 1, EASTW, SIGN, IHMSF15)
      
       IF(N1 == 1 .or. N2 == 1 .or. N3 == 1) THEN               ! Twilight parameter
          commentP = '---- 00:00:00.0 =   CONTIUOUS TWILIGHT --- '
        ELSE IF(N1 == 0 .and. N2 == 0 .and. N3 == 0) THEN
          commentP = ""
       END IF
      ELSE IF (NTARG /= 11) THEN
          commentP = ""
      END IF 
           
        IF (LX == 1 ) THEN                                      ! Rising,Setting Transit parameter
         commentM = '    ---- Body circumpolar -----   '
        ELSE IF (LX == 0) THEN
         commentM = ""
        END IF
        IF(LM == 1) THEN
         commentN = '--- CONTIUOUS TWILIGHT --- '
         ELSE IF (LM == 2) THEN
         commentN = ' .... NO SET - SUN always up the horizon '
         commentP = ""
         ELSE IF (LM == 0) THEN
         commentN = " "
        END IF
        IF(LN == 1) THEN
         commentA = ' .... DARK .... '
         ELSE IF (LN ==2) THEN
          commentA =  '... NO RISE - body always below the horizon '
         ELSE IF(LN == 0) THEN
          commentA =  " "
        END IF 

        IF ( NR == 1 )THEN                                       
           AS1 = "* "
           commentR = '** Rise - event occurs the next day'
         ELSEIF ( NR == 2) THEN
           AS1 = "**"
           commentR = '** Rise - event occurs the previus day'
         ELSEIF ( NR == 0) THEN
           AS1 = "  " 
           commentR = ""
        END IF
        IF ( MS == 1 )THEN
           commentS = '** Set - event occurs the next day'
           AS2 = "* "
         ELSEIF ( MS == 2 )THEN
           AS2 = "**"
           commentS = '** Set- event occurs the previus day'
         ELSEIF ( MS == 0) THEN
           AS2 = "  " 
           commentS = "" 
        END IF
        IF ( PT == 1 )THEN
           AS3 = "* "
           commentT = '** Transit- event occurs the next day'
         ELSEIF ( PT == 2) THEN
           AS3 = "**" 
           commentT = '** Transit- event occurs the previus day'
         ELSEIF ( PT == 0) THEN
           AS3 = " " 
           commentT = ""
        END IF

!-----------------------------------------------------------------------------------       
       IF (NTARG == 6 ) THEN
        commentU = "In calculating the Saturn's magnitude is also taking into account of the rings." 
        ELSE
        commentU = " "
       END IF 

       IDEM = "******"       
      
     IF (NTARG == 10) THEN
       DISTKm = DIST3 * AUM
      
       CALL moon ( ET,IY,MO,ID,DPSI,RA6,DE6,DISTKm,R_sun,D_sun,RA7,DE7,OBL,GD1,PL,HLON,HLAT, &
                OMEGA,OME0,AGE,JH0,CHI,JAP,APE,DLR,TLong,TLat,LTL,LTB,PA,THL,THB,JH1,JHP,JHU)
     
      MJ = 2400000.5d0
      DJ2 = JH0 - MJ                     ! 2a part of date in  MJD for compute New Moon
      DJ3 = JAP - MJ                     ! 2a part of date in  MJD for compute Apogee - Perigee
      DJ4 = JH1 - MJ                     ! 2a part of date in  MJD for compute First Quarter Moon
      DJ5 = JHP - MJ                     ! 2a part of date in  MJD for compute Full Moon 
      DJ6 = JHU - MJ                     ! 2a part of date in  MJD for compute Last Quarter Moon
      PH = ASIN(6378.137d0 / DLR)        ! APOGEE- PERIGEE Horizontal lunar parallax. (radians)
     END IF

!----------- Geocentric longitude of the Sun for compute the aberration. 

      a0 = (1.396971d0 * T_0 + 0.0003086d0 * T_0**2d0) * DD2R
      b0 = (0.013056d0 * T_0 - 0.0000092d0 * T_0**2d0) * DD2R
      c0 = (5.12362d0 + 0.241614d0*T_0 + 0.0001122d0*T_0**2d0) *DD2R
      
      CALL GEPV(PVH,OBL,1,SLO,SLA) 

      Bsun = (SLA + b0 * sin(SLO + c0)) *(-1d0)
      Lsun = SLO + a0 - b0 * cos(SLO + c0) * tan(SLA) + PI
   
!---------- Convert Geocentric Equatorial (GMB) to Geocentric Ecliptic Coordinate 

      CALL GEPV(GMB,OBL,1,LOGE,BGE1)    ! LOGE = Geocentric - Ecliptic longitude. (radians)
      CALL GEPV(u4,OBL6,1,GELO,GELA)    ! BGE1 = Geocentric - Ecliptic latitude. (radians)

      A_ = cos(eta) * cos(BGE1) * sin(PI_ - LOGE) - sin(eta)*sin(BGE1)
      B_ = cos(BGE1) * cos(PI_ - LOGE)
      C_ = cos(eta) * sin(BGE1) + sin(eta) * cos(BGE1) * sin(PI_ - LOGE)

      atlon = atan2(A_,B_)              ! tan A_/B_ = (p_ + PI_ - HLON)   
      GELO_ = p_ + PI_ - atlon          ! HLON = ELIOCENTRIC ECLIPTIC LONGITUDE (mean of the date)
      GELA_ = asin(C_)                  ! HLAT = ELIOCENTRIC ECLIPTIC LATITUDE  (mean of the date)
    
      IF ( GELO < 0D0) THEN
       GELO = GELO + 2D0 * PI
      END IF
      
      LOGE = LOGE * DR2D                ! LOGE = Geocentric - Ecliptic longitude J2000. (degrees)
      BGE1 = BGE1 * DR2D                ! BGE1  = Geocentric - Ecliptic latitude.J2000 (degrees)

      DLO = (-KA*COS(Lsun-GELO)+ECC*KA*COS(PHE-GELO)) / COS(GELA)        ! Difference in longitude (arcsec)

      DLA = -KA*SIN(GELA)*(SIN(Lsun-GELO)-ECC*SIN(PHE-GELO))             ! Difference in Latitude (arcsec)

      IF (NTARG == 10)THEN
       GELOAP = GELO   
       GELAAP = GELA 

      ELSE IF(NTARG < 10) THEN
       GELOAP = GELO                     ! GELOAP = Geocentric Apparent Ecliptic longitude mean the date(rad).             
       GELAAP = GELA                     ! GELAAP = Geocentric Apparent Ecliptic latitude mean the date rad).  

      ELSE IF(NTARG == 11) THEN
       GELOAP = GELO 
       GELAAP = GELA                          
      END IF

      GELOD = GELO * DR2D                ! GELOD = Geocentric - Ecliptic longitude mean the date. (degrees) 
      GELAD = GELA * DR2D                ! GELAD = Geocentric - Ecliptic latitude mean the date. (degrees)
      GELOAPD = GELOAP * DR2D
      GELAAPD = GELAAP * DR2D 

      Tcy= (ET-DJ00)/36525d0
      LEarth = iau_FAE03 ( Tcy )
      PEarth = iau_FAPA03(Tcy)
      LSole = (LEarth + PEarth + PI)* DR2D

      WRITE(*,*)"=========================================================================================="

      WRITE(*,*)" --------------------------------DATA OUTPUT ----------------------------------"
      WRITE(*,*)
      WRITE(*,999) Ytoday,Mtoday,Dtoday,MJDtoday,FINAL_,EOP_final
999   FORMAT(2x,"TODAY Y/M/D ",i4,"/",i2,"/",i2,4x,"MJDtoday ",i5,4x,"IERS FINAL_(MJD) ",i5,4x,"JUNP(MJD) EOP14 C04 to IERS A ",i5)     
      WRITE(*,*)
      WRITE(*,*) "Observer Geographic coordinate  "
      WRITE(*,5) longit,lS,latit,bS,height
5     FORMAT (1x,"Longitude",1x,"deg",F10.5,1x,a,3x,"Latitude",1x,"deg",F9.5,1x,a,3x,"Height",1x,F6.1,1x,"meter") 
      WRITE(*,*)
      WRITE(*,*) "Observer Geocentric coordinate X,Y,Z "     
      WRITE(*,*) " X = ",XYZ(1)
      WRITE(*,*) " Y = ",XYZ(2)
      WRITE(*,*) " Z = ",XYZ(3)
      WRITE(*,*) " "
      WRITE(*,6) Dday
6     FORMAT (1x,"Day of week",1x,a10)
      CALL iau_D2DTF ( "UTC", 6, UTC1, UTC2, IY, MO, ID, IHMSF, J )
      WRITE ( *, 1 ) "UTC", IY, MO, ID, IHMSF, UTC
      CALL iau_D2DTF ( "UT1", 6, UT11, UT12, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "UT1", IY, MO, ID, IHMSF, UT1
      CALL iau_D2DTF ( "TAI", 6, TAI1, TAI2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "TAI", IY, MO, ID, IHMSF, TAI
      CALL iau_D2DTF ( "TT", 6, TT1, TT2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "TT ", IY, MO, ID, IHMSF, TT
      CALL iau_D2DTF ( "TCG", 6, TCG1, TCG2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "TCG", IY, MO, ID, IHMSF, TCG
      CALL iau_D2DTF ( "TDB", 6, TDB1, TDB2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "TDB", IY, MO, ID, IHMSF, TDB
      CALL iau_D2DTF ( "TCB", 6, TCB1, TCB2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 1 ) "TCB", IY, MO, ID, IHMSF, TCB
1     FORMAT ( 1X,A,I5,2("/",I2.2),1X,I3,2(":",I2.2),".",I6.6,"   MJD ",F18.12 )

      WRITE(*,13)" UT1-UTC",DUT,data_b,LDUT
13    FORMAT(a8,"  ",F10.7," sec.",4x,"DATABASE ",x,a6,a30)
      WRITE(*,1313),ut1_notid
1313  FORMAT(x,"ut1_int , without tidal corrections ",f11.8," sec.")
      WRITE(*,14)" DELTA_T",DET,LDET
14    FORMAT(a8,"  ",F9.6," sec.",a30)
      WRITE(*,15)" LOD ",LOD,Slod
      
15    FORMAT(a4,6x,F7.4,2x," msec.",1x,a21)
      WRITE(*,150)" TT - TDB at geocenter: ",PVV(1)," sec.  / ",PVV(2),"sec/day"
150   FORMAT(a23,x,f15.13,a5,2x,f15.13,a10)
      WRITE(*,*)
 
      WRITE(*,16) str,EQTM,EQTS,EQ_label
16    FORMAT('     Equation of time -----------    =',2x,a,i2," min ",f5.2,' sec.',a20)  
      WRITE(*,17) EQEQ 
17    FORMAT('     Equation of equinox --------    =',2x,f13.9,4x,' degs')  
      CALL iau_A2TF ( 6, GAST, SIGN, IHMSF )      
      WRITE(*,18) GAST*DR2D,IHMSF
18    FORMAT('     Greenwich App. Sid. Time (GAST) =',f19.13,' degs =',1X,I3,2(":",I2.2),".",I6.6," h,m,s")
      CALL iau_A2TF ( 6, LAST, SIGN, IHMSF )      
      WRITE(*,180) LAST*DR2D,IHMSF
180   FORMAT('     Local App. Sid. Time (LAST)     =',f19.13,' degs =',1X,I3,2(":",I2.2),".",I6.6," h,m,s")
      CALL iau_A2TF ( 6, ERA, SIGN, IHMSF )
      WRITE(*,19) ERA*DR2D,IHMSF
19    FORMAT('     Earth Rotation Angle ( ERA )    =',f19.13,' degs =',1X,I3,2(":",I2.2),".",I6.6," h,m,s")
      
      WRITE(*,*)
      WRITE(*,*)" ============ Use IERS STANDARD EOP DATA (iau 2000) ==============="
      WRITE(*,*)
      WRITE (*,20) XP * DR2AS,LDUT
20    FORMAT('     Polar motion XP  ',f9.6,' arcsec.',a30)   
      WRITE (*,21) YP * DR2AS,LDUT
21    FORMAT('     Polar motion YP  ',f9.6,' arcsec.',a30)   
      WRITE(*,*)
      WRITE (*,22) DX0,LDUT
22    FORMAT('     Celestial pole offset dX  ',f8.4,' mas',a30)   
      WRITE (*,23) DY0,LDUT
23    FORMAT('     Celestial pole offset dY  ',f8.4,' mas',a30)   
      WRITE(*,*)
      WRITE (*,24)DPSI * DR2AS
24    FORMAT("     Nutation in longitude ",3x,F10.6,1x,"arcsec.") 
      WRITE (*,25)DEPS * DR2AS
25    FORMAT("     Nutation in obliquity ",3x,F10.6,1x,"arcsec.")  
      WRITE (*,*)
      WRITE (*,*)"=================================================================="
      WRITE (*,*)
      WRITE (*,36)" Body barycenter : ",planet
36    FORMAT (a17,a24 )      
      WRITE (*,48)"=================================================================================================="
      WRITE(*,*)
      WRITE(*,*) "      Body's geocentric  equatorial position and velocity     "
      WRITE(*,*) "         mean equator and equinox of J2000.           "
      WRITE (*,37)"             X                   Y                     Z "
      WRITE (*, 38)" POS ", BGV(1,1), BGV(2,1), BGV(3,1)," AU"
      WRITE (*, 38)" VEL ", BGV(1,2), BGV(2,2), BGV(3,2)," AU/day"
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,48)"=================================================================================================="
      WRITE(*,*)
      WRITE(*,*) "        SUN geocentric equatorial position and velocity  "
      WRITE(*,*) "             mean equator and equinox of J2000.       "
      WRITE (*,37)"             X                   Y                     Z "
      WRITE (*, 38)" POS ", SGE(1,1), SGE(2,1), SGE(3,1)," AU"
      WRITE (*, 38)" VEL ", SGE(1,2), SGE(2,2), SGE(3,2)," AU/day"
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,48)"=================================================================================================="
      WRITE(*,*)
      WRITE(*,*) "        Barycentric equatorial position and velocity  "
      WRITE(*,*) "             mean equator and equinox of J2000.       "
      WRITE (*,37)"             X                   Y                     Z "
      WRITE (*, 38)" POS ", PVB(1,1), PVB(2,1), PVB(3,1)," AU"
      WRITE (*, 38)" VEL ", PVB(1,2), PVB(2,2), PVB(3,2)," AU/day"
      WRITE (*,*)
      WRITE (*,48)"=================================================================================================="
      WRITE(*,*)
      WRITE(*,49)" Heliocentric rectangular position and velocity mean equator and equinox J2000.(λ,B)"
      WRITE (*,37)"             X                   Y                     Z "
      WRITE (*, 38)" POS ", BEV(1,1), BEV(2,1), BEV(3,1)," AU"
      WRITE (*, 38)" VEL ", BEV(1,2), BEV(2,2), BEV(3,2)," AU/day"
      WRITE (*,*)
      WRITE(*,*) " Heliocentric ecliptic coordinate mean equator and equinox J2000 ( λ, B, RV )"
      WRITE (*,31) "     o , ' , ''               o , ' , ''               RV AU          Rdot (Km/s)              "
      CALL iau_A2AF ( 5 ,LON0 ,sign , IDMSF0 )
      CALL iau_A2AF ( 4 ,LAT0 ,sign , IDMSF )
      WRITE ( *, 40 ) IDMSF0 ,sign,IDMSF,RVET0,Rdot
      WRITE (*,*)
      WRITE(*,*) " Heliocentric ecliptic coordinate mean equator and equinox of the date ( λ, B, RV )"
      WRITE (*,31) "     o , ' , ''               o , ' , ''               RV AU          Rdot (Km/s)              "
      CALL iau_A2AF ( 5 ,HLON ,sign , IDMSF0 )
      CALL iau_A2AF ( 4 ,HLAT ,sign , IDMSF )
      WRITE ( *, 40 ) IDMSF0 ,sign,IDMSF,RVET0,Rdot

       WRITE (*,48)"=================================================================================================="
      WRITE (*,*) "Geocentric Astrometric  J2000 coordinate with relativistic corrections"
      WRITE (*,32) "              R.A                    Dec.            Distance       Light time       deldot (Km/s)"
      IF ( NTARG == 10 ) THEN
       WRITE (*,64) 
      ELSE IF( NTARG /= 10) THEN
       WRITE (*,65)
      END IF  
      CALL iau_A2TF ( 5 ,RA3 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE3 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,54 ) IHMSF ,sign,IDMSF,DIST1 * AUM,Mtime,deldot
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *, 53 ) IHMSF ,sign,IDMSF,DIST1,Ltime,deldot
      END IF
            
       WRITE (*,48)"=================================================================================================="
      WRITE (*,*) "Topocentric Astrometric  J2000 coordinate with relativistic corrections"
      WRITE (*,32) "              R.A                    Dec.            Distance       Light time       deldot (Km/s)"
      IF ( NTARG == 10 ) THEN
       WRITE (*,64) 
      ELSE IF( NTARG /= 10) THEN
       WRITE (*,65)
      END IF  
      CALL iau_A2TF ( 5 ,RA0 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE0 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,54 ) IHMSF ,sign,IDMSF,DIST4 * AUM,Mtime1,deldot1
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *, 53 ) IHMSF ,sign,IDMSF,DIST,Ltime1,deldot1
      END IF
      
      WRITE (*,48)"=================================================================================================="

      WRITE (*,*) "Geocenter mean J2000  equatorial and  ................................ ecliptic coordinate        "
      WRITE (*,32) "             R.A                     Dec.            Distance         Longit.        Latit.      "
      IF ( NTARG == 10 ) THEN
      WRITE (*,62) 
      ELSE IF( NTARG /= 10) THEN
      WRITE (*,63)
      END IF        
      CALL iau_A2TF ( 5 ,RA4 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE4 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,52 ) IHMSF ,sign,IDMSF,DIST2 * AUM,LOGE,BGE1
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *, 30 ) IHMSF ,sign,IDMSF,DIST2,LOGE,BGE1
      END IF

      WRITE (*,48)"=================================================================================================="
      WRITE (*,*)equinox
      WRITE (*,32) "            R.A                      Dec.          Distance          Distance                    "
      IF ( NTARG == 10 ) THEN
       WRITE (*,64) 
      ELSE IF( NTARG /= 10) THEN
       WRITE (*,65)
      END IF  
       IF (mode == 1) THEN
      CALL iau_A2TF ( 5 ,RA6 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE6 ,sign , IDMSF )
       ELSE IF (mode == 2) THEN
      CALL iau_A2TF ( 5 ,RA6 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE6 ,sign , IDMSF )
       END IF

      IF ( NTARG == 10 ) THEN
        WRITE ( *,58 ) IHMSF ,sign,IDMSF,DIST3 * AUM,Distmoon
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *,59 ) IHMSF ,sign,IDMSF,DIST3
      END IF

      WRITE (*,48)"==================================================================================================" 
      WRITE (*,*)equinox2 
      WRITE (*,61)                                      
      WRITE (*,60)                         
      CALL iau_A2TF ( 5 ,RA8 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE8 ,sign , IDMSF )
      WRITE ( *, 35 )IHMSF ,sign,IDMSF,AZ,ZD0,HOT
      WRITE (*,48)"==================================================================================================" 
      WRITE(*,*)equinox3
      WRITE (*,32) "           R.A                       Dec.             Dist.      (*)Elevation   (*)Hour Angle  "
      IF ( NTARG == 10 ) THEN
      WRITE (*,62) 
      ELSE IF( NTARG /= 10) THEN
      WRITE (*,63)
      END IF        
      CALL iau_A2TF ( 5 ,RA9 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE9 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,58 ) IHMSF ,sign,IDMSF,DIST4 * AUM,ZOB,HOB
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *,59 ) IHMSF ,sign,IDMSF,DIST,ZOB,HOB
      END IF

      IF( NTARG /= 10) THEN
      WRITE (*,48)"=====================================================================================================" 
      WRITE(*,*)
       CALL iau_A2AF ( 4 ,GELOAP ,sign , IDMSF )
        WRITE (*,310)sign,IDMSF,GELOAPD !, GELOD 
        WRITE(*,*)
       CALL iau_A2AF ( 4 ,GELAAP ,sign , IDMSF )
        WRITE (*,320)sign,IDMSF,GELAAPD !, GELAD      !!<<<<<<<<<<<<<<<< controllare
      END IF 
      WRITE (*,48)"=====================================================================================================" 
      WRITE(*,*)
      WRITE (*,32) "       ILL. %       MAGN.        Angular Diameter        Phase Angle (o)     Sun's Elong.         "
      IF ( NTARG == 10 .OR. NTARG == 11 ) THEN
        WRITE ( *,56 ) ILL,MAGN,DIAM,PANGLE,ELONG,EL
       ELSE IF( NTARG < 10 ) THEN
        WRITE ( *,55 ) ILL,MAGN,DIAM,PANGLE,ELONG,EL 
       IF ( NTARG == 5 .OR. NTARG == 6 .OR. NTARG == 7 .OR. NTARG == 8) THEN
        WRITE (*,50) PDIAM
       ELSE IF ( NTARG /= 5 .OR. NTARG /= 6 .OR. NTARG /= 7 .OR. NTARG /= 8) THEN
        WRITE ( *,51 )IDEM  
       END IF
      END IF
      WRITE(*,48)"======================================================================================================"
      WRITE (*,*)
      IF (NTARG == 10)THEN
      WRITE (*,*)"                  ------------------ OTHER MOON'S  PARAMETER ------------------------                         "
     CALL iau_D2DTF ( "TT", 6, MJ, DJ2, IY, IM, ID, IHMSF4, J ) !Transform Julian day of date,tine New Moon
     CALL iau_D2DTF ( "TT", 6, MJ, DJ4, IY1, IM1, ID1, IHMSF6, J ) !Transform Julian day of date,tine First Quarter Moon
     CALL iau_D2DTF ( "TT", 6, MJ, DJ5, IY2, IM2, ID2, IHMSF7, J ) !Transform Julian day of date,tine Full Moon
     CALL iau_D2DTF ( "TT", 6, MJ, DJ6, IY3, IM3, ID3, IHMSF8, J ) !Transform Julian day of date,tine Last Quarter Moon
       CALL iau_A2AF ( 4 ,GELO ,sign , IDMSF )
      WRITE (*,301)sign,IDMSF,GELOD
       CALL iau_A2AF ( 4 ,GELOAP ,sign , IDMSF )
      WRITE (*,311)sign,IDMSF,GELOAPD
       CALL iau_A2AF ( 4 ,GELA ,sign , IDMSF )
      WRITE (*,321)sign,IDMSF,GELAD 
      WRITE(*,*)
      WRITE(*,*)"-------------------------------------------------------------------------------------"     
      WRITE (*,*)
      WRITE (*,299),LEA(1),LEA(2),LEA(3)
299   FORMAT( " Lunar Euler angles: phi, theta, psi (rad)........... ",x,f12.10,2x,f12.10,2x,f14.8)
      WRITE (*,300) LEA(4),LEA(5),LEA(6)
300   FORMAT(" Lunar librations: dphi/dt, theta/dt, dpsi/dt (rad/day)",x,f12.10,2x,f12.10,2x,f12.10)      
      WRITE(*,*)
      WRITE (*,302) "           NEW MOON (TT)", IY, IM, ID, IHMSF4,AGE               !New Moon in Gregorian date
      WRITE (*,*)
      WRITE (*,3020)" FIRST QUARTER MOON (TT)", IY1, IM1, ID1, IHMSF6                !Firs Quarter Moon in Gregorian date
      WRITE (*,*)
      WRITE (*,3021)"          FULL MOON (TT)", IY2, IM2, ID2, IHMSF7                !Firs Quarter Moon in Gregorian date
      WRITE (*,*)
      WRITE (*,3022)"  LAST QUARTER MOON (TT)", IY3, IM3, ID3, IHMSF8                !Last Quarter Moon in Gregorian date
      WRITE (*,*)
      WRITE (*,303)OMEGA,OME0
      WRITE (*,*)
      CALL iau_D2DTF ( "UTC", 6, MJ, DJ3, IY4, IM4, ID4, IHMSF5, J ) !Transform Julian day of 
      CALL iau_A2AF ( 4, PH, SIGN, IDMSF1 )
      WRITE (*,304)APE,IY4,IM4,ID4,IHMSF5,DLR,IDMSF1                 ! Apogee-Perigee
      WRITE (*,*) 
      WRITE (*,305) CHI
      WRITE (*,*)
      WRITE (*,306)LTL,LTB
      WRITE (*,*)
      WRITE (*,307)TLong,TLat
      WRITE (*,*)
      WRITE (*,308)PA
      WRITE (*,*)
      WRITE (*,309)THL,THB
      WRITE (*,*)  "-----------------------------------------------------------------------------------------------------"      
      END IF
      WRITE (*,*)
      WRITE (*,41) "          RISE (UT)                TRANSIT (UT)               SET (UT)                               "
      WRITE (*,*)
      WRITE (*,42) IHMSF1,AS1, IHMSF3,AS3,IHMSF2,AS2
      WRITE (*,*)
      WRITE (*,41) "        Azimuth                    Altitude                 Azimuth                              " 
      WRITE (*,43) AZM1,HAL,As,AZM2 
      WRITE (*,*)
             IF (NTARG == 11) THEN
      WRITE (*,41)"               MORNING                                          EVENING                           " 
      WRITE (*,70)IHMSF10,IHMSF11
      WRITE (*,71)IHMSF12,IHMSF13
      WRITE (*,72)IHMSF14,IHMSF15
             END IF
      WRITE(*,*)
      WRITE (*,47) 
      WRITE (*,44)planet                             
      WRITE (*,239) commentP
      WRITE (*,39) commentR
      WRITE (*,39) commentS
      WRITE (*,39) commentT 
      WRITE (*,45) commentM
      WRITE (*,45) commentN
      WRITE (*,45) commentA
      WRITE (*,48)"======================================================================================================"
      WRITE (*,*)" (*) SOFA NOTES for option (IAU 2006 CIO BASED)                         "
      WRITE (*,*)"                                                                        "
      WRITE (*,*)"       The accuracy of the result is limited by the corrections for     "      
      WRITE (*,*)"       refraction, which use a simple A*tan(z) + B*tan^3(z) model.      "
      WRITE (*,*)"       Providing the meteorological parameters are known accurately and " 
      WRITE (*,*)"       there are no gross local effects, the predicted observed         "
      WRITE (*,*)"       coordinates should be within 0.05 arcsec (optical) or 1 arcsec   "                                  
      WRITE (*,*)"       (radio) for a zenith distance of less than 70 degrees, better    "                           
      WRITE (*,*)"       than 30 arcsec (optical or radio) at 85 degrees and better than  "
      WRITE (*,*)"       20 arcmin (optical) or 30 arcmin (radio) at the horizon.         "
      WRITE (*,*)"       Relative to the Sun and the Moon rise-setting they are taken into account "
      WRITE (*,*)"       respectively (- 0 °. 50 ') and (-0 ° .34 '; lunar semid. and parallax)"
      WRITE (*,*)
      WRITE (*,*)"       The option 'IUA 2000 Equinox based'          "
      WRITE (*,*)"       compute refraction only for body elevation   "
      WRITE (*,*)"       between -1.0 and 89.9 degrees.               "
      WRITE (*,*)
      WRITE (*,*)"  (*)  The calculated value of the magnitude of Mercury and Venus is    "
      WRITE (*,*)"       not real when the phase angle of the two celestial bodies is too "
      WRITE (*,*)"       close to the Sun to be observed photometrically, so the reported "
      WRITE (*,*)"       value is invalid.                                                "
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,48)"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
      
30    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,5x,F11.6,5x,F11.6)       
31    FORMAT (A98)
32    FORMAT (A98)
35    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F10.6,5x,F11.6,4x,F11.6)      
37    FORMAT (A60)
38    FORMAT (A5,' ',F22.16,' ', F22.16,' ', F22.16,' ',A7)    
40    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,5x,F10.6) 
41    FORMAT (a98)
42    FORMAT ( 6X," h ",I3,2(":",I2.2),".",I1.1,a2,11X," h ",I3,2(":",I2.2),".",I2.2,a2,10X," h ",I3,2(":",I2.2),".",I1.1,a2) 
43    FORMAT ( 9X,F5.1,"°",22X,F4.1,"°",a,19X,F5.1,"°")
44    FORMAT (10x," The",1x,A10,1x,A90)
39    FORMAT (25x,A45) 
239   FORMAT (25x,A45)
45    FORMAT (25x,A45)
46    FORMAT (25x,A45) 
47    FORMAT (5x," *** NOTES *** ")
48    FORMAT (A100)
49    FORMAT (A84) 
50    FORMAT (30x,"Polar.",F10.3," arcsec.") 
51    FORMAT (30x,"Polar.",3x,a6,1x," arcsec.")!(15x,a14,x,"Polar.",xxx," *****",x," arcsec.")      
52    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F11.3,5x,F11.6,4x,F11.6)
53    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,4x,F10.6," min.",5x,F10.6)
54    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F11.3,4x,F9.6," sec.",5x,F10.6)        
55    FORMAT (7x,F6.2,5x,F6.2,6x,"Equat.",F10.3," arcsec.",6x,F6.2,12x,F5.1,1x,5a)  
56    FORMAT (7x,F6.2,5x,F6.2,9x,F8.3," arcmin.",10x,F6.2,12x,f5.1,1x,5a) 
58    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,4x,F11.3,4x,F11.6,5x,F11.6) 
59    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,4x,F13.10,4x,F11.6,5x,F11.6) 
60    FORMAT (9x,"h , m ,  s",14x,"o , ' , ''",13x,"degs",12x,"degs",11x,"degs")
61    FORMAT (14x,"R.A",21x,"Dec.",13x,"Azimut",5x,"(*)Elevation .",2x,"(*)Hour Angle")
62    FORMAT (9x,"h , m ,  s",14x,"o , ' , ''",13x," Km ",12x,"degs",11x,"degs")
63    FORMAT (9x,"h , m ,  s",14x,"o , ' , ''",13x," AU ",12x,"degs",11x,"degs")
64    FORMAT (9x,"h , m ,  s",14x,"o , ' , ''",13x," Km ",9x,"Earth radius")
65    FORMAT (9x,"h , m ,  s",14x,"o , ' , ''",13x," AU ")
70    FORMAT (5x,"Beginnig civil twilight        "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending civil twilight        ", &
              " h ",I3,2(":",I2.2),".",I1.1)
71    FORMAT (5x,"Beginnig nautical twilight     "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending nautical twilight     ", &
              " h ",I3,2(":",I2.2),".",I1.1)
72    FORMAT (5x,"Beginnig astronomical twilight "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending astronomical twilight ", &
              " h ",I3,2(":",I2.2),".",I1.1)
301   FORMAT (2x," Moon longitude mean the date  ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
311   FORMAT (2x," Moon apparent longitude       ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
321   FORMAT (2x," Moon's latitude mean the date ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
302   FORMAT (3X,A,I5,2("/",I2.2),3X,"h",1x,I3,2(":",I2.2),".",I6.6,5x,"Days after the new moon",4x,f5.2)
3020  FORMAT (3X,A,I5,2("/",I2.2),3X,"h",1x,I3,2(":",I2.2),".",I6.6)
3021  FORMAT (3X,A,I5,2("/",I2.2),3X,"h",1x,I3,2(":",I2.2),".",I6.6)
3022  FORMAT (3X,A,I5,2("/",I2.2),3X,"h",1x,I3,2(":",I2.2),".",I6.6)
303   FORMAT (3x,"Mean Longit. ascend. node",1x,f10.6," degs",3x," True Longit. ascend. node ",f11.6," degs")
304   FORMAT (3x,a7,I5,2("/",I2.2),4X,"h",1x,I3,2(":",I2.2),".",I6.6,4x," Distance ",f10.3," Km = Parallax " & 
              ,1x,I3,'°',1x,I2,"'",1x,I2,"''.",I4.4)  
305   FORMAT (3x,"Position angle of the Moon's bright limb  ",f8.4," degs")
306   FORMAT (3x,"Earth's Selenographic Longitude     ",f8.4," degs","    Earth Selenographic Latitude       ",f7.4," degs")
307   FORMAT (3x,"Topocentric libration in Longitude  ",f8.4," degs","    Topocentric Libration in Latitude  ",f7.4," degs")
308   FORMAT (3x,"Position angle of axis              ",f8.4," degs")
309   FORMAT (3x,"Sun's Selenographic Longitude       ",f8.4," degs","    Sun's Selenographic Latitude       ",f7.4," degs")      
310   FORMAT (5x,"Apparent Ecliptic Longit.mean the date",2x,a,I3,"°",x,I2,"'",x,I2,"''",".",I4.4,4x,"degs",x,f10.6 )
320   FORMAT (5x,"Apparent Ecliptic Latit. mean the date",2x,a,I3,"°",x,I2,"'",x,I2,"''",".",I4.4,4x,"degs",x,f10.6 )

 GOTO 3000 
      
!=====================================================================================================================

     
1640  WRITE(*,*)"*********************************************************************"              
      WRITE(*,*)""
      WRITE(*,*)"         THE NUMBERING CONVENTION FOR 'NTARG' IS:"
      WRITE(*,*)""
      WRITE(*,*)"     1 = MERCURY          8 = NEPTUNE"
      WRITE(*,*)"     2 = VENUS            9 = PLUTO"
      WRITE(*,*)"     3 = EARTH           10 = MOON"
      WRITE(*,*)"     4 = MARS            11 = SUN"
      WRITE(*,*)"     5 = JUPITER         12 = SOLAR-SYSTEM BARYCENTER"
      WRITE(*,*)"     6 = SATURN          13 = EARTH-MOON BARYCENTER"
      WRITE(*,*)"     7 = URANUS          14 = ************** "
      WRITE(*,*)""
      WRITE(*,*)
      WRITE(*,*)"*********************************************************************"

1650      WRITE(*,*) 'Enter Target Body ' 
      READ *, NTARG
                

      DO I = 1,15
       IF ( NTARG == I) THEN
        planet = body(I)
        GOTO 1660
       ELSE IF ( NTARG > 15) THEN
        PRINT *,"Uncorrect target body!"
        GOTO 1650 
       END IF
      END DO  

1660 continue        

!---------------------------------------------------------------------

      IF(IY == 1900 .OR. IY == 2100) THEN
        numday = 365
      END IF
  
      BIS = IY/4d0 - int(IY/4d0)
      IF(BIS == 0) THEN
        numday = 366
        step = 18.3d0
      ELSE IF (BIS /= 0) THEN
        numday = 365
        step = 18.25d0
      END IF
   
      GN = (100D0*IY+MO-190002.5)                    !I find the sign of the numerical expression
      IF (GN < 0D0) THEN
       SGN = -1
      ELSE IF (GN > 0D0) THEN
       SGN = +1
      END IF

!      JD0 = 367*IY - (7*(IY+((MO+9)/12)) / 4)+(275*MO/9) + ID + 1721013.5  & 
!             + TIME - 0.5 * SGN +0.5
             
      JD0 = 367*IY - (7*(IY+((1+9)/12)) / 4)+(275*1/9) + 1 + 1721013.5  &     ! JD 1900.0
             + 0D0 - 0.5 * SGN +0.5

      day0 = JD -JD0  

      fract_year = day0/numday
      
      IYdate = IY + fract_year

!-------- Compute TT and UT1 
      
      IF (modetime == 1) THEN
        CALL iau_CAL2JD ( IY, MO, ID, DJM0, DATE_, J )            ! input UT1 (IY,MO,ID)
        UT1 = DATE_ + TIME
        GOTO 1750
      END IF 

      IF (modetime == 2) GOTO 1700 
      
!--------TT (MJD).
1700  CALL iau_CAL2JD ( IY, MO, ID, DJM0, DATE_, J )              ! input TT (IY,MO,ID)
      TDB = DATE_ + TIME 

1750  continue                            
     
!----------Call the interpolation routine for date above JD 2437665.5 (1962/01/01) 

!  input  JD    :  date 
!         
!  output XP0   : polar motion (arcsec)      
!         YP0   : polar motion (arcsec)
!         DUT   : UT1 - TAI since 1962.00 (sec.)
!         dpsi_ : celestial pole offset dpsi / UAI 1980 (arcsec)
!         deps_ : --------------------- deps ------------------
!         DELTA : DeltaT (TDB-UT sec.)
!    
!  Coded by A. Nicola - October 2011
      
      CALL iers_1900 (JD,IYdate,step,numday,XPO,YPO,DUT,dpsi_,deps_,DELTA)

      DELTAT = DELTA / 86400D0       ! Value at day
      TUT = TIME                             
      
      XP = XPO * AS2R                !  Polar motion (arcsec---> radians)fva<vfaaev
      YP = YPO * AS2R 
      DDP80 = dpsi_ * AS2R           !  CIP correction ( mas-----> radians)
      DDE80 = deps_ *AS2R 
   

     IF (modetime == 1) THEN                 
      TDB = UT1 + DELTAT             ! TDB in MJD 
      TDB1 = INT(TDB + DJMJD0)       ! integer part of TDB (JD) 
      TDB2 = (TDB + DJMJD0) - TDB1   ! fractionale part of TDB (JD)
      UT11 = UT1 + DJMJD0            ! UT1 in JD  
      UT = UT11 - INT(UT11)          ! UT = fract. part of UT1 (JD)
     
      DTR = iau_DTDB ( TDB1, TDB2, UT, LONG, U, V )
     
!----------- TDB---> TT
      CALL iau_TDBTT(TDB1,TDB2,DTR,TT1,TT2,J)
      IF ( J.NE.0 ) STOP
    
!----------- TT---> TCG
      CALL iau_TTTCG ( TT1, TT2, TCG1, TCG2, J )
      IF ( J.NE.0 ) STOP
!----------- TDB---> TCB
      CALL iau_TDBTCB ( TDB1, TDB2, TCB1, TCB2, J )
      IF ( J.NE.0 ) STOP 
      GOTO 1900 
     ELSE IF (modetime == 2) THEN
       UT1 = TDB - DELTAT             ! UT1 in MJD  Equalizing TT at TDB
       UT11 = UT1 + DJMJD0           ! UT1 in JD  
       UT = UT11 - INT(UT11)
      TDB1 = INT(TDB + DJMJD0)       ! integer part of TDB (JD) 
      TDB2 = (TDB + DJMJD0) - TDB1   ! fractionale part of TDB (JD)
       DTR = iau_DTDB ( TDB1, TDB2, UT, LONG, U, V )  
!----------- TDB---> TT
      CALL iau_TDBTT(TDB1,TDB2,DTR,TT1,TT2,J)
      IF ( J.NE.0 ) STOP
!----------- TT---> TCG
      CALL iau_TTTCG ( TT1, TT2, TCG1, TCG2, J )
      IF ( J.NE.0 ) STOP
!----------- TDB---> TCB
      CALL iau_TDBTCB ( TDB1, TDB2, TCB1, TCB2, J )
      IF ( J.NE.0 ) STOP 
      GOTO 1900 
     END IF
      
1900  continue
      TT= TT1+TT2-DJMJD0                  
      TDB = TDB1+TDB2         !
      TCG = TCG1+TCG2         ! 
      TCB = TCB1+TCB2         !
     
      UT_0 = INT(MJD) + DJMJD0        ! JD Time to 0.0 h UT1.      
      TT_0 = UT_0 + DELTAT            ! JD Time to 0.0 h TT the day.
!
! =============
! IAU 1976/1980
! =============
! ----IAU 1976 precession matrix, J2000.0 to date.

     CALL iau_PMAT76 ( DJMJD0, TT, RP )
!---- IAU 1980 nutation.
     CALL iau_NUT80 ( DJMJD0, TT, DP80, DE80 )

!---- Add adjustments: frame bias, precession-rates, geophysical.
     DPSI = DP80 + DDP80
     DEPS = DE80 + DDE80

!---- Mean obliquity.
     EPSA = iau_OBL80 ( DJMJD0, TT )

!---- Build the rotation matrix.
     CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )

!---- Combine the matrices: PN = N x P.
     CALL iau_RXR ( RN, RP, RNPB )

!---- Equation of the equinoxes, including nutation correction.
     EE = iau_EQEQ94 ( DJMJD0, TT ) + DDP80 * COS ( EPSA )

!---- Greenwich apparent sidereal time (IAU 1982/1994).
     GST = iau_ANP ( iau_GMST82 ( DJMJD0+DATE_, TUT ) + EE )
 
!=========================================================================================

!-----Nutation, IAU 2000A. to 0.0h TT (for compute body Transit )
      CALL iau_NUT00A ( TT_0, 0D0, DP00, DE00 )
      
!---- Precession-nutation quantities, IAU 2000. to 0.0h
      CALL iau_PN00 ( TT_0, 0D0, DP00, DE00, &
           EPSA, RB, RP, RPB, RN, RNPB0 )

!-----Greenwich apparent sidereal time (IAU 2006).
      GAST0 = iau_ANP ( iau_GMST06 ( UT_0, 0D0, TT_0, 0D0 ) &
          + iau_EE06A ( TT_0, 0D0, EPSA, DPSI ) )
     
      GAST = iau_ANP ( iau_GMST06 ( DJMJD0, UT1, TT1, TT2 )&
          + iau_EE06A ( TT1, TT2, EPSA, DPSI ) )
      
      CALL GEOCPV ( LONG,LAT,HM,GAST,POS,VEL,J)     ! Geocenter position and velocity observer 
      IF ( J.NE.0 ) STOP
   
!-----Earth  rotation  angle. 
       ERA  =  iau_ERA00  ( DJMJD0+DATE_,TUT)  
       
!==========================================================================================
     
      ET = TDB 
      T0 =2451545.0D0
      print*," ET,TDB",ET,TDB
      CALL NUTATION (-ET,DJMJD0,TT,POS,POS2)
      CALL NUTATION (-ET,DJMJD0,TT,VEL,VEL2)

!     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   GP )

!     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   GV )

!--------------------------------------------------------------------------------------------
      ISTEP=1

      NSTEP=1

!       ET0(1) = TDB1
!       ET0(2) = TDB2 
            
!      CALL  Const (NAMS, VALS, SS, NVS)

      CALL PLEPH ( ET, 17, 0, PVV)              ! PVV(1) = TT-TDB in sec./ PVV(2)=rate of change of TT-TDB in sec/day
      
      CALL PLEPH ( ET, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
      CALL PLEPH ( ET, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( ET, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector
 
      CALL PLEPH ( ET, NTARG, 11, BEV)          ! BEV Body's heliocentric state vector 

      call PLEPH ( ET, 11, 3, SGE )             ! SGE = Sun geocentric state vector J2000.0

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

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP Body's heliocentric position

      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus

      CALL iau_PMP(PVB,PVE,GMB)                  ! GMB  Body's geocentric mean state vector

      CALL iau_PVUP (DT,GMB,GPB)                 ! GPB  Body's geocentric mean J2000 position    

      CALL  iau_C2S   ( GPB, RA4, DE4 )       
      RA4 = iau_ANP(RA4)                         ! range angle (0,2!pi)

      DIST2 = DSQRT(GPB(1)*GPB(1)+GPB(2)*GPB(2)+GPB(3)*GPB(3))     ! Geocentric Distance 

      EARTHVEC = DSQRT(HPE(1)*HPE(1)+HPE(2)*HPE(2)+HPE(3)*HPE(3))  ! Distance Sun - Earth
      
!---------- Astrometric coordinate      
      
      CALL L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,HE0,THD,EPV)
      
      BEV1 = PV - PVE                              ! Geocentric state vector

      BODY_G = BODY_G
      DIST1 = GD
      RVETT = HD
!-----------Form  astrometric position RA3,DE3
 
      CALL iau_C2S ( BODY_G, RA3, DE3 )   ! Call for direction cosines to spherical coordinates     
    
      RA3 = iau_ANP(RA3)                  ! range angle (0,2!pi)
   
!-----------Form RVET (vector ray mean J2000 )and  Rdot.
      
      RVET = DSQRT(BODY_H(1)**2D0+BODY_H(2)**2D0+BODY_H(3)**2D0)

      Rdot = (HE0(1,2)*HE0(1,1)+HE0(2,2)*HE0(2,1)+HE0(3,2)*HE0(3,1))/ HD  

      Rdot = Rdot*AU/86400D0/1000D0                 ! Rdot Km/sec.              
   
!-------------Form " deldot "            
       
      DOT = (BEV1(1,2)*BEV1(1,1)+BEV1(2,2)*BEV1(2,1)+BEV1(3,2)*BEV1(3,1))/GD 

      deldot = DOT*AU/86400D0/1000D0 
      
!-----------------------------------------------------------------------------------------------      

      Ltime = DIST1 * 8.316746395045D0           ! Light time (min.)
      Mtime = DIST1 * 8.316746395045D0 * 60D0    ! Moon Light time (sec)

!----------------------------------------------------------------------------------------------
      
      COSI = (RVET**2D0+DIST1**2D0-EARTHVEC**2D0) / (2D0*RVET*DIST1)

      ANGLE = acos(COSI) * DR2D                     ! ANGLE = Phase angle
     
      ILL = ((RVET + DIST1)**2D0 - EARTHVEC**2D0) / (4D0 * RVET * DIST1) *100D0 ! % Illuminated fraction disk   
     
      OBL = 0.409092804022D0                         !new value OBL = 0.4090926006005829D0     
      
!-----------HLON = ELIOCENTRIC ECLIPTIC LONGITUDE ; HLAT = ELIOCENTRIC ECLIPTIC LATITUDE  

      RAR=atan2(HE0(2,1),HE0(1,1))
      DER=asin(HE0(3,1)/HD)

      HLON = atan2((sin(RAR)*cos(OBL)+tan(DER)*sin(OBL)),(cos(RAR)))
      HLAT = asin(sin(DER)*cos(OBL)-cos(DER)*sin(OBL)*sin(RAR))
      
      IF ( HLON < 0D0) THEN
       HLON = HLON + 2D0 * PI
      END IF

!----------- Form Magnitude and angular diameter 

      PANGLE = ANGLE
      CALL MAGNIT (ET,NTARG,ANGLE,RVET,DIST1,PVH,EARTHVEC,OBL,HLON,HLAT,ILL,MAGN,comment) 

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
       ILL = 1.0D0 *100D0
      END IF

!-----------Form  topocentric  position RA8,DE8
       
      CALL L_TIME(ET,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST,THD,XET,PV,HE0,THD,EPV)
      POS_T = POS_T        
      DIST4 = DIST
      CALL  iau_C2S ( POS_T, RA8, DE8 )       
      RA8 = iau_ANP(RA8)                  ! range angle (0,2!pi)

      IF(NTARG == 10) THEN
      DIAM = (1737.53D0 * 2D0 * DR2AS /( DIST * AU)*1000d0) ! arcmin-Topoc.Moon Diam. 
      END IF
      IF ( NTARG == 11 ) THEN
       DIAM = 1919.29D0 / DIST                              ! arcmin-Topoc. Sun Diam.
      END IF
             
           
!----------------------------------------------------------------------------
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun ( for apparent geocentric coordinate)

      CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
 
      CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

      CALL  iau_C2S   ( p, RA5, DE5 )       
      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
      
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent geocentric coordinate)

      CALL ABERR (AUDAY,u1,BVE,u0)

      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 
      
   
! ......... FINAL_ GCRS (RA5,DE5) of body 
! ......... Call for direction cosines to spherical coordinates
 
      CALL  iau_C2S   ( u3, RA10, DE10 )       
      RA10 = iau_ANP(RA10)                 ! range angle (0,2!pi)

!=====================================================================================
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun (for Topocentric coordinate)

      CALL DEFLIGHT_0(POS_T,BEP,HPE,u01)
 
!--------- Add in deflection due to Earth

      CALL DEFLIGHT_T(u01,BEP,HPE,u01)

      CALL iau_PN(u01,u02,p)          ! p = unit vector in order MU/C^2, u02 = modulus

!---------------------------------------------------------------------------       
!         RELATIVISTIC ANNUAL ABERRATION  (for Topocentric coordinate)

      CALL ABERR (AUDAY,u01,OV,u03)

      CALL iau_PN(u03,MO2,u5)         !Form unit vector of u03 = u5 ) 
      
   
!=====================================================================================
        
!---------- Form apparent geocentric coordinate (RA6,DE6)       

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

!      IF (mode == 1) then 
       CALL   iau_RXP   ( RNPB, u3, u4 )
       equinox = "Apparent geoc. posit. true equinox and ecliptic of the date.(IAU 2000 EQUINOX BASED)" 
!      ELSE IF (mode == 2) then 
!       CALL   iau_RXP   ( RC2I, u3, u4 )
!       equinox = "Apparent geoc. posit. true equinox and ecliptic of the date.(IAU 2006 CIO BASED)"  
!      END IF 

! ......... Call for direction cosines to spherical coordinates

      JDT = 2437663.5D0 - JD            ! Tempo trascorso tra JD e JDT = 1962/01/01
      TIME62 = 22645.5D0                ! Tempo in JD tra 1900/01/01 e 1962/01/01 
      CALL   iau_C2S   ( u4, RA6, DE6 )
      RA6 = RA6 - ((0.213D0 / TIME62 * JDT) * AS2R)    ! 50mas = Diff. RA origin
      DE6 = DE6 
         RA6= iau_ANP(RA6)                       ! range angle (0,2!pi)
 
      DIST3 = DIST1                              ! Distance 
!-------------------------------------------------------------------------------------------------
 
!---------- Form apparent topocentric coordinate RA7,DE7)

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

!      IF (mode == 1) then
       CALL   iau_RXP  ( RNPB, u5, u6 )
       equinox2 = "Apparent Topocentric position  ' no refraction'.(IAU 2000 EQUINOX BASED)"


! ......... Call for direction cosines to spherical coordinates

5000    CALL   iau_C2S   ( u6, RA7, DE7 )
        RA7 = RA7 - (0.0492D0 * AS2R)             ! 50mas = Diff. RA origin 
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)
         RA7H = RA7 * DR2D / 15D0                ! RA7 in hour
         DE7D = DE7 * DR2D
!--------------------------------------------------------------------------------------------

!--------Transform topocentric right ascension and declination to zenith distance and azimuth.
!        and compute atmosferic refraction.
      
      UJD = UT11                      ! UJD = UT1  Julian day       
                                      ! IREFR = 1 whit atmoswferic refraction
                                      ! IREFR = 0 no atmosferic refraction 
                        !  INPUT -----! XPO ,YPO polar motion (arcsec)
                                      ! longit,latit (degrees)  
                                      ! HM = heigh (meter)
                                      ! RA7H, DE7D (hour,degrees)

                                      ! ZD Zenith distance (degrees) 
                        !  OUTPUT-----! AZ Azimuth (degrees)
                                      ! RA8  RA corrected for refraction (degrees)
                                      ! DE8  DEC corrected for refraction (degrees)
     
      CALL SETDT ( DELTA )
     
      CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 0, &
                        ZD, AZ, RA8, DE8 )
      
      RA8 = RA8 * 15D0* DD2R
      DE8 = DE8 * DD2R 
      ZD0 = 90d0-ZD
      HOT = (LONG + GMST0 - RA7)* DR2D          ! HOT = Hour angle 

      ! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

      equinox3 = "Apparent Topocentric position with refraction.(IAU 2000 EQUINOX BASED)"
       
       CALL REFDAT (PHPA,TC,0.0,WL)
       CALL SETDT ( DElTA )
       
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 1, &
                        ZD, AZ, RA9, DE9 )
       RA9 = RA9 * 15D0* DD2R
       DE9 = DE9 * DD2R 
       HOB = (LONG + GMST0 - RA7)* DR2D
       ZOB = 90d0-ZD 
       
!=======================================================================================================
   
      IF ( NTARG == 10 ) THEN
        distance = " Km "
      ELSE
        distance = " AU "
      END IF  
      
      AUM = AU/1000D0

!-------- Compute Astrometric Sun Right ascension J2000
  
      Rsun = sqrt( SGE(1,1)**2d0 + SGE(2,1)**2d0 + SGE(3,1)**2d0)   ! Distance Earth - Sun 

      xs = SGE(1,1);  ys = SGE(2,1)

      RAsun = datan2( ys,xs)                  !Astrometric Sun Right ascension J2000(rad)
      RAsun =iau_ANP(RAsun)
      RAsun =RAsun * DR2D                     !Astrometric Sun Right ascension J2000(degrees)

      
!------------ Compute body's elongation to Sun

!-------Compute the Elongation " Elong " to the Sun in degrees

      ELONG = acos((Rsun**2d0 + DIST1**2d0 - RVET**2d0) / (2d0 * Rsun * DIST1))
      ELONG = ELONG * DR2D


      
      EVET = DSQRT( EPV(1,1)**2D0 + EPV(2,1)**2D0 + EPV(3,1)**2D0)   ! Dist. Sun-Earth at XET time
      EPOS = ATAN2(EPV(2,1),EPV(1,1)) 
      EDEC = ASIN (EPV(3,1) / EVET)
      
      SUNLONG = ATAN2((sin(EPOS)*cos(OBL)+tan(EDEC)*sin(OBL)),(cos(EPOS))) + PI  !Geocentric Longit. Sun  
      
      
      GLON = ATAN2((sin(RA6)*cos(OBL)+tan(DE6)*sin(OBL)),(cos(RA6)))
      GLAT = ASIN(sin(DE6)*cos(OBL) - cos(DE6)*sin(OBL)*sin(RA6))
      COSELON = cos(GLAT) * cos(GLON - SUNLONG)
      
      
!      ELONG = DACOS ( COSELON) * DR2D
      GLON = GLON * DR2D
      SUNLONG = SUNLONG * DR2D


!------- Search relative position body to Sun ( E = East  or W = West) 

      if(RA >= RAsun) then
       if((RA-RAsun)>=180)then
        EL="West"
        else if((RA-RAsun)<180)then
        EL="East"
       end if
      else if(RA<RAsun) then
        if((RAsun-RA)>=180)then
        EL="East"
        else if((RAsun-RA)<180)then
        EL="West"
       end if
      end if
 
      
      
      IF (NTARG == 11) ELONG = 0.0D0    
!------------------------------------------------------------------------------


!--------- Compute  -TDB0- JD Time to 0.0 h TDB the day (for routine TRANSIT and TWILIGHT)

        TDB0 = INT(MJD) + DJMJD0  
       
       PL = ASIN(6378.137/(DIST3 * AUM)) * DR2D  ! Lunar horizontal parallax (degs)
       SD = (DIAM / 3600D0 ) / 2D0               ! Lunar semidiameter in degs
          
       CALL TRANSIT(NTARG,TDB0,GAST0,AUDAY,DELTAT,PL,SD,LAT,LONG,DX00,DY00,RIS,SET,TRS,  &
            AZM1,AZM2,HAL,As,NR,MS,PT,LX,LM,LN) 
            
       RIS = RIS / 24D0
       SET = SET / 24D0
       TRS = TRS / 24D0 
       
       CALL iau_D2TF ( 1, RIS, SIGN, IHMSF1)     
       CALL iau_D2TF ( 1, SET, SIGN, IHMSF2)
       CALL iau_D2TF ( 1, TRS, SIGN, IHMSF3)
       
       IF (LX == 1 ) THEN
         commentM = '      ------ Body circumpolar ----- '
        ELSE IF (LX == 0) THEN
         commentM = ""
        END IF
       IF(LM == 1) THEN
         commentN = ' .... BRIGHT .... '
         ELSE IF (LM == 2) THEN
         commentN = ' .... NO SET - body always up the horizon '
         ELSE IF (LM == 0) THEN
         commentN = " "
       END IF
       IF(LN == 1) THEN
         commentA = ' .... DARK .... '
         ELSE IF (LN ==2) THEN
          commentA =  '... NO RISE - body always below the horizon '
         ELSE IF(LN == 0) THEN
          commentA =  " "
       END IF 
   
       IF ( NR == 1 )THEN
           AS1 = "* "
           commentR = '** Rise - event occurs the next day'
         ELSEIF ( NR == 2) THEN
           AS1 = "**"
           commentR = '** Rise - event occurs the previus day'
         ELSEIF ( NR == 0) THEN
           AS1 = "  " 
           commentR = ""
       END IF
       IF ( MS == 1 )THEN
           commentS = '** Set - event occurs the next day'
           AS2 = "* "
        ELSEIF ( MS == 2 )THEN
           AS2 = "**"
           commentS = '** Set- event occurs the previus day'
        ELSEIF ( MS == 0) THEN
           AS2 = "  " 
           commentS = "" 
       END IF
       IF ( PT == 1 )THEN
           AS3 = "* "
           commentT = '** Transit- event occurs the next day'
         ELSEIF ( PT == 2) THEN
           AS3 = "**" 
           commentT = '** Transit- event occurs the previus day'
         ELSEIF ( PT == 0) THEN
           AS3 = " " 
           commentT = ""
       END IF
      
     IF(NTARG == 11) THEN
      CALL TWILIGHT(NTARG,TDB0,GAST0,AUDAY,DELTAT,LAT,LONG,DX00,DY00,Ha, &
           BCTW,ECTW,BNTW,ENTW,BASTW,EASTW,N1,N2,N3)
      
      CALL iau_D2TF ( 1, BCTW, SIGN, IHMSF10)     
      CALL iau_D2TF ( 1, ECTW, SIGN, IHMSF11)
      CALL iau_D2TF ( 1, BNTW, SIGN, IHMSF12)
      CALL iau_D2TF ( 1, ENTW, SIGN, IHMSF13)     
      CALL iau_D2TF ( 1, BASTW, SIGN, IHMSF14)
      CALL iau_D2TF ( 1, EASTW, SIGN, IHMSF15)  
      
        IF(N1 == 1 .or. N2 == 1 .or. N3 == 1) THEN               ! Twilight parameter
          commentP = '---- 00:00:00.0 =   CONTIUOUS TWILIGHT --- '
        ELSE IF(N1 == 0 .and. N2 == 0 .and. N3 == 0) THEN
          commentP = ""
        END IF
       ELSE IF (NTARG /= 11) THEN
          commentP = ""
     END IF 

        IF(LM == 2) THEN
          commentP = " " 
        END IF

!====================================================================================================

!---------- Output.

      WRITE(*,*)
      WRITE(*,115) longit,lS,latit,bS,height
115   FORMAT ( 1x,"Longitude",1x,"deg",F10.5,1x,a,3x,"Latitude",1x,"deg",F9.5,1x,a,3x,"Height",1x,F6.1,1x,"meter")       
      WRITE(*,*)
      WRITE(*,116) Dday
116   FORMAT ( 1x,"Day of week",1x,a10)
      WRITE(*,*)
      CALL iau_D2DTF ( "UT1", 6, UT1, DJMJD0, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 11 ) "UT1", IY, MO, ID, IHMSF, UT1
      CALL iau_D2DTF ( "TT", 6, TT, DJMJD0, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 11 ) "TT ", IY, MO, ID, IHMSF, TT
      CALL iau_D2DTF ( "TCG", 6, TCG1, TCG2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 11 ) "TCG", IY, MO, ID, IHMSF, TCG - DJMJD0
      CALL iau_D2DTF ( "TDB", 6, TDB1, TDB2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 11 ) "TDB", IY, MO, ID, IHMSF, TDB - DJMJD0
      CALL iau_D2DTF ( "TCB", 6, TCB1, TCB2, IY, MO, ID, IHMSF, J )
      IF ( J.NE.0 ) STOP
      WRITE ( *, 11 ) "TCB", IY, MO, ID, IHMSF, TCB - DJMJD0
11    FORMAT ( 1X,A,I5,2("/",I2.2),1X,I3,2(":",I2.2),".",I6.6,"   MJD",F20.10 )

      WRITE(*,117) " DELTA_T (TDB-UT) ",DELTA
117   FORMAT(a17,"  ",F9.6," sec.")
      WRITE(*,1117) DUT
1117  FORMAT(x,"UT1-TAI value (1900-1962) ",f11.7," sec")
      WRITE(*,150)" TT - TDB at geocenter: ",PVV(1),"sec",PVV(2),"sec/day"
      WRITE(*,*)
      WRITE(*,118) GST*DR2D
118   FORMAT('     Greenwich Sideral Time (GAST) =',f19.13,' deg', " (IAU 1976-1980 )")
      WRITE(*,1118) GAST*DR2D
1118   FORMAT('     Greenwich Sideral Time (GAST) =',f19.13,' deg', " (IAU 2006 )")
      WRITE(*,119) ERA*DR2D
119   FORMAT('     Earth Rotation Angle ( ERA )  =',f19.13,' deg'," (aproximate prior the year 1973 )")
      WRITE(*,*)
      WRITE(*,*)" ============= Use IERS EOP C01 data(1900-now) ==============="
      WRITE(*,*)
      WRITE(*,*)
      WRITE (*,136) XP * DR2AS
136   FORMAT('     Polar motion XP = ',f9.6,' arcsec.')   
      WRITE (*,137) YP * DR2AS
137   FORMAT('     Polar motion YP = ',f9.6,' arcsec.')   
      WRITE(*,*)
      WRITE (*,138) dpsi_*1000D0
138   FORMAT('     Celestial pole offset dpsi = ',f9.4,' mas')   
      WRITE (*,139) deps_*1000D0
139   FORMAT('     Celestial pole offset deps = ',f9.4,' mas') 
      WRITE(*,*)
      WRITE (*,122)DPSI * DR2AS
122   FORMAT("     Nutation in longitude ",3x,F10.6,1x,"arcsec.")             

      WRITE (*,123)DEPS * DR2AS
123   FORMAT("     Nutation in obliquity ",3x,F10.6,1x,"arcsec.")        
      WRITE (*,*)
      WRITE (*,*)"================================================================"
      WRITE (*,*)
 

      WRITE (*,125)" Body barycenter : ", planet

125   FORMAT (a17,a24 )
      WRITE(*,*)
      WRITE(*,*) "        Barycentric equatorial position and velocity  "
      WRITE(*,*) "             mean equator and equinox of J2000.       "
      WRITE (*,121)"             X                   Y                     Z "
      WRITE (*, 129)" POS ", PVB(1,1), PVB(2,1), PVB(3,1)," AU"
      WRITE (*, 129)" VEL ", PVB(1,2), PVB(2,2), PVB(3,2)," AU/day"
      WRITE (*,*)
      WRITE (*,*)
      WRITE(*,*) " Heliocentric ecliptic coordinate mean equator and equinox J2000 ( λ, B, RV )"
      WRITE (*,131) "        o , ' , ''               o , ' , ''            RV AU          Rdot (Km/s)"
      CALL iau_A2AF ( 5 ,HLON ,sign , IDMSF0 )
      CALL iau_A2AF ( 4 ,HLAT ,sign , IDMSF )
      WRITE ( *,140 ) IDMSF0 ,sign,IDMSF,RVET,Rdot
      WRITE (*,*)"========================================================================================"
      WRITE (*,*) "Astrometric geocentric J2000 coordinate with relativistic corrections"
      WRITE (*,132) "            R.A                      Dec.            Distance         deldot (Km/s)               "
      WRITE (*,131)"  h , m ,  s              o , ' , ''                                     "         
      CALL iau_A2TF ( 5 ,RA3 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE3 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,1302 ) IHMSF ,sign,IDMSF,DIST1 * AUM,deldot
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *,1301 ) IHMSF ,sign,IDMSF,DIST1,deldot
      END IF
      WRITE (*,1200)"=================================================================================================="
      WRITE (*,*) "Geocenter mean J2000  equatorial coordinate "
      WRITE (*,132) "            R.A                      Dec.            Distance          Light time                 "
      WRITE (*,131) "  h , m ,  s              o , ' , ''                                     "              
      CALL iau_A2TF ( 5 ,RA4 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE4 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,1300 ) IHMSF ,sign,IDMSF,DIST2 * AUM,Ltime*60D0
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *, 130 ) IHMSF ,sign,IDMSF,DIST2,Ltime
      END IF
      WRITE (*,1200)"=================================================================================================="
      WRITE (*,*)equinox
      WRITE (*,132) "            R.A                      Dec.            Distance                                     "
      WRITE (*,131) " h , m ,  s              o , ' , ''                                     "
      CALL iau_A2TF ( 5 ,RA6 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE6 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,1306 ) IHMSF ,sign,IDMSF,DIST3 * AUM
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *,1307 ) IHMSF ,sign,IDMSF,DIST3
      END IF
      WRITE (*,1200)"==================================================================================================" 
      WRITE (*,*)equinox2 
      WRITE (*,161)                                      
      WRITE (*,160)                         
      CALL iau_A2TF ( 5 ,RA8 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE8 ,sign , IDMSF )
      WRITE ( *,135 )IHMSF ,sign,IDMSF,AZ,ZD0,HOT
      WRITE (*,1200)"==================================================================================================" 
      WRITE(*,*)equinox3
      WRITE (*,132) "           R.A                       Dec.             Dist.        (*)Elevation   (*)Hour Angle  "
      WRITE (*,131) "         h , m ,  s              o , ' , ''                            degs.           degs.     "
      CALL iau_A2TF ( 5 ,RA9 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE9 ,sign , IDMSF )
      IF ( NTARG == 10 ) THEN
        WRITE ( *,1306 ) IHMSF ,sign,IDMSF,DIST4 * AUM,ZOB,HOB
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *,1307 ) IHMSF ,sign,IDMSF,DIST,ZOB,HOB
      END IF

      WRITE (*,*)"========================================================================================" 
      WRITE(*,*)
      WRITE(*,*)
      WRITE (*,132) "       ILL. %       MAGN.        Angular Diameter        Phase Angle (°)       Sun's Elong.         "
      IF ( NTARG == 10 .OR. NTARG == 11 ) THEN
        WRITE ( *,1304 ) ILL,MAGN,DIAM/60D0,PANGLE,ELONG,EL
       ELSE IF( NTARG < 10 ) THEN
        WRITE ( *,1303 ) ILL,MAGN,DIAM,PANGLE,ELONG,EL 
       IF ( NTARG == 5 .OR. NTARG == 6 .OR. NTARG == 7 .OR. NTARG == 8) THEN
        WRITE (*,1205) PDIAM
       ELSE IF ( NTARG /= 5 .OR. NTARG /= 6 .OR. NTARG /= 7 .OR. NTARG /= 8) THEN
        WRITE ( *,1206 )  
       END IF
       IF ( NTARG == 1 .OR. NTARG == 2) THEN
        WRITE ( *,1206 )comment 
       END IF 
      END IF
      WRITE(*,1200)"======================================================================================================"
      WRITE (*,141) "          RISE (UT)                TRANSIT (UT)               SET (UT)                               "
      WRITE (*,*)
      WRITE (*,142) IHMSF1,AS1, IHMSF3,AS3, IHMSF2,AS2
      WRITE (*,*)
      WRITE (*,141) "        Azimuth                    Altitude                 Azimuth                              " 
      WRITE (*,143) AZM1,HAL,As,AZM2 
      WRITE (*,*)
            IF(NTARG == 11)THEN 
      WRITE (*,141)"               MORNING                                          EVENING                           "         
      WRITE (*,270)IHMSF10,IHMSF11
      WRITE (*,271)IHMSF12,IHMSF13
      WRITE (*,272)IHMSF14,IHMSF15
            END IF
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,147) 
      WRITE (*,144)planet,commentR
      WRITE (*,145)commentP 
      WRITE (*,145)commentS
      WRITE (*,145)commentT
      WRITE (*,146)commentM
      WRITE (*,146)commentN
      WRITE (*,146)commentA
      WRITE (*,*)
      WRITE (*,*)"   (*) The option 'IUA 2000 Equinox based'          "
      WRITE (*,*)"       compute refraction only for zenith distance "
      WRITE (*,*)"       between 0.1 and 91 degrees.                 "
      WRITE (*,*)"   (*) The calculated value of the magnitude of Mercury and Venus is    "
      WRITE (*,*)"       not real when the phase angle of the two celestial bodies is too "
      WRITE (*,*)"       close to the Sun to be observed photometrically, so the reported "
      WRITE (*,*)"       value is invalid.                                                "
      WRITE (*,*)      
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,*)
      WRITE (*,1200)"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
      

121    FORMAT (A60)
129    FORMAT (A5,' ',F22.16,' ', F22.16,' ', F22.16,' ',A7)    
130    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,5x,F10.6," min.")       
131    FORMAT (A98)
132    FORMAT (A98)
135    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F10.6,7x,F11.6,5x,F11.6)
140    FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,5x,F10.6) 
141    FORMAT (a98)
142    FORMAT ( 6X," h ",I3,2(":",I2.2),".",I1.1,a2,11X," h ",I3,2(":",I2.2),".",I1.1,a2,13X," h ",I3,2(":",I2.2),".",I1.1,a2) 
143    FORMAT ( 9X,F5.1,"°",22X,F4.1,"°",a,19X,F5.1,"°")
144    FORMAT (10x," The",1x,A10,1x,A90)
145    FORMAT (25x,A40)
146    FORMAT (25x,A40) 
147    FORMAT (5x," *** NOTES *** ")
160    FORMAT (10x,"h , m ,  s",14x,"o , ' , ''",12x,"degs",12x,"degs",11x,"degs")
161    FORMAT (14x,"R.A",21x,"Dec.",13x,"Azimut",6x,"(*)Elevation",2x,"(*)Hour Angle")
270    FORMAT (5x,"Beginnig civil twilight        "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending civil twilight        ", &
              " h ",I3,2(":",I2.2),".",I1.1)
271    FORMAT (5x,"Beginnig nautical twilight     "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending nautical twilight     ", &
              " h ",I3,2(":",I2.2),".",I1.1)
272    FORMAT (5x,"Beginnig astronomical twilight "," h ",I3,2(":",I2.2),".",I1.1,5x," Ending astronomical twilight ", &
              " h ",I3,2(":",I2.2),".",I1.1)
1200   FORMAT (A100) 
1205   FORMAT (30x,"Polar.",F10.3," arcsec.") 
1206   FORMAT (15x,a14,1x,"Polar.",3x," *****",1x," arcsec.")      
1300   FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F11.3,7x,F9.6," sec.")
1301   FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F13.10,5x,F10.6)
1302   FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,5x,F11.3,5x,F10.6)        
1303   FORMAT (7x,F6.2,5x,F6.2,6x,"Equat.",F10.3," arcsec.",6x,F6.2,12x,F5.1,1x,5a)  
1304   FORMAT (7x,F6.2,5x,F6.2,9x,F8.3," arcmin.",10x,F6.2,12x,f5.1,1x,5a) 
1306   FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,3x,F11.3,6x,F11.6,5x,F11.6) 
1307   FORMAT (7x,I3,2x,I2,2x,I2,".",I5.5,6x,a,I3,2x,I2,2x,I2,".",I4.4,3x,F13.10,6x,F11.6,5x,F11.6) 
  
     
3000  continue 
      WRITE(*,*)
3001  WRITE(*,*)" *****************************************************************" 
      WRITE(*,*)" want to search for a new ephemeris ?  'Y' or 'N'                 "
      WRITE(*,*)" or or choose different times ( '+' forwards or '-' for backwards "
      WRITE(*,*)" *****************************************************************" 
      WRITE(*,*)""
      READ *, NEW
      IF(NEW == "y")NEW ="Y"
      IF(NEW == "n")NEW ="N" 
             
      IF ( NEW == "Y" .or. NEW == "y") THEN
        PRINT*, "  CHOICE ... ",NEW," ===> CONTINUE " 
        GOTO 2
       ELSE IF (NEW == "N" .or. NEW == "n" ) THEN
        PRINT*, "  CHOICE ... ",NEW," ===> EXIT " 
        GOTO 6000
       ELSE IF (NEW == "+" ) THEN
        PRINT*, "  CHOICE ... ",NEW," ===> CONTINUE " 
        GOTO 4000    
       ELSE IF (NEW == "-" ) THEN
        PRINT*, "  CHOICE ... ",NEW," ===> CONTINUE " 
        GOTO 4000    
      END IF
      IF ( NEW /= 'Y' .OR. NEW /= "N" .OR. NEW /= '+' .OR. NEW /= '-' )GOTO 3001

4000  print*," Input time interval: ' DAY, HOURS, MINUTES, SECONDS '. Example:  0,10,0,0  (for + 10 hours)"
      read *,d_d,h_h,m_m,s_s
      print*,d_d,"day",h_h,"hours",m_m,"minutes  ",s_s,"seconds" 
   
      pas = d_d+h_h/24d0+m_m/1440d0+s_s/86400d0
      
      IF ( NEW == "+") THEN
       Julian = MJD+pas
      ELSE IF( NEW == "-") THEN
       Julian = MJD-pas
      END IF
      newtime = Julian 
      
      CALL iau_JD2CAL ( DJMJD0, newtime, IY, MO, ID, FD, J )
      FD=FD
      CALL time_hms ( FD, IH, IM, SEC)
      WRITE(*,*)" "
      WRITE(*,4500) IY, MO, ID, IH, IM, SEC
4500  FORMAT( " NEW TIME ",x,i4,"/",i2,"/",i2,4x,i2,":",i2,":",f9.6)
      WRITE(*,*)" "
   
      GOTO 5100 
       

6000  continue    
        
 END program planets
      
