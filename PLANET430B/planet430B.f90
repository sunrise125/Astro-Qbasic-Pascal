 program planets

!     gfortran -w planet430B.f90 1818.for jpl_430.f sky6.for

!     THIS PROGRAN COMPUTES THE APPARENT DIRECTION OF SOLAR SYSTEM
!     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM
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



      INTEGER IY, MO, ID, IH, IM, J, I,mode,IREFR,SGN,ios,modetime,NVS,NTARG,NCENT,ISTEP,NSTEP
      INTEGER IHMSF1(4),IHMSF2(4),IHMSF3(4),IHMSF(4),IDMSF(4), IDMSF0(4),IHMSF4(4),IHMSF5(4)
      INTEGER FINAL,JW,EQTM,IY1, IM1, ID1,IDMSF1(4),d_d,h_h,m_m,conta
      DOUBLE PRECISION SEC,XP,YP,DUT,TT,TCB1,TCB2,T,AU,PS,DEL,HA,TDB,Ltime,TDB1,TDB2,TDB0,OBL,s_s
      DOUBLE PRECISION UTC1, UTC2, UT11, UT12, UT,TAI1, TAI2, TT1, TT2, TCG1, TCG2,RVETT,xs,ys,pas
      DOUBLE PRECISION DDP80, DDE80, DX00, DY00, DX06, DY06,DELTA_T,DE_T,u01(3),Rdot,GN,Mtime,newtime
      DOUBLE PRECISION DATE, TIME, UTC, DAT,TCB,TCG,POS2(3),DOT,TAI, TUT, UT1, DP80, DE80, MJD0,MJDprova,DJM
      DOUBLE PRECISION RA1,DE1,RA2,DE2,RA3,DE3,RA4,DE4,RA5,DE5,RA6,DE6,RA7,DE7,RA8,DE8,RA9,DE9,ddec,DJM0
      DOUBLE PRECISION RA10,DE10,RMASSE,SUNLONG,COSELON,EDEC,EVET,GLON,GLAT,XEL,DIST,DIST4,ECC,dra         
      DOUBLE PRECISION DPSI, DEPS, EPSA, ERA, DP00, DE00, DDP00, DDE00 ,ELON, DELTAT, EARTHVEC,RA_3 
      DOUBLE PRECISION GAST,LONG ,LAT,HM,DTR,MO1,FD,julian,DSQRT,HLON,HLAT,radius,T0,TDBJD,ET,XET,DE_3
      DOUBLE PRECISION  X, Y, S , DT, ELONG, u02, MO2, longit, DELTA, DS, AUDAY, XPO, YPO, DIAM,AE      
      DOUBLE PRECISION latit,UJD,MJD,ANGLE,PANGLE,PDIAM,AZ,height,ZD,DJMJD0,MAGN,AUM,Dtime,IMJD,PVV(2)
      DOUBLE PRECISION RA7H, DE7D,JD,RVET,RVET0,EPOS,EPSA0,RB0,RP0,RPB0,RN0,DP000,DE000,DPSI0,DEPS0,MM
      DOUBLE PRECISION DIST5,DX0,DY0,UT_0,TT_0,GAST0,PL,SD,RIS,SET,TRS,AZM1,AZM2,HAL,MUC,BP,de_S
      DOUBLE PRECISION C,MU,u2,RAR,DER,K,TAU,THD,TGD,ILL,COSI,DET,HD,GD,HMOD,DIST1,DIST2,DIST3,cosDE,sinRA
      DOUBLE PRECISION RC2IT(3,3),RC2I(3,3),RPB(3,3),RNPB0(3,3),RN(3,3),RNPB(3,3), RP(3,3),RBP(3,3)
      DOUBLE PRECISION PV(3,2),PVH(3,2),PVB(3,2),BEV1(3,2),BEV(3,2),RB(3,3),GMB(3,2),PVE(3,2),HSV(3,2)
      DOUBLE PRECISION HE0(3,2),SUNB(3,2),EPV(3,2),GPB(3),HV(3),POS(3),VEL(3),XYZ(3),V1(3),RA66,RA,RA0
      DOUBLE PRECISION EB(3),u0(3),u03(3),u1(3),u3(3),u4(3),u5(3),u6(3),U(3),V(3),V2(3),p(3),RAsun,DE0 
      DOUBLE PRECISION VEL2(3),VEL3(3),POS3(3),POS_H(3),GE(3),HPE(3),BVE(3),BEP(3),TSL,EO ,BGV(3,2)
      DOUBLE PRECISION BODY_G(3),BODY_H(3),HE(3),RS(3),OB(3),OV(3),GP(3),GV(3),POS_T(3),EO0,GEB,eta,PI_
      DOUBLE PRECISION iau_ANP,iau_GMST00,iau_S06,iau_DTDB,iau_GMST06,iau_GST06A,iau_GST00A,CHI,p_,A_,B_    
      DOUBLE PRECISION iau_ERA00, iau_EE00,iau_EE06A,iau_EORS,iau_EE00A,off_dpsi,off_deps,deldot,C_
      DOUBLE PRECISION AOB,ZOB,HOB,DOB,ROB,PHPA,TC,RH,WL,ZOT,AOT, HOT, DOOT,ROT,EQEQ,GMST0,ZD0,GEL,atlon  
      DOUBLE PRECISION Sun_lon,Earth_lon,TDBJ,iau_FAE03,CORR,EQT,EQTS,EQT0,LOD,EQTD,LatLu,AGE,DJ2,Rsun 
      DOUBLE PRECISION TT01,TT02,DE_GMST,EQTMG,UTC0,UTC01,UTC02,RA06,EQTX,BODY_sun(3),R_sun,JHE,MJ,tanAX,AX 
      DOUBLE PRECISION BODY_H1,GD1,HD1,XET1,PV1,HE01,THD1,E_P(3,2),D_sun,DISTKm,LoLu,LoLu_ap,OMEGA,OME0
      DOUBLE PRECISION JAP,DLR,DJ3,PH,TLong,TLat,LTL,LTB,PA,iau_OBL06,LON0,LAT0,T_0,a0,b0,c0,THL,THB,LAV(6) 
      DOUBLE PRECISION LOGE,BGE,Distmoon,R_sun1,GELO,GELA,GELOAP,GELOD,GELAD,GELOAPD,VC,CVE(3),BM1,BCV,er
      DOUBLE PRECISION Tcy,PHE,KA,Lsun,Bsun,SLO,SLA,SGP(3),DLO,DLA,GELAAPD,GELAAP,AZ1,prova,SGE(3,2),DETDB
      DOUBLE PRECISION RVET1,Rdot1,Ltime1,deldot1,Mtime1,HD0,BEV0(3,2),SUNBA(3,2),HSVT(3,2),PVBT(3,2),LEA(6)
      DOUBLE PRECISION  R(6)
      DOUBLE PRECISION  SS(3)
      DOUBLE PRECISION  VALS(400)
      DIMENSION body(17),radius(10),dayweek(7)
      
      
      CHARACTER(len=6):: NAMS(400),TX,choose,yes,no,distance,EL
      CHARACTER(len=16):: XDIAM,dayweek,Dday,IDEM 
      CHARACTER(len=26):: body,planet,LDUT,LDET,EQ_label,Slod
      CHARACTER(len=1):: sign,aS,lS,bS,update,NEW,str
      CHARACTER(len=95):: equinox,equinox2,equinox3,commentR,commentS,commentT,commentU,labelL
      CHARACTER(len=95):: labelR,labelS,labelT,labelM,labelN,labelO,labelX,labelY,commentM,commentN
      CHARACTER(len=7):: APE 
      CHARACTER(len=14):: comment 
      common SGE,RVET0,HSV                   

      DATA body / "MERCURY","VENUS ","EARTH ","MARS ","JUPITER","SATURN ","URANUS ","NEPTUNE","PLUTO ", &
                  "MOON ","SUN ","SOLAR-SYSTEM BARYCENTER","EARTH-MOON BARYCENTER","NUTATIONS", &
                  "Lunar Euler angles"," Lunar libration","TT - TDB" / 
                  
      DATA radius / 2440D0, 6051.8D0, 6378.137D0, 3397D0, 71492D0, &
                    60268D0, 25559D0, 24766D0, 1151D0, 1737.53D0 /

      DATA dayweek / "SUNDAY","MONDAY","TUESDAY","WEDNESDAY","THURSDAY","FRIDAY", &
                     "SATURDAY" /

      WRITE (*,*)"  **************************************************************************"
      WRITE (*,*)"  *     Last release of Planet, renamed planet430B, to avoid mismatches    *"
      WRITE (*,*)"  *     New SOFA routines added into 1818.for , JPL DE430t database,       *" 
      WRITE (*,*)"  *     IERS Final data (IAU2000), and EOP C01 1900-now modified           *"           
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *      --------------------planet430B.f90 -------- (Ago.30,2018) --------*"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     THIS PROGRAN COMPUTES THE APPARENT DIRECTION OF SOLAR SYSTEM       *"
      WRITE (*,*)"  *     BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE SYSTEM      *"
      WRITE (*,*)"  *            for this purpose we use JPL DE430t database                 *"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     Compiling command:                                                 *"
      WRITE (*,*)"  *     gfortran planet430B.f90 sky6.for 1818.for jpl_430.f                *"
      WRITE (*,*)"  *     which generates the executable file: a.out                         *"
      WRITE (*,*)"  *     You can also use the program with INTEL_FORTRAN in quadruple       *"
      WRITE (*,*)"  *     precision with the command:                                        *"
      WRITE (*,*)"  *     ifort -r16 planet430B.f90 1818.for jpl16.for sky6.for              *"
      WRITE (*,*)"  *     which generates the executable file: a.out                         *"
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  *     Range of calendar dates:                                           *"
      WRITE (*,*)"  *      a) Not less 1900-01-01                                            *"
      WRITE (*,*)"  *      b) Not over MJD of line 128 (or nearby), i.e. FINAL = 58492 of    *"
      WRITE (*,*)"  *         main file planet430B.f90.                                      *"
      WRITE (*,*)"  *      c) Ephemerides computations may be forced up to 2100/01/01 upper  *"                 
      WRITE (*,*)"  *         limit of JPL database,but with less precision.                 *"
      WRITE (*,*)"  *         ------------------------------------------------------------   *" 
      WRITE (*,*)"  *         For final-iers.txt updating, please refer to : www.iers.org/   *"
      WRITE (*,*)"  *       - Data/Products - Standard EOP data files - Finals all IAU2000 - *"  
      WRITE (*,*)"  *         and change the FINAL value accordingly, say MJD shown 3 lines  *"
      WRITE (*,*)"  *         prior the last full row.                                       *" 
      WRITE (*,*)"  *                                                                        *"
      WRITE (*,*)"  **************************************************************************"
      WRITE (*,*)
                   

2      FINAL = 58730       ! 2019/09/04   reported in MJD IERS final database (for update) 
!      PRINT (*,*) " MJD currently final number, FINAL = ",FINAL,"2019/09/04  "


            
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
      100   WRITE (*,*)" Choose Time input!  'UTC' = 1 ( 'UT1' prior the 1973.0) or 'TT' = 2 " 
      READ *, modetime

      IF(modetime > 2)THEN
       GOTO 100
      ELSE IF(modetime == 1) THEN
       TX = "UTC"
       GOTO 200
      ELSE IF (modetime == 2) THEN
       TX = "TT"
       GOTO 200 
      END IF 
 
200   WRITE (*,*)
      WRITE (*,*) "       Date and Time " ,TX   
      WRITE (*,*) "   =============================="
210   WRITE (*,*)
      WRITE (*,*) " Input : Year , Month, Day (yyyy,mm,dd) "
      READ *, IY, MO, ID
       IF(IY < 1900  .OR. IY == 2040 .AND. MO == 2) THEN              
        WRITE(*,*)"*****  YEAR out of range ! *****"
        GOTO 210
       ELSE  
        GOTO 215
       END IF 
                
215   WRITE (*,*)
      WRITE (*,*) " Input : Hour, Minutes, Seconds (hh,mm,ss.xxx..) "
      READ *, IH, IM, SEC

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
      TSL = 20.0                           ! Is the approximate sea-level air temperature in K"
      RH = 0.5                             ! relative humidity at the observer (range 0-1)
      WL = 0.55                            ! wavelength (micrometers)
      KA = 20.49552d0                      ! Costant of aberration 
!---------------------------------------------------------------------------------------  

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
       
!---------- Compute JD (Julian day) date .--------------------------

      GN = (100D0*IY+MO-190002.5)                    !I find the sign of the numerical expression
      IF (GN < 0D0) THEN
       SGN = -1
      ELSE IF (GN > 0D0) THEN
       SGN = +1
      END IF

!      JD = 367*IY - (7*(IY+((MO+9)/12)) / 4)+(275*MO/9) + ID + 1721013.5d0  & 
!             + TIME - 0.5d0 * SGN +0.5d0

      call iau_CAL2JD( IY, MO, ID, DJM0, DJM, J)
      MJD = DJM + TIME                          ! UTC or TT in MJD (Julian Day modified)
      DJMJD0 = 2400000.5D0
      JD = DJMJD0 + MJD                         ! UTC or TT in Julian Day
             
!----------Compute the day of week---------------------------------

      JW = INT((JD + 1.5D0) - 7D0 * Int((JD + 1.5D0) / 7D0)) 
      
      DO I = 1,7
       IF ( JW == I-1) THEN
        Dday = dayweek(I)
       END IF
      END DO  
!-------------------------------------------------------------------    
 
      IF (MJD .GE. 41687) then        ! MJD 41687 = JD 2441687.5 = 1973/01/05.0
       GOTO 178 
      ELSE IF (MJD .LT. 41687)THEN
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
      
      IF( MJD > FINAL ) THEN
       de_S = (MJD - 51544.5D0) * DS                ! Time in sec. from J2000 (JD 2451545.0 TT)
       MM = 6.239996D0 + de_S * 1.99096871D-7
       DELTA = 32.184 + DAT + 0.001657D0 *sin(MM + 0.01671D0 * sin(MM))
      END IF                                                 
       
!----------------------------------------------------------------------iau_GST00A
! Compute TT and trasform TT to UTC
      
      IF (modetime == 1) GOTO 175
      IF (modetime == 2) GOTO 170 
!----------TT (MJD).

170   CALL iau_CAL2JD ( IY, MO, ID, DJMJD0, DATE, J )              ! input TT (IY,MO,ID)
      TT = DATE + TIME                                             ! TT value (MJD)
                                                                                   
!----------UTC,TAI (MJD)      
       
      UTC = TT - (32.184D0 + DAT)/DS
      CALL iau_JD2CAL ( DJMJD0, UTC, IY, MO, ID, FD, J )
      CALL time_hms ( FD, IH, IM, SEC)

!---------Transform into internal format.
175   CALL iau_DTF2D ( "UTC", IY, MO, ID, IH, IM, SEC, UTC1, UTC2, J )   !UTC1 + UTC2 = Julian Day
      IF ( J.NE.0 ) STOP                                                 !UTC1 is the Julian Day number and
      CALL iau_DTF2D ( "UTC", IY, 1, 1, 0, 0, 0D0, UTC01, UTC02, J )     !UTC2 is the fraction of a day.   
    
      UTC = UTC1 + UTC2 - DJMJD0                                         ! UTC in MJD (modified JD)

      DATE = UTC1 - DJMJD0                                               ! DATE is the MJD number
      UTC0 = UTC01 + UTC02                                               ! UTC JD at 0.0h YYYY/01/01  
!----------Call the interpolation routine for per XP,YP,DUT,dX,dY
       
      CALL iers_calc(UTC,DATE,FINAL,XPO,YPO,DUT,LOD,DX0,DY0)
      
      DUT = DUT                          !  sec.
      TUT = UTC2 + DUT/86400D0           !  in fraction of day          
      LOD = LOD
      XP = XPO * AS2R                    ! arcsec---> radians
      YP = YPO * AS2R   
      
      DX00 = DX0/1000D0 * AS2R           ! mas---> radians 
      DY00 = DY0/1000D0 * AS2R 

      IF( MJD > FINAL ) THEN             !
       DUT = 0D0
       LOD = 0D0                         !
       TUT = UTC2                        !
       XP = 0D0                          !
       YP = 0D0                          !
       DX0 = 0D0                         !
       DY0 = 0D0                         ! IF MJD > FINAL
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

      TT  = DATE+TT2          ! 
      TDB = DATE+TDB2         !
      TAI = DATE+TAI2         ! Time show in MJD
      TCG = DATE+TCG2         ! 
      TCB = DATE+TCB2         !
      UT1 = DATE+UT12         !

      IF( MJD > FINAL ) THEN
       DET = DELTA
       DELTAT = DELTA/86400D0
       GOTO 195
      ELSE IF ( MJD <= FINAL ) THEN
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

!====================================================================================


             
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

      prova = (GAST-RA4) * DR2D

!----------- Transform RA3, DE3 in astrometric-topocentric 

      
!--------------Form RVET0 (vector ray mean J2000 )and  Rdot.
      
      RVET0 = DSQRT(BODY_H(1)**2D0+BODY_H(2)**2D0+BODY_H(3)**2D0)
      
      Rdot = (HE0(1,2)*HE0(1,1)+HE0(2,2)*HE0(2,1)+HE0(3,2)*HE0(3,1))/ HD  

      Rdot = Rdot*AU/86400D0/1000D0                 ! Rdot Km/sec.              
   
!-------------Form " deldot "            
       
      DOT = (BEV1(1,2)*BEV1(1,1)+BEV1(2,2)*BEV1(2,1)+BEV1(3,2)*BEV1(3,1))/GD 

      deldot = DOT*AU/86400D0/1000D0                ! deldot Km/sec. 

      BODY_G = BODY_G
      DIST1 = GD                                 ! true distance
      HD = HD                                    ! True vector ray 
      Ltime = DIST1 * 8.316746395045D0           ! Light time (min.)
      Mtime = DIST1 * 8.316746395045D0 * 60D0    ! Moon Light time (sec)
      Distmoon = GD*AU/1000d0/6378.137d0         ! Moon's distance in Earth-radius


!------------- Compute the phase angle of body

      RVET = RVET0
      
      COSI = (RVET**2D0+DIST1**2D0-EARTHVEC**2D0) / (2D0*RVET*DIST1)

      ANGLE = acos(COSI) * DR2D                     ! ANGLE = Phase angle

            
      ILL = ((RVET + DIST1)**2D0 - EARTHVEC**2D0) / (4D0 * RVET * DIST1) *100D0 ! % Illuminated fraction disk   

!-----------LON0 = ELIOCENTRIC ECLIPTIC LONGITUDE ; LAT0 = ELIOCENTRIC ECLIPTIC LATITUDE  (mean J2000)
!      OBL0 = iau_OBL06(TT1,TT2)
      CALL GEPV(HSV,OBL,1,LON0,LAT0)
     
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
!-----------------------------   no rigorous method  ---------------------------------------------------------------
!      a0 = (1.396971d0 * T_0 + 0.0003086d0 * T_0**2d0) * DD2R
!      b0 = (0.013056d0 * T_0 - 0.0000092d0 * T_0**2d0) * DD2R
!      c0 = (5.12362d0 + 0.241614d0*T_0 + 0.0001122d0*T_0**2d0) *DD2R
!      HLAT = LAT0 + b0 * sin(LON0 + c0)                     ! HLAT = ELIOCENTRIC ECLIPTIC LATITUDE  (mean of the date)
!      HLON = LON0 + a0 - b0 * cos(LON0 + c0) * tan(HLAT)    ! HLON = ELIOCENTRIC ECLIPTIC LONGITUDE (mean of the date)      
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
!===================================================================================================
!--------- Compute  -TDB0- JD Time to 0.0 h TDB the day (for routine TRANSIT)


       TDB0 = INT(MJD) + DJMJD0   
       AUM = AU/1000D0
       PL = ASIN(6378.137/(DIST1 * AUM)) * DR2D  ! Lunar horizontal parallax 
       SD = (DIAM / 3600D0 ) / 2D0               ! Lunar semidiameter in arcmin
          
       
       
       CALL TRANSIT(NTARG,TDB0,GAST0,AUDAY,DELTAT,PL,SD,LAT,LONG,RNPB0,RIS,SET,TRS,  &
            AZM1,AZM2,HAL,labelR,labelS,labelT,labelM,labelN,labelY,labelX,As) 
      
              
      
!----------------------------------------------------------------------------------------------
      CALL L_TIME(ET,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST,HD0,XET,PV,HE0,THD,EPV) ! RECALL for Topocentric data
   
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

      
!--------------Form RVET1 (vector ray mean J2000 )and  Rdot.
      
      RVET1 = DSQRT(POS_H(1)**2D0+POS_H(2)**2D0+POS_H(3)**2D0)
      
      Rdot1 = (HE0(1,2)*HE0(1,1)+HE0(2,2)*HE0(2,1)+HE0(3,2)*HE0(3,1))/ HD0  

      Rdot1 = Rdot1*AU/86400D0/1000D0                 ! Rdot Km/sec.              
   
!-------------Form " deldot "            
       
      DOT = (BEV0(1,2)*BEV0(1,1)+BEV0(2,2)*BEV0(2,1)+BEV0(3,2)*BEV0(3,1))/DIST   

      deldot1 = DOT*AU/86400D0/1000D0                ! deldot Km/sec. 

!----------------------------------------------------------------------------
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun ( for apparent geocentric coordinate)

      CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
 
      CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

      CALL  iau_C2S   ( p, RA5, DE5 )       
      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
      
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent geocentric coordinate)

         

      CALL iau_SXP ( VC,BVE,CVE)    ! CVE body velocity vector in unit of C 
      
      BCV = SQRT(CVE(1)**2d0+CVE(2)**2d0+CVE(3)**2d0) !form module of vector CVE

      BM1 = SQRT(1d0-abs(BCV)**2d0)            ! for use with routine SOFA iau_AB()  

      CALL ABERR (AUDAY,u1,BVE,u0)

      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 

   
! ......... Final GCRS (RA5,DE5) of body 
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

      CALL  iau_C2S   ( u5, RA10, DE10 )       
      RA10 = iau_ANP(RA10)                 ! range angle (0,2!pi)
     

         
!=====================================================================================

!---------- Form apparent geocentric coordinate (RA6,DE6)       

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

      IF (mode == 1) then 
       CALL   iau_RXP   ( RNPB, u3, u4 )
       equinox = "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2000 EQUINOX BASED)" 
      ELSE IF (mode == 2) then 
       CALL   iau_RXP   ( RC2I, u3, u4 )
       equinox = "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2006 CIO BASED)"  
      END IF 
! ......... Call for direction cosines to spherical coordinates

      CALL   iau_C2S   ( u4, RA6, DE6 )

         RA6= iau_ANP(RA6)                       ! range angle (0,2!pi)
         RA06 = RA6
         RA66= RA6-EO
         DIST3 = DIST1                              ! Distance 
   
!-------------------------------------------------------------------------------------------------

       IF( MJD > FINAL ) THEN
        DET = DELTA
        LDUT =" value not determinated"
        LDET =" value aproximate "
        XPO = 0D0
        YPO = 0D0
       ELSE IF (MJD <= FINAL) THEN
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
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)

         RA7H = RA7 * DR2D / 15D0                ! RA7 in hour
         DE7D = DE7 * DR2D

       CALL SETDT ( DET )
       
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 0, &
                        ZD, AZ, RA8, DE8 )
     
        AZ = AZ
       ZD0 = ZD
       HOT = (LONG + GMST0 - RA7)* DR2D          ! HOT = Hour angle 
       RA8 = RA8 * 15D0* DD2R
       DE8 = DE8 * DD2R 
       
      ELSE IF (mode == 2) then

       equinox2 = "Apparent Topocentric position no refraction (IAU 2006 CIO BASED)"
       CALL   iau_RXP  ( RC2I, u5, u6 )
       CALL   iau_C2S   ( u6, RA7, DE7 )
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)

       CALL iau_ATIO13 ( RA7, DE7, UTC1, UTC2, DUT, LONG, LAT, HM, &      ! AOT = Azimut (CIO Based)
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
       CALL REFDAT (PHPA,TC,0.0,WL)
       CALL SETDT ( DET )
       
       CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 1, &
                        ZD, AZ, RA9, DE9 )
      
       RA9 = RA9 * 15D0* DD2R
       DE9 = DE9 * DD2R 
       HOB = (LONG + GMST0 - RA9)* DR2D
       ZOB = ZD 
       AZ1 = AZ
      ELSE IF (mode == 2) then
              equinox3 = "Topocentric observed position with refraction.(IAU 2006 CIO BASED)"

       CALL iau_ATIO13 ( RA7, DE7, UTC1, UTC2, DUT,LONG, LAT, HM, XP, &
                        YP,PHPA,TC,RH,WL, AOB, ZOB, HOB, DOB, ROB, J ) 
       IF (J .NE. 0 ) STOP                                             
      AZ1 = AOB *DR2D; ZOB = ZOB * DR2D; HOB = HOB *DR2D        ! AOB = Azimut (CIO Based)                     
      RA9 = ROB                                                 ! ZOB = Dist.Zenith(CIO Based)
      DE9 = DOB                                                 ! HOB = Hour angle (CIO Based)
                                                                !PHPA = pressure at the observe(hPa=mB)
      END IF                                                    !  TC = ambient temperature(deg C)
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
      ELONG = Elong * DR2D  

      RA = RA3 * DR2D          !Trasform "rad " to degrees astrometri position body                          
      
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

!------Search  RISE, SET, TRANSIT
       
       RIS = RIS / 24D0
       SET = SET / 24D0
       TRS = TRS / 24D0 
       
       CALL iau_D2TF ( 1, RIS, SIGN, IHMSF1)     
       CALL iau_D2TF ( 1, SET, SIGN, IHMSF2)
       CALL iau_D2TF ( 1, TRS, SIGN, IHMSF3)

        labelL = labelN
        labelO = " "
       
        IF ( RIS == 0D0 .AND. SET == 0D0) THEN
         commentM = labelM 
         IF(commentM == labelM) THEN
          labelN = labelO
         END IF
         commentN = labelL 
         IF(commentN == labelL) THEN
          labelM = labelO
         END IF
        END IF

        IF( RIS /= 0D0 .AND. SET /= 0D0) THEN
        commentR = labelR
        commentS = labelS
        commentT = labelT 
        commentM = labelO
        commentN = labelO
        ELSE
        commentR = labelO
        commentS = labelO
        commentT = labelO  
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
                OMEGA,OME0,AGE,JHE,CHI,JAP,APE,DLR,TLong,TLat,LTL,LTB,PA,THL,THB)
     
      MJ = 2400000.5d0
      DJ2 = JHE - MJ                     ! 2a part of date in  MJD for compute New Moon
      DJ3 = JAP - MJ                     ! 2a part of date in  MJD for compute Apogee - Perigee
      PH = ASIN(6378.137d0 / DLR)        ! APOGEE- PERIGEE Horizontal lunar parallax. (radians)
     END IF

!----------- Geocentric longitude of the Sun for compute the aberration. 
      
      CALL GEPV(SGP,OBL,1,SLO,SLA) 

      Bsun = SLA + b0 * sin(SLO + c0)
      Lsun = SLO + a0 - b0 * cos(SLO + c0) * tan(SLA)

      Lsun = iau_ANP(Lsun)

!---------- Convert Geocentric Equatorial (GMB) to Geocentric Ecliptic Coordinate 

      CALL GEPV(GMB,OBL,1,LOGE,BGE)       ! LOGE = Geocentric - Ecliptic longitude. (radians)
                                          ! BGE  = Geocentric - Ecliptic latitude. (radians)

      IF ( LOGE < 0D0) THEN
       LOGE = LOGE + 2D0 * PI
      END IF
      
      
      GELA = BGE + b0 * sin(LOGE + c0)                     ! GELA = GEOCENTRIC ECLIPTIC LATITUDE  (mean of the date,radians)
      GELO = LOGE + a0 - b0 * cos(LOGE + c0) * tan(BGE)    ! GELO = GEOCENTRIC ECLIPTIC LONGITUDE (mean of the date,radians)
      
      IF ( GELO < 0D0) THEN
       GELO = GELO + 2D0 * PI
      END IF
      
      LOGE = LOGE * DR2D                 ! LOGE = Geocentric - Ecliptic longitude J2000. (degrees)
      BGE = BGE * DR2D                   ! BGE  = Geocentric - Ecliptic latitude.J2000 (degrees)

      DLO = (-KA*COS(Lsun-GELO)+ECC*KA*COS(PHE-GELO)) / COS(GELA)        ! Difference in longitude (arcsec)

      DLA = -KA*SIN(GELA)*(SIN(Lsun-GELO)-ECC*SIN(PHE-GELO))             ! Difference in Latitude (arcsec)
      
      IF(NTARG /= 10) THEN
       GELOAP = GELO + DPSI+(DLO/DIST1) * AS2R   ! GELOAP = Geocentric Apparent Ecliptic longitude mean the date(rad).             
       GELAAP = GELA + DEPS + DLA * AS2R      ! GELAAP = Geocentric Apparent Ecliptic latitude mean the daterad).                         
      ELSE IF (NTARG == 10)THEN
       GELOAP = GELO + DPSI  
       GELAAP = GELA + DEPS
      END IF

      GELOD = GELO * DR2D                 ! GELOD = Geocentric - Ecliptic longitude mean the date. (degrees) 
      GELAD = GELA * DR2D                 ! GELAD = Geocentric - Ecliptic latitude mean the date. (degrees)
      GELOAPD = GELOAP * DR2D
      GELAAPD = GELAAP * DR2D 


      WRITE(*,*)"=========================================================================================="

      WRITE(*,*)" --------------------------------DATA OUTPUT ----------------------------------"

      WRITE(*,*)
      WRITE(*,5) longit,lS,latit,bS,height
5     FORMAT (1x,"Longitude",1x,"deg",F10.5,1x,a,3x,"Latitude",1x,"deg",F9.5,1x,a,3x,"Height",1x,F6.1,1x,"meter")       
      WRITE(*,*)
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
1     FORMAT ( 1X,A,I5,2("/",I2.2),1X,I3,2(":",I2.2),".",I6.6,"   MJD",F21.11 )

      WRITE(*,13)" UT1-UTC",DUT,LDUT
13    FORMAT(a8,"  ",F9.6," sec.",a30)
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
        WRITE ( *,52 ) IHMSF ,sign,IDMSF,DIST2 * AUM,LOGE,BGE
        ELSE IF( NTARG /= 10) THEN
        WRITE ( *, 30 ) IHMSF ,sign,IDMSF,DIST2,LOGE,BGE
      END IF

      WRITE (*,48)"=================================================================================================="
      WRITE (*,*)equinox
      WRITE (*,32) "            R.A                      Dec.          Distance          Distance                    "
      IF ( NTARG == 10 ) THEN
       WRITE (*,64) 
      ELSE IF( NTARG /= 10) THEN
       WRITE (*,65)
      END IF  
      CALL iau_A2TF ( 5 ,RA6 ,sign , IHMSF )
      CALL iau_A2AF ( 4 ,DE6 ,sign , IDMSF )
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
      WRITE (*,310)sign,IDMSF,GELOAPD ! GELOD
       CALL iau_A2AF ( 4 ,GELAAP ,sign , IDMSF )
      WRITE (*,320)sign,IDMSF,GELAAPD !GELAD
      END IF      
      WRITE (*,48)"=====================================================================================================" 
      WRITE(*,*)
      WRITE (*,32) "       ILL. %       MAGN.        Angular Diameter        Phase Angle (o)     Sun's Elong.         "
      IF ( NTARG == 10 .OR. NTARG == 11 ) THEN
        WRITE ( *,56 ) ILL,MAGN,DIAM/60D0,PANGLE,ELONG,EL
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
     CALL iau_D2DTF ( "UTC", 6, MJ, DJ2, IY, IM, ID, IHMSF4, J ) !Transform Julian day of 
       CALL iau_A2AF ( 4 ,GELO ,sign , IDMSF )
      WRITE (*,301)sign,IDMSF,GELOD
       CALL iau_A2AF ( 4 ,GELOAP ,sign , IDMSF )
      WRITE (*,311)sign,IDMSF,GELOAPD
       CALL iau_A2AF ( 4 ,GELA ,sign , IDMSF )
      WRITE (*,321)sign,IDMSF,GELAD 
      WRITE (*,*)
      WRITE (*,299),LEA(1),LEA(2),LEA(3)
299   FORMAT( " Lunar Euler angles: phi, theta, psi (rad)........... ",x,f12.10,2x,f12.10,2x,f14.8)
      WRITE (*,300) LEA(4),LEA(5),LEA(6)
300   FORMAT(" Lunar librations: dphi/dt, theta/dt, dpsi/dt (rad/day)",x,f12.10,2x,f12.10,2x,f12.10)      
      WRITE(*,*)
      WRITE (*,302)"NEW MOON", IY, IM, ID, IHMSF4 ,AGE            !New Moon in Gregorian date
      WRITE (*,*)
      WRITE (*,303)OMEGA,OME0
      WRITE (*,*)
      CALL iau_D2DTF ( "UTC", 6, MJ, DJ3, IY1, IM1, ID1, IHMSF5, J ) !Transform Julian day of 
      CALL iau_A2AF ( 4, PH, SIGN, IDMSF1 )
      WRITE (*,304)APE,IY1,IM1,ID1,IHMSF5,DLR,IDMSF1                 ! Apogee-Perigee
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
      WRITE (*,42) IHMSF1, IHMSF3, IHMSF2
      WRITE (*,*)
      WRITE (*,41) "        Azimuth                    Altitude                 Azimuth                              " 
      WRITE (*,43) AZM1,HAL,As,AZM2 
      WRITE (*,*)
      WRITE (*,47) 
      WRITE (*,44)planet
      WRITE (*,39) commentM,commentN
      WRITE (*,45)commentR 
      WRITE (*,45)commentS
      WRITE (*,46)commentT
      WRITE (*,48)commentU
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
42    FORMAT ( 6X," h ",I3,2(":",I2.2),".",I1.1,11X," h ",I3,2(":",I2.2),".",I1.1,13X," h ",I3,2(":",I2.2),".",I1.1) 
43    FORMAT ( 9X,F5.1,"°",22X,F4.1,"°",a,19X,F5.1,"°")
44    FORMAT (10x," The",1x,A10,1x,A90)
39    FORMAT (25x,A45) 
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
301   FORMAT (2x," Moon longitude mean the date  ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
311   FORMAT (2x," Moon apparent longitude       ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
321   FORMAT (2x," Moon's latitude mean the date ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
302   FORMAT (3X,A,I5,2("/",I2.2),3X,"h",1x,I3,2(":",I2.2),".",I6.6,5x,"Days after the new moon",4x,f9.6)
303   FORMAT (3x,"Mean Longit. ascend. node",1x,f10.6," degs",3x," True Longit. ascend. node ",f11.6," degs")
304   FORMAT (3x,a7,I5,2("/",I2.2),4X,"h",1x,I3,2(":",I2.2),".",I6.6,4x," Distance ",f10.3," Km = Parallax " & 
              ,1x,I3,'°',1x,I2,"'",1x,I2,"''.",I4.4)  
305   FORMAT (3x,"Position angle of the Moon's bright limb  ",f8.4," degs")
306   FORMAT (3x,"Earth's Selenographic Longitude     ",f8.4," degs","    Earth Selenographic Latitude       ",f7.4," degs")
307   FORMAT (3x,"Topocentric libration in Longitude  ",f8.4," degs","    Topocentric Libration in Latitude  ",f7.4," degs")
308   FORMAT (3x,"Position angle of axis              ",f8.4," degs")
309   FORMAT (3x,"Sun's Selenographic Longitude       ",f8.4," degs","    Sun's Selenographic Latitude       ",f7.4," degs")      
!310   FORMAT (11x,"Longitude mean the date      ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
!320   FORMAT (11x,"Latitude mean the date       ",2x,a,I3,"°",1x,I2,"'",1x,I2,"''",".",I4.4,4x,"degs",1x,f10.6 )
310   FORMAT (5x,"Apparent Ecliptic Longit.mean the date",2x,a,I3,"°",x,I2,"'",x,I2,"''",".",I4.4,4x,"degs",x,f10.6 )
320   FORMAT (5x,"Apparent Ecliptic Latit. mean the date",2x,a,I3,"°",x,I2,"'",x,I2,"''",".",I4.4,4x,"degs",x,f10.6 )

 GOTO 3000 
      
!=====================================================================================================================
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

!-------- Compute TT and UT1 
      
      IF (modetime == 1) THEN
        CALL iau_CAL2JD ( IY, MO, ID, DJMJD0, DATE, J )            ! input UT1 (IY,MO,ID)
        UT1 = DATE + TIME
        GOTO 1750
      END IF 

      IF (modetime == 2) GOTO 1700 
      
!--------TT (MJD).
1700  CALL iau_CAL2JD ( IY, MO, ID, DJMJD0, DATE, J )              ! input TT (IY,MO,ID)
      TT = DATE + TIME 
1750  continue                            

!----------Call the interpolation routine for date above MJD 41684 

!  input  MJD   : modified julian date 
!         
!  output XP0   : polar motion (arcsec)  UT11    
!         YP0   : polar motion (arcsec)
!         DUT   : UT1 - UTC since 1962.00 (sec.)
!         dpsi  : celestial pole offset dpsi / UAI 2000 en arcsec
!         deps  : --------------------- deps ------------------
!         DELTA : DeltaT (sec.)
!    
!  Coded by A. Nicola - October 2011
      
      CALL iers_1900 (MJD,XPO,YPO,DUT,dpsi,deps,DELTA)
      DELTAT = DELTA / 86400D0
      
!--------- tranform "dpsi ,deps" en mas

      off_dpsi = dpsi * 1000D0
      off_deps = deps * 1000D0
     
!---------  Transformation of the celestial pole offsets (dpsi,deps)_2000 into  
!           the celestial pole offsets (dX,dY)_2000  by using SOFA 
!           matrix transformation recommanded by UAI 2000 (Wallace, 2006).

!  input  MJD   : modified julian date 
!         dpsi  : celestial pole offset dpsi / UAI 2000 en mas
!         deps  : --------------------- deps ------------------
!
!  output dX    : celestial pole offset dX / UAI 2000 en mas
!         dY    : --------------------- dY   ---------------
!
! 
!  Coded by Ch. Bizouard - May 2009
   
      CALL DXDY2000_DPSIDEPS1980 (MJD,dpsi,deps,DX0,DY0)
      
      
      DUT = DUT                           !  DUT = UT1-TAI sec.
!      IF ( MJD < 41687D0) DUT = 0D0       
      IF ( MJD < 35473D0) DUT = 0D0       !  before MJD 35473 UT1-TAI (DUT) = 0D0
      TUT = TIME                          !         

      XP = XPO * AS2R                     !  Polar motion (arcsec---> radians)
      YP = YPO * AS2R                     ! 
      
      DX00 = DX0  * AS2R/1000D0           !  CIP correction ( mas-----> radians)
      DY00 = DY0 * AS2R/1000D0            !
    
      IF (modetime == 2) THEN 
       UT1 = TT - DELTAT                   ! UT1 in MJD
       UT11 = UT1 + DJMJD0                 ! UT1 in JD  
       UT = UT11 - INT(UT11)
       TT1 = INT(TT + DJMJD0)
       TT2 = (TT + DJMJD0) - TT1 
       GOTO 1900 
      ELSE IF (modetime == 1) THEN
       TT = UT1 + DELTAT
       TT1 = INT(TT + DJMJD0)
       TT2 = (TT + DJMJD0) - TT1
       UT11 = UT1 + DJMJD0                 ! UT1 in JD  
       UT = UT11 - INT(UT11) 
       GOTO 1900 
      END IF
           
!----------- TT --> TDB -> TCB.
1900  DTR = iau_DTDB ( TT1, TT2, UT, LONG, U, V )
      CALL iau_TTTDB ( TT1, TT2, DTR, TDB1, TDB2, J )
      IF ( J.NE.0 ) STOP
      
       
!----------- TT---> TCG
      CALL iau_TTTCG ( TT1, TT2, TCG1, TCG2, J )
      IF ( J.NE.0 ) STOP
      
      CALL iau_TDBTCB ( TDB1, TDB2, TCB1, TCB2, J )
      IF ( J.NE.0 ) STOP
      
          
      TDB = TDB1+TDB2         !
      TCG = TCG1+TCG2         ! 
      TCB = TCB1+TCB2         !

      UT_0 = INT(MJD) + DJMJD0        ! JD Time to 0.0 h UT1.      
      TT_0 = UT_0 + DELTAT          ! JD Time to 0.0 h TT the day.
!--------- Greenwich apparent sidereal time (IAU 2006).

      GAST0 = iau_ANP ( iau_GMST06 ( UT_0, 0D0, TT_0, 0D0 ) &
          + iau_EE06A ( TT_0, 0D0, EPSA, DPSI ) )
     
      
!==============================================================================================
      

      GAST = iau_ANP ( iau_GMST06 ( DJMJD0, UT1, TT1, TT2 )&
          + iau_EE06A ( TT1, TT2, EPSA, DPSI ) )
      
      CALL GEOCPV ( LONG,LAT,HM,GAST,POS,VEL,J)     ! Geocenter position and velocity observer 
      IF ( J.NE.0 ) STOP
   
!------- Earth  rotation  angle. 
       ERA  =  iau_ERA00  ( DJMJD0+DATE,TUT)  !(   UTC1,  TUT    ) 
      
! ========================
! IAU 2000A, equinox based
! ========================
!---------Nutation, IAU 2000A.
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
!-------- Corrected nutation.
      DPSI = DP00 + DDP00
      DEPS = DE00 + DDE00
  
!-------- Build the rotation matrix.
      CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )
!-------- Combine the matrices: PN = N x P.
      CALL iau_RXR ( RN, RPB, RNPB )

!=========================================================================================

!---------Nutation, IAU 2000A. to 0.0h TT (for compute body Transit )
      CALL iau_NUT00A ( TT_0, 0D0, DP00, DE00 )
      
!-------- Precession-nutation quantities, IAU 2000. to 0.0h
      CALL iau_PN00 ( TT_0, 0D0, DP00, DE00, &
           EPSA, RB, RP, RPB, RN, RNPB0 )
 
!==========================================================================================
     
      ET = TDB 
      T0 =2451545.0D0
      
      CALL NUTATION (-ET,DJMJD0,TT,POS,POS2)
      CALL NUTATION (-ET,DJMJD0,TT,VEL,VEL2)

!     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   GP )

!     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL PRECES ( ET, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   GV )

     
!---------------------------------------------------------------------------------------------------
      ISTEP=1
      NSTEP=1
            
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

! ......... Call for direction cosines to spherical coordinates
 
      CALL  iau_C2S ( BODY_G, RA3, DE3 )       
    
      RA3 = iau_ANP(RA3)                  ! range angle (0,2!pi)
   
!--------------Form RVET (vector ray mean J2000 )and  Rdot.
      
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
       DIAM = 1919.29D0 / DIST1
       ILL = 1.0D0 *100D0
      END IF

!-----------Form  topocentric  position RA8,DE8
       
      CALL L_TIME(ET,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST,THD,XET,PV,HE0,THD,EPV)
      POS_T = POS_T        
      DIST4 = DIST
      CALL  iau_C2S ( POS_T, RA8, DE8 )       
      RA8 = iau_ANP(RA8)                  ! range angle (0,2!pi)
             
           
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
      
   
! ......... Final GCRS (RA5,DE5) of body 
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

      CALL   iau_C2S   ( u4, RA6, DE6 )
         RA6= iau_ANP(RA6)                       ! range angle (0,2!pi)
 
      DIST3 = DIST1                              ! Distance 
!-------------------------------------------------------------------------------------------------
 
!---------- Form apparent topocentric coordinate RA7,DE7)

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006 CIO based using X,Y series" 

!      IF (mode == 1) then
       CALL   iau_RXP  ( RNPB, u5, u6 )
       equinox2 = "Apparent Topocentric position  ' no refraction'.(IAU 2000 EQUINOX BASED)"
!       goto 5000
!      ELSE IF (mode == 2) then
!       CALL   iau_RXP   ( RC2I, u5, u6 )
!       equinox2 = "Apparent Topocentric position  ' no refraction .(IAU 2006 CIO BASED)"  
!       goto 6000 
!      END IF 


! ......... Call for direction cosines to spherical coordinates

5000    CALL   iau_C2S   ( u6, RA7, DE7 )
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
!----------------------------------------------------------------

!--------- Compute  -TDB0- JD Time to 0.0 h TDB the day (for routine TRANSIT)

        TDB0 = INT(MJD) + DJMJD0  
       
       PL = ASIN(6378.137/(DIST3 * AUM)) * DR2D  ! Lunar horizontal parallax 
       SD = (DIAM / 3600D0 ) / 2D0               ! Lunar semidiameter in arcmin
          
      
       CALL TRANSIT(NTARG,TDB0,GAST0,AUDAY,DELTAT,PL,SD,LAT,LONG,RNPB0,RIS,SET,TRS,  &
            AZM1,AZM2,HAL,labelR,labelS,labelT,labelM,labelN,labelO,labelX,As)  
     
       
       RIS = RIS / 24D0
       SET = SET / 24D0
       TRS = TRS / 24D0 
       
       CALL iau_D2TF ( 1, RIS, SIGN, IHMSF1)     
       CALL iau_D2TF ( 1, SET, SIGN, IHMSF2)
       CALL iau_D2TF ( 1, TRS, SIGN, IHMSF3)
       
       IF ( RIS == 0D0 .AND. SET == 0D0) THEN
        commentR = labelM
        commentS = labelN
        commentT = labelO
       ELSE
        commentR = labelR
        commentS = labelS
        commentT = labelT   
       END IF
       
       IF (NTARG == 6 ) THEN
        commentU = "In calculating the Saturn's magnitude is also taking into account of the rings."
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

      WRITE(*,117) " DELTA_T",DELTA
117   FORMAT(a8,"  ",F7.4," sec.")
      WRITE(*,1117) DUT
1117  FORMAT(x,"UT1-TAI value (1956-1972) ",f11.7," sec")
      WRITE(*,150)" TT - TDB at geocenter: ",PVV(1),"sec",PVV(2),"sec/day"
      WRITE(*,*)
      WRITE(*,118) GAST*DR2D
118   FORMAT('     Greenwich Sideral Time (GAST) =',f19.13,' deg')
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
      WRITE (*,138) off_dpsi
138   FORMAT('     Celestial pole offset dPsi = ',f9.4,' mas')   
      WRITE (*,139) off_deps
139   FORMAT('     Celestial pole offset dEps = ',f9.4,' mas') 
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
      WRITE (*,142) IHMSF1, IHMSF3, IHMSF2
      WRITE (*,*)
      WRITE (*,141) "        Azimuth                    Altitude                 Azimuth                              " 
      WRITE (*,143) AZM1,HAL,As,AZM2 
      WRITE (*,*)
      WRITE (*,147) 
      WRITE (*,144)planet,commentR
      WRITE (*,145)commentS
      WRITE (*,146)commentT
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
142    FORMAT ( 6X," h ",I3,2(":",I2.2),".",I1.1,11X," h ",I3,2(":",I2.2),".",I1.1,13X," h ",I3,2(":",I2.2),".",I1.1) 
143    FORMAT ( 9X,F5.1,"°",22X,F4.1,"°",a,19X,F5.1,"°")
144    FORMAT (10x," The",1x,A10,1x,A90)
145    FORMAT (25x,A40)
146    FORMAT (25x,A40) 
147    FORMAT (5x," *** NOTES *** ")
160    FORMAT (10x,"h , m ,  s",14x,"o , ' , ''",12x,"degs",12x,"degs",11x,"degs")
161    FORMAT (14x,"R.A",21x,"Dec.",13x,"Azimut",6x,"(*)Elevation",2x,"(*)Hour Angle")
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
      
