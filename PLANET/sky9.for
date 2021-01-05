      subroutine iers_1962(MJD_62,DATE,EOP_final,UTC,RJD,X,Y,UT1,DDT,
     . DDX,DDY,N,rjd_int)
                 

C
C This subroutine detects POLAR MOTION (XP, YP, UT1-UTC, dX and DY from
C IERS database "iers2000A all" (1973-2012) and passes the values
C for the subsequent interpolation and computation of diurnal gravitational  
C and tidal ocanic effects.                
C                    coded by Aldo Nicola (revised Set.2012)  

C choose between A or B bulletin IERS final_all (IAU2000) Standard EOP Data
C Bulletin A , between lines 38 - 46


      implicit none

      DOUBLE PRECISION  XP,YP,DUT1,X1,Y1,RJDINT,UTC,RJD(7)        
      DOUBLE PRECISION YEAR,X(7),Y(7),UT1(7),rjd_int
      DOUBLE PRECISION JD,DATE,Xx,Yy,DUT,x_int,y_int,ut1_int
      DOUBLE PRECISION dX,dY,DDX(7),DDY(7),DEX,DEY,lo_d
      DOUBLE PRECISION Ra,Rb,Rc,Rd,Re,Rf,DDT(7),h,k,l,m,o,p
      INTEGER A,I,A1,N,J,FINAL,EOP_final,MJD                        
      INTEGER aa,bb,cc,ios,loop,MJD_62,date1,date2
      character(len=1):: B,B1,dd,ee,ff  
      INTEGER E(7)
      INTEGER F(7)
      INTEGER G(7)

      A1 = int(4)         
      A=int(DATE - A1)
C      FINAL = FINAL + 3
C      Mtoday = FINAL -350 !ricalcola il giorno per il cambio database (EOP14 C04 == B )verso A.
      
C    Open database  "iers1962_now" 
      open (unit=15, file="iers1962_now.txt", status="old",         
     .  action="read",form="formatted",position="rewind")
      
      do loop = 37665,EOP_final          
      
C------ Use iers1962_now (IAU2000) EOP14 C04 Data- 
C      read (unit=14,fmt="(3(i4),i7,2(f11.6),2(f12.7),2(f11.6))",
      read (unit=15,fmt="(3(I4),I7,2(F11.6),2(F12.7),2(F11.6),2(F11.6),
     .  2(F11.7),2(F12.6))",iostat=ios)aa,bb,cc,MJD,Xx,Yy,DUT,lo_d,
     .  dX,dY,h,k,l,m,o,p
    
      if (ios==0) then 
       continue
      else
       exit
      end if
      
      do I = 1,7              !  reads  7  values of each parameter           
       if(loop == A + I) then 
        
        E(I)= aa
        F(I)= bb  
        G(I)= cc
        RJD(I)= MJD
        X(I)= Xx 
        Y(I)= Yy
        UT1(I)= DUT
        DDT(I) = lo_d *1000D0
        DDX(I)= dX *1000D0
        DDY(I)= dY *1000D0 

       end if
       end do
      end do
      
      rewind (unit=15)
C    " N ,NUMBER OF PAIRS FOR INTERPOLATION ( X, Y)"

      N = 7
      
C    "INTERPOLATION FOR UTC (MJD)
      rjd_int = UTC

C    Calling Sub INTERP,LAGINT,PMUT1_OCEANS,PM_GRAVI
      
      END
C----------------------------------------------------------------------------------

      subroutine ast_search1(NTARG,NAMES,num,body,H,G,Mean,Peri,Node, 
     . Incl,Ecc,a,code)

      implicit none

      DOUBLE PRECISION  Mean,Peri,Node,Incl,Ecc,No,a 
      DOUBLE PRECISION  Sa,E,I,N,P,F 
      REAL  H, G, C, D
      INTEGER NTARG, J, ios 
      CHARACTER *17 body,B,NAMES
      CHARACTER *5 num,nu
      CHARACTER *5 code,co
      
C     Open asteroids database  NX=1 "20000ast.dat for orbital elements" 
       
      open (unit=14, file="20000ast.dat", status="old",         
     .  action="read",form="formatted",position="rewind") 
       
      do J = 1,20000
 
      read (unit=14,fmt="(a5,3x,f5.2,2x,f4.2,x,a5,x,f9.5,2x,f9.5,2x,  
     .  f9.5,3x,f8.5,2x,f9.7,2x,f10.8,3x,f9.7,72x,a17)",iostat=ios) 
     .  nu,C,D,co,F,P,No,I,E,N,Sa,B

      if (ios.eq.0) then                         
       continue
      else
       exit
      end if
       
      if (NTARG .eq. J)then         ! NTARG   = data number of asteroid
        num = nu
          H = C                      ! H and G = parameters for for     
          G = D                      !           calculating the magnitude
       code = co
       Mean = F                      ! Mean    = Mean anomaly
       Peri = P                      ! Peri    = Perielio of asteroid
       Node = No                     ! Node    = Ascending Node
       Incl = I                      ! Incl    = Inclination
        Ecc = E                      ! Ecc     = Eccentricity  
          a = Sa                     ! a       = Semimajor axis  
       body = B                      ! body    = name of asteroid 
       goto 10 
      end if 
      end do
             
10    continue      
      rewind (unit=14)
      end
      
*-----------------------------------------------------------------------

      subroutine ast_search2(NTARG,NAMES,num,body,H,G,Mean,Peri,Node, 
     . Incl,Ecc,a,code)

      implicit none

      DOUBLE PRECISION  Mean,Peri,Node,Incl,Ecc,No,a 
      DOUBLE PRECISION  Sa,E,I,N,P,F 
      REAL  H, G, C, D
      INTEGER NTARG, J, ios
      CHARACTER *17 body,B,NAMES
      CHARACTER *5 num,nu
      CHARACTER *5 code,co
      
C     Open asteroids database  NX=2 "NEA203.txt for orbital elements" 
       
      open (unit=14, file="NEA203.txt", status="old",        
     .  action="read",form="formatted",position="rewind") 
       
       do J = 1,203
 
      read (unit=14,fmt="(a5,3x,f5.2,2x,f4.2,x,a5,x,f9.5,2x,f9.5,2x,  
     .  f9.5,3x,f8.5,2x,f9.7,2x,f10.8,3x,f9.7,72x,a17)",iostat=ios) 
     .  nu,C,D,co,F,P,No,I,E,N,Sa,B

      if (ios.eq.0) then                         
       continue
      else
       exit
      end if
       
      if (NTARG .eq. J)then         ! NTARG   = data number of asteroid
        num = nu
          H = C                      ! H and G = parameters for for     
          G = D                      !           calculating the magnitude
       code = co
       Mean = F                      ! Mean    = Mean anomaly
       Peri = P                      ! Peri    = Perielio of asteroid
       Node = No                     ! Node    = Ascending Node
       Incl = I                      ! Incl    = Inclination
        Ecc = E                      ! Ecc     = Eccentricity  
          a = Sa                     ! a       = Semimajor axis  
       body = B                      ! body    = name of asteroid 
       goto 10 
      end if 
      end do
             
10    continue      
      rewind (unit=14)
      end
      


*-----------------------------------------------------------------------

      subroutine iers_calc(UTC,DATE,FINAL_,EOP_final,XP,YP,DUT1,LOD,
     .  dX,dY,ut1_notid)
                 

C
C This subroutine detects POLAR MOTION (XP, YP, UT1-UTC, dX and DY from
C IERS database "iers2000A all" (1973-2012) and passes the values
C for the subsequent interpolation and computation of diurnal gravitational  
C and tidal ocanic effects.                
C                    coded by Aldo Nicola (revised Set.2012)  

C choose between A or B bulletin IERS final_all (IAU2000) Standard EOP Data
C Bulletin A , between lines 38 - 46


      implicit none

      DOUBLE PRECISION  XP,YP,DUT1,X1,Y1,RJDINT,UTC,MJD,RJD(7)        
      DOUBLE PRECISION YEAR,X(7),Y(7),UT1(7),rjd_int,LD(7)
      DOUBLE PRECISION JD,DATE,Xx,Yy,DUT,x_int,y_int,ut1_int
      DOUBLE PRECISION dX,dY,DDX(7),DDY(7),DEX,DEY,lo_d,ut1_out
      DOUBLE PRECISION Ra,Rb,Rc,Rd,Re,Rf,DDT(7),LOD,U,ut1_notid
      
      INTEGER :: A,I,A1,N,J,FINAL_,EOP_final,Mtoday ,MJD62                       
      INTEGER :: aa,bb,cc,ios,loop,MJD_62,FIN
      character(len=1):: B,B1,dd,ee,ff  
      INTEGER E(7)
      INTEGER F(7)
      INTEGER G(7)

      
      A1 = int(3)         
      A=int(DATE - A1)
      FIN = FINAL_ + 3

      MJD62 = 37665          !   initial time of EOP14 C04 (MJD = 1962/01/01) 
      
      IF (UTC .LE. (EOP_final-4) ) THEN
       CALL iers_1962(MJD_62,DATE,EOP_final,UTC,RJD,X,Y,UT1,DDT,DDX,DDY,
     .      N,rjd_int)

        GOTO 100
      ELSE 
        GOTO 80
      END IF

C     Open database  "iers_finals_2000A" 
80    open (unit=13, file="finals-iers.txt", status="old",         
     .  action="read",form="formatted",position="rewind")
      
      do loop = 48622,FIN !FINAL_  ! INI = start MJD (final data'IAU2000') 
                   
C------ Use IERS final_all (IAU2000) Standard EOP Data- from Bulletin B                                     
C      read (unit=13,fmt="(i6,x,f8.2,65x,f7.4,48x,f9.6,x,f9.6,x,f10.7,  
C     .  3x,f7.3,4x,f6.3)",iostat=ios)aa,MJD,lo_d,Xx,Yy,DUT,dX,dY
             
C------ Use IERS final_all (IAU2000) Standard EOP Data- from Bulletin A
      read (unit=13,fmt="(i2,i2,i2,x,f8.2,x,a1,x,f9.6,f9.6,x,f9.6,f9.6,  
     .  2x,a1,f10.7,a11,f7.4,a14,f6.3,a13,f6.3)",iostat=ios)aa,bb,cc,
     .  MJD,B,Xx,X1,Yy,Y1,B1,DUT,ff,lo_d,dd,dX,ee,dY     
             
      if (ios==0) then 
       continue
      else
       exit
      end if
      
      do I = 1,7              !  reads  7  values of each parameter           
       if(loop == A + I) then 
       
        E(I)= aa
        F(I)= bb  
        G(I)= cc
        RJD(I)= MJD
        X(I)= Xx 
        Y(I)= Yy
        UT1(I)= DUT
        DDT(I) = lo_d
        DDX(I)= dX
        DDY(I)= dY  
       end if
      
      end do
      end do
      
      rewind (unit=13)

100   continue      
C    " N ,NUMBER OF PAIRS FOR INTERPOLATION ( X, Y)"
      N = 7
C    "INTERPOLATION FOR UTC (MJD)
      rjd_int = UTC
C    Calling Sub INTERP,LAGINT,PMUT1_OCEANS,PM_GRAVI
C      print*,"RJD(7) ",RJD,UT1  
      
      CALL INTERP(RJD,X,Y,UT1,DDX,DDY,N,rjd_int,x_int,y_int,
     .     ut1_int,DEX,DEY,ut1_notid)
     
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                                                                           !
C  Recoding values of Polar motion (XP , YP), (UT1-UTC), LOD and            ! C(dX ,dY)  !                                                                !
C                                                                           ! 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C      print*,"YPO",y_int
      CALL DE_DUT(rjd_int,UT1,ut1_int,ut1_out)      
      XP = x_int
      YP = y_int
      DUT1 = ut1_out 
      dX = DEX
      dY = DEY
      LOD = DDT(3)
      ut1_notid = ut1_notid
      
!      close (13)            ! added Mar.30,2013  (Giuseppe)
      end 
C-----------------------------------------------------------------------------------------------------------
      subroutine DE_DUT (rjd_int,UT1,ut1_int,ut1_out)
     
C---- Questa routine serve per eliminare l'ambiguità che si crea
C     quando il Leap-second trova UT1-UTC negativo.

      implicit none

      INTEGER NDAT
      PARAMETER ( NDAT = 182 )
      DOUBLE PRECISION UT1(7),ut1_int,rjd_int,ut1_out
      REAL MD(NDAT)
      REAL DAY
      INTEGER N,flag, N1
      

       DATA MD /
     . 41680, 41681, 41682, 41683, 41684, 41685, 41686,
     . 42045, 42046, 42047, 42048, 42049, 42050, 42051,
     . 42410, 42411, 42412, 42413, 42414, 42415, 42416,
     . 42775, 42776, 42777, 42778, 42779, 42780, 42781,
     . 43141, 43142, 43143, 43144, 43145, 43146, 43147,
     . 43506, 43507, 43508, 43509, 43510, 43511, 43512,
     . 43871, 43872, 43873, 43874, 43875, 43876, 43877,
     . 44236, 44237, 44238, 44239, 44240, 44241, 44242,
     . 44783, 44784, 44785, 44786, 44787, 44788, 44789,
     . 45148, 45149, 45150, 45151, 45152, 45153, 45154,
     . 45513, 45514, 45515, 45516, 45517, 45518, 45519,
     . 46244, 46245, 46246, 46247, 46248, 46249, 46250,
     . 47158, 47159, 47160, 47161, 47162, 47163, 47164,
     . 47889, 47890, 47891, 47892, 47893, 47894, 47895,
     . 48254, 48255, 48256, 48257, 48258, 48259, 48260,
     . 48801, 48802, 48803, 48804, 48805, 48806, 48807,
     . 49166, 49167, 49168, 49169, 49170, 49171, 49172,
     . 49531, 49532, 49533, 49534, 49535, 49536, 49537,
     . 50080, 50081, 50082, 50083, 50084, 50085, 50086,
     . 50627, 50628, 50629, 50630, 50631, 50632, 50633,
     . 51176, 51177, 51178, 51179, 51180, 51181, 51182,
     . 53733, 53734, 53735, 53736, 53737, 53738, 53739,
     . 54829, 54830, 54831, 54832, 54833, 54834, 54835,
     . 56106, 56107, 56108, 56109, 56110, 56111, 56112,
     . 57201, 57202, 57203, 57204, 57205, 57206, 57207, 
     . 57751, 57752, 57753, 57754, 57755, 57756, 57757 /

      DAY = int(rjd_int)
      flag = 0
      
       DO N=1,NDAT 
        DO N1 = 1,NDAT
         IF (MD(N1) == DAY) THEN
         ut1_out = UT1(4)
         flag=1
         goto 10
         ELSE     
         ut1_out = ut1_int
         END IF
        END DO
       END DO
10    continue      
      
      END
C------------------------------------------------------------------------------------------------------


      SUBROUTINE INTERP(RJD,X,Y,UT1,DDX,DDY,N,rjd_int,  
     .  x_int,y_int,ut1_int,DEX,DEY,ut1_notid)
C
C     This subroutine takes a series of x, y, and UT1-UTC values
C     and interpolates them to an epoch of choice. This routine
C     assumes that the values of x and y are in seconds of
C    arc and that UT1-UTC is in seconds of time. At least
C     one point before and one point after the epoch of the
C     interpolation point are necessary in order for the
C     interpolation scheme to work. 
C
C     parameters are :
C     RJD     - array of the epochs of data (given in mjd)
C     X       - array of x polar motion (arcsec)
C     Y       - array of y polar motion (arcsec)
C     UT1     - array of UT1-UTC (sec)
C     n       - number of points in arrays
C     rjd_int - epoch for the interpolated value
C     x_int   - interpolated value of x
C     y_int   - interpolated value of y
C     ut1_int - interpolated value of ut1-utc
C
C     CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
C                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
C                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects)
C
C      coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
C                                          Corrected : September 2007   
  
      implicit none
      integer I,N
      DOUBLE PRECISION RJD(N),X(N), Y(N), UT1(N),DDX(N),DDY(N),DDT(N),
     . rjd_int, x_int, y_int, ut1_int, DEX, DEY,DET,lod,
     . cor_x, cor_y, cor_ut1, cor_lod,ut1_notid
      
       
      CALL LAGINT (RJD,X,n,rjd_int,x_int)
  
      CALL LAGINT (RJD,Y,n,rjd_int,y_int)
      
      CALL LAGINT (RJD,UT1,n,rjd_int,ut1_int)

      CALL LAGINT (RJD,DDX,n,rjd_int,DEX)

      CALL LAGINT (RJD,DDY,n,rjd_int,DEY)

C      CALL LAGINT (RJD,DDT,n,rjd_int,DET)
       ut1_notid =  ut1_int
C --------------
C Oceanic effect      
C --------------
   
      CALL PMUT1_OCEANS (rjd_int,cor_x,cor_y,cor_ut1,cor_lod)

      x_int = x_int + cor_x
      y_int = y_int + cor_y
      ut1_int = ut1_int + cor_ut1
C      lod = lod + cor_lod
      
C Lunisolar effect 
      CALL PM_GRAVI (rjd_int,cor_x,cor_y)
      
      x_int   = x_int + cor_x
      y_int   = y_int + cor_y
     
      RETURN

      END
C
C----------------------------------------------------------------
C----------------------------------------------------------------
C



      SUBROUTINE LAGINT (X,Y,n,xint,yout)
C 
C     This subroutine performs lagrangian interpolation
C     within a set of (X,Y) pairs to give the y
C     value corresponding to xint. This program uses a
C     window of 4 data points to perform the interpolation.
C     if the window size needs to be changed, this can be
C     done by changing the indices in the do loops for
C     variables m and j.
C
C     PARAMETERS ARE :
C     X     - array of values of the independent variable
C     Y     - array of function values corresponding to x
C     n     - number of points
C     xint  - the x-value for which estimate of y is desired
C     yout  - the y value returned to caller

      implicit none
 
      DOUBLE PRECISION X(n),Y(n),xint,yout,term
      INTEGER i,j,k,m,n

      yout = 0.d0
      do  i = 1,n-1
        if ( xint .ge. X(i) .and. xint .lt. X(i+1) ) k = i
      enddo
    
      if ( k .lt. 2 ) k = 2
      if ( k .gt. n-2 ) k = n-2
   
      do m = k-2,k+4
        term = y(m)
        do  j = k-2,k+4
          if ( m .ne. j ) then
            term = term * (xint - X(j))/(X(m) - X(j))
          end if
        enddo 
        yout = yout + term
      enddo
   
      return
      end
 
C----------------------------------------------------------------
      SUBROUTINE PMUT1_OCEANS (rjd,cor_x,cor_y,cor_ut1,cor_lod)
C
C    This subroutine provides, in time domain, the diurnal/subdiurnal
C    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
C    listed in the program above, have been extracted from the procedure   
C    ortho_eop.f coed by Eanes in 1997.
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
C     cor_lod  - tidal correction in length of day (sec. of time)
C
C     coded by Ch. Bizouard (2002), initially coded by McCarthy and 
C     D.Gambis(1997) for the 8 prominent tidal waves.  
      
      IMPLICIT NONE
      
      INTEGER nlines
      PARAMETER(nlines=71)
      DOUBLE PRECISION ARG(6),    ! Array of the tidal arguments   
     .                 DARG(6)    ! Array of their time derivative 
      
      REAL*4 XCOS(nlines),XSIN(nlines),
     .YCOS(nlines),YSIN(nlines),UTCOS(nlines),UTSIN(nlines)
     
      DOUBLE PRECISION t,ag,dag,rjd,halfpi,secrad,
     .       cor_x,cor_y,cor_ut1,cor_lod
      INTEGER NARG(nlines,6),i,j
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),UTSIN(j),UTCOS(j),j=1,nlines)/
     &1,-1, 0,-2,-2,-2,  -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078,
     &1,-2, 0,-2, 0,-1,   0.06,   0.64,  -0.64,   0.06,  0.195, -0.059,
     &1,-2, 0,-2, 0,-2,   0.30,   3.42,  -3.42,   0.30,  1.034, -0.314,
     &1, 0, 0,-2,-2,-1,   0.08,   0.78,  -0.78,   0.08,  0.224, -0.073,
     &1, 0, 0,-2,-2,-2,   0.46,   4.15,  -4.15,   0.45,  1.187, -0.387,
     &1,-1, 0,-2, 0,-1,   1.19,   4.96,  -4.96,   1.19,  0.966, -0.474,
     &1,-1, 0,-2, 0,-2,   6.24,  26.31, -26.31,   6.23,  5.118, -2.499,
     &1, 1, 0,-2,-2,-1,   0.24,   0.94,  -0.94,   0.24,  0.172, -0.090,
     &1, 1, 0,-2,-2,-2,   1.28,   4.99,  -4.99,   1.28,  0.911, -0.475,
     &1, 0, 0,-2, 0, 0,  -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070,
     &1, 0, 0,-2, 0,-1,   9.22,  25.06, -25.06,   9.22,  3.025, -2.280,
     &1, 0, 0,-2, 0,-2,  48.82, 132.91,-132.90,  48.82, 16.020,-12.069,
     &1,-2, 0, 0, 0, 0,  -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078,
     &1, 0, 0, 0,-2, 0,  -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154,
     &1,-1, 0,-2, 2,-2,  -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074,
     &1, 1, 0,-2, 0,-1,  -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050,
     &1, 1, 0,-2, 0,-2,  -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271,
     &1,-1, 0, 0, 0, 0,  -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751,
     &1,-1, 0, 0, 0,-1,  -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151,
     &1, 1, 0, 0,-2, 0,  -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137,
     &1, 0,-1,-2, 2,-2,   1.54,   3.03,  -3.03,   1.54,  0.315, -0.189,
     &1, 0, 0,-2, 2,-1,  -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035,
     &1, 0, 0,-2, 2,-2,  26.13,  51.25, -51.25,  26.13,  5.512, -3.095,
     &1, 0, 1,-2, 2,-2,  -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025,
     &1, 0,-1, 0, 0, 0,  -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070,
     &1, 0, 0, 0, 0, 1,   1.54,   3.00,  -3.00,   1.54,  0.348, -0.171,
     &1, 0, 0, 0, 0, 0, -77.48,-151.74, 151.74, -77.48,-17.620,  8.548,
     &1, 0, 0, 0, 0,-1, -10.52, -20.56,  20.56, -10.52, -2.392,  1.159,
     &1, 0, 0, 0, 0,-2,   0.23,   0.44,  -0.44,   0.23,  0.052, -0.025,
     &1, 0, 1, 0, 0, 0,  -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065,
     &1, 0, 0, 2,-2, 2,  -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111,
     &1,-1, 0, 0, 2, 0,  -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043,
     &1, 1, 0, 0, 0, 0,  -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187,
     &1, 1, 0, 0, 0,-1,  -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037,
     &1, 0, 0, 0, 2, 0,  -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005,
     &1, 2, 0, 0, 0, 0,  -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005,
     &1, 0, 0, 2, 0, 2,  -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037,
     &1, 0, 0, 2, 0, 1,  -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023,
     &1, 0, 0, 2, 0, 0,  -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005,
     &1, 1, 0, 2, 0, 2,  -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024,
     &1, 1, 0, 2, 0, 1,  -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015,
     &2,-3, 0,-2, 0,-2,  -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011,
     &2,-1, 0,-2,-2,-2,  -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032,
     &2,-2, 0,-2, 0,-2,  -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177,
     &2, 0, 0,-2,-2,-2,  -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222,
     &2, 0, 1,-2,-2,-2,  -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015,
     &2,-1,-1,-2, 0,-2,   0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013,
     &2,-1, 0,-2, 0,-1,   2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058,
     &2,-1, 0,-2, 0,-2, -56.87, -12.93,  11.15,  32.88, -3.795, -1.556,
     &2,-1, 1,-2, 0,-2,  -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015,
     &2, 1, 0,-2,-2,-2, -11.01,  -2.40,   1.89,   6.41, -0.698, -0.298,
     &2, 1, 1,-2,-2,-2,  -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014,
     &2,-2, 0,-2, 2,-2,   0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022,
     &2, 0,-1,-2, 0,-2,   1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025,
     &2, 0, 0,-2, 0,-1,  12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266,
     &2, 0, 0,-2, 0,-2,-330.15, -26.96,  37.58, 195.92,-16.195, -7.140,
     &2, 0, 1,-2, 0,-2,  -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021,
     &2,-1, 0,-2, 2,-2,   2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034,
     &2, 1, 0,-2, 0,-2,   9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117,
     &2,-1, 0, 0, 0, 0,  -2.35,   0.37,   0.47,   1.41, -0.106, -0.029,
     &2,-1, 0, 0, 0,-1,  -1.04,   0.17,   0.21,   0.62, -0.047, -0.013,
     &2, 0,-1,-2, 2,-2,  -8.51,   3.50,   3.29,   5.11, -0.437, -0.019,
     &2, 0, 0,-2, 2,-2,-144.13,  63.56,  59.23,  86.56, -7.547, -0.159,
     &2, 0, 1,-2, 2,-2,   1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000,
     &2, 0, 0, 0, 0, 1,   0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001,
     &2, 0, 0, 0, 0, 0, -38.48,  19.14,  17.72,  23.11, -2.104,  0.041,
     &2, 0, 0, 0, 0,-1, -11.44,   5.75,   5.32,   6.87, -0.627,  0.015,
     &2, 0, 0, 0, 0,-2,  -1.24,   0.63,   0.58,   0.75, -0.068,  0.002,
     &2, 1, 0, 0, 0, 0,  -1.77,   1.79,   1.71,   1.04, -0.146,  0.037,
     &2, 1, 0, 0, 0,-1,  -0.77,   0.78,   0.75,   0.45, -0.064,  0.017,
     &2, 0, 0, 2, 0, 2,  -0.33,   0.62,   0.65,   0.19, -0.049,  0.018/

      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   
      DARG(1) = (876600d0*3600d0 + 8640184.812866d0 
     .         + 2.d0 * 0.093104d0 * T - 3.d0 * 6.2d-6*T**2)*15.d0
      DARG(1) = DARG(1)* secrad / 36525.0D0   ! rad/day


      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      
      DARG(2) = -4.d0*0.00024470d0*T**3 + 3.d0*0.051635d0*T**2 
     .  + 2.d0*31.8792d0*T + 1717915923.2178d0 
      DARG(2) = DARG(2)* secrad / 36525.0D0   ! rad/day

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

      DARG(3) = -4.D0*0.00001149d0*T**3 - 3.d0*0.000136d0*T**2
     .  -  2.D0*0.5532d0*T + 129596581.0481d0
      DARG(3) = DARG(3)* secrad / 36525.0D0   ! rad/day
          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

      DARG(4) = 4.d0*0.00000417d0*T**3 - 3.d0*0.001037d0*T**2 
     .- 2.d0 * 12.7512d0*T + 1739527262.8478d0 
      DARG(4) = DARG(4)* secrad / 36525.0D0   ! rad/day
    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

      DARG(5) = -4.d0*0.00003169d0*T**3 + 3.d0*0.006593d0*T**2
     . - 2.d0 * 6.3706d0*T + 1602961601.2090d0
      DARG(5) = DARG(5)* secrad / 36525.0D0   ! rad/day

      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad

      DARG(6) = -4.d0*0.00005939d0*T**3 + 3.d0 * 0.007702d0*T**2
     .  + 2.d0 * 7.4722d0*T - 6962890.2665d0
      DARG(6) = DARG(6)* secrad / 36525.0D0   ! rad/day

C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0
	cor_ut1= 0.d0
	cor_lod= 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
 	dag = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
 		dag = dag + dble(narg(j,i))*DARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x= cor_x + dble(XCOS(j))*dcos(ag) + dble(XSIN(j)) * dsin(ag)
        cor_y= cor_y + dble(YCOS(j))*dcos(ag) + dble(YSIN(j)) * dsin(ag)
        cor_ut1= cor_ut1+dble(UTCOS(j))*dcos(ag)+dble(UTSIN(j))*dsin(ag)
        cor_lod= cor_lod -(-dble(UTCOS(j)) * dsin(ag) 
     &                    + dble(UTSIN(j)) * dcos(ag) ) * dag   	 

        enddo
  
       cor_x   = cor_x * 1.0d-6   ! arcsecond (")
       cor_y   = cor_y * 1.0d-6   ! arcsecond (")
       cor_ut1 = cor_ut1 * 1.0d-6 ! second (s)
       cor_lod = cor_lod * 1.0d-6 ! second (s)
 
      RETURN
      END
      	
C----------------------------------------------------------------
      SUBROUTINE PM_GRAVI (rjd,cor_x,cor_y)
C
C    This subroutine provides, in time domain, the diurnal
C    lunisolar effet on polar motion (")
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C
C     coded by Ch. Bizouard (2002)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER nlines
      PARAMETER(nlines=10)
      DOUBLE PRECISION ARG(6)    ! Array of the tidal arguments   
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines)
      INTEGER NARG(nlines,6)
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),j=1,nlines)/    
     & 1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44,
     & 1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31,
     & 1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44,
     & 1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14,
     & 1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36,
     & 1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84,
     & 1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76,
     & 1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27,
     & 1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93,
     & 1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76/
 
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   

      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

  
      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad


C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x+dble(XCOS(j))*dcos(ag)+dble(XSIN(j))*dsin(ag)
        cor_y =cor_y+dble(YCOS(j))*dcos(ag)+dble(YSIN(j))*dsin(ag) 

        enddo
  
      cor_x = cor_x * 1.0d-6   ! arcsecond (")
      cor_y = cor_y * 1.0d-6   ! arcsecond (")
 
      RETURN

      END
!-------------------------------------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION sla_PA (HA, DEC, PHI)
*+
*     - - -
*      P A
*     - - -
*
*  HA, Dec to Parallactic Angle (double precision)
*
*  Given:
*     HA     d     hour angle in radians (geocentric apparent)
*     DEC    d     declination in radians (geocentric apparent)
*     PHI    d     observatory latitude in radians (geodetic)
*
*  The result is in the range -pi to +pi
*
*  Notes:
*
*  1)  The parallactic angle at a point in the sky is the position
*      angle of the vertical, i.e. the angle between the direction to
*      the pole and to the zenith.  In precise applications care must
*      be taken only to use geocentric apparent HA,Dec and to consider
*      separately the effects of atmospheric refraction and telescope
*      mount errors.
*
*  2)  At the pole a zero result is returned.
*
*  P.T.Wallace   Starlink   16 August 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*
*  License:
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program (see SLA_CONDITIONS); if not, write to the 
*    Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
*    Boston, MA  02111-1307  USA
*
*-

      IMPLICIT NONE

      DOUBLE PRECISION HA,DEC,PHI

      DOUBLE PRECISION CP,SQSZ,CQSZ



      CP=COS(PHI)
      SQSZ=CP*SIN(HA)
      CQSZ=SIN(PHI)*COS(DEC)-CP*SIN(DEC)*COS(HA)
      IF (SQSZ.EQ.0D0.AND.CQSZ.EQ.0D0) CQSZ=1D0
      sla_PA=ATAN2(SQSZ,CQSZ)

      END
!---------------------------------------------------------------------------------------------------


      SUBROUTINE L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD, 
     .            XET,PV,UB,THD,EPV)
!            
!     This subroutine compute the geocentric and heliocentric true vector  
!          position , distance and radius vector of the body.       
      
!       Input:  ET     = TDB time
!               NTARG  = The numbering convention for the body 
!               EB     = Earth Barycentric position vector (ET time)
!               BVE    = Earth Barycentric velocity vector (ET time)
!               HMOD   = Earth heliocentric position modulus (ET time)  
!       Output  BODY_G = Geocentric position true vector 
!               BODY_H = Heliocentric position true vector
!               GD     = true distance (AU) 
!               HD     = Radius vector (AU)
!               XET    = TDB time corrected to light-time                            
!               PV     = Body's barycentric state vector (XET time)
!               UB     = Body's Heliocentric state vector
!               EPV    = Earth's Heliocentric state vector 

      IMPLICIT NONE
      INTEGER NTARG
      DOUBLE PRECISION C,AU,DS,AUDAY,MU,MUC,DET,DT,XET,HD,GD,THD,TGD
      DOUBLE PRECISION ET,HMOD,R0,R1,RTERM,TAU,RV,DST,DOT,HV(3)
      DOUBLE PRECISION EBT(3,2),SB(3,2),PV(3,2),EB(3),BODY_H(3),BGV(3)
      DOUBLE PRECISION PV1(3),EBP(3),HE(3),GE(3),RQ(3),RS(3),BODY_G(3)
      DOUBLE PRECISION BVE(3),GV(3),BV1(3),EBV(3),UB(3,2),UV(3,2)
      DOUBLE PRECISION EPV(3,2),GK,HEP(3)

      C  = 299792458D0                     ! speed of light  m/s 
      AU = 149597870691D0                  ! Astronomical Unit (meter)(TDB)
      DS = 86400D0                         ! n° seconds  Day
      AUDAY = (C * DS / AU)                ! Speed of light (AU per day  "TDB")
      MU = 1.32712440041D20                ! heliocentric gravitational costant (2010 system)   
      GK = 0.01720209895D0                 ! Gaussian gravitational constant
      MUC = 2D0 * (MU / (C * C))/AU

!---------Body's astrometric position (recursive determination)
      DT = 0.0D0
      DET = 0.0D0
  
      DO
       XET = ET - DET                          ! XET = new reference time
      
      CALL PLEPH ( XET, NTARG, 12, PV)         ! PV = body's barycentric state vector 
      
      CALL PLEPH ( XET, 11, 12, SB)            ! SB = sun's barycentric state vector  
      
      CALL PLEPH ( XET, 3, 11, EPV)            ! EPV = Earth heliocentric state vector

      CALL iau_PVUP (DT,EPV,HEP)               ! HEP = Earth Heliocentric position         

      CALL iau_PVMPV(PV,SB,UB)                 ! UB = Body's heliocentric state vector 
      
      CALL iau_PVUP (DT,PV,PV1)                ! PV1 = body's barycentric  position  

      CALL iau_PVUV (DT,PV,BV1)                ! BV1 = body's barycentric  velocity       

      CALL iau_PVUV (DT,UB,HV)                 ! HV Body's Heliocentric velocity        
!---------------------------------------------------------------------------------------------------
!------Approximation to light time 

      CALL iau_PVUP (DT,UB,HE)                 ! HE Body's Heliocentric position     
        
      CALL iau_PMP (PV1,EB,GE)                 ! GE Body's Geocentric  position 

      CALL iau_PMP (BV1,BVE,BGV)               ! BGV Body's Geocentric velocity  

!-----Form modulus of HE , GE, HEP
      
      CALL iau_PM(HE,HD)         ! HD modulus Heliocentric position body

      CALL iau_PM(GE,GD)         ! GD modulus Geocentric position Earth

      CALL iau_PM(HEP,HMOD)
      
      IF (DET .EQ. 0D0) THEN
        THD = HD                 ! THD True Heliocentric distance
        TGD = GD                 ! TGD True Geocentric distance
      END IF

      IF (HD .EQ. 0D0) THEN
        RTERM = 0D0
      ELSE
        RTERM = MUC * DLOG ((HMOD+GD+HD) /DABS (HMOD-GD+HD))
      END IF

      TAU = (GD + RTERM) / AUDAY
      
      IF ( DABS( DET - TAU) .LT. 1D-15 )  EXIT 
       
      DET = TAU

      END DO
      
!----------Form unit vector of " RQ , RS "     

      CALL iau_PN(HE,R0,RQ)            ! RQ Unic vector body_Helioc
        
      CALL iau_PN(GE,R1,RS)            ! RS Unic vector body_Geocentric

        
!----------------------------------------------------------------------------

      BODY_G = RS * TGD                ! True Vector Geocentric  Position
       
      BODY_H = RQ * THD                ! True Vector Heliocentric  Position
      PV = PV
      XET = XET 
      END
!============================================================================================================


      SUBROUTINE RA_SUN (mode,AUDAY,BODY_sun,RNPB,RC2I,BEP,HPE,BVE, 
     .                   MO1,R_sun,D_sun)

!     -------------------------------------------------------------------
!     This routine search the RA of Sun for compute the Equation of Time
!     -------------------------------------------------------------------      

      IMPLICIT NONE

!     Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

      DOUBLE PRECISION iau_ANP,BODY_sun(3),RNPB(3,3),RC2I(3,3)
      DOUBLE PRECISION BEP(3),HPE(3),BVE(3),MO1,AUDAY,p(3)
      DOUBLE PRECISION u0(3),u1(3),u2(3),u3(3),u4(3),R_sun,DE6,D_sun
      INTEGER  mode 
!----------------------------------------------------------------------------
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun ( for apparent geocentric coordinate)

      CALL DEFLIGHT_0(BODY_sun,BEP,HPE,u1)
 
!      CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

!      CALL  iau_C2S   ( p, RA5, DE5 )       
!      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
      
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent geocentric coordinate)

      CALL ABERR (AUDAY,u1,BVE,u0)

      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 
      
   
! ......... Final GCRS (RA5,DE5) of body 
! ......... Call for direction cosines to spherical coordinates
 
!      CALL  iau_C2S   ( u3, RA10, DE10 )       
!      RA10 = iau_ANP(RA10)                 ! range angle (0,2!pi)
      
!-----------------------------------------------------------------------------

!---------- Form apparent geocentric coordinate (RA6,DE6)       

!--------- Option from the method of calculating "IUA 2000 Equinox based" or 
!          "IAU 2006 CIO based using X,Y series" 

      IF (mode == 1) then 
       CALL   iau_RXP   ( RNPB, u3, u4 )
      ELSE IF (mode == 2) then 
       CALL   iau_RXP   ( RC2I, u3, u4 )
      END IF 
! ......... Call for direction cosines to spherical coordinates

      CALL   iau_C2S   ( u4, R_sun, DE6 )            ! RA of Sun in radians  
         R_sun= iau_ANP(R_sun)                       ! range angle (0,2!pi)
         R_sun = R_sun * DR2D                        ! RA of Sun in degrees
         D_sun = DE6 *DR2D                           ! Decl. in degrees
      END

!============================================================================================================


      SUBROUTINE ABERR (AUDAY,u1,BVE,u0)
   
!-------  RELATIVISTIC ANNUAL ABERRATION 
!         Input  AUDAY(constant, velocity of light " c " x 86400 sec. in  AU)
!                u1    Body's geocentric position state vector (AU)
!                BVE   Barycentric Earth velocity (AU/Day)
!         Output u0    Body's geocentric position state vector corrected 
!                       for aberration (AU)  


      IMPLICIT NONE

      DOUBLE PRECISION u2,V0,B1,F1,F2,AUDAY
   
      DOUBLE PRECISION p(3),V(3),BU(3),XV(3),FXV(3),V4(3) 

      DOUBLE PRECISION u1(3),BVE(3),u0(3)
      
      CALL iau_PN(u1,u2,p)          ! p = unit vector of geocentric position(u1)

      V = BVE / AUDAY                !AUDAY = 173.14463268465693 
 
      CALL iau_PM(V,V0)              ! V0 modulus V

      B1 = DSQRT(1.0D0-V0*V0)
     
      CALL iau_PDP(p,V,F1)           ! F1 product scalar p x V

      F2 = 1.0D0 +F1 /(1.0D0 + B1)

      CALL iau_SXP(B1,u1,BU)         ! multiply scalar B1 x vector u1

      CALL iau_SXP(u2,V,XV)          ! multiply scalar u2 x vector V  ( call iau_PDP(u1,V,UU))

      CALL iau_SXP(F2,XV,FXV)        ! multiply scalar F2 x vector XV

      CALL iau_PPP(BU,FXV,V4)        ! add vector (BU + FXV)

      u0 =  V4 / (1.0D0 + F1)        ! u0 = Position vector corrected for aberration

      END
      
!-----------------------------------------------------------------------------------------------------

      SUBROUTINE DEFLIGHT_0(BODY,BEP,HPE,u1)

!-----RELATIVISTIC DEFLECTION OF LIGHT DUE THE SUN
!     Input:  BODY  Geocentric vector of body position 
!             BEP   Heliocentric vector position of body
!             HPE   Heliocentric vector Earth position
!     Output  u1    Geocentric vector position corrected for 
!                   relativistic deflection of light due the Sun

      IMPLICIT NONE

      DOUBLE PRECISION  C                  ! speed of light  m/s 
      PARAMETER ( C = 299792458D0 )
     
      DOUBLE PRECISION AU
      PARAMETER (AU = 149597870691D0)      ! Astronomical Unit (meter)(TDB)

      DOUBLE PRECISION MU
      PARAMETER (MU = 1.32712440017987D20) ! Heliocentric gravitational constant (m3 s-2 : TDB)

      DOUBLE PRECISION R0,R1,R2,QES,UQS,EUS,g1,g2,UM
      DOUBLE PRECISION BODY(3),BEP(3),HPE(3),U(3),Q(3),E(3),u1(3) 
           
!     Form unit vector of " U , Q , E "

      CALL iau_PN(BODY,R0,U)       

      CALL iau_PN(BEP,R1,Q)           

      CALL iau_PN(HPE,R2,E)        ! R2 = modulus of HPE 
      
      UM = R0                      ! UM = R0 modulus of vector BODY

      CALL iau_PDP(Q,E,QES)        ! QES product scalar Q x E 

      CALL iau_PDP(U,Q,UQS)        ! UQS product scalar Q x U

      CALL iau_PDP(E,U,EUS)        ! EUS product scalar E x U

            
      
      g1 = 2D0 * ( MU / (C * C * R2))/AU
      g2 = 1D0 + QES     
      
      u1 = UM * (U + (g1 / g2) * (UQS * E - EUS * Q))
       
      END   
!----------------------------------------------------------------------------------------------------

      SUBROUTINE DEFLIGHT_T(BODY,BEP,HPE,u2)

!-----RELATIVISTIC DEFLECTION OF LIGHT DUE THE EARTH
!             ( use for geocentric topocentric position)
!     Input:  BODY  Geocentric vector of body position 
!             BEP   Heliocentric vector position of body
!             HPE   Heliocentric vector Earth position
!     Output  u2    Geocentric vector position corrected for 
!                   relativistic deflection of light due the Earth

      IMPLICIT NONE

      DOUBLE PRECISION  C                  ! speed of light  m/s 
      PARAMETER ( C = 299792458D0 )
     
      DOUBLE PRECISION AU
      PARAMETER (AU = 149597870691D0)      ! Astronomical Unit (meter)(TDB)

      DOUBLE PRECISION MU
      PARAMETER (MU = 3.98600433D14)        ! Geocentric gravitational constant (m3 s-2 : TDB)) 

      DOUBLE PRECISION R0,R1,R2,QES,UQS,EUS,g1,g2,UM
      DOUBLE PRECISION BODY(3),BEP(3),HPE(3),U(3),Q(3),E(3),u2(3) 
           
!     Form unit vector of " U , Q , E "

      CALL iau_PN(BODY,R0,U)       

      CALL iau_PN(BEP,R1,Q)           

      CALL iau_PN(HPE,R2,E)        ! R2 = modulus of HPE 
      
      UM = R0                      ! UM = R0 modulus of vector BODY

      CALL iau_PDP(Q,E,QES)        ! QES product scalar Q x E 

      CALL iau_PDP(U,Q,UQS)        ! UQS product scalar Q x U

      CALL iau_PDP(E,U,EUS)        ! EUS product scalar E x U

            
      
      g1 = 2D0 * ( MU / (C * C * R2))/AU
      g2 = 1D0 + QES     
      
      u2 = UM * (U + (g1 / g2) * (UQS * E - EUS * Q))
       
      END   
 
!--------------------------------------------------------------------------------------------------------
      SUBROUTINE GEOCPV ( ELONG, PHI, HEIGHT, GST, POS, VEL, J )
!     This subroutine compute the velocity vector of a terrestrial observer 
!          with respect to the geocenter.       
      
!       Input:  ELONG  = Longitude of geodetic observer in radians ( + EAST)
!               PHI    = Latitude  of geodetic observer in radians ( + NORD)
!               GST    = Greenwich apparent sideral time in radians
!               HEIGHT = HEIGHT of observer in meter 
!       Output  VEL    = Velocity vector of observer (Geocenter equatorial 
!                          rectangular coordinate in AU/Day)
!               POS    = Position vector of observer (Geocenter equatorial 
!                          rectangular coordinate in AU)
!               J      = status:  0 = OK
!                                -1 = illegal case 
!
      IMPLICIT NONE
      
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )  ! Radians to degrees

      DOUBLE PRECISION AU
      PARAMETER ( AU = 149597870.691D0)                    ! Astronomical Unit (Km)(TDB) 
   
!     NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF EARTH
!     RADIANS/SECOND, FROM IERS CONVENTIONS (2003)      
      DOUBLE PRECISION OMEGA
      PARAMETER (OMEGA = 7.2921150D-5)                   

!     EQUATORIAL RADIUS OF EARTH IN Km, FROM IERS CONVENTIONS (2003)
      DOUBLE PRECISION A
      PARAMETER ( A = 6378.1366D0 )

!     FLATTENING FACTOR OF EARTH, FROM IERS CONVENTIONS (2003)
      DOUBLE PRECISION F
      PARAMETER ( F = 1.D0 / 298.25642D0 )

      DOUBLE PRECISION ELONG, PHI, HEIGHT, GST, VEL(3),POS(3)
      INTEGER J,I

      DOUBLE PRECISION SP, CP, W, D, AC, AS, R, GP(3),GV(3)


!  Functions of geodetic latitude.
      SP = SIN(PHI)
      CP = COS(PHI)
      W = 1D0-F
      W = W*W
      D = CP*CP + W*SP*SP
      IF ( D .GT. 0D0 ) THEN
         AC = A / SQRT(D)
         AS = W * AC

*     Geocentric position vector.
      R = ( AC + HEIGHT/1000D0 ) * CP     
      POS(1) = R * COS(ELONG+GST)
      POS(2) = R * SIN(ELONG+GST)
      POS(3) = ( AS + HEIGHT/1000D0 ) * SP
      

*     Geocentric velocity vector in KM/SEC
      VEL(1) = -OMEGA * R * sin(ELONG+GST)
      VEL(2) =  OMEGA * R * cos(ELONG+GST)
      VEL(3) =  0.D0

*     Success.
         J = 0
      ELSE

*     Fail.
         J = -1
      END IF

*     CONVERT POSITION AND VELOCITY COMPONENTS TO AU AND AU/DAY
      DO 10 I=1,3
      POS(I) = POS(I) / AU
      VEL(I) = VEL(I) / AU * 86400.D0
   10 CONTINUE 

      END


      SUBROUTINE PRECES (TJD1,POS1,TJD2,POS2)
*
*     THIS SUBROUTINE PRECESSES EQUATORIAL RECTANGULAR COORDINATES FROM
*     ONE EPOCH TO ANOTHER.  THE COORDINATES ARE REFERRED TO THE MEAN
*     DYNAMICAL EQUATOR AND EQUINOX OF THE TWO RESPECTIVE EPOCHS.  SEE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, PP. 103-104,
*     AND CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412,
*     567-586.
*
*          TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
*          POS1 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF FIRST EPOCH (IN)
*          TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
*          POS2 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF SECOND EPOCH (OUT)
*
*     NOTE:  EITHER TJD1 OR TJD2 MUST BE 2451545.0 (J2000.0) TDB.
*
*
      DOUBLE PRECISION TJD1,TJD2,POS1,POS2,PI,SECCON,T0,TLAST,T,
     .     EPS0,PSIA,OMEGAA,CHIA,SA,CA,SB,CB,SC,CC,SD,CD,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /
      
*     INITIALIZE PRECESSION ROTATION MATRIX AS IDENTITY MATRIX
      DATA XX, XY, XZ / 1.D0, 0.D0, 0.D0 /
      DATA YX, YY, YZ / 0.D0, 1.D0, 0.D0 /
      DATA ZX, ZY, ZZ / 0.D0, 0.D0, 1.D0 /

   3  FORMAT ( ' PRECES ERROR: PRECESSION FROM JD ', F10.1, ' TO ',
     .     F10.1, ' NOT TO/FROM J2000' )

      IF ( TJD1 .NE. T0 .AND. TJD2 .NE. T0 ) THEN
          WRITE ( *, 3 ) TJD1, TJD2
          GO TO 50
      END IF

*     T IS TIME IN TDB CENTURIES BETWEEN THE TWO EPOCHS
      T = ( TJD2 - TJD1 ) / 36525.D0
      IF ( TJD2 .EQ. T0 ) T = -T
      IF ( DABS ( T - TLAST ) .LT. 1.D-15 ) GO TO 20

*     NUMERICAL COEFFICIENTS OF PSI_A, OMEGA_A, AND CHI_A, ALONG WITH
*     EPSILON_0, THE OBLIQUITY AT J2000.0, ARE 4-ANGLE FORMULATION
*     FROM CAPITAINE ET AL. (2003), EQS. (4), (37), & (39)
      EPS0   = 84381.406D0
      PSIA   = ( ( ( ( -    0.0000000951D0   * T
     .                 +    0.000132851D0  ) * T
     .                 -    0.00114045D0   ) * T
     .                 -    1.0790069D0    ) * T
     .                 + 5038.481507D0     ) * T
      OMEGAA = ( ( ( ( +    0.0000003337D0   * T
     .                 -    0.000000467D0  ) * T
     .                 -    0.00772503D0   ) * T
     .                 +    0.0512623D0    ) * T
     .                 -    0.025754D0     ) * T + EPS0
      CHIA   = ( ( ( ( -    0.0000000560D0   * T
     .                 +    0.000170663D0  ) * T
     .                 -    0.00121197D0   ) * T
     .                 -    2.3814292D0    ) * T
     .                 +   10.556403D0     ) * T
      EPS0 = EPS0 / SECCON
      PSIA = PSIA / SECCON
      OMEGAA = OMEGAA / SECCON
      CHIA = CHIA / SECCON
      SA = DSIN ( EPS0 )
      CA = DCOS ( EPS0 )
      SB = DSIN ( -PSIA )
      CB = DCOS ( -PSIA )
      SC = DSIN ( -OMEGAA )
      CC = DCOS ( -OMEGAA )
      SD = DSIN ( CHIA )
      CD = DCOS ( CHIA )

*     COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
*     EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
      XX =  CD * CB - SB * SD * CC
      YX =  CD * SB * CA + SD * CC * CB * CA - SA * SD * SC
      ZX =  CD * SB * SA + SD * CC * CB * SA + CA * SD * SC
      XY = -SD * CB - SB * CD * CC
      YY = -SD * SB * CA + CD * CC * CB * CA - SA * CD * SC
      ZY = -SD * SB * SA + CD * CC * CB * SA + CA * CD * SC
      XZ =  SB * SC
      YZ = -SC * CB * CA - SA * CC
      ZZ = -SC * CB * SA + CC * CA

      TLAST = T

   20 IF ( TJD2 .EQ. T0 ) GO TO 30

*     PERFORM ROTATION FROM J2000.0 TO EPOCH
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM EPOCH TO J2000.0
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END
!----------------------------------------------------------------------------------------------

      SUBROUTINE FRAME (POS1,K,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM THE DYNAMICAL REFERENCE
*     SYSTEM TO THE INTERNATIONAL CELESTIAL REFERENCE SYSTEM (ICRS),
*     OR VICE VERSA.  THE DYNAMICAL REFERENCE SYSTEM IS BASED ON THE
*     DYNAMICAL MEAN EQUATOR AND EQUINOX OF J2000.0.  THE ICRS IS
*     BASED ON THE SPACE-FIXED ICRS AXES DEFINED BY THE RADIO CATALOG
*     POSITIONS OF SEVERAL HUNDRED EXTRAGALACTIC OBJECTS.  THE ROTATION
*     MATRIX USED HERE IS EQUIVALENT TO THAT GIVEN BY HILTON AND
*     HOHENKERK (2004), ASTRONOMY AND ASTROPHYSICS 413, 765-770,
*     EQ. (6) AND (8).
*
*          POS1   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (IN)
*          K      = DIRECTION OF ROTATION (IN)
*                   SET K < 0 FOR DYNAMICAL TO ICRS
*                   SET K > 0 FOR ICRS TO DYNAMICAL
*          POS2   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (OUT)
*
*     NOTE:  FOR GEOCENTRIC COORDINATES, THE SAME TRANSFORMATION IS
*     USED BETWEEN THE DYNAMICAL REFERENCE SYSTEM AND THE GCRS.
*
*
      DOUBLE PRECISION POS1,POS2,PI,SECCON,XI0,ETA0,DA0,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     XI0, ETA0, AND DA0 ARE ICRS FRAME BIASES IN ARCSECONDS TAKEN
*     FROM IERS CONVENTIONS (2003), CHAPTER 5
      DATA XI0, ETA0, DA0 / -0.0166170D0, -0.0068192D0, -0.01460D0 /

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1

*     COMPUTE ELEMENTS OF ROTATION MATRIX (TO FIRST ORDER)
      IF ( NTIMES .GT. 1 ) GO TO 20
      XX =  1.D0
      YX = -DA0  / SECCON
      ZX =  XI0  / SECCON
      XY =  DA0  / SECCON
      YY =  1.D0
      ZY =  ETA0 / SECCON
      XZ = -XI0  / SECCON
      YZ = -ETA0 / SECCON
      ZZ =  1.D0
*     INCLUDE SECOND-ORDER CORRECTIONS TO DIAGONAL ELEMENTS
      XX = 1.D0 - 0.5D0 * ( YX**2 + ZX**2 )
      YY = 1.D0 - 0.5D0 * ( YX**2 + ZY**2 )
      ZZ = 1.D0 - 0.5D0 * ( ZY**2 + ZX**2 )
   20 IF ( K .GE. 0 ) GO TO 30

*     PERFORM ROTATION FROM DYNAMICAL SYSTEM TO ICRS
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM ICRS TO DYNAMICAL SYSTEM
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END

!---------------------------------------------------------------------------------------------



*      SUBROUTINE GETMOD ( MODE )
*
*     THIS SUBROUTINE RETURNS THE 'MODE' VALUE, WHICH
*     DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (OUT)
*                   MODE=0 MEANS CIO-BASED METHOD, FULL ACCURACY
*                   MODE=1 MEANS CIO-BASED METHOD, REDUCED ACCURACY
*                   MODE=2 MEANS EQUINOX-BASED METHOD, FULL ACCURACY
*                   MODE=3 MEANS EQUINOX-BASED METHOD, REDUCED ACCURACY
*
*
*      MODE = IMODE

*      RETURN

*      END
!-----------------------------------------------------------------------------------------

      SUBROUTINE NOD (T,DPSI,DEPS)
*
*     THIS SUBROUTINE RETURNS THE VALUES FOR NUTATION IN LONGITUDE AND
*     NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
*
*          T     = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          DPSI  = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
*          DEPS  = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
*
*
      DOUBLE PRECISION T,DPSI,DEPS,PI,SECCON,T0,T1,DP,DE
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      T1 = T * 36525.D0

*     =================================================================
*     EVALUATE NUTATION SERIES
*     RESULTING NUTATION IN LONGITUDE AND OBLIQUITY IN ARCSECONDS

*     CALL SUBROUTINE TO EVALUATE NUTATION SERIES
      IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*         HIGH ACCURACY MODE -- IERS SUBROUTINE
          CALL iau_NUT00A ( T0, T1, DP, DE )
*      ELSE
*         LOW ACCURACY MODE -- MODIFICATION OF IERS SUBROUTINE
*          CALL NU2000K ( T0, T1, DP, DE )
      END IF
      DPSI = DP * SECCON
      DEPS = DE * SECCON

*     =================================================================

      RETURN

      END

!---------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION ANMP ( A )

*  Normalize angle into the range -pi <= A < +pi.

      IMPLICIT NONE

      DOUBLE PRECISION A

      DOUBLE PRECISION DPI, D2PI
      PARAMETER ( DPI = 3.141592653589793238462643D0,
     :            D2PI = 6.283185307179586476925287D0 )

      DOUBLE PRECISION W

      W = MOD(A,D2PI)
      IF ( ABS(W) .GE. DPI ) W = W - SIGN(D2PI,A)
      ANMP = W

      END
!--------------------------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION  EECT2000 ( DATE1, DATE2 )
*+
*  - - - - - - - - -
*   E E C T 2 0 0 0
*  - - - - - - - - -
*
*  Equation of the equinoxes complementary terms, consistent with
*  IAU 2000 resolutions.
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Capitaine, N., Wallace, P.T., & McCarthy, D.D. (2003). Astron. &
*    Astrophys. 406, pp. 1135-1149, Table 3.
*  IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
*    (Table 5.2e presented in the printed publication is a truncated
*    series. The full series, which is used in NOVAS, is available on
*    the IERS Conventions Center website in file tab5.2e.txt.)
*    ftp://tai.bipm.org/iers/conv2010/chapter5/
*
*  Given:
*     DATE1,DATE2   d    TT date (JD = DATE1+DATE2)
*
*  Returned:
*     EECT00        d    complementary terms (radians)
*
*  This revision:  2002 November 13
*                  References updated 2010 November 26
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Time since J2000, in Julian centuries
      DOUBLE PRECISION T

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION A, S0, S1
      DOUBLE PRECISION ANMP

*  Fundamental arguments
      DOUBLE PRECISION FA(14)

*  -----------------------------------------
*  The series for the EE complementary terms
*  -----------------------------------------

*  Number of terms in the series
      INTEGER NE0, NE1
      PARAMETER ( NE0=  33, NE1=  1 )

*  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
      INTEGER KE0 ( 14, NE0 ),
     :        KE1 ( 14, NE1 )

*  Sine and cosine coefficients
      DOUBLE PRECISION SE0 ( 2, NE0 ),
     :                 SE1 ( 2, NE1 )

*  Argument coefficients for t^0
      DATA ( ( KE0(I,J), I=1,14), J =    1,   10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   11,   20 ) /
     :  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   21,   30 ) /
     :  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1,
     :  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   31,  NE0 ) /
     :  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

*  Argument coefficients for t^1
      DATA ( ( KE1(I,J), I=1,14), J =    1,  NE1 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

*  Sine and cosine coefficients for t^0
      DATA ( ( SE0(I,J), I=1,2), J =    1,   10 ) /
     :            +2640.96D-6,          -0.39D-6,
     :              +63.52D-6,          -0.02D-6,
     :              +11.75D-6,          +0.01D-6,
     :              +11.21D-6,          +0.01D-6,
     :               -4.55D-6,          +0.00D-6,
     :               +2.02D-6,          +0.00D-6,
     :               +1.98D-6,          +0.00D-6,
     :               -1.72D-6,          +0.00D-6,
     :               -1.41D-6,          -0.01D-6,
     :               -1.26D-6,          -0.01D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   11,   20 ) /
     :               -0.63D-6,          +0.00D-6,
     :               -0.63D-6,          +0.00D-6,
     :               +0.46D-6,          +0.00D-6,
     :               +0.45D-6,          +0.00D-6,
     :               +0.36D-6,          +0.00D-6,
     :               -0.24D-6,          -0.12D-6,
     :               +0.32D-6,          +0.00D-6,
     :               +0.28D-6,          +0.00D-6,
     :               +0.27D-6,          +0.00D-6,
     :               +0.26D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   21,   30 ) /
     :               -0.21D-6,          +0.00D-6,
     :               +0.19D-6,          +0.00D-6,
     :               +0.18D-6,          +0.00D-6,
     :               -0.10D-6,          +0.05D-6,
     :               +0.15D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               +0.13D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   31,  NE0 ) /
     :               -0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6 /

*  Sine and cosine coefficients for t^1
      DATA ( ( SE1(I,J), I=1,2), J =    1,  NE1 ) /
     :               -0.87D-6,          +0.00D-6 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  Fundamental Arguments (from IERS Conventions 2000)

*  Mean Anomaly of the Moon.
      FA(1) = ANMP ( ( 485868.249036D0 +
     :               ( 715923.2178D0 +
     :               (     31.8792D0 +
     :               (      0.051635D0 +
     :               (     -0.00024470D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1325D0*T, 1D0 ) * D2PI )

*  Mean Anomaly of the Sun.
      FA(2) = ANMP ( ( 1287104.793048D0 +
     :               ( 1292581.0481D0 +
     :               (      -0.5532D0 +
     :               (      +0.000136D0 +
     :               (      -0.00001149D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 99D0*T, 1D0 ) * D2PI )

*  Mean Longitude of the Moon minus Mean Longitude of the Ascending
*  Node of the Moon.
      FA(3) = ANMP ( (  335779.526232D0 +
     :               (  295262.8478D0 +
     :               (     -12.7512D0 +
     :               (      -0.001037D0 +
     :               (       0.00000417D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1342D0*T, 1D0 ) * D2PI )

*  Mean Elongation of the Moon from the Sun.
      FA(4) = ANMP ( ( 1072260.703692D0 +
     :               ( 1105601.2090D0 +
     :               (      -6.3706D0 +
     :               (       0.006593D0 +
     :               (      -0.00003169D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1236D0*T, 1D0 ) * D2PI )

*  Mean Longitude of the Ascending Node of the Moon.
      FA(5) = ANMP ( (  450160.398036D0 +
     :               ( -482890.5431D0 +
     :               (       7.4722D0 +
     :               (       0.007702D0 +
     :               (      -0.00005939D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( -5D0*T, 1D0 ) * D2PI )

      FA( 6) = ANMP ( 4.402608842D0 + 2608.7903141574D0 * T )
      FA( 7) = ANMP ( 3.176146697D0 + 1021.3285546211D0 * T )
      FA( 8) = ANMP ( 1.753470314D0 +  628.3075849991D0 * T )
      FA( 9) = ANMP ( 6.203480913D0 +  334.0612426700D0 * T )
      FA(10) = ANMP ( 0.599546497D0 +   52.9690962641D0 * T )
      FA(11) = ANMP ( 0.874016757D0 +   21.3299104960D0 * T )
      FA(12) = ANMP ( 5.481293872D0 +    7.4781598567D0 * T )
      FA(13) = ANMP ( 5.311886287D0 +    3.8133035638D0 * T )
      FA(14) =      ( 0.024381750D0 +    0.00000538691D0 * T ) * T

*  Evaluate the EE complementary terms.
      S0 = 0D0
      S1 = 0D0

      DO I = NE0,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE0(J,I))*FA(J)
         END DO
         S0 = S0 + ( SE0(1,I)*SIN(A) + SE0(2,I)*COS(A) )
      END DO
      DO I = NE1,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE1(J,I))*FA(J)
         END DO
         S1 = S1 + ( SE1(1,I)*SIN(A) + SE1(2,I)*COS(A) )
      END DO
      EECT2000 = ( S0 + S1 * T ) * DAS2R

*  Finished.

      END

!----------------------------------------------------------------------------------------

      SUBROUTINE FUNARG (T,EL,ELP,F,D,OMEGA)
*
*     THIS SUBROUTINE COMPUTES FUNDAMENTAL ARGUMENTS (MEAN ELEMENTS)
*     OF THE SUN AND MOON.  SEE SIMON ET AL. (1994) ASTRONOMY AND
*     ASTROPHYSICS 282, 663-683, ESPECIALLY SECTIONS 3.4-3.5.
*
*          T      = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          EL     = MEAN ANOMALY OF THE MOON IN RADIANS
*                   AT DATE TJD (OUT)
*          ELP    = MEAN ANOMALY OF THE SUN IN RADIANS
*                   AT DATE TJD (OUT)
*          F      = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
*                   OF THE MOON'S ASCENDING NODE IN RADIANS
*                   AT DATE TJD (OUT)
*          D      = MEAN ELONGATION OF THE MOON FROM THE SUN IN
*                   RADIANS AT DATE TJD (OUT)
*          OMEGA  = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
*                   IN RADIANS AT DATE TJD (OUT)
*
*
      DOUBLE PRECISION T,EL,ELP,F,D,OMEGA,PI,SECCON,REV,DMOD

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
      PARAMETER ( REV    = 360.D0 * 3600.D0      )

*     FUNDAMENTAL (DELAUNAY) ARGUMENTS FROM SIMON ET AL. (1994)

*     MEAN ANOMALY OF THE MOON
      EL    = DMOD (         485868.249036D0 +
     .               T*( 1717915923.2178D0 +
     .               T*(         31.8792D0 +
     .               T*(          0.051635D0 +
     .               T*(        - 0.00024470D0 )))), REV ) / SECCON

*     MEAN ANOMALY OF THE SUN
      ELP   = DMOD (        1287104.79305D0 +
     .               T*(  129596581.0481D0 +
     .               T*(        - 0.5532D0 +
     .               T*(          0.000136D0 +
     .               T*(        - 0.00001149D0 )))), REV ) / SECCON

*     MEAN ARGUMENT OF THE LATITUDE OF THE MOON
      F     = DMOD (         335779.526232D0 +
     .               T*( 1739527262.8478D0 +
     .               T*(       - 12.7512D0 +
     .               T*(       -  0.001037D0 +
     .               T*(          0.00000417D0 )))), REV ) / SECCON

*     MEAN ELONGATION OF THE MOON FROM THE SUN
      D     = DMOD (        1072260.70369D0 +
     .               T*( 1602961601.2090D0 +
     .               T*(        - 6.3706D0 +
     .               T*(          0.006593D0 +
     .               T*(        - 0.00003169D0 )))), REV ) / SECCON

*     MEAN LONGITUDE OF THE ASCENDING NODE OF THE MOON (FROM SIMON
*     SECTION 3.4(b.3), PRECESSION=5028.8200 ARCSEC/CY)
      OMEGA = DMOD (         450160.398036D0 +
     .               T*(  - 6962890.5431D0 +
     .               T*(          7.4722D0 +
     .               T*(          0.007702D0 +
     .               T*(        - 0.00005939D0 )))), REV ) / SECCON

      RETURN

      END

!------------------------------------------------------------------------------------------

      SUBROUTINE ETILT (TJD,OBLM,OBLT,EQEQ,DPSI,DEPS)
*
*     THIS SUBROUTINE COMPUTES QUANTITIES RELATED TO THE ORIENTATION
*     OF THE EARTH'S ROTATION AXIS AT JULIAN DATE TJD.
*
*          TJD    = TDB JULIAN DATE FOR ORIENTATION PARAMETERS (IN)
*          OBLM   = MEAN OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          OBLT   = TRUE OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          EQEQ   = EQUATION OF THE EQUINOXES IN TIME SECONDS AT
*                   DATE TJD (OUT)
*          DPSI   = NUTATION IN LONGITUDE IN ARCSECONDS AT
*                   DATE TJD (OUT)
*          DEPS   = NUTATION IN OBLIQUITY IN ARCSECONDS AT
*                   DATE TJD (OUT)
*
*     NOTE:  THE EQUATION OF THE EQUINOXES INCLUDES THE COMPLEMENTARY
*     TERMS.
*
*
      DOUBLE PRECISION TJD,OBLM,OBLT,EQEQ,DPSI,DEPS,PI,SECCON,
     .     T0,TLAST,T,PSI,EPS,PSICOR,EPSCOR,CTERMS,DELPSI,DELEPS,
     .     EL,ELP,F,D,OMEGA,OBM,OBT,EE,
     .     DPOLE1,DPOLE2,DX,DY,DZ,SINE,X,DP1,DP2,DP3,
     .     OBLIQ,DABS,DSIN,DCOS,EECT2000,OCC
      INTEGER ACCDIF
      DIMENSION DP1(3), DP2(3), DP3(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   MLAST / 0 /
      DATA DELPSI, DELEPS, CTERMS, PSICOR, EPSCOR / 5 * 0.D0 /

*     FUNCTION TO COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
*     CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412, 567-586,
*     EXPRESSION FROM EQ. (39) WITH OBLIQUITY AT J2000.0 TAKEN FROM
*     EQ. (37) OR TABLE 8
      OBLIQ(T) = ( ( ( ( -  0.0000000434D0   * T
     .                   -  0.000000576D0  ) * T
     .                   +  0.00200340D0   ) * T
     .                   -  0.0001831D0    ) * T
     .                   - 46.836769D0     ) * T + 84381.406D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

*     CHECK FOR DIFFERENCE IN ACCURACY MODE FROM LAST CALL
      ACCDIF = MOD ( MODE, 2 ) - MOD ( MLAST, 2 )
      
      T = ( TJD - T0 ) / 36525.D0

      IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 .OR. ACCDIF .NE. 0 ) THEN

*         OBTAIN NUTATION PARAMETERS IN ARCSECONDS
          CALL NOD ( T,   PSI, EPS )

*         OBTAIN COMPLEMENTARY TERMS FOR EQUATION OF THE EQUINOXES
*         IN ARCSECONDS
          IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*             HIGH-ACCURACY MODE
              CTERMS = EECT2000 ( TJD, 0.D0 ) * SECCON
          ELSE
*             LOW-ACCURACY MODE
              CALL FUNARG ( T,   EL, ELP, F, D, OMEGA )
*             SERIES FROM IERS CONVENTIONS (2003), CHAPTER 5,
*             TABLE 5.2C, WITH SOME ADJUSTMENTS TO COEFFICIENT VALUES
*             COPIED FROM IERS FUNCTION EECT2000, WHICH HAS A MORE
*             COMPLETE SERIES
              CTERMS =
     .          2640.96D-6 * DSIN ( OMEGA )
     .        +   63.52D-6 * DSIN ( 2.D0 * OMEGA )
     .        +   11.75D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 3.D0 * OMEGA )
     .        +   11.21D-6 * DSIN ( 2.D0 * F - 2.D0 * D +        OMEGA )
     .        -    4.55D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 2.D0 * OMEGA )
     .        +    2.02D-6 * DSIN ( 2.D0 * F            + 3.D0 * OMEGA )
     .        +    1.98D-6 * DSIN ( 2.D0 * F            +        OMEGA )
     .        -    1.72D-6 * DSIN ( 3.D0 * OMEGA )
     .        -    0.87D-6 * T * DSIN ( OMEGA )
*             (TERMS SMALLER THAN 2 MICROARCSECONDS OMITTED)
          END IF
          TLAST = TJD
          MLAST = MODE

      END IF

      DELPSI = PSI + PSICOR
      DELEPS = EPS + EPSCOR

*     COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBM = OBLIQ(T)

*     COMPUTE TRUE OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBT = OBM + DELEPS

*     COMPUTE EQUATION OF THE EQUINOXES IN ARCSECONDS
      EE = DELPSI * DCOS ( OBM / SECCON ) + CTERMS

*     CONVERT TO OUTPUT UNITS
      OBLM = OBM / 3600.D0
      OBLT = OBT / 3600.D0
      EQEQ = EE  / 15.D0
      DPSI = DELPSI
      DEPS = DELEPS

      RETURN

      END
*--------------------------------------------------------------------------------------

      SUBROUTINE NUTATE (TJD,POS1,POS2)
*
*     THIS SUBROUTINE NUTATES EQUATORIAL RECTANGULAR COORDINATES FROM
*     THE MEAN DYNAMICAL EQUATOR AND EQUINOX OF EPOCH TO THE TRUE
*     EQUATOR AND EQUINOX OF EPOCH.  SEE EXPLANATORY SUPPLEMENT TO THE
*     ASTRONOMICAL ALMANAC, PP. 114-115.
*
*          TJD    = TDB JULIAN DATE OF EPOCH (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                   EQUINOX OF EPOCH (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF EPOCH (OUT)
*
*     NOTE:  IF TJD IS NEGATIVE, INVERSE NUTATION (TRUE TO MEAN)
*     IS APPLIED.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,TJD1,PI,SECCON,OBLM,OBLT,EQEQ,
     .     DPSI,DEPS,COBM,SOBM,COBT,SOBT,CPSI,SPSI,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      TJD1 = DABS(TJD)

      CALL ETILT ( TJD1,   OBLM, OBLT, EQEQ, DPSI, DEPS )
      OBLM = OBLM * 3600.D0 / SECCON
      OBLT = OBLT * 3600.D0 / SECCON
      DPSI = DPSI / SECCON
      DEPS = DEPS / SECCON
      COBM = DCOS ( OBLM )
      SOBM = DSIN ( OBLM )
      COBT = DCOS ( OBLT )
      SOBT = DSIN ( OBLT )
      CPSI = DCOS ( DPSI )
      SPSI = DSIN ( DPSI )

*     COMPUTE ELEMENTS OF NUTATION ROTATION MATRIX
      XX =  CPSI
      YX = -SPSI * COBM
      ZX = -SPSI * SOBM
      XY =  SPSI * COBT
      YY =  CPSI * COBM * COBT + SOBM * SOBT
      ZY =  CPSI * SOBM * COBT - COBM * SOBT
      XZ =  SPSI * SOBT
      YZ =  CPSI * COBM * SOBT - SOBM * COBT
      ZZ =  CPSI * SOBM * SOBT + COBM * COBT
   10 IF ( TJD .LT. 0.D0 ) GO TO 30

*     PERFORM ROTATION FROM MEAN TO TRUE
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM TRUE TO MEAN
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END

*-------------------------------------------------------------------------------------------------
      SUBROUTINE NUTATION (TJD,DJMJD0,TT,POS1,POS2)
*
*     THIS SUBROUTINE NUTATES EQUATORIAL RECTANGULAR COORDINATES FROM
*     THE MEAN DYNAMICAL EQUATOR AND EQUINOX OF EPOCH TO THE TRUE
*     EQUATOR AND EQUINOX OF EPOCH.  SEE EXPLANATORY SUPPLEMENT TO THE
*     ASTRONOMICAL ALMANAC, PP. 114-115.
*
*          TJD    = TDB JULIAN DATE OF EPOCH (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                   EQUINOX OF EPOCH (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF EPOCH (OUT)
*
*     NOTE:  IF TJD IS NEGATIVE, INVERSE NUTATION (TRUE TO MEAN)
*     IS APPLIED.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,OBLM,OBLT,
     .     DPSI,DEPS,COBM,SOBM,COBT,SOBT,CPSI,SPSI,DP00,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DE00,DABS,DCOS,DSIN
*      Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )
      DOUBLE PRECISION iau_OBL06
      DIMENSION POS1(3), POS2(3)
      
       
      OBLM = iau_OBL06(DJMJD0,TT) 

*         OBTAIN NUTATION PARAMETERS IN RADIANS
      CALL iau_NUT00A (DJMJD0, TT, DP00, DE00 ) 
     

      DPSI = DP00  
      DEPS = DE00  
      OBLT = OBLM + DEPS  
      COBM = DCOS ( OBLM )
      SOBM = DSIN ( OBLM )
      COBT = DCOS ( OBLT )
      SOBT = DSIN ( OBLT )
      CPSI = DCOS ( DPSI )      
      SPSI = DSIN ( DPSI )
       
      
*     COMPUTE ELEMENTS OF NUTATION ROTATION MATRIX
      XX =  CPSI
      YX = -SPSI * COBM
      ZX = -SPSI * SOBM
      XY =  SPSI * COBT
      YY =  CPSI * COBM * COBT + SOBM * SOBT
      ZY =  CPSI * SOBM * COBT - COBM * SOBT
      XZ =  SPSI * SOBT
      YZ =  CPSI * COBM * SOBT - SOBM * COBT
      ZZ =  CPSI * SOBM * SOBT + COBM * COBT
   10 IF ( TJD .LT. 0.D0 ) GO TO 30

*     PERFORM ROTATION FROM MEAN TO TRUE
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM TRUE TO MEAN
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END
*-------------------------------------------------------------------------------------------------
      SUBROUTINE iau_PN03 ( DATE1, DATE2, DPSI, DEPS,
     :                      EPSA, RB, RP, RBP, RN, RBPN )

      IMPLICIT NONE


      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS,
     :                 EPSA, RB(3,3), RP(3,3), RBP(3,3),
     :                 RN(3,3), RBPN(3,3)

      DOUBLE PRECISION DPSIPR, DEPSPR

      DOUBLE PRECISION  iau_OBL80, iau_OBL06

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  IAU 2000 precession-rate adjustments.
      CALL iau_PR00 ( DATE1, DATE2, DPSIPR, DEPSPR )

*  Mean obliquity, consistent with IAU 2000 precession-nutation.
      EPSA = iau_OBL06 ( DATE1, DATE2 ) 

*  Frame bias and precession matrices and their product.
      CALL iau_BP03 ( DATE1, DATE2, RB, RP, RBP )

*  Nutation matrix.
      CALL iau_NUMAT ( EPSA, DPSI, DEPS, RN )

*  Bias-precession-nutation matrix (classical).
      CALL iau_RXR ( RN, RBP, RBPN )

*  Finished.
      END
*---------------------------------------------------------------------------------------------

      SUBROUTINE IAU_BP03 ( DATE1, DATE2, RB, RP, RBP)

      IMPLICIT NONE

      

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000.0), JD
      DOUBLE PRECISION DJ00
      PARAMETER ( DJ00 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

!  J2000.0 obliquity (Lieske et al. 1977)
!      DOUBLE PRECISION EPS0
!      PARAMETER ( EPS0 = 84381.406D0 *DAS2R )  

      DOUBLE PRECISION DATE1, DATE2, RB(3,3), RP(3,3), RBP(3,3)
      DOUBLE PRECISION T, DPSIBI, DEPSBI, DRA0,EPS0,
     :                  CHIA, DPSIPR, DEPSPR, PSIA, OMA


*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ00 ) + DATE2 ) / DJC

*  Frame bias.
      CALL iau_BI00 ( DPSIBI, DEPSBI, DRA0 )
*========================================================================================
*     NUMERICAL COEFFICIENTS OF PSI_A, OMEGA_A, AND CHI_A, ALONG WITH
*     EPSILON_0, THE OBLIQUITY AT J2000.0, ARE 4-ANGLE FORMULATION
*     FROM CAPITAINE ET AL. (2003), EQS. (4), (37), & (39)
      EPS0   = 84381.406D0
      PSIA   =      ( ( ( ( - 0.0000000951D0 * T
     .                 +    0.000132851D0  ) * T
     .                 -    0.00114045D0   ) * T
     .                 -    1.0790069D0    ) * T
     .                 + 5038.481507D0     ) * T
      OMA =         ( ( ( ( + 0.0000003337D0 * T
     .                 -    0.000000467D0  ) * T
     .                 -    0.00772503D0   ) * T
     .                 +    0.0512623D0    ) * T
     .                 -    0.025754D0     ) * T + EPS0
      CHIA =        ( ( ( ( - 0.0000000560D0 * T
     .                 +    0.000170663D0  ) * T
     .                 -    0.00121197D0   ) * T
     .                 -    2.3814292D0    ) * T
     .                 +   10.556403D0     ) * T

      PSIA = PSIA * DAS2R
      OMA  = OMA * DAS2R 
      CHIA = CHIA * DAS2R
      EPS0 = EPS0 * DAS2R
*  Frame bias matrix: GCRS to J2000.0.
      CALL iau_IR ( RB )
      CALL iau_RZ ( DRA0, RB )
      CALL iau_RY ( DPSIBI*SIN(EPS0), RB )
      CALL iau_RX ( -DEPSBI, RB )

*  Precession matrix: J2000.0 to mean of date.
      CALL iau_IR ( RP )
      CALL iau_RX ( EPS0, RP )
      CALL iau_RZ ( -PSIA, RP )
      CALL iau_RX ( -OMA, RP )
      CALL iau_RZ ( CHIA, RP )

*  Bias-precession matrix: GCRS to mean of date.
      CALL iau_RXR ( RP, RB, RBP )

*  Finished.

      END
*======================================================================================================

      SUBROUTINE time_hms ( FD, IH, IM, SEC) 

      IMPLICIT NONE
      DOUBLE PRECISION FD,DH,DM,SEC
      INTEGER IH,IM

      IH = INT( FD * 24D0)        ! transform part of a day ( FD) 
      DH = ( FD *24D0) - IH       ! in hours, minuts and seconds
      IM = INT( DH * 60D0)
      DM = ( DH * 60D0) - IM
      SEC = DM * 60D0

      IF(SEC >= 59.999999d0) THEN
       SEC = 0d0
       IM = IM + 1
      END IF

      IF (IM == 60) THEN
       IM = 0
       IH = IH +1
      END IF 
     
      END 


*==============================================================================================

      SUBROUTINE TERCEL ( TJDH, TJDL, XP, YP, VEC1, VEC2 )
*
*     THIS SUBROUTINE ROTATES A VECTOR FROM THE TERRESTRIAL TO THE
*     CELESTIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
*     ITRS (A ROTATING EARTH-FIXED SYSTEM) TO THE GCRS (A LOCAL
*     SPACE-FIXED SYSTEM) BY APPLYING ROTATIONS FOR POLAR MOTION,
*     EARTH ROTATION, NUTATION, PRECESSION, AND THE DYNAMICAL-TO-GCRS
*     FRAME TIE.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
*                   SYSTEM) (IN)
*          VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
*                   SYSTEM) (OUT)
*
*     NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
*
*     NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*     NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
*     (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
*     SPECIAL OPTION WHERE THE OUTPUT VECTOR IS PRODUCED WITH RESPECT
*     TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH THE
*     TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
*     WORKS ONLY IN 'EQUINOX' MODE.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE XPOLE, YPOLE IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*

! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

      DOUBLE PRECISION TJDH,TJDL,XP,YP,VEC1,VEC2,
     .     T0,DELTAT,UTJDH,UTJDL,UTJD,TTJD,TDBJD,T,SECDIF,
     .     GAST,RCIO,THETA,V1,V2,V3,V4,X,Y,Z
      DIMENSION VEC1(3), VEC2(3), V1(3), V2(3), V3(3), V4(3),
     .     X(3), Y(3), Z(3)

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

      CALL GETDT ( DELTAT )
       
      IF ( TJDH .GE. 0.D0 ) THEN
          UTJDH = TJDH
          UTJDL = TJDL
      ELSE
          UTJDH = TJDL
          UTJDL = 0.D0
      END IF
      UTJD = UTJDH + UTJDL

*     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0
      
*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MODE .GE. 2 ) THEN
*         'EQUINOX' MODE

*         APPLY POLAR MOTION
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF
          
*         APPLY EARTH ROTATION
          CALL SIDTIM ( UTJDH, UTJDL, 1,   GAST )
          CALL SPIN ( -GAST *15D0, V1,   V2 )

*         SPECIAL OPTION SKIPS REMAINING TRANSFORMATIONS
          IF ( TJDH .LT. 0.D0 ) THEN
              VEC2(1) = V2(1)
              VEC2(2) = V2(2)
              VEC2(3) = V2(3)
          ELSE

*         APPLY NUTATION AND PRECESSION
          CALL NUTATE ( -TDBJD, V2,   V3 )

          CALL PRECES ( TDBJD, V3, T0,   V4 )

*         APPLY FRAME-TIE MATRIX
          CALL FRAME ( V4, -1, VEC2 )
          
          END IF

      ELSE
*         'CIO-TIO-THETA' MODE
*         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
*         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
*         EQ. (3) AND (4).

*         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE TERRESTRIAL
*         INTERMEDIATE SYSTEM
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF

*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )

*         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
*         THE VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
          CALL EROT ( UTJDH, UTJDL,   THETA )
          CALL SPIN ( -THETA, V1,   V2 )

*         TRANSFORM THE VECTOR FROM THE CELESTIAL INTERMEDIATE SYSTEM
*         TO THE GCRS
          VEC2(1) = X(1) * V2(1) + Y(1) * V2(2) + Z(1) * V2(3)
          VEC2(2) = X(2) * V2(1) + Y(2) * V2(2) + Z(2) * V2(3)
          VEC2(3) = X(3) * V2(1) + Y(3) * V2(2) + Z(3) * V2(3)

      END IF

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
  50  CALL SETVEC ( VEC2 )

      RETURN

      END
*=================================================================================================

      SUBROUTINE SETDT ( DELT )
*
*     THIS SUBROUTINE ALLOWS FOR THE SPECIFICATION OF THE DELTA-T VALUE
*     (DELTA-T = TT - UT1) TO BE USED IN THE CALCULATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS
*     THESE CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE
*     TIME ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME
*     ARGUMENT FOR THE PRECESSION AND NUTATION COMPONENTS.  THIS
*     SUBROUTINE, IF USED, SHOULD BE CALLED BEFORE ANY SUBROUTINE
*     RELATED TO EARTH ROTATION (E.G., SIDTIM OR TERCEL) FOR A GIVEN
*     DATE.  THE VALUE OF DELTA-T SPECIFIED HERE WILL BE USED UNTIL
*     EXPLICITLY CHANGED.
*
*          DELT   = VALUE OF DELTA-T (TT-UT1) IN SECONDS (IN)
*
*     NOTE 1:  THE COMPUTED VALUE OF SIDEREAL TIME, AND THE EQUIVALENT
*     EARTH ORIENTATION ANGLES, ARE RELATIVELY INSENSITIVE TO THE VALUE
*     OF DELTA-T: UP TO ONLY ~3 MICROARCSECONDS PER SECOND OF DELTA-T.
*     THEREFORE, FOR MANY APPLICATIONS, THIS SUBROUTINE EITHER NEED NOT
*     BE CALLED AT ALL, OR CAN BE CALLED JUST ONCE FOR A WIDE RANGE OF
*     DATES (E.G., A YEAR).  IF THIS CALL IS NOT USED, A DEFAULT
*     VALUE OF DELTA-T OF 64 SECONDS IS USED, WHICH IS APPROPRIATE TO
*     2000.0.
*
*     NOTE 2:  THE INPUT TIME ARGUMENTS TO SIDTIM AND TERCEL (TJDH AND
*     TJDL) ARE EXPRESSED IN UT1 REGARDLESS OF WHETHER THIS CALL IS
*     USED.
*
*
      DOUBLE PRECISION DELTAT, DT, DELT
      SAVE DT

*     DEFAULT VALUE OF DELTA-T IN DAYS, EQUIVALENT TO 64 SECONDS,
*     THE APPROXIMATE VALUE AT 2000.0
*      DATA DT / 0.00074074D0 /

      DT = DELT / 86400.D0

      RETURN


      ENTRY GETDT ( DELTAT )

*     THIS ENTRY RETURNS THE CURRENT VALUE OF DELTA-T
*     (DELTA-T = TT - UT1), PREVIOUSLY SET BY THE USER.  THE VALUE
*     RETURNED IS TO BE USED IN THE CALCULATION OF SIDEREAL TIME AND
*     THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS THESE
*     CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE TIME
*     ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME ARGUMENT
*     FOR THE PRECESSION AND NUTATION COMPONENTS.
*
*          DELTAT = VALUE OF DELTA-T (TT-UT1) IN DAYS (OUT)


      DELTAT = DT

      RETURN

      END

*=====================================================================================================

      SUBROUTINE ZDAZ ( UJD, XP, YP, GLON, GLAT, HT, RA, DEC, IREFR,
     .                   ZD, AZ, RAR, DECR )
*
*     THIS SUBROUTINE TRANSFORMS TOPOCENTRIC RIGHT ASCENSION AND
*     DECLINATION TO ZENITH DISTANCE AND AZIMUTH.  THIS ROUTINE USES
*     A METHOD THAT PROPERLY ACCOUNTS FOR POLAR MOTION, WHICH IS
*     SIGNIFICANT AT THE SUB-ARCSECOND LEVEL.  THIS SUBROUTINE
*     CAN ALSO ADJUST COORDINATES FOR ATMOSPHERIC REFRACTION.
*
*          UJD    = UT1 JULIAN DATE(IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE (IN)
*          DEC    = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE (IN)
*          IREFR  = ATMOSPHERIC REFRACTION OPTION (IN)
*                   SET IREFR=0 FOR NO REFRACTION
*                   SET IREFR=1 TO INCLUDE REFRACTION
*          ZD     = TOPOCENTRIC ZENITH DISTANCE IN DEGREES,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          AZ     = TOPOCENTRIC AZIMUTH (MEASURED EAST FROM NORTH)
*                   IN DEGREES (OUT)
*          RAR    = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE, AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          DECR   = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*
*     NOTE 1:  XP AND YP CAN BE SET TO ZERO IF SUB-ARCSECOND ACCURACY IS
*     NOT NEEDED.  HT IS USED ONLY FOR REFRACTION, IF IREFR=1.  RA AND
*     DEC CAN BE OBTAINED FROM TPSTAR, TPPLAN, OR PLACE.
*
*     NOTE 2:  THE DIRECTONS ZD=0 (ZENITH) AND AZ=0 (NORTH) ARE
*     HERE CONSIDERED FIXED IN THE TERRESTRIAL SYSTEM.  SPECIFICALLY,
*     THE ZENITH IS ALONG THE GEODETIC NORMAL, AND NORTH IS TOWARD
*     THE ITRS POLE.
*
*     NOTE 3:  IF IREFR=0, THEN RAR=RA AND DECR=DEC.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION UJD,XP,YP,GLON,GLAT,HT,RA,DEC,ZD,AZ,RAR,DECR,
     .     PI,DEGRAD,RADDEG,
     .     SINLAT,COSLAT,SINLON,COSLON,SINDC,COSDC,SINRA,COSRA,
     .     UZE,UNE,UWE,UZ,UN,UW,P,PR,PZ,PN,PW,PROJ,
     .     ZD0,ZD1,REFR,SINZD,COSZD,SINZD0,COSZD0,
     .     DSIN,DCOS,DSQRT,DATAN2
      DIMENSION UZE(3), UNE(3), UWE(3), UZ(3), UN(3), UW(3),
     .     P(3), PR(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )
      PARAMETER ( RADDEG = 180.D0 / PI           )
      
      RAR    = RA 
      DECR   = DEC 
      SINLAT = DSIN ( GLAT * DEGRAD )
      COSLAT = DCOS ( GLAT * DEGRAD )
      SINLON = DSIN ( GLON * DEGRAD )
      COSLON = DCOS ( GLON * DEGRAD )
      SINDC  = DSIN ( DEC * DEGRAD )
      COSDC  = DCOS ( DEC * DEGRAD )
      SINRA  = DSIN ( RA * 15.D0 * DEGRAD )
      COSRA  = DCOS ( RA * 15.D0 * DEGRAD )
      
* --- SET UP ORTHONORMAL BASIS VECTORS IN LOCAL EARTH-FIXED SYSTEM ----

*     DEFINE VECTOR TOWARD LOCAL ZENITH IN EARTH-FIXED SYSTEM (Z AXIS)
      UZE(1) =  COSLAT * COSLON
      UZE(2) =  COSLAT * SINLON
      UZE(3) =  SINLAT
      
*     DEFINE VECTOR TOWARD LOCAL NORTH IN EARTH-FIXED SYSTEM (X AXIS)
      UNE(1) = -SINLAT * COSLON
      UNE(2) = -SINLAT * SINLON
      UNE(3) =  COSLAT
      
*     DEFINE VECTOR TOWARD LOCAL WEST IN EARTH-FIXED SYSTEM (Y AXIS)
      UWE(1) =  SINLON
      UWE(2) = -COSLON
      UWE(3) =  0.D0
      
* --- OBTAIN VECTORS IN CELESTIAL SYSTEM ------------------------------

*     ROTATE EARTH-FIXED ORTHONORMAL BASIS VECTORS TO CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)

      CALL EQINOX
!      CALL CIOTIO
      CALL HIACC
      CALL TERCEL ( -1D0, UJD, XP, YP, UZE,   UZ )
      CALL TERCEL ( -1D0, UJD, XP, YP, UNE,   UN )
      CALL TERCEL ( -1D0, UJD, XP, YP, UWE,   UW )
      CALL RESUME
      
     
*     DEFINE UNIT VECTOR P TOWARD OBJECT IN CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)
      P(1) = COSDC * COSRA
      P(2) = COSDC * SINRA
      P(3) = SINDC
* --- COMPUTE COORDINATES OF OBJECT WRT ORTHONORMAL BASIS -------------

*     COMPUTE COMPONENTS OF P -- PROJECTIONS OF P ONTO ROTATED
*     EARTH-FIXED BASIS VECTORS
      PZ = P(1) * UZ(1) + P(2) * UZ(2) + P(3) * UZ(3)
      PN = P(1) * UN(1) + P(2) * UN(2) + P(3) * UN(3)
      PW = P(1) * UW(1) + P(2) * UW(2) + P(3) * UW(3)
      
*     COMPUTE AZIMUTH AND ZENITH DISTANCE
      PROJ = DSQRT ( PN**2 + PW**2 )
      AZ = 0.D0
      IF ( PROJ .GT. 0.D0 ) AZ = -DATAN2 ( PW, PN ) * RADDEG
      
      IF ( AZ .LT.   0.D0 ) AZ = AZ + 360.D0
      IF ( AZ .GE. 360.D0 ) AZ = AZ - 360.D0
      ZD = DATAN2 ( PROJ, PZ ) * RADDEG
      
* --- APPLY ATMOSPHERIC REFRACTION IF REQUESTED -----------------------

      IF ( IREFR .EQ. 1 ) THEN

*         GET REFRACTION IN ZENITH DISTANCE
*         ITERATIVE PROCESS REQUIRED BECAUSE REFRACTION ALGORITHMS ARE
*         ALWAYS A FUNCTION OF OBSERVED (NOT COMPUTED) ZENITH DISTANCE
           
          ZD0 = ZD
  40      ZD1 = ZD
          CALL REFRAC ( HT, ZD,   REFR )
          ZD = ZD0 - REFR
          
*         REQUIRE CONVERGENCE TO 0.1 ARCSEC (ACTUAL ACCURACY LESS)
          IF ( DABS ( ZD - ZD1 ) .GT. 3.D-5 ) GO TO 40

*         APPLY REFRACTION TO CELESTIAL COORDINATES OF OBJECT
          IF ( REFR .GT. 0.D0 .AND. ZD .GT. 3.D-4 ) THEN
             
*             SHIFT POSITION VECTOR OF OBJECT IN CELESTIAL SYSTEM
*             TO ACCOUNT FOR FOR REFRACTION (SEE USNO/AA TECHNICAL
*             NOTE 1998-09)
              SINZD  = DSIN ( ZD * DEGRAD )
              COSZD  = DCOS ( ZD * DEGRAD ) 
              SINZD0 = DSIN ( ZD0 * DEGRAD )
              COSZD0 = DCOS ( ZD0 * DEGRAD )
*             COMPUTE REFRACTED POSITION VECTOR
              DO  J = 1, 3
              PR(J) = ( ( P(J) - COSZD0 * UZ(J) ) / SINZD0 ) * SINZD
     .                 +                  UZ(J)              * COSZD   
              END DO
*             COMPUTE REFRACTED RIGHT ASCENSION AND DECLINATION
              PROJ = DSQRT ( PR(1)**2 + PR(2)**2 )
              RAR = 0.D0
              IF ( PROJ .GT. 0.D0 ) RAR = DATAN2 ( PR(2), PR(1) )
     .                                   *  RADDEG / 15.D0
       
              IF ( RAR .LT.  0.D0 ) RAR = RAR + 24.D0
              IF ( RAR .GE. 24.D0 ) RAR = RAR - 24.D0

              DECR = DATAN2 ( PR(3), PROJ )  * RADDEG
               
          END IF

      END IF

      
      END


*===================================================================================================

      SUBROUTINE SIDTIM ( TJDH, TJDL, K,   GST )
*
*     THIS SUBROUTINE COMPUTES THE GREENWICH SIDEREAL TIME
*     (EITHER MEAN OR APPARENT) AT JULIAN DATE TJDH + TJDL.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          K      = TIME SELECTION CODE (IN)
*                   SET K=0 FOR GREENWICH MEAN SIDEREAL TIME
*                   SET K=1 FOR GREENWICH APPARENT SIDEREAL TIME
*          GST    = GREENWICH (MEAN OR APPARENT) SIDEREAL TIME
*                   IN HOURS (OUT)
*
*     NOTE:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*
      DOUBLE PRECISION TJDH,TJDL,GST,PI,DEGCON,DELTAT,
     .     T0,UTJD,TTJD,TDBJD,SECDIF,A,THETA,RCIO,
     .     UNITX,W1,W2,X,Y,Z,EQ,HAEQ,EE,DMOD,DATAN2
      DIMENSION UNITX(3), W1(3), W2(3), X(3), Y(3), Z(3), EQ(3)
      SAVE 

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGCON = 180.D0 / PI           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA UNITX / 1.D0, 0.D0, 0.D0 /

   3  FORMAT ( ' SIDTIM ERROR: CANNOT RETURN SIDEREAL TIME FOR ',
     .     'JD ', F10.1 )

      CALL GETDT ( DELTAT )
      
*     TIME ARGUMENT FOR PRECESSION AND NUTATION COMPONENTS OF SIDEREAL
*     TIME IS TDB
      UTJD = TJDH + TJDL
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, A,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MODE .GE. 2 ) THEN
*          EQUINOX-BASED MODE
*          SEE USNO CIRCULAR 179, SECTION 2.6.2   
     
*         GET -1 TIMES THE MEAN OR TRUE RIGHT ASCENSION OF THE CIO     
          CALL EQXRA ( TDBJD, K,   RCIO )
*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         COMBINE TO OBTAIN SIDEREAL TIME       
          GST = DMOD ( THETA / 15.D0 - RCIO, 24.D0 )
          IF ( GST .LT. 0.D0 ) GST = GST + 24.D0

      ELSE
*         CIO-BASED MODE
*         SEE USNO CIRCULAR 179, SECTION 6.5.4

*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          IF ( RCIO .EQ. 99.D0 ) THEN
              WRITE ( *, 3 ) TDBJD
              GST = 99.D0
              GO TO 50
          END IF
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
*         COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
          CALL NUTATE ( -TDBJD, UNITX,   W1 )
          CALL PRECES ( TDBJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,    EQ )
*         COMPUTE THE HOUR ANGLE OF THE EQUINOX WRT THE TIO MERIDIAN
*         (NEAR GREENWICH, BUT PASSES THROUGH THE CIP AND TIO)
          HAEQ = THETA - DATAN2 ( EQ(1)*Y(1) + EQ(2)*Y(2) + EQ(3)*Y(3),
     .                            EQ(1)*X(1) + EQ(2)*X(2) + EQ(3)*X(3) )
     .                   * DEGCON

*         FOR MEAN SIDEREAL TIME, OBTAIN THE EQUATION OF THE EQUINOXES
*         AND SUBTRACT IT
          IF ( K .EQ. 0 ) THEN
              CALL ETILT ( TDBJD,   A, A, EE, A, A )
              HAEQ = HAEQ - EE / 240.D0
          END IF

          HAEQ = DMOD ( HAEQ, 360.D0 ) / 15.D0
          IF ( HAEQ .LT. 0.D0 ) HAEQ = HAEQ + 24.D0
          GST = HAEQ

      END IF

  50  RETURN

      END

*===================================================================================================
      SUBROUTINE CIOLOC ( TJD,   RACIO, K )
*
*     THIS SUBROUTINE RETURNS THE LOCATION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO) FOR A GIVEN JULIAN DATE, AS A
*     RIGHT ASCENSION WITH RESPECT TO EITHER THE GCRS (GEOCENTRIC ICRS)
*     ORIGIN OR THE TRUE EQUINOX OF DATE.  THE CIO IS ALWAYS LOCATED ON
*     THE TRUE EQUATOR (=INTERMEDIATE EQUATOR) OF DATE.
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (OUT)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   GIVEN (OUT)
*                   K=1 MEANS GCRS
*                   K=2 MEANS TRUE EQUATOR AND EQUINOX OF DATE
*
*     NOTE:  IF AN EXTERNAL FILE OF CIO RIGHT ASCENSIONS IS AVAILABLE,
*     IT WILL BE USED AND K WILL BE SET TO 1.  OTHERWISE AN INTERNAL
*     COMPUTATION WILL BE USED AND K WILL BE SET TO 2.
*
*
      DOUBLE PRECISION TJD,RACIO,A,TLAST,RLAST,JD,RA,P,EQOR,DABS
      LOGICAL USEFIL
      DIMENSION JD(8), RA(8), A(1)
      SAVE

*     NUMBER OF VALUES IN ARRAYS FOR LAGRANGIAN INTERPOLATION 
      DATA M / 6 /
      
      DATA TLAST, RLAST, KLAST / 0.D0, 0.D0, 0 /
      
   3  FORMAT ( ' CIOLOC ERROR: CANNOT RETURN CIO RA VALUE FOR JD ',
     .     F10.1 )

*     CHECK IF EXTERNAL FILE EXISTS
      CALL CIORD ( 0.D0, 1,   A, A, IERR )
      USEFIL = IERR .EQ. 0

*     CHECK IF PREVIOUSLY COMPUTED RA VALUE CAN BE USED
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 ) THEN
          RACIO = RLAST
          K = KLAST
          GO TO 77
      END IF

* --- IF EXTERNAL FILE EXISTS, INTERPOLATE RA VALUES FROM FILE --------

      IF ( USEFIL ) THEN
      
          K = 1

*         GET ARRAYS OF VALUES TO INTERPOLATE
          CALL CIORD ( TJD, M,   JD, RA, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) TJD
              RACIO = 99.D0
              GO TO 77
          END IF

*         PERFORM LAGRANGIAN INTERPOLATION FOR RA AT TJD
          RACIO = 0.D0
          DO 40 J = 1, M
              P = 1.D0
              DO 30 I = 1, M
                  IF ( I .EQ. J ) GO TO 30
                  P = P * ( TJD - JD(I) ) / ( JD(J) - JD(I) )
  30          CONTINUE
              RACIO = RACIO + P * RA(J)
  40      CONTINUE

          RACIO = RACIO / 54000.D0
 
* --- OTHERWISE, USE INTERNAL COMPUTATION ----------------------------

      ELSE

          K = 2

*         GET EQUATION OF THE ORIGINS
          CALL EQXRA ( TJD, 1,   EQOR )
          
          RACIO = -EQOR
 
      END IF

* ---------------------------------------------------------------------

      TLAST = TJD
      RLAST = RACIO
      KLAST = K

  77  RETURN

      END

*============================================================================================
      SUBROUTINE CIOBAS ( TJD, RACIO, K,   X, Y, Z )
*
*     THIS SUBROUTINE RETURNS THE ORTHONORMAL BASIS VECTORS, WITH
*     RESPECT TO THE GCRS (GEOCENTRIC ICRS), OF THE CELESTIAL
*     INTERMEDIATE SYSTEM DEFINED BY THE CELESTIAL INTERMEDIATE POLE
*     (CIP) (IN THE Z DIRECTION) AND THE CELESTIAL INTERMEDIATE ORIGIN
*     (CIO) (IN THE X DIRECTION).  A TDB JULIAN DATE AND THE RIGHT
*     ASCENSION OF THE CIO AT THAT DATE IS REQUIRED AS INPUT.  THE
*     RIGHT ASCENSION OF THE CIO CAN BE WITH RESPECT TO EITHER THE
*     GCRS ORIGIN OR THE TRUE EQUINOX OF DATE -- DIFFERENT ALGORITHMS
*     ARE USED IN THE TWO CASES. 
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (IN)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   EXPRESSED (IN)
*                   SET K=1 FOR GCRS
*                   SET K=2 FOR TRUE EQUATOR AND EQUINOX OF DATE
*          X      = UNIT VECTOR TOWARD THE CIO, EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO THE GCRS (OUT)
*          Y      = UNIT VECTOR TOWARD THE Y-DIRECTION, EQUATORIAL
*                   RECTANGULAR COORDINATES, REFERRED TO THE GCRS (OUT)
*          Z      = UNIT VECTOR TOWARD NORTH CELESTIAL POLE (CIP),
*                   EQUATORIAL RECTANGULAR COORDINATES, REFERRED TO
*                   THE GCRS (OUT)
*
*
      DOUBLE PRECISION TJD,RACIO,X,Y,Z,XX,YY,ZZ,W0,W1,W2,Z0,
     .     PI,RADCON,T0,TLAST,SINRA,COSRA,XMAG,DABS,DSIN,DCOS,DSQRT
      DIMENSION X(3), Y(3), Z(3), XX(3), YY(3), ZZ(3), Z0(3),
     .     W0(3), W1(3), W2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA Z0 / 0.D0, 0.D0, 1.D0 /,   TLAST / 0.D0 /,   KLAST / 0 /
      
   3  FORMAT ( ' CIOBAS ERROR: INVALID VALUE FOR K FOR JD ',
     .     F10.1 )      

*     USE LAST-COMPUTED BASIS VECTORS IF POSSIBLE
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 .AND. K .EQ. KLAST )
     .   GO TO 60

*     COMPUTE UNIT VECTOR Z TOWARD CELESTIAL POLE (CIP)      
      CALL NUTATE ( -TJD, Z0,   W1 )
      CALL PRECES ( TJD, W1, T0,   W2 )
      CALL FRAME ( W2, -1,   ZZ )
      
* --- RA OF CIO EXPRESSED IN GCRS -------------------------------------      
      
      IF ( K .EQ. 1 ) THEN

*         COMPUTE VECTOR X TOWARD CIO IN GCRS
          SINRA = DSIN ( RACIO * 15.D0 * RADCON )
          COSRA = DCOS ( RACIO * 15.D0 * RADCON )
          XX(1) =  ZZ(3) * COSRA
          XX(2) =  ZZ(3) * SINRA
          XX(3) = -ZZ(1) * COSRA - ZZ(2) * SINRA

*         NORMALIZE VECTOR X 
          XMAG = DSQRT ( XX(1)**2 + XX(2)**2 + XX(3)**2 )
          XX(1) = XX(1) / XMAG
          XX(2) = XX(2) / XMAG
          XX(3) = XX(3) / XMAG

*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
          
* --- RA OF CIO EXPRESSED IN EQUATOR-AND-EQUINOX OF DATE SYSTEM -------          
      
      ELSE IF ( K .EQ. 2 ) THEN
      
*         CONSTRUCT UNIT VECTOR TOWARD CIO
*         IN EQUATOR-AND-EQUINOX-OF-DATE SYSTEM
          W0(1) = DCOS ( RACIO * 15.D0 * RADCON )
          W0(2) = DSIN ( RACIO * 15.D0 * RADCON )
          W0(3) = 0.D0
       
*         ROTATE THE VECTOR INTO THE GCRS TO FORM UNIT VECTOR X      
          CALL NUTATE ( -TJD, W0,   W1 )
          CALL PRECES ( TJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,   XX )
          
*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
      
* ---------------------------------------------------------------------      
      
      ELSE
          
          WRITE ( *, 3 ) TJD
          GO TO 77
      
      END IF   

* ---------------------------------------------------------------------      

      TLAST = TJD
      KLAST = K
      
  60  DO 66 J = 1, 3
          X(J) = XX(J)
          Y(J) = YY(J)
          Z(J) = ZZ(J)
  66  CONTINUE    

  77  RETURN

      END

*==================================================================================================

      SUBROUTINE EROT (DATE1,DATE2,THETA)
*
*     THIS SUBROUTINE RETURNS THE VALUE OF THE EARTH ROTATION ANGLE
*     (THETA) FOR A GIVEN UT1 JULIAN DATE.  THE EXPRESSION USED IS
*     TAKEN FROM THE NOTE TO IAU RESOLUTION B1.8 OF 2000.
*
*         DATE1  = HIGH-ORDER PART OF UT1 JULIAN DATE (IN)
*         DATE2  = LOW-ORDER PART OF UT1 JULIAN DATE (IN)
*         THETA  = EARTH ROTATION ANGLE IN DEGREES (OUT)
*
*
      DOUBLE PRECISION DATE1, DATE2, THETA, T0, THET1, THET2, THET3,
     .     DMOD

      DATA T0 / 2451545.0D0 /

*     THE ALGORITHM USED BELOW IS EQUIVALENT TO THE CANNONICAL
*     THETA = 0.7790572732640D0 + 1.00273781191135448D0 * T,
*     WHERE T IS THE TIME IN UT1 DAYS FROM 2000.0 (T=DATE1+DATE2-T0),
*     BUT IT AVOIDS MANY TWO-PI 'WRAPS' THAT DECREASE PRECISION
*     (ADOPTED FROM SOFA ROUTINE IAU_ERA00 BY PAT WALLACE; SEE ALSO
*     EXPRESSION AT TOP OF PAGE 35 OF IERS CONVENTIONS (1996))

      THET1 = 0.7790572732640D0 + 0.00273781191135448D0 * ( DATE1 - T0 )
      THET2 =                     0.00273781191135448D0 *   DATE2
      THET3 = DMOD ( DATE1, 1.D0 ) + DMOD ( DATE2, 1.D0 )
      THETA = DMOD ( THET1 + THET2 + THET3, 1.D0 ) * 360.D0
      IF ( THETA .LT. 0.D0 ) THETA = THETA + 360.D0

      RETURN

      END
* ===========================================================================================

      SUBROUTINE WOBBLE (TJD,XP,YP,POS1,POS2)
*
*     THIS SUBROUTINE CORRECTS A VECTOR IN THE ITRS (A ROTATING EARTH-
*     FIXED SYSTEM) FOR POLAR MOTION, AND ALSO CORRECTS THE LONGITUDE
*     ORIGIN (BY A TINY AMOUNT) TO THE TERRESTRIAL INTERMEDIATE ORIGIN
*     (TIO).  THE ITRS VECTOR IS THEREBY TRANSFORMED TO THE TERRESTRIAL
*     INTERMEDIATE SYSTEM, BASED ON THE TRUE (ROTATIONAL) EQUATOR AND
*     THE TERRESTRIAL INTERMEDIATE ORIGIN (TIO).  SINCE THE TRUE EQUATOR
*     IS THE PLANE ORTHOGONAL TO THE DIRECTION OF THE CELESTIAL
*     INTERMEDIATE POLE (CIP), THE COMPONENTS OF THE OUTPUT VECTOR ARE
*     REFERRED TO Z AND X AXES TOWARD THE CIP AND TIO, RESPECTIVELY.
*
*          TJD    = TT OR UT1 JULIAN DATE (IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND TIO (OUT)
*
*     NOTE 1:  IF TJD IS NEGATIVE, THE INVERSE TRANSFORMATION (TERRESTRIAL 
*     INTERMEDIATE SYSTEM TO ITRS) IS APPLIED.
*
*     NOTE 2: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION TJD,XP,YP,POS1,POS2,PI,SECCON,T0,T,XPOLE,YPOLE,
     .     SPRIME,TIOLON,SINX,COSX,SINY,COSY,SINL,COSL,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DSIN,DCOS
      DIMENSION POS1(3), POS2(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TT JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.0D0 /

      XPOLE = XP / SECCON
      YPOLE = YP / SECCON

      T = ( DABS(TJD) - T0 ) / 36525.D0

*     COMPUTE APPROXIMATE LONGITUDE OF TIO, USING EQ. (10) OF
*     LAMBERT & BIZOUARD (2002), ASTRONOMY AND ASTROPHYSICS 394,
*     317-321
      SPRIME = -47.0D-6 * T
      TIOLON = -SPRIME / SECCON
*     NOTE THAT TIOLON, THE LONGITUDE CORRECTION, IS NEGLIGIBLE FOR
*     MOST ASTRONOMICAL PURPOSES

*     COMPUTE ELEMENTS OF ROTATION MATRIX
*     EQUIVALENT TO R3(-S')R2(X)R1(Y) AS PER IERS CONVENTIONS (2003)
      SINX = DSIN ( XPOLE )
      COSX = DCOS ( XPOLE )
      SINY = DSIN ( YPOLE )
      COSY = DCOS ( YPOLE )
      SINL = DSIN ( TIOLON )
      COSL = DCOS ( TIOLON )
      XX =  COSX * COSL
      YX =  SINX * SINY * COSL + COSY * SINL
      ZX = -SINX * COSY * COSL + SINY * SINL
      XY = -COSX * SINL
      YY =  SINX * SINY * SINL + COSY * COSL
      ZY =  SINX * COSY * SINL + SINY * COSL
      XZ =  SINX
      YZ = -COSX * SINY
      ZZ =  COSX * COSY
   10 IF ( TJD .LT. 0.D0 ) GO TO 30

*     PERFORM ROTATION FROM ITRS TO TERRESTRIAL INTERMEDIATE SYSTEM
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM TERRESTRIAL INTERMEDIATE SYSTEM TO ITRS
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END
*================================================================================================

      SUBROUTINE SPIN (ANGL,POS1,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM ONE COORDINATE SYSTEM
*     TO ANOTHER WITH SAME ORIGIN AND AXES ROTATED ABOUT THE
*     Z AXIS.
*
*          ANGL   = ANGLE OF COORDINATE SYSTEM ROTATION, POSITIVE
*                   COUNTERCLOCKWISE WHEN VIEWED FROM +Z,
*                   IN DEGREES (IN)
*          POS1   = POSITION VECTOR (IN)
*          POS2   = POSITION VECTOR EXPRESSED IN NEW COORDINATE
*                   SYSTEM ROTATED ABOUT Z BY ANGLE ANG (OUT)
*
*
      DOUBLE PRECISION ANGL,POS1,POS2,PI,ALAST,ANG,COSANG,SINANG,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )

      DATA ALAST / -999.D0 /

      IF ( DABS ( ANGL - ALAST ) .GT. 1.D-12 ) THEN

          ANG = ANGL / 180.D0 * PI
          COSANG = DCOS ( ANG )
          SINANG = DSIN ( ANG )

*         ROTATION MATRIX FOLLOWS
          XX =  COSANG
          YX =  SINANG
          ZX =  0.D0
          XY = -SINANG
          YY =  COSANG
          ZY =  0.D0
          XZ =  0.D0
          YZ =  0.D0
          ZZ =  1.D0

          ALAST = ANGL

      END IF

*     PERFORM ROTATION
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)

      RETURN

      END
*===================================================================================================

      SUBROUTINE TIMES (TDBJD,TTJD,SECDIF)
*
*     THIS SUBROUTINE COMPUTES THE TERRESTRIAL TIME (TT) JULIAN DATE
*     CORRESPONDING TO A BARYCENTRIC DYNAMICAL TIME (TDB) JULIAN DATE.
*     THE EXPRESSION USED IN THIS VERSION IS A TRUNCATED FORM OF A
*     LONGER AND MORE PRECISE SERIES GIVEN BY FAIRHEAD & BRETAGNON
*     (1990) A&A 229, 240.  THE RESULT IS GOOD TO ABOUT 10 MICROSECONDS.
*
*          TDBJD  = TDB JULIAN DATE (IN)
*          TTJD   = TT JULIAN DATE (OUT)
*          SECDIF = DIFFERENCE TDBJD-TTJD, IN SECONDS (OUT)
*
*
      DOUBLE PRECISION TDBJD,TTJD,SECDIF,T,T0,DSIN

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

      T = ( TDBJD - T0 ) / 36525.D0

*     EXPRESSION GIVEN IN USNO CIRCULAR 179, EQ. 2.6 
      SECDIF = 0.001657D0 * DSIN (  628.3076D0 * T + 6.2401D0)                     
     .       + 0.000022D0 * DSIN (  575.3385D0 * T + 4.2970D0)   
     .       + 0.000014D0 * DSIN ( 1256.6152D0 * T + 6.1969D0)            
     .       + 0.000005D0 * DSIN (  606.9777D0 * T + 4.0212D0)                                     
     .       + 0.000005D0 * DSIN (   52.9691D0 * T + 0.4444D0)     
     .       + 0.000002D0 * DSIN (   21.3299D0 * T + 5.5431D0) 
     .       + 0.000010D0 * T * DSIN ( 628.3076D0 * T + 4.2490D0)

      TTJD = TDBJD - SECDIF / 86400.D0

      RETURN

      END
*=======================================================================================================

      SUBROUTINE SETMOD ( MODE )
*
*     THIS SUBROUTINE ALLOWS THE USER TO SPECIFY THE 'MODE' VALUE,
*     WHICH DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (IN)
*                   SET MODE=0 FOR CIO-BASED METHOD, FULL ACCURACY
*                   SET MODE=1 FOR CIO-BASED METHOD, REDUCED ACCURACY
*                   SET MODE=2 FOR EQUINOX-BASED METHOD, FULL ACCURACY
*                   SET MODE=3 FOR EQUINOX-BASED METHOD, REDUCED
*                                  ACCURACY
*
*     NOTE: OTHER ENTRY POINTS ARE PROVIDED TO ALLOW THE METHOD AND
*     ACCURACY TO BE SPECIFIED IN A MORE OBVIOUS WAY:
*     MODE=0 CAN BE SET BY CALL CIOTIO AND CALL HIACC
*     MODE=1 CAN BE SET BY CALL CIOTIO AND CALL LOACC
*     MODE=2 CAN BE SET BY CALL EQINOX AND CALL HIACC
*     MODE=3 CAN BE SET BY CALL EQINOX AND CALL LOACC
*
*
      SAVE IMODE, LMODE

      DATA IMODE, LMODE / 0, 0 /

      LMODE = IMODE
      IMODE = MODE

      RETURN


      ENTRY CIOTIO
      LMODE = IMODE
      IF ( IMODE .GE. 2 ) IMODE = IMODE - 2
      RETURN


      ENTRY EQINOX
      LMODE = IMODE
      IF ( IMODE .LE. 1 ) IMODE = IMODE + 2
      RETURN


      ENTRY LOACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 0 ) IMODE = IMODE + 1
      RETURN


      ENTRY HIACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 1 ) IMODE = IMODE - 1
      RETURN


      ENTRY RESUME
      IMODE = LMODE
      RETURN


      ENTRY GETMOD ( MODE )
*
*     THIS SUBROUTINE RETURNS THE 'MODE' VALUE, WHICH
*     DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (OUT)
*                   MODE=0 MEANS CIO-BASED METHOD, FULL ACCURACY
*                   MODE=1 MEANS CIO-BASED METHOD, REDUCED ACCURACY
*                   MODE=2 MEANS EQUINOX-BASED METHOD, FULL ACCURACY
*                   MODE=3 MEANS EQUINOX-BASED METHOD, REDUCED ACCURACY
*
*
      MODE = IMODE
      
      RETURN

      END

*===========================================================================================

      SUBROUTINE GETVEC ( UNITV )
*
*     THIS SUBROUTINE ALLOWS THE USER TO RETRIEVE THE LAST COMPUTED
*     POSITION ON THE SKY AS A UNIT VECTOR.
*
*          UNITV  = UNIT VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (OUT)
*
*
      DOUBLE PRECISION UNITV, P, POS, R, DSQRT
      DIMENSION UNITV(3), P(3), POS(3)
      SAVE P

      R = DSQRT ( P(1)**2 + P(2)**2 + P(3)**2 )

      DO 20 J = 1, 3
          UNITV(J) = P(J) / R
  20  CONTINUE

      RETURN


      ENTRY SETVEC ( POS )
*
*     THIS ENTRY STORES THE LAST COMPUTED POSITION ON THE SKY.
*
*          POS    = VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (IN)
*
*
      DO 30 J = 1, 3
          P(J) = POS(J)
  30  CONTINUE

      RETURN

      END
*=================================================================================================

      SUBROUTINE REFRAC (HEIGHT,ZDOBS,REFR)
*
*     THIS SUBROUTINE COMPUTES ATMOSPHERIC REFRACTION IN ZENITH
*     DISTANCE.  THIS VERSION COMPUTES APPROXIMATE REFRACTION FOR
*     OPTICAL WAVELENGTHS.  IT CAN BE USED FOR PLANNING OBSERVATIONS
*     OR TELESCOPE POINTING, BUT SHOULD NOT BE USED FOR THE REDUCTION
*     OF PRECISE OBSERVATIONS.  BASIC ALGORITHM IS DESCRIBED IN THE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, P. 144,
*     AND IS AN ADAPTATION OF A FORMULA IN BENNETT (1982), JOURNAL
*     OF NAVIGATION (ROYAL INSTITUTE) 35, 255-259.
*
*          HEIGHT = HEIGHT OF OBSERVER IN METERS (IN)
*          ZDOBS  = OBSERVED ZENITH DISTANCE IN DEGREES (IN)
*          REFR   = ATMOSPHERIC REFRACTION IN DEGREES (OUT)
*
*     NOTE:  HEIGHT IS NOT USED IF ENTRY REFDAT HAS BEEN CALLED
*     TO SPECIFY ATMOSPHERIC PRESSURE.
*
*
      DOUBLE PRECISION HEIGHT,ZDOBS,REFR,PI,DEGRAD,S,
     .     POBS,TOBS,DOBS,WLOBS,OBSP,OBST,OBSD,OBSWL,P,T,D,WL,H,R,
     .     DEXP,DTAN
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )

*     S IS APPROXIMATE SCALE HEIGHT OF ATMOSPHERE IN METERS
      DATA S / 9.1D3 /
      DATA POBS, TOBS, DOBS, WLOBS / 4 * -999.D0 /

*     COMPUTE REFRACTION ONLY FOR ZENITH DISTANCES
*     BETWEEN 0.1 AND 91 DEGREES
      IF ( ZDOBS .LT. 0.1D0 .OR. ZDOBS .GT. 91.D0 ) THEN
          REFR = 0.D0
          GO TO 77
      END IF

*     IF OBSERVED WEATHER DATA ARE AVAILABLE, USE THEM
*     OTHERWISE, USE CRUDE ESTIMATES OF AVERAGE CONDITIONS
      IF ( POBS .GE. 1.D0 .AND. TOBS .GT. -100.D0 ) THEN
          P  = POBS
          T  = TOBS
          D  = DOBS
          WL = WLOBS
      ELSE
          P  = 1010.D0 * DEXP ( -HEIGHT / S )
          T  = 10.D0
          D  =  0.D0
          WL =  0.5D0
      
      END IF
*     D AND WL NOT USED IN THIS VERSION

      H = 90.D0 - ZDOBS
      R = 0.016667D0 / DTAN ( ( H +  7.31D0 / ( H + 4.4D0 ) ) * DEGRAD )
      REFR = R * ( 0.28D0 * P / ( T + 273.D0 ) )
      
  77  RETURN


      ENTRY REFDAT (OBSP,OBST,OBSD,OBSWL)
*
*     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF WEATHER OBSERVATIONS
*     AND OTHER DATA TO BE USED IN THE ATMOSPHERIC REFRACTION
*     CALCULATION.  THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE
*     SUBROUTINE REFRAC OR ZDAZ FOR A GIVEN DATE/TIME.  DATA SPECIFIED
*     VIA A CALL TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
*
*          OBSP   = OBSERVED ATMOSPHERIC PRESSURE IN MILLIBARS (IN)
*          OBST   = OBSERVED TEMPERATURE IN DEGREES CELSIUS (IN)
*          OBSD   = OBSERVED DEW POINT IN DEGREES CELSIUS (IN)
*          OBSWL  = OBSERVING WAVELENGTH IN MICRONS (IN)
*
*     NOTE:  OBSD AND OBSWL ARE NOT USED IN THIS VERSION'S REFRACTION
*     ALGORITHM, AND CAN BE SET TO ANY VALUE.
*
*
      POBS  = OBSP
      TOBS  = OBST
      DOBS  = OBSD
      WLOBS = OBSWL
      RETURN

      END

*=================================================================================================

      SUBROUTINE EQXRA ( TJD, K,    RAEQ )
*
*     THIS SUBROUTINE COMPUTES THE INTERMEDIATE RIGHT ASCENSION
*     OF THE EQUINOX AT JULIAN DATE TJD, USING AN ANALYTICAL EXPRESSION
*     FOR THE ACCUMULATED PRECESSION IN RIGHT ASCENSION.  FOR THE
*     TRUE EQUINOX THE RESULT IS THE EQUATION OF THE ORIGINS. 
*
*          TJD    = TDB JULIAN DATE (IN)
*          K      = EQUINOX SELECTION CODE (IN)
*                   SET K=0 FOR MEAN EQUINOX
*                   SET K=1 FOR TRUE EQUINOX (EQUATION OF THE ORIGINS)
*          RADIF  = INTERMEDIATE RIGHT ASCENSION OF THE EQUINOX,
*                   IN HOURS (+ OR -) (OUT)
*
*
      DOUBLE PRECISION TJD,RAEQ,T0,TLAST,EE,EQEQ,T,A,PRECRA,DABS
      SAVE 

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   EE / 0.D0 / 
      
      T = ( TJD - T0 ) / 36525.D0	

*     FOR THE TRUE EQUINOX, OBTAIN THE EQUATION OF THE EQUINOXES IN
*     TIME SECONDS, WHICH INCLUDES THE 'COMPLIMENTARY TERMS'
      IF ( K .EQ. 1 ) THEN
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( TJD,   A, A, EE, A, A )
              TLAST = TJD
          END IF
          EQEQ = EE 
      ELSE
          EQEQ = 0.D0
      END IF
      
*     PRECESSION IN RA IN ARCSECONDS TAKEN FROM CAPITAINE ET AL. (2003),
*     ASTRONOMY AND ASTROPHYSICS 412, 567-586, EQ. (42)
      PRECRA = 0.014506D0 +
     .         ( ( ( ( -    0.0000000368D0   * T
     .                 -    0.000029956D0  ) * T
     .                 -    0.00000044D0   ) * T
     .                 +    1.3915817D0    ) * T
     .                 + 4612.156534D0     ) * T

      RAEQ = - ( PRECRA / 15.D0 + EQEQ ) / 3600.D0

      RETURN
      
      END

*=================================================================================================

      SUBROUTINE CIORD ( TJD, NVALS,   TLIST, RALIST, IERR )
*
*     GIVEN AN INPUT TDB JULIAN DATE AND THE NUMBER OF DATA POINTS
*     DESIRED, THIS SUBROUTINE RETURNS A SET OF JULIAN DATES AND
*     CORRESPONDING VALUES OF THE GCRS RIGHT ASCENSION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO).  THE RANGE OF DATES IS CENTERED (AT LEAST
*     APPROXIMATELY) ON THE REQUESTED DATE.  THE SUBROUTINE OBTAINS THE
*     DATA FROM AN EXTERNAL DATA FILE.
*
*         TJD    = TDB JULIAN DATE (IN)
*         NVALS  = NUMBER OF JULIAN DATES AND RIGHT ASCENSION VALUES
*                  REQUESTED (NOT LESS THAN 2 OR MORE THAN 20) (IN)
*         TLIST  = ARRAY OF TDB JULIAN DATES (OUT)
*         RALIST = ARRAY OF GCRS RIGHT ASCENSIONS OF THE CIO, FOR THE
*                  JULIAN DATES IN TLIST, IN ARCSECONDS (OUT)
*         IERR   = ERROR INDICATOR (OUT)
*                  IERR=0 MEANS EVERYTHING OK
*                  IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
*                  IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
*                  IERR=3 MEANS BAD VALUE OF NVALS
*                  IERR=4 MEANS EXTERNAL FILE CANNOT BE FOUND
*
*     NOTE:  TJD=0.D0 WITH NVALS=1 INDICATES A SPECIAL CALL JUST TO
*     DETERMINE IF EXTERNAL FILE EXISTS.
*
*
      DOUBLE PRECISION TJD,TLIST,RALIST,T,T1,R,TBEG,TEND,TINT,DIF
      CHARACTER FILNAM*40, FILEID*(*)
      LOGICAL FILEOK
      DIMENSION TLIST(NVALS), RALIST(NVALS), T(20), R(20)
      SAVE

*     LOGICAL UNIT NUMBER AND FILE ID OF TIME SERIES OF CIO RA VALUES
      DATA LU, ITYP / 24, 1 /
      DATA FILNAM / 'CIO_RA.TXT                              ' /

      DATA NTIMES, TBEG, TEND, FILEOK / 0, 0.D0, 1.D10, .FALSE. /

   1  FORMAT ( A )
   2  FORMAT ( ' CIORD ERROR: CANNOT FIND FILE ', A )
   3  FORMAT ( ' CIORD ERROR: REQUESTED JD ', F10.1, 1X, A,
     .      ' ALLOWED JD ', F10.1 )
   4  FORMAT ( F16.6, F24.14 )

*     SPECIAL CALL JUST TO DETERMINE IF FILE EXITS
*     (NO PRINTED ERROR MESSAGE IF NOT)
      IF ( TJD .EQ. 0.D0 .AND. NVALS .EQ. 1 ) THEN
          IERR = 4
          IF ( NTIMES .EQ. 0 ) INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
          IF ( FILEOK ) IERR = 0
          GO TO 77
      END IF

*     IF EXTERNAL FILE IS ALREADY KNOWN NOT TO EXIST, IMMEDIATELY
*     EXIT WITH IERR=4 
      IF ( NTIMES .GT. 0 .AND. .NOT. FILEOK ) THEN  
          WRITE ( *, 2 ) FILNAM
          IERR = 4
          GO TO 77
      END IF 

*     CHECK FOR VALID VALUE OF NVALS
      IF ( NVALS .LT. 2 .OR. NVALS .GT. 20 ) THEN
          IERR = 3
          GO TO 77
      END IF

      MIDDL = NVALS / 2
   
*     CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF FILE
  10  IF ( TJD .LT. TBEG ) THEN
          WRITE ( *, 3 )  TJD, 'BEFORE FIRST', TBEG
          IERR = 1
          GO TO 77
      ELSE IF ( TJD .GT. TEND ) THEN
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          IERR = 2
          GO TO 77
      END IF

      IF ( ITYP .EQ. 1 ) THEN

*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM FORMATTED
*         SEQUENTIAL INPUT FILE

*         EACH RECORD OF THE FILE MUST CONTAIN A TDB JULIAN DATE
*         AND A CORRESPONDING CIO RA VALUE (WRT GCRS) IN ARCSECONDS

*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL
*         -------------------------------------------------------------

*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='FORMATTED',
     .               ACCESS='SEQUENTIAL', STATUS='OLD' )
              READ ( UNIT=LU, FMT=1 )
              DO 19 I = 1, NVALS
                  READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  19          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              IF ( TJD .LT. TBEG ) GO TO 10
          END IF

*         -------------------------------------------------------------

*         FILE READ SEQUENCE

  20      DIF = ( TJD - T(MIDDL) ) / TINT
          NDIF = DIF

*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( DIF .LT. 0.D0 ) THEN
*             REQUESTED JULIAN DATE BEFORE DATA PREVIOUSLY READ
              IREC = ( T(NVALS) - TBEG ) / TINT
              NBACK = 3d0 * NVALS
              IF ( -DIF .LE. 2 * NVALS .AND. IREC .GT. NBACK ) GO TO 34
              GO TO 25
          ELSE IF ( NDIF .GT. ( NVALS + 1 ) ) THEN
*             REQUESTED JULIAN DATE FAR AHEAD OF DATA PREVIOUSLY READ
              NSKIP = NDIF - NVALS - 1
              GO TO 30
          ELSE
*             REQUESTED JULIAN DATE A BIT AHEAD OF DATA PREVIOUSLY READ
              GO TO 40
          END IF

*         REPOSITION FILE AT BEGINNING
  25      REWIND ( UNIT=LU )
          READ ( UNIT=LU, FMT=1 )
          GO TO 36

*         FAST SKIP FORWARD
  30      DO 32 I = 1, NSKIP
              READ ( UNIT=LU, FMT=1, END=50 )
  32      CONTINUE
          GO TO 36

*         BACKSPACE FILE
  34      DO 35 I = 1, NBACK
              BACKSPACE ( UNIT=LU )
  35      CONTINUE

*         FILL UP ARRAYS
  36      DO 38 I = 1, NVALS
              READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  38      CONTINUE
          GO TO 20

*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  40      DO 44 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  44      CONTINUE
          READ ( UNIT=LU, FMT=4, END=50 ) T(NVALS), R(NVALS)

          GO TO 20

*         -------------------------------------------------------------

*         END OF FILE ENCOUNTERED
  50      BACKSPACE ( UNIT=LU )
          BACKSPACE ( UNIT=LU )
          READ ( UNIT=LU, FMT=4 ) TEND
          TEND = TEND - ( NVALS - MIDDL - 1 ) * TINT
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          T(MIDDL) = TEND + TINT
          IERR = 2
          GO TO 77

      ELSE IF ( ITYP .EQ. 2 ) THEN

*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM UNFORMATTED
*         DIRECT ACCESS INPUT FILE

*         EACH RECORD OF THE FILE (EXCEPT THE FIRST) MUST CONTAIN A
*         TDB JULIAN DATE AND A CORRESPONDING CIO RA VALUE (WRT GCRS)
*         IN ARCSECONDS

*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL

*         THE FIRST RECORD OF THE FILE IS SPECIAL AND MUST CONTAIN THE
*         TOTAL NUMBER OF RECORDS IN THE FILE
*         -------------------------------------------------------------

*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='UNFORMATTED',
     .               ACCESS='DIRECT', RECL=16, STATUS='OLD' )
              READ ( UNIT=LU, REC=1 ) NRECS
              DO 59 I = 1, NVALS
                  READ ( UNIT=LU, REC=I+1 ) T(I), R(I)
  59          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              TEND = T(1) + ( NRECS - NVALS + MIDDL ) * TINT
              T1 = T(1)
              LREC = 1
              MAXREC = NRECS - NVALS + 1
              IF ( TJD .LT. TBEG .OR. TJD .GT. TEND ) GO TO 10
          END IF

*         -------------------------------------------------------------

*         FILE READ SEQUENCE

  60      DIF = ( TJD - T(MIDDL) ) / TINT
*         IREC IS THE DATA RECORD NUMBER (PHYSICAL RECORD NUMBER - 1)
*         OF THE FIRST RECORD IN THE SEQUENCE OF NVALS RECORDS WITH
*         THE RELEVANT DATA TO BE RETURNED
          IREC = ( ( TJD - T1 ) / TINT ) - MIDDL + 1.9D0
          IF ( IREC .LT. 1      ) IREC = 1
          IF ( IREC .GT. MAXREC ) IREC = MAXREC

*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( IREC .GT. LREC .AND. IREC - LREC .LE. MIDDL ) THEN
*             REQUESTED JULIAN DATE JUST AHEAD OF DATA PREVIOUSLY READ
              GO TO 62
          ELSE
*             REQUESTED JULIAN DATE IN DIFFERENT PART OF FILE
              GO TO 66
          END IF

*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  62      DO 64 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  64      CONTINUE
          READ ( UNIT=LU, REC=LREC+NVALS+1 ) T(NVALS), R(NVALS)
          LREC = LREC + 1
          GO TO 60

*         GO DIRECTLY TO PROPER RECORD RANGE AND FILL UP ARRAYS
  66      DO 68 I = 1, NVALS
              READ ( UNIT=LU, REC=IREC+I ) T(I), R(I)
  68      CONTINUE
          LREC = IREC

*         -------------------------------------------------------------

      END IF

*     GOT DATA, SO FILL OUTPUT ARRAYS
  70  DO 75 I = 1, NVALS
          TLIST(I) = T(I)
          RALIST(I) = R(I)
  75  CONTINUE
      IERR = 0

  77  RETURN


      ENTRY CIOFIL ( LUNIT, FILEID, ITYPE )
*
*     THIS ENTRY ALLOWS SPECIFICATION OF THE LOGICAL UNIT NUMBER AND
*     FILE IDENTIFIER OF THE INPUT FILE CONTAINING A TIME SERIES OF CIO
*     RA VALUES.
*
*          LUNIT  = LOGIAL UNIT NUMBER TO BE USED FOR FILE (IN)
*          FILEID = FILE ID (IN)
*          ITYPE  = TYPE OF FILE (IN)
*                   SET ITYPE=1 FOR FORMATTED SEQUENTIAL FILE
*                   SET ITYPE=2 FOR UNFORMATTED BINARY FILE
*                   SET ITYPE=0 OR ANYTHING ELSE TO CLOSE THE CURRENT
*                               FILE (LUNIT AND FILEID IGNORED)
*
*     NOTE:  AFTER A CALL TO CIOFIL WITH ITYPE=0, THE ORIGINAL OR A
*     DIFFERENT FILE OF CIO RA VALUES CAN BE ACCESSED BY SUBSEQUENT
*     CALLS TO CIORD, BUT ONLY AFTER ANOTHER CALL TO CIOFIL WITH
*     ITYPE=1 OR 2.
*
*
      IF ( ITYPE .EQ. 1 .OR. ITYPE .EQ. 2 ) THEN
          LU = LUNIT
          FILNAM = FILEID
          ITYP = ITYPE
      ELSE
          CLOSE ( UNIT=LU )
          NTIMES = 0
          TBEG = 0.D0
          TEND = 1.D10
          FILEOK = .FALSE.  
      END IF

      RETURN

      END
  
*===================================================================================================


      subroutine iers_1900 (JD,UTCD,step,numday,XP,YP,DUT1,dpsi,deps,
     . DET)

      implicit none

      
      DOUBLE PRECISION  XP,YP,DUT1,X1,Y1,RJDINT,DET,B        
      DOUBLE PRECISION YEAR,X(7),Y(7),UT1(7),rjd_int,RJD(7)
      DOUBLE PRECISION JD,Xx,Yy,DUT,x_int,y_int,ut1_int
      DOUBLE PRECISION dX,dY,DDX(7),DDY(7),DEX,DEY,DELTA_T,dt
      DOUBLE PRECISION dpsi,deps,DDT(7),UTCD,UTC0,step,JDX,de,MJ0
      INTEGER :: A,I,A1,N,J,DATE,num,numday,MJD,MJ                        
      INTEGER :: ios,loop

C        INPUT :     " UTCD " UTC date (es: 1952.023,1942.11)
C                    " JD " Julian day.
C                    " step " number of day for step.
C                    " numday " number of day for year       
C        OUTPUT:     " XP " Polar motion , arcsec. 
C                    " YP " Polar motion , arcsec
C                    " DUT1 " UT1-TAI time ,sec.
C                    " dpsi " Celestial Pole Offset (longit.), as.                   
C                    " deps " Celestial Pole Offset (obliq.), as.  
C                    " DET " Value of DeltaT ,sec.
C
C        Coded by A. Nicola - October 2010 - revised October 2020

C   calculate B as compared to 2001 intervals of 18.262 MJD (UTC 1900.00-2000.00)  


!      UTC0 = 15020.0000D0             ! (MJD) UTC0 =  1900/00/00  correzione (06.04.2018)         
!      B = 1D0 + (UTCD - UTC0 ) / 18.262D0    ! 18.262 = number of JD of 2001 intervals I_N
!      DATE = int(B)                   !    to MJD 15020.0 - 51544.0

      UTC0 =1900.0d0
      MJ0 = 2400000.5d0
        B = 1+((UTCD - UTC0) * numday) / step     ! step varies from 18.25 for normal year to 18.3 for leap years 
       
      DATE = int(B)
      A1 = int(4)                 ! 4 intervals before calculating new date       
      A=int(DATE - A1)            ! as a starting point for the interpolation          
                                           
C     Open database  "Iers Delta1900B" 

      open (unit=13, file="delta1900B.txt", status="old",         
     .  action="read",form="formatted",position="rewind")

      open (unit=14, file="TDB-UT.txt", status="old",         
     .  action="read",form="formatted",position="rewind")
    
      do loop = 1,1251 ! correzione (2020/09/28   loop =2000)
      
C------ Use IERS EOP C01 data modified by A. Nicola .

C      read (unit=13,fmt="(f10.4,8x,f9.6,15x,f9.6,13x,f11.7,15x,f9.6,15x,  ! SOSTITUITO IL 2020/09/28
C     .  f9.6,19x,f6.3,8x,i4)",iostat=ios)MJD,Xx,Yy,DUT,dX,dY,dt,num         

      read (unit=13,fmt="(f10.2,2x,f10.6,14x,f10.6,14x,f11.7,13x,f10.6,  
     . 14x,f10.6,16x,i4)",iostat=ios)YEAR,Xx,Yy,DUT,dX,dY,num
          
      if (ios==0) then 
       continue
      else
       exit
      end if

      do I = 1,7              !  reads  6  values of each parameter           
       if(loop == A + I) then 
       
        RJD(I)= YEAR           ! numbers of YEAR 
        X(I)= Xx 
        Y(I)= Yy
        UT1(I)= DUT
        DDX(I)= dX
        DDY(I)= dY   
           
       end if
      end do
      end do
      rewind (unit=13)

C----- Read value dt=TDB-UT (1900/01/01 to 1962/06/01)  
      MJD = int(JD - MJ0)
      DO J = 1,23010
       read (unit=14,fmt="(f17.9,8x,f9.6)",iostat=ios)JDX,dt 
        MJ = int(JDX - MJ0) 
        IF( MJD == MJ) THEN
          DET = dt
        END IF
      END DO
      DET = DET 
      REWIND (UNIT=14)

C    " N ,NUMBER OF PAIRS FOR INTERPOLATION ( X, Y)"
      N = 7

C    "INTERPOLATION FOR UTCD (YEAR 1900.00-1962.50)
      rjd_int = UTCD
C     
C    Calling Sub INTERP LAGINT
C
      CALL INTERP2 (RJD,X,Y,UT1,DDX,DDY,N,rjd_int,x_int,y_int, 
     .  ut1_int,DEX,DEY)   !,DET) 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                                                            !        
C     Recoding values of Polar motion (XP , YP), (UT1-UTC)   !
C                 ,(dX , dY ) and  Delta T                   !                                             
C                                                            !                                                                               
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      XP = x_int
      YP = y_int
      DUT1 = ut1_int
      dpsi = DEX
      deps = DEY

      end 
C----------------------------------------------------------------------------------------------------------

      SUBROUTINE INTERP2(RJD,X,Y,UT1,DDX,DDY,N,rjd_int,  
     .  x_int,y_int,ut1_int,DEX,DEY)  !,DET)
C
C     This subroutine takes a series of x, y, and UT1-UTC values
C     and interpolates them to an epoch of choice. This routine
C     assumes that the values of x and y are in seconds of
C    arc and that UT1-UTC is in seconds of time. At least
C     one point before and one point after the epoch of the
C     interpolation point are necessary in order for the
C     interpolation scheme to work. 
C
C     parameters are :
C     RJD     - array of the epochs of data (given in mjd or step day)
C     X       - array of x polar motion (arcsec)
C     Y       - array of y polar motion (arcsec)
C     UT1     - array of UT1-UTC (sec)
C     DDX     - array of dpsi  : celestial pole offset dpsi / UAI 2000 en arcsec
C     DDY     - array of deps  : --------------------- deps ------------------ 
C     DDT     - array of Delta T (sec.) ---------- Not used
C     n       - number of points in arrays
C     rjd_int - epoch for the interpolated value
C     x_int   - interpolated value of x
C     y_int   - interpolated value of y
C     ut1_int - interpolated value of ut1-utc
C     DEX     - interpolated value of dpsi or dX
C     DEY     - interpolated value of deps or dY
C     DET     - interpolated value of Delta T  ------- Not used
C
C     CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
C                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
C                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects)
C
C      coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
C                                          Corrected : September 2007   
  
      implicit none
      integer I,N
      DOUBLE PRECISION RJD(N), X(N), Y(N), UT1(N),DDX(N),DDY(N),DDT(N),
     . rjd_int, x_int, y_int, ut1_int, DEX, DEY, DET
     
      
      
      CALL LAGINT (RJD,X,N,rjd_int,x_int)
  
      CALL LAGINT (RJD,Y,N,rjd_int,y_int)
      
      CALL LAGINT (RJD,UT1,N,rjd_int,ut1_int)
      
      CALL LAGINT (RJD,DDX,N,rjd_int,DEX)
      
      CALL LAGINT (RJD,DDY,N,rjd_int,DEY)

!      CALL LAGINT (RJD,DDT,N,rjd_int,DET)
      

      x_int = x_int 
      y_int = y_int 
      ut1_int = ut1_int 
      DEX = DEX
      DEY = DEY
!      DET = DET
     
      RETURN
      
      END
C
C----------------------------------------------------------------
C
c#################################################################
        subroutine DXDY2000_DPSIDEPS1980(dmjd,dpsi,deps,dX,dY)
c
c Transformation of the celestial pole offsets (dpsi,deps)_1980 into  
c the celestial pole offsets (dX,dY)_2000  by using SOFA 
c matrix transformation recommanded by UAI 2000 (Wallace, 2006).
c 
c 
c  input dmjd : modified julian date 
c        dpsi : celestial pole offset dpsi / UAI 1980 en mas
c        deps : --------------------- deps -----------------
c
c  output dX  : celestial pole offset dX / UAI 2000  en mas
c         dY  : --------------------- dY ------------------

c 
c  Coded by Ch. Bizouard - December 2006 
c#####################################################################


	implicit none
	double precision dmjd,dpsi,deps,dX,dY
	double precision psi,eps,epsa,RP(3,3),RN(3,3),R(3,3),
     :  X1, Y1, X2, Y2,r2mas,dj0
        double precision iau_OBL80
		
        parameter ( dj0 = 2400000.5D0,
     :              r2mas =206264.80624709635515647335733D3)
		

*  IAU 1976 precession matrix.
        CALL iau_PMAT76 ( dj0, dmjd, RP )

*  IAU 1980 nutation components.
        CALL iau_NUT80 ( dj0, dmjd , psi, eps )

*  Add nutation corrections.
        psi = psi + dpsi / r2mas
        eps = eps + deps / r2mas

*  Obliquity (IAU 1980)
        epsa = iau_OBL80 ( dj0, dmjd )

*  Nutation matrix, with respect to IAU 1976 precession.
        CALL iau_NUMAT ( epsa, psi, eps, RN )

*  Form observed precession-nutation matrix.
        CALL iau_RXR ( RN, RP, R )

*  Extract observed CIP X,Y.
        CALL iau_BPN2XY ( R, X1, Y1 )

*  -------------------------
*  GCRS X,Y of IAU 2000A CIP
*  -------------------------

*  IAU 2000A GCRS to CIRS matrix.
        CALL iau_C2I00A ( dj0, dmjd, R )

*  Extract CIP X,Y.
        CALL iau_BPN2XY ( R, X2, Y2 )

*  --------------
*  Observed dX,dY
*  --------------

*  Differences (mas).
        dX = r2mas * ( X1 - X2 )
        dY = r2mas * ( Y1 - Y2 )

        END
!-------------------------------------------------------------------------------------------


      subroutine MAGNIT (ET,NTARG,ANGLE,RVET,DIST1,PVH,EARTHVEC,OBL,
     .  HLON,HLAT,ILL,MAGN) 

!       This routine compute the magnitude of major planets      

!     ---Input        ET : TDB time. ( unit JD )
!                  NTARG : Number relative to planet.
!                  ANGLE : Phase angle. (degrees) 
!               PVH(3,2) : Earth ecliptic state vector. ( unit AU, AU/day )
!               EARTHVEC : Sun - Earth vector. ( unit AU )
!                    OBL : Earth orbit obliquity respect
!                          ecliptic plane. ( unit Rad )
!                  DIST1 : Distance body to Earth. ( unit AU )
!                   HLON : Ecliptic longitude of Saturn. ( unit Rad )
!                   HLAT : Ecliptic latitude of Saturn. ( unit Rad )
!                   RVET : Sun  vector. ( unit  AU )
!     ---Output     MAGN : Visual magnitude of planet. 
!       Coded by A. Nicola - December 2008 
      IMPLICIT NONE

! Pi
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )
      DOUBLE PRECISION DR2AS
      PARAMETER ( DR2AS = 206264.8062470963551564734D0 )
!
      INTEGER NTARG 
      DOUBLE PRECISION RVET,DIST1,MAGN,delta,EARTHVEC,deM,ANGLE
      DOUBLE PRECISION PVH(3,2),Bmagn,B0,OBL,HLON,HLAT,ET,ILL
      
      delta = DIST1

      IF (NTARG == 1) THEN
      MAGN = -0.42D0 + 5D0 * (Log10(RVET * delta)) + 0.038D0 * ANGLE - 
     . 0.000273D0 * ANGLE**2D0 + 0.000002D0 * ANGLE**3D0         ! Magnitude of Mercury to TT 0.0h 
      END IF

      IF (NTARG == 2) THEN
      MAGN = -4.40D0 + 5D0 * (Log10(RVET * delta)) + 0.0009D0 * ANGLE + 
     .     0.000239D0 * ANGLE**2D0 - 0.00000065D0 * ANGLE**3D0   ! Magnitude of Venus to TT 0.0h 
      END IF

      IF (NTARG == 4) THEN
      MAGN = -1.52D0 + 5D0 * (Log10(RVET * delta)) + 0.016D0 * ANGLE    
      END IF                                                     ! Magnitude of Mars to TT 0.0h 

      IF (NTARG == 5) THEN
      MAGN = -9.4D0 + 5D0 * (Log10(RVET * delta)) + 0.005D0 * ANGLE                                 
      END IF                                                     ! Magnitude of Jupiter to TT 0.0h 

      IF (NTARG == 6) THEN
      CALL SATURN_RING (ET,PVH,EARTHVEC,OBL,DIST1,HLON,HLAT,RVET,
     .  Bmagn,B0) 
      MAGN = -8.88D0 + 5D0 * (Log10(RVET * delta)) + 0.044D0 * ANGLE        
     .  - 2.6D0 * sin(Bmagn) + 1.25D0 * (sin(B0))**2D0      
      END IF                                                     ! Magnitude of Saturn to TT 0.0h

      IF (NTARG == 7) THEN
      MAGN = -7.19D0 + 5D0 * (Log10(RVET * delta)) + 0.0028D0 * ANGLE
      END IF                                                      ! Magnitude of Uranus to TT 0.0h

      IF (NTARG == 8) THEN
      MAGN = -6.87D0 + 5D0 * (Log10(RVET * delta))
      END IF                                                      ! Magnitude of Neptune to TT 0.0h

      IF (NTARG == 9) THEN
      MAGN = -1.01D0 + 5D0 * (Log10(RVET * delta)) + 0.041D0 * ANGLE 
      END IF                                                      ! Magnitude of Pluto to TT 0.0h
      
      IF (NTARG == 10) THEN
      
      CALL MOONLIGHT( ANGLE, deM)                                 ! deM variation in magnit. to the
      MAGN = 0.21D0 + 5D0 * (Log10(RVET * delta)) - 2.5*log10(deM)  ! law of B. Hapke.
      END IF                                                      ! Magnitude of Moon to TT 0.0h

      IF (NTARG == 11) THEN
      MAGN = 4.83D0 - 5D0 + 5D0 * LOG10(delta/DR2AS) 
      END IF                                                      ! Magnitude of Sun to TT 0.0h

      END
!====================================================================================================================================

      subroutine SATURN_RING (ET,PVH,EARTHVEC,OBL,DIST1,HLON,HLAT,
     .  RVET,Bmagn,B0)

!     This routine compute the saturnicentric latitude of the Earth
!     respect to Saturn rigs. 
!     Reference: Meeus. " Astronomical formulae for calculators"
!               " 2th edition.   pag.317"            
     
!
!     ---Input        ET : TDB time. ( unit JD )
!               PVH(3,2) : Earth ecliptic state vector. ( unit AU, AU/day )
!               EARTHVEC : Sun - Earth vector. ( unit AU )
!                    OBL : Earth orbit obliquity respect
!                          ecliptic plane. ( unit Rad )
!                  DIST1 : Distance Saturn to Earth. ( unit AU )
!                   HLON : Ecliptic longitude of Saturn. ( unit Rad )
!                   HLAT : Ecliptic latitude of Saturn. ( unit Rad )
!                   RVET : Sun - Saturn vector. ( unit  AU )
!     ---Output       B0 : Saturnicentric latitude of Earth respect to
!                          plane of rings. (unit Rad)
!                  Bmagn : Absolute value of B0 |B0|       
!        Coded by A. Nicola - December 2008 
      IMPLICIT NONE

! Pi
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )
! Degrees to radians      
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

       
      DOUBLE PRECISION RVET,DIST1,incl,delta,T,ET,R,DE,sinB1
      DOUBLE PRECISION ELON,ELAT,sinB,B0,OBL,a_axis,b_axis,N,DeU
      DOUBLE PRECISION B1,L,B,Bmagn,HLON,HLAT,U1,U2,EARTHVEC,ome
      DOUBLE PRECISION PVH(3,2)
      
      delta = DIST1                          ! distance in AU of Saturn to Earth
      T = (ET - 2451545.0d0) / 36525d0       ! T = Time in centuries from Epoch J2000

!     compute the saturnicentric latitude. 
 
      incl = (28.075216d0 - 0.012998d0 * T + 0.000004d0 * T** 2d0)*DD2R  ! inclinaz. of rings respect to
                                                                           !  ecliptic plane
      ome = (169.50847d0 + 1.394681d0 * T + 0.000412d0 * T** 2d0)*DD2R   ! longit. of ascending node respect
                                                                           !  to ecliptic plane 
 
      R=atan2(PVH(2,1),PVH(1,1))                                           
      DE=asin(PVH(3,1)/EARTHVEC)

      ELON = atan2((sin(R)*cos(OBL)+tan(DE)*sin(OBL)),(cos(R)))           ! Ecliptic Earth longit.
      ELAT = asin(sin(DE)*cos(OBL)-cos(DE)*sin(OBL)*sin(R))               ! Ecliptic Earth latit.
      
      IF ( ELON < 0D0) THEN
       ELON = ELON + 2D0 * PI
      END IF
                                                                           

      sinB = dsin(incl)*dcos(HLAT)*dsin(HLON-ome)-dcos(incl) *
     . dsin(HLAT)
      B0 = dasin (sinB)                                    ! saturnicentric latitude of Earth respect to
                                                           ! plane of rings 
 
      a_axis = 375.35d0 / delta                            ! major axis of the Saturn rings in  arcsec.
      b_axis = a_axis * sin(abs(B0))                               ! minor axis of the Saturn rings in  arcsec.
 
      N = 113.66552d0 * DD2R + 0.8771d0 * DD2R * T         ! compute the long. ascending node of Saturn orbit.
                                                           !B0 = datan(sinB / dsqrt(-sinB * sinB + 1D0))  
 
!------ Correction to long. e Latit. for Saturn diurnal aberraz.  
      L = HLON - (0.01759d0 * DD2R / RVET)
      B = HLAT - (0.000764d0 * DD2R) * (Cos(HLON - N ) / RVET)

      sinB1 = dsin(incl)*dcos(B) * dsin(L-ome) - dcos(incl) * dsin(B)
      B1 = dasin(sinB1)                                    ! saturnicentric latitudine del Sole  riferita
                                                           ! al piano degli anelli( positiva verso Nord
 

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!  Now we calculate " deU " , it's the difference between the Saturnicentric longitudes of the Sun
!  and the Earth , measured in the plane of the rings.  
!   Saturnicentric long. Earth = U1+180°,  Saturnicentric long. Sun = U2+180°

 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 
      U1 = datan((dsin(incl)*dsin(B) + dcos(incl)*dcos(B)*dsin(L-ome))
     .  / (dcos(B)  * dcos(L - ome)))                            ! Saturn geocentric longit measured in to ring plane
  
      U2 = datan((dsin(incl) * dsin(HLAT) + dcos(incl) * dcos(HLAT) *
     . dsin(HLON - ome)) / (dcos(HLAT) * dcos(HLON - ome)))      ! Saturn Eliocentric longit. measured in to ring plane
  

      U1 = U1 * DR2D                       ! transorm U1 and U2 to Degrees
      U2 = U2 * DR2D      
  
      deU = Abs(U1 - U2)                   ! (deU) never more than 7 ° angle
 
      Bmagn = B0 * DR2D
      Bmagn = Abs(Bmagn) * DD2R


      END
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      SUBROUTINE R_vect(ET,Dtime,NTARG,HMOD,BODY_G,BODY_H,HD,THD) 
    
!            
!     This subroutine compute the geocentric and heliocentric true vector  
!          position , distance and radius vector of the body.       
      
!       Input:  ET     = TDB time
!               NTARG  = The numbering convention for the body 
!               PVE    = Earth Barycentric state  vector(ET time)
!               HMOD   = Earth heliocentric position modulus (ET time)  
!       Output  BODY_G = Geocentric position true vector 
!               BODY_H = Heliocentric position true vector
!               GD     = true distance (AU) 
!               HD     = Radius vector (AU)
!               XET    = TDB time corrected to light-time                            
!               PV     = Body's barycentric state vector (XET time)

      IMPLICIT NONE
      DOUBLE PRECISION C,AU,DS,AUDAY,MU,MUC,Dtime,DT,XET,HD,THD
      DOUBLE PRECISION ET,NTARG,HMOD,R0,R1,RTERM,TAU,RV,DST,DOT,HV(3)
      DOUBLE PRECISION EBT(3,2),SB(3,2),PV(3,2),BODY_H(3)
      DOUBLE PRECISION PV1(3),EBP(3),HE(3),GE(3),RQ(3),RS(3),BODY_G(3)
      DOUBLE PRECISION GV(3),BV1(3),EBV(3),UB(3,2),UV(3,2)
      
      C  = 299792458D0                     ! speed of light  m/s 
      AU = 149597870691D0                  ! Astronomical Unit (meter)(TDB)
      DS = 86400D0                         ! n° seconds  Day
      AUDAY = (C * DS / AU)                ! Speed of light (AU per day  "TDB")
      MU = 1.32712440041D20                ! heliocentric gravitational costant (2010 system)   
      MUC = 2D0 * (MU / (C * C))/AU

!---------Body's astrometric position (recursive determination)

      XET = ET - Dtime                          ! XET = new reference time

      CALL PLEPH ( XET, NTARG, 12, PV)         ! PV = body's barycentric state vector 

      CALL PLEPH ( XET, 11, 12, SB)            ! SB = sun's barycentric state vector  

      UB = PV - SB                             ! UB = Body's eliocentric state vector 

      CALL iau_PVUP (UB,DT,HV)                 ! HV Body's Heliocentric velocity 

!------Approximation to light time 

      CALL iau_PVUP (DT,UB,HE)                 ! HE Body's Heliocentric position     
        

!-----Form modulus of HE 
      
      CALL iau_PM(HE,HD)         ! HD modulus Heliocentric position body


      
      IF (Dtime == 0D0) THEN
        THD = HD                 ! THD True Heliocentric distance

      
      END IF


      
!----------Form unit vector of " RQ , RS "     

      CALL iau_PN(HE,R0,RQ)            ! RQ Unic vector body_Helioc
        
!      CALL iau_PN(GE,R1,RS)            ! RS Unic vector body_Geocentric

        
!----------------------------------------------------------------------------
       
!      BODY_G = RS * TGD                ! True Vector Geocentric  Position
       
      BODY_H = RQ * THD                ! True Vector Heliocentric  Position
     
      END
!==============================================================================================================


      subroutine MOONLIGHT (ANGLE,deM)

!      This routine use the law of B. Hapke for compute deM (function to phase)
!      for estimate the apparent Lunar magnitude. 
!      Modified algorithm by Giuseppe Matarazzo (2015/01/26)
 
        
      IMPLICIT NONE

! Pi
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )
! Degrees to radians      
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

      DOUBLE PRECISION DeM, theta0, ANGLE, A, y0
      DOUBLE PRECISION  w, h, S0, b, c, TH, B0, A1,A2,K, coef
      DOUBLE PRECISION B_ih, P_i, P0, r0, p, P_HI, PHI,C0,x
      INTEGER n,m
      
      A = ANGLE * DD2R              !  the phase angle (radians)
      theta0 = 20.0D0               !  parameter of roughness (angle of reflectance)

      TH = theta0 * DD2R

      w = 0.21D0                    !  albedo at ph = 0°

      h = 0.07D0                    ! half angular width of the opposition effect

      S0 = 0.71D0                   ! terms of amplitude effect of the opposition

      b = 0.29D0                    ! cofficients of the Legendre polynomial

      C0 = 0.39D0                   ! phase function of single scattering particle P(i)
      
      B0 = 2.01D0                   ! total amplitude the opposition effect

      
      B_ih = B0 / (1D0 + 1D0/h * tan(A/2D0)) 

      P_i = 1D0 + b * cos(A) + C0 * (1.5D0 * (cos(A))**2D0 - 0.5D0)

      P0 = S0 / ( w * B0 )

      y0 = SQRT( 1D0 - w )

      r0 = (1D0 - y0) / (1D0 + y0)

      p = (w/8D0) * ((1D0+B0)*P0-1D0) + ( r0/2D0 + (r0**2D0)/6)

!--------- The term PHI is the phase function for smooth surfaces 

      
      PHI = (1D0/p)* (w/8D0 * ((1D0 + B_ih) * P_i - 1D0) 
     .  + r0/2D0 * (1D0 - r0) * (1D0 - sin(A/2D0) * tan(A/2D0) 
     .  * LOG(1/tan(A/4)))+(2D0*r0**2D0)/(3D0*PI)*(sin(A)
     .  + (PI-A) * cos(A))) 

!      write (*,*) ' Moon phase angle= ', i
! ---------- finalK: Jan.26,2015 ---------------------
      IF (ANGLE.GE.0d0 .AND. ANGLE.LE.58.520d0) THEN
!     as per formula (15)
!     theta = TH is in radians; A= ANGLE in radians, too         
      coef= -TH*(0.32d0*sqrt(tan(TH)*tan(A/2.0))
     .      +0.52d0*tan(TH)*tan(A/2.0)) 
      K= exp(coef) 
!     curve #1 (cubic) 
      ELSE IF (ANGLE.GT.58.520d0 .AND. ANGLE.LE.106.644d0) THEN  
      K= 4.227D-6*ANGLE**3 - 7.048513D-4*ANGLE**2
     .   +3.87961041D-2*ANGLE + 9.30748447D-2  
!     curve #2 (polynomial 5th degree)            
      ELSE IF (ANGLE.GT.106.644d0 .AND. ANGLE.LE.120.589d0) THEN   
      K=  4.5974D-6*ANGLE**5 -2.6170774D-3*ANGLE**4
     .   +0.5955197293D0*ANGLE**3 -6.77127016892D+1*ANGLE**2
     .   +3.8472367031696D+3*ANGLE -8.73831821270379D+4        
!     curve #3 (polynomial 5th degree)             
      ELSE IF (ANGLE.GT.120.589d0 .AND. ANGLE.LE.128.887d0) THEN   
      K= 1.98652D-5*ANGLE**5 - 1.24476522D-2*ANGLE**4
     .  +3.119376959D0*ANGLE**3 -3.907937048388D+2*ANGLE**2 
     .  +2.44753934296309D+4*ANGLE -6.130635997016D+5   
!     curve #4 (straight line)             
      ELSE IF (ANGLE.GT.128.887d0 .AND. ANGLE.LE.139.800d0) THEN  
      K= -6.9710897D-2*ANGLE + 10.2783976288D0      
!     curve #5 (parabola)             
      ELSE IF (ANGLE.GT.139.800d0 .AND. ANGLE.LE.161.48d0) THEN  
      K= 1.1439325D-3*ANGLE**2 - 0.3679413977D0*ANGLE
     .   + 2.96122777626D+1 
!     curve #6 (exponential "mathmagic")    
      ELSE IF (ANGLE.GT.161.48d0 .AND. ANGLE.LE.172.332d0) THEN
      K= 2.044900838707796D+15*exp(-0.2418415532*ANGLE)
     .  -1.15707D-5*ANGLE**2 + 3.7865713D-3*ANGLE - 0.3097395215D0
!     curve #7 (my parabola, 172.33 to 180)
      ELSE
      K= 2.01766D-5*ANGLE**2 - 7.2132306D-3*ANGLE + 0.6446609563D0      
      END IF   
! ---------- end finalK: Jan.26,2015 ---------------------
      
! ----skipped ----         
! ----------------- NEW interpolation: cubic spline --------       
!  call spline to calculate spline coefficients
!      call spline (xi, yi, bs, cs, ds,ns) 
!      K=ispline(ANGLE, xi, yi, bs, cs, ds, ns)    
! ----------------------------------------------------------
!      call bilinear (ANGLE, theta0, K)
! --- end skip -------

      deM =  PHI * K
!      write (*,*) ' Coefficient K= ', K      
      END

!==========================================================================================================================


      subroutine COORD (NTARG,TDB0,AUDAY,DX00,DY00,RA3,R2,R3,R6,R8,R9,       
     .  R10,DE3,D2,D3,D6,D8,D9,D10)
!===========================================================================
!
!     This routine compute the apparent geocentric coordinate
!     for the routine "Transit
!=========================================================================== 
!      Parameter
!      INPUT:
!           NTARG:     is Solar System Body number     
!           TDB0 :     Barycentric Dinamical Time to 0.0h UT1
!           GST0 :     Greenwich Sideral Time 
!           RNPB0:      Precession-Nutation matrix to 0.0h UT1. 
!               
!      OUTPUT
!           RA2/6/8/9: RA value 
!           DE2/4/5/6: Decl. value
!=========================================================================== 
 

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

      INTEGER NTARG,I ,J,II,C
      DOUBLE PRECISION TDB0,T0,DT,HMOD,GD,HD,DJMJD0,TT,DX00,DY00
      DOUBLE PRECISION PVB(3,2),BEV(3,2),PVH(3,2),PVE(3,2),HPE(3)
      DOUBLE PRECISION BEP(3),BVE(3),EB(3),BODY_G(3),BODY_H(3),PVV(2)
      DOUBLE PRECISION XET,PV(3,2),HE0,THD,EPV(3,2),u1(3),u0(3),u3(3)
      DOUBLE PRECISION MO1,RNPB0(3,3),u4(3),RA,DE,RA1,RA2,RA3,DEPS
      DOUBLE PRECISION RA4,RA5,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,DE5
      DOUBLE PRECISION D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,DE1,DE2,DE3,DE4
      DOUBLE PRECISION DE6,T,iau_ANP,AUDAY,UB(3,2),T1,T2,T3,RR,DD,DPSI
      DOUBLE PRECISION DP00,DE00,EPSA,RB(3,3),DDP00,DDE00,V1(3),V2(3)
      DOUBLE PRECISION RPB(3,3),RN(3,3),RNPB(3,3), RP(3,3)

      DJMJD0 = 2400000.5D0      
      T = TDB0 - 2D0
!      print*, "T----",T
20    T = T     
      IF ( T == TDB0 + 4) GOTO 30 

      CALL PLEPH ( T, 17, 0, PVV)              ! PVV(1) = TT-TDB in sec./ PVV(2)=rate of change of TT-TDB in sec/day

      CALL PLEPH (T,NTARG,12,PVB)

      CALL PLEPH ( T, NTARG, 11, BEV)          ! BEV = Body's heliocentric state vector

      CALL PLEPH ( T, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( T, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
!------Update a PVH-vector, discarding the velocity component.

      DT = 0.0D0     
      CALL iau_PVUP (DT,PVH,HPE)               ! HPE Heliocentric position Earth vector 

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP Body's heliocentric position

      CALL iau_PVUP (DT,PVE,EB)                  ! EB Barycentric  position Earth vector
 
      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus
      
      CALL iau_PVUV (DT,PVE,BVE)                 ! BVE Barycentric velocity Earth vector      
      
      CALL L_TIME(T,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,
     . UB,THD,EPV)

!-------- Nutation, IAU 2000A.------------------------------------------------------------

      TT = (T - (PVV(1)/86400D0))- DJMJD0

      CALL iau_NUT00A ( DJMJD0, TT, DP00, DE00 )
   
!-------- Precession-nutation quantities, IAU 2000.
      CALL iau_PN00 ( DJMJD0, TT, DP00, DE00, 
     .      EPSA, RB, RP, RPB, RN, RNPB )
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
      CALL iau_RXR ( RN, RPB, RNPB0 )
           
!----------------------------------------------------------------------------
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun ( for apparent geocentric coordinate)

      CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
 
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent geocentric coordinate)

      CALL ABERR (AUDAY,u1,BVE,u0)

      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 
      
!---------- Form apparent geocentric coordinate (RA,DECL)       

! Option from the method of calculating "IUA 2000 Equinox based" or "IAU 2006
!                   CIO based using X,Y series" 
      
      CALL   iau_RXP ( RNPB0, u3, u4 )
      
! ......... Call for direction cosines to spherical coordinates

      CALL   iau_C2S   ( u4, RR, DD )
      
      IF (T == TDB0-2D0)THEN  
        RA1 = RR 
        DE1 = DD
       ELSE IF (T == TDB0-1D0)THEN
        RA2 = RR
        DE2 = DD
       ELSE IF (T == TDB0)THEN
        RA3 = RR
        DE3 = DD
       ELSE IF (T == TDB0+1D0)THEN
        RA4 = RR
        DE4 = DD
       ELSE IF (T == TDB0+2D0)THEN
        RA5 = RR
        DE5 = DD
      END IF
        T = T + 1D0
      GOTO 20
      
       
30    If ((RA2 - RA1) < -1.7) RA2 = RA2 + D2PI   ! routine that eliminates the step of ambigità equinox in RA
      If ((RA3 - RA2) < -1.7) RA3 = RA3 + D2PI
      If ((RA4 - RA3) < -1.7) RA4 = RA4 + D2PI
      If ((RA5 - RA4) < -1.7) RA5 = RA5 + D2PI
            
      R1 = RA2 - RA1                                 ! differences Bessel interpolation "RA" (1 st order)
      R2 = RA3 - RA2
      R3 = RA4 - RA3
      R4 = RA5 - RA4                                 ! differences Bessel interpolation "RA" (2 st order)
      R5 = R2 - R1
      R6 = R3 - R2
      R7 = R4 - R3
      R8 = R6 - R5
      R9 = R7 - R6
      R10 = R9 - R8

      D1 = DE2 - DE1                                 ! differences Bessel interpolation "Decl"(1 st order)
      D2 = DE3 - DE2
      D3 = DE4 - DE3
      D4 = DE5 - DE4                                 ! differences Bessel interpolation "RA" (2 st order)
      D5 = D2 - D1
      D6 = D3 - D2
      D7 = D4 - D3
      D8 = D6 - D5
      D9 = D7 - D6
      D10 = D9 - D8
      
      END
!-------------------------------------------------------------------------------------------------------------------------------


      SUBROUTINE iau_PVUV ( DT, PV, V )
*+
*  - - - - - - - - -
*   i a u _ P V U P
*  - - - - - - - - -
*
*  Update a pv-vector, discarding the position component.
*
*  This routine is part of the International Astronomical Union's
*  SOFA (Standards of Fundamental Astronomy) software collection.
*
*  Status:  vector/matrix support routine.
*
*  Given:
*     DT       d           time interval
*     PV       d(3,2)      pv-vector
*
*  Returned:
*     V        d(3)        v-vector
*
*  Notes:
*
*  1) "Update" means "refer the velocity component of the vector to a
*     new date DT time units from the existing date".
*
*  2) The time units of DT must match those of the velocity.
*
*  This revision:  2012 June 14
*
*  NICOLA ALDO release 2012-06-14
*
*  Copyright (C) 2012.
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DT, PV(3,2), V(3)

      INTEGER I

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO 1 I=1,3
         V(I) = DT * PV(I,1) + PV(I,2)
 1    CONTINUE

*  Finished.
      END
*+----------------------------------------------------------------------
*
*  Copyright (C) 2012
*  Standards Of Fundamental Astronomy Board
*  of the International Astronomical Union.
*  Software modified
*  =====================
*  SOFA Software License
*  =====================
*
*--------------------------------------------------------------------------------------------------------------------------

       subroutine Interpol(GST0,oldm,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,H0,newm,Ha,delta,azm)

     
      IMPLICIT NONE
      
* 2PI
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI = 6.283185307179586476925287D0 )
* PI
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )
      
      DOUBLE PRECISION gast, N, alpha, H, dt, denom, B2, HH, GST0 
      DOUBLE PRECISION oldm,DELTAT,RA3,R2,R3,R6,R8,R9,R10,LONG, H0
      DOUBLE PRECISION newm, Ha, delta, azm, LAT, DE3,D2,D3,D6,D8
      DOUBLE PRECISION D9,D10
      INTEGER C
      
      newm = oldm
      
      DO
      gast = GST0 + 6.300388093D0 * newm

      gast = MODULO ( gast,D2PI) 

      N = newm + DELTAT 

!      B2 = N * (N - 1D0) / 4D0                                            ! B2 is parameter relate a day fraction 

!      alpha = RA2 + N * RA6 + B2 * (RA8 + RA9)

      alpha = RA3 + N *((R2+R3) / 2D0 - (R8+R9) / 12D0) + N * N *
     . (R6 / 2D0 - R10 / 24D0) + N**3D0 * ((R8+R9) / 12D0) + N**4D0 *
     . (R10 / 24D0) 

      alpha = MODULO (alpha, D2PI)                                        ! Mod (alpha,D2PI)

      delta = DE3 + N *((D2+D3) / 2D0 - (D8+D9) / 12D0) + N * N *
     . (D6 / 2D0 - D10 / 24D0) + N**3D0 * ((D8+D9) / 12D0) + N**4D0 *
     . (D10 / 24D0) 

!      delta = DE2 + 0.5D0 * N * (DE4 + DE5 + N * DE6)                     ! other Bessel 2° order formula (Decl)

      H = gast + LONG - alpha                                             ! local hour angle of body

      Ha = ASIN (SIN(LAT) * SIN(delta) + COS(LAT) * COS(delta) * COS(H))  ! Ha = altitude of body
      
      azm = ATAN (SIN(H) / (COS(H) * SIN(LAT) - TAN(delta) * COS(LAT)))   ! azimuth of body
      
      IF (azm < 0D0) azm = azm + D2PI
      
      IF ( C == 0 ) THEN
         H = MODULO ( H, D2PI )                                           ! = Mod ( H ,2d0 * PI) 
       IF (H > PI ) H = H - D2PI                                          
         dt = -H / D2PI  
       ELSE 
         dt = (Ha - PI/2D0 + H0)/(D2PI * COS(delta) * COS(LAT) * SIN(H))
      END IF
      
      newm = newm + dt
      Ha = Ha
      delta = delta
      azm = azm
      
      IF ( ABS(dt) .LT. 1D-15) EXIT
      
      END DO
      
      END                         
*-------------------------------------------------------------------------------------------------------------------------------

* --------------- Interpolazione BILINEAR ------------------
*
* Lettura dei parametri tabulati, scelta della griglia di
* contorno di un punto dato e calcolo dell'interpolazione   
* -----------------------------------------------------------
* Risultati dell'esempio di calcolo
*   Parametri della griglia
*       0.984       0.974
*       0.965       0.944
*        d1          d2
*      10.000       3.000
*
*    Valore di interpol. BILINEARE:  0.964304
* ---------------------------------------------------------------
      subroutine bilinear (i0, theta0, Y)

      implicit none
      INTEGER N,M
      DOUBLE PRECISION  x1a,x1b,x2a,x2b, x1,x2, d1,d2, y1,y2,y3,y4
      DOUBLE PRECISION  t,u, Y
      DOUBLE PRECISION  K(21,7), i(21), theta(7), i0, theta0
      INTEGER  nup,ndown, mup,mdown
      DATA i /0,2,5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,
     :        150,160,170,180/
      DATA theta /0,10,20,30,40,50,60/

*  Parametri K(i,theta)
      DATA (K(N,1),K(N,2),K(N,3),K(N,4),K(N,5),K(N,6),K(N,7), N=1,21)
     :  / 1D0,     1D0,     1D0,     1D0,     1D0,     1D0,     1D0,  
     :    1D0,  997D-3,  991D-3,  984D-3,  974D-3,  961D-3,  943D-3,
     :    1D0,  994D-3,  981D-3,  965D-3,  944D-3,  918D-3,  881D-3,
     :    1D0,  991D-3,  970D-3,  943D-3,  909D-3,  866D-3,  809D-3,
     :    1D0,  988D-3,  957D-3,  914D-3,  861D-3,  797D-3,  715D-3,    
     :    1D0,  986D-3,  947D-3,  892D-3,  825D-3,  744D-3,  644D-3,      
     :    1D0,  984D-3,  938D-3,  871D-3,  789D-3,  692D-3,  577D-3,  
     :    1D0,  982D-3,  926D-3,  846D-3,  748D-3,  635D-3,  509D-3,
     :    1D0,  979D-3,  911D-3,  814D-3,  698D-3,  570D-3,  438D-3,        
     :    1D0,  974D-3,  891D-3,  772D-3,  637D-3,  499D-3,  366D-3,    
     :    1D0,  968D-3,  864D-3,  719D-3,  566D-3,  423D-3,  296D-3,
     :    1D0,  959D-3,  827D-3,  654D-3,  487D-3,  346D-3,  231D-3,
     :    1D0,  946D-3,  777D-3,  575D-3,  403D-3,  273D-3,  175D-3,
     :    1D0,  926D-3,  708D-3,  484D-3,  320D-3,  208D-3,  130D-3,
     :    1D0,  894D-3,  617D-3,  386D-3,  243D-3,  153D-3,   94D-3,
     :    1D0,  840D-3,  503D-3,  290D-3,  175D-3,  107D-3,   64D-3,
     :    1D0,  747D-3,  374D-3,  201D-3,  117D-3,   70D-3,   41D-3,    
     :    1D0,  590D-3,  244D-3,  123D-3,   69D-3,   40D-3,   23D-3,      
     :    1D0,  366D-3,  127D-3,   60D-3,   32D-3,   18D-3,   10D-3,
     :    1D0,  128D-3,   37D-3,   10D-3,   85D-4,   47D-4,   26D-4,     
     :    1D0,     0D0,     0D0,     0D0,     0D0,     0D0,     0D0 /     
         

         
* Coordinate Punto dato (i0,theta0)   
!      i0=176.64D0    !41.100991317d0    !3.956d0
!      theta0=20d0   !34.256d0

* ---- Inizio scelta della griglia ---- 
      N=1
      DO WHILE (i(N)<=i0)
       N=N+1 
	  END DO 
	  nup=N
	  ndown=N-1
      
      M=1
      DO WHILE (theta(M)<=theta0)
       M=M+1          
      END DO
	  mup=M
	  mdown=M-1
* ---- Fine scelta della griglia ---- 
   
* Valori della matrice agli angoli della griglia selezionata   
      y1= K(ndown,mdown)
      y2= K(ndown,mup)
      y3= K(nup,mup)
      y4= K(nup,mdown)
* ---------------------------
     
* Estremità griglia         
      x1a=theta(mdown)
      x1b=theta(mup)
      x2a=i(ndown)
      x2b=i(nup)

* -------------- Inizio Calcoli di interpolazione  ------------
*
* base griglia
      d1= theta(mup)-theta(mdown) 
* alt. griglia      
      d2= i(nup)-i(ndown)         
      
*
* fraz. asc. griglia [0,1]
      t= (theta0-x1a)/d1  
* fraz. ord. griglia [0,1]      
      u= (i0-x2a)/d2      
*      
*     Interpolazione bilineare                                        
*
      Y= (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4
*
!      write(*,*)
!      write (*,48),' Valore di interpol. BILINEARE: ', Y
!48    format(a35,f9.6," ")

 
      end 
!--------------------------------------------------------------------------------------


      SUBROUTINE moon(TDB,YA,MES,ID,DPSI,RA6,DE6,DISTKm,R_sun,D_sun,
     .        RA7,DE7,OBL,GD1,PL,HLON,HLAT,OMEGA,
     .        OME0,AGE,JH0,CHI,JAP,APE,DLR,TLong,TLat,LTL,LTB,PA, 
     .        THL,THB,JH1,JHP,JHU)

      implicit none

! Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )! Pi
      DOUBLE PRECISION PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )

      double precision TDB,EL,ELP,F,D,OMEGA,L1,L1R,PRG,A01,A02,A03
      double precision PRGR,E,T,T0,LoLu,LoLu_ap,DPSI,RA6,DE6,DISTKm
      double precision R_sun,D_sun,PL,GD1,Day,HI,HL,CHI,KA,KO,JDE0,E0
      double precision Tn,MMR,MLR,FFR,OMR,LatLu,BL,K0,SL,OME0,A1,A2
      double precision A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,AP,PE
      double precision JH,JH0,JE,AGE,KB,FKB,KC,TAP,DAP,MPA,FAP,JAP
      double precision APS,SR,TAPOG,DLR,DR,FRX,OM,M1R,MR,JDA,OME1,K10
      double precision INC,WL,XA,XAA,ATXA,OMEGAR,E1,RO,SIGMA,TAU,HB
      double precision LOL,LOB,SIN_B,LFL,LFB,LTL,LTB,Beta,K1,K2,PA,K20
      double precision RA7,DE7,OBL,TLong,TLat,HLON,HLAT,ome,VV,X,Y,NV
      double precision WS,XS,XSS,ATXS,LOH,SIN_H,LBH,LFH,HFB,THL,THB
      double precision JE1,JH1,JHP,W,JHU,K3  !JEP,JEU,JE0,K3
      integer MES,YA,ID
      character*8 APE
      
!      T0 = (TDB - 2451545d0) / 365250d0    !Time in fract/millennium a by 12h  TDB of 2000 Gen 01
      T = (TDB - 2451545d0) / 36525d0      !Time in fract/century a by 12h  TDB of 2000 Gen 01 

      CALL FUNARG (T,EL,ELP,F,D,OMEGAR)   ! This routine computes fundamental argument
                                           !(MEAN ELEMENTS)of the SUN and MOON.
      
*                                            T = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*                                            EL  = MEAN ANOMALY OF THE MOON IN RADIANS
*                                                  AT DATE TDB (OUT)
*                                            ELP = MEAN ANOMALY OF THE SUN IN RADIANS
*                                                  AT DATE TDB (OUT)
*                                            F   = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
*                                                  OF THE MOON'S ASCENDING NODE IN RADIANS
*                                                  AT DATE TDB (OUT)
*                                            D   = MEAN ELONGATION OF THE MOON FROM THE SUN IN
*                                                  RADIANS AT DATE TDB (OUT)
*                                          OMEGAR = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
*                                                  IN RADIANS AT DATE TDB (OUT)
*
!  Compute the mean lunar long. L1 in degrees taking into account the effect of the time light 
      L1 = 218.3164477d0 + 481267.88123421d0 * T - 0.0015786d0 * 
     . T**2d0 + T**3d0 / 538841d0 - T**4d0 / 65194000d0
      L1R = L1 * DD2R           ! in Radians
!  Compute the mean lunar longit. of Perigee (PRG)
      PRG = 83.3532465d0 + 4069.0137287d0 * T - 0.01032d0 * T**2d0 
     .      - T**3d0 / 80053d0 + T**4d0 / 18999000d0
      PRGR = PRG * DD2R         ! in Radians
 
      A01 = (119.75d0 + 131.849d0 * T) * DD2R         ! Argom. due to the action of Venus. 
      A02 = (53.09d0 + 479264.29d0 * T) * DD2R        ! Argom. due to the action of Jupiter
      A03 = (313.45d0 + 481266.484d0 * T) * DD2R      ! Additional Argom. 
   
!  Compute the mean lunar longit.(Coeff. in 0.000001 degrees)
   
      E = 1d0 - 0.002516d0 * T - 0.0000074d0 * T**2d0  ! Additional Coeff. due at eccentric 
                                                           ! terrestrial orbit ( MR ).


!------Compute the additional term for lunar longit.- latit.- long. ascending node.----------

      CALL MEANLON(E,L1R,EL,ELP,F,D,OMEGAR,A01,A02,A03,OME0,SL,BL,SR )

      LoLu = L1 + SL / 1000000d0

      LoLu = MOD(LoLu,360d0)                    ! Mean lunar longitude (degrees)

      IF(LoLu < 0d0) LoLu = LoLu + 360d0

      LoLu_ap = LoLu + DPSI * DR2D 

      LoLu_ap = MOD(LoLu_ap,360d0)              ! Apparent lunar longitude (degrees) 
  
      IF(LoLu_ap < 0d0) LoLu_ap = LoLu_ap + 360d0

      LatLu = BL / 1000000d0                    ! Moon's latitude (degrees)
!----------------------------------------------------------------------------------------------
   
!   ' Mean Longit. of ascending node between 0.0° and 360°

      OMEGA = OMEGAR * DR2D
      OMEGA = MOD(OMEGA,360D0) 
      IF(OMEGA < 0d0) OMEGA = OMEGA + 360d0
     
!   ' True Longit. of ascending  node between 0.0° e 360°

      OME0 = MOD(OME0,360D0) 
      IF(OME0 < 0d0) OME0 = OME0 + 360d0
     
!----------------------------------------------------------------------------------------------
  
      R_sun = R_sun * DD2R
      D_sun = D_sun * DD2R
      HI = Cos(D_sun) * Sin(R_sun - RA6)                                     ! numerator  
      HL = Sin(D_sun) * Cos(DE6) - Cos(D_sun) * Sin(DE6) *                   !denominator  
     .     Cos(R_sun - RA6)       
      
      CHI = ATAN2(HI,HL)
      CHI = CHI * DR2D                    ! Position angle for the Moon's bright limb (degrees)
      if (CHI < 0d0) CHI = CHI+360d0
!----------------------------------------------------------------------------------------------


!------ calculating of the lunar phase of New Moon
  
      If (MES ==  1) Day = 15.5
      If (MES ==  3) Day = 15.5
      If (MES ==  4) Day = 15
      If (MES ==  5) Day = 15.5
      If (MES ==  6) Day = 15
      If (MES ==  7) Day = 15.5
      If (MES ==  8) Day = 15.5
      If (MES ==  9) Day = 15
      If (MES ==  10) Day = 15.5
      If (MES ==  11) Day = 15
      If (MES ==  12) Day = 15.5
      If (MES ==  2 .And. YA == 1700 .Or. YA == 1800 .Or. YA == 1900
     .   .Or. YA == 2100) Day = 14
      If (MES ==  2 .And. Mod (YA, 4 ) .NE. 0 )Day = 14
      If (MES ==  2 .And. Mod (YA, 4 ) == 0 ) Day = 14.5

      KA = ((YA + MES / 12D0 - Day / 365.25d0) - 2000d0) * 12.3685d0    ! K0 = aproximate parameter 
      K0 = FLOOR(KA)
      K1 = K0 + 0.25d0 
      K2 = K0 + 0.5d0
      K3 = K0 + 0.75d0
      
      CALL TIMEPHASE (K0,0,JH0)      ! JE0  time for New Moon      
                       
      CALL TIMEPHASE (K1,1,JH1)      ! JE1  time for First Quarter
                   
      CALL TIMEPHASE (K2,2,JHP)      ! JEP  time for Full Moon
                       
      CALL TIMEPHASE (K3,3,JHU)      ! JEU  time for Last Quarter
                      
       

      AGE = TDB - JH0               !AGE ,age of last New Moon (day)
!      print*,"TDB,JH0",TDB,JH0    
      If (AGE < 0) Then
       AGE = 29.4 + AGE
      ElseIf (AGE > 29.530589) Then
       AGE = AGE - 29.530589
      End if
  
!--------Compute Moon's Apogee and Perigee.
  
      KB = (YA + (MES - 1d0) / 12d0 + ID / 365.25d0 - 1999.97d0)     ! KB = aproximate parameter
     .     * 13.2555d0 
      FKB = KB - Int(KB)                                             ! FKB = fract. part of KB
   
      If (KB > 0) Then
        If (Abs(FKB) < 0.75d0 .And. Abs(FKB) > 0.25d0) FKB = 0.5d0
        If (FKB < 0.25d0) FKB = 0d0
        If (FKB > 0.75d0) FKB = 1d0
        KC = Int(KB) + FKB
      ElseIf (KB < 0) Then
        If (Abs(FKB) < 0.75d0 .And. Abs(FKB) > 0.25d0) FKB = -0.5d0
        If (Abs(FKB) < 0.25d0) FKB = 0d0
        If (Abs(FKB) > 0.75d0) FKB = -1d0
        KC = Int(KB) + FKB                 ! FIX = INT
      End If
  
      TAP = KC / 1325.55d0                                        ! TAP = time in julian centuries
  
!-------- JDA  Julian date 
      JDA = 2451534.6698d0+27.55454989d0*KC-0.0006691d0*TAP**2d0 - 
     .      0.000001098d0*TAP**3d0+0.0000000052d0*TAP**4d0
    
!-------- Moon's mean elongation at time TAP in Rad.
      DAP = (171.9179d0+335.9106046d0*KC-0.0100383d0*TAP**2d0 -
     .      0.00001156d0*TAP**3d0+0.000000055d0*TAP**4d0) * DD2R
  
!-------- Moon's mean anomaly at tyme TAP in Rad.
      MPA = (347.3477d0+27.1577721d0*KC-0.000813d0*TAP**2d0 - 
     .      0.000001d0*TAP**3d0) * DD2R
  
!-------- Moon's argument Latitude at time TAP in Rad.
      FAP = (316.6109d0+364.5287911d0*KC-0.0125053d0*TAP**2d0 - 
     .      0.0000148d0*TAP**3d0) * DD2R

      CALL PER_APOG(TAP,DAP,MPA,FAP,PE,AP,APS )
  
      If (Abs(FKB)== 0.5)Then
       JAP = JDA + AP +(APS/(86400*15))
       APE = "APOGEE"
      Else If (Abs(FKB) == 1 .Or. FKB == 0) Then
       JAP = JDA + PE 
       APE = "PERIGEE"                    
      End If
    
      TAPOG = (JAP - 2451545d0) / 36525d0                ! Compute the fraction of century at the date of APOGEE

      CALL FUNARG (TAPOG,M1R,MR,FRX,DR,OM)   ! This routine computes fundamental argument
                                                !(MEAN ELEMENTS)of the SUN and MOON.
      
!                                            TAPOG = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
!                                            M1R  = MEAN ANOMALY OF THE MOON IN RADIANS
!                                                  AT DATE TDB (OUT)
!                                            MR   = MEAN ANOMALY OF THE SUN IN RADIANS
!                                                  AT DATE TDB (OUT)
!                                            FRX  = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
!                                                  OF THE MOON'S ASCENDING NODE IN RADIANS
!                                                  AT DATE TDB (OUT)
!                                            DR   = MEAN ELONGATION OF THE MOON FROM THE SUN IN
!                                                  RADIANS AT DATE TDB (OUT)
!                                            OM   = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
!                                                  IN RADIANS AT DATE TDB (OUT) !!! NOT USED !!!
!

  
!------- Compute the lunar mean long. L1 in degrees taking into account the effect of the time light

      L1=218.3164477d0+481267.88123421d0*TAPOG-0.0015786d0*TAPOG**2d0
     .   + TAPOG**3d0 / 538841d0-TAPOG**4d0 / 65194000d0
      L1R = L1 * DD2R           !'  Radians
  
   
      E1 = 1d0 - 0.002516d0 * TAPOG - 0.0000074d0 * TAPOG**2d0  ! Additional Coeff. due at eccentric 
                                                               ! terrestrial orbit ( MR ).

      CALL MEANLON(E1,L1R,M1R,MR,FRX,DR,OM,A01,A02,A03,OME1,SL,BL,SR )
      
      DLR = 385000.56d0 + SR / 1000d0                    ! Lunar distance in Apogee and Perigee(Km)


!---------------------------------------------------------------------------------------------------------
 
!------- Compute the optical libration

      INC = 1.54242d0 * DD2R        ! Inclinaction of mean lunar equator over the ecliptic(1° 32' 32".7)

      WL = (LoLu - OMEGA)* DD2R; Beta = LatLu * DD2R  
  
      XA = (Sin(WL) * Cos(Beta) * Cos(INC) - Sin(Beta) * Sin(INC))    ! XA = numerat.tan Sin (WL) ecc....
      XAA = (Cos(WL) * Cos(Beta))                                     ! XAA = denominat.
      ATXA = Atan2(XA , XAA)                                          ! Arctang. in Rad.
     
      LOL = (ATXA - F) * DR2D                                         ! Optics Libration by Long. (degrees)
                                                                        
      SIN_B = -Sin(WL) * Cos(Beta) * Sin(INC) - Sin(Beta) * Cos(INC)  ! Optics Libration by Latit. (degrees)
      LOB = asin(SIN_B)* DR2D                                          
    
      K10 = (119.75d0 + 131.849 * T) * DD2R
      K20 = (72.56 + 20.186 * T) * DD2R

!------- Compute the physical libration , RO,SIGMA,TAU

      CALL LIBRATION(E,ELP,EL,F,D,K10,K20,OMEGAR,RO,SIGMA,TAU)

       
      LFL = -TAU + (RO * Cos(ATXA) + SIGMA*Sin(ATXA))*Tan(LOB*DD2R)        ! LFL = Physical Librat. in Longit.
      LFB = SIGMA * Cos(ATXA) - RO * Sin(ATXA)                             ! LFB = Physical Librat. in Latit.
   
!------- Compute total libration
      LTL = (LOL + LFL)                                                    ! LTL = Longit.total Libration.(degrees)
      LTL = MOD(LTL,360d0)                                                 ! Selenographic Terrestrial Longit. 
                                                                           
      LTB = (LOB + LFB)                                                    ! LTB = Latit. total libration.(degrees)
                                                                         ! Selenographic Terrestrial Latit.       
      if( LTL <= -270d0) then
       LTL = (LTL + 360d0)
      else if (LTL >= 270d0) then
       LTL = LTL - 360d0 
      else   
       LTL = LTL 
      end if
                                                                  
!-------- Correction for Topocentric Libration

      TLong = atan((sin(RA7)*cos(OBL) + tan(DE7)*sin(OBL)) / cos(RA7))     ! TLong = Topocentric longitude

      TLat = sin(DE7)*cos(OBL) - cos(DE7)*sin(OBL)*sin(RA7)                ! TLat = Topocentric latitude   

!-------- Compute the Position Angle of Axis.
  

      VV = OMEGAR + DPSI + (SIGMA * DD2R /sin(INC))

      X = sin(INC + RO * DD2R) * sin(VV)
      
      Y = sin(INC+RO*DD2R)*cos(VV)*cos(OBL)-cos(INC+RO*DD2R)*sin(OBL)

      ome = atan2(X,Y)

      PA = asin((sqrt(X*X+Y*Y)*cos(RA6-ome)) / cos(LTB*DD2R))              ! PA = Position angle of axis

      PA = PA * DR2D

      if( PA < 0d0) then 
        PA = PA + 360d0
      else if (PA <= 90) then
        PA = PA
      end if


!-------- Compute the Sun's selenographic position.
!-------- HLON = Heliocentric ecliptic longit.  HLAT = Heliocentric ecliptic latit.
      
      WS = HLON - OMEGAR - DPSI  
      XS = (Sin(WS) * Cos(HLAT) * Cos(INC) - Sin(HLAT) * Sin(INC))         ! XS = numerat.tan Sin (WL) ecc....
      XSS = (Cos(WS) * Cos(HLAT))                                          ! XSS = denominat.
      ATXS = Atan2(XS , XSS)                                               ! Arctang. in Rad.
      
      ATXS= mod(ATXS,2*PI)
      
      LOH = (ATXS - F) * DR2D                                              ! Optics Libration by Helioc. Long. (degrees)
                                                                        
      SIN_H = -Sin(WS) * Cos(HLAT) * Sin(INC) - Sin(HLAT) * Cos(INC)       ! Optics Libration by Helioc. Latit. (degrees)
      LBH = asin(SIN_H)* DR2D                                          

      LFH = -TAU + (RO * Cos(ATXS) + SIGMA*Sin(ATXS))*Tan(LBH*DD2R)        ! LFH = Helioc. physical Librat. in Longit.
      HFB = SIGMA * Cos(ATXS) - RO * Sin(ATXS)                             ! HFB = Helioc. physical Librat. in Latit.
   
!------- Compute total libration
      THL = (LOH + LFH)                                                    ! THL = Longit.total Libration.(degrees)
                                                                           ! Selenographic Sun's Longit. 
      THB = (LBH + HFB)                                                    ! THB = Latit. total libration.(degrees)
                                                                           ! Selenographic Sun's Latit. 
      THL = 450d0 - THL 

      THL=mod(THL,360d0)  
                                                               
      END
!-------------------------------------------------------------------------------------------------------

      SUBROUTINE TIMEPHASE(K,N,JEX)
    


      implicit none

! Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
      integer N
      double precision K,JDE,E0,MMR,MLR,FFR,OMR,A1,Tn,JE,KX
      double precision A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14
      double precision JE0,JE1,JEP,JEU,W,K0,JEX,JH00,JH10,JH20,JH30 
      double precision JDE0,JDE1,JDEP,JDEU,JH,JH0,JH1,JHP,JHU
      
10    Tn = K / 1236.85d0                                                ! time related to the lunar phases
!      print*," K---N ",K,N
      JDE = 2451550.09766d0 + 29.530588861d0 * K + 0.00015437d0 * 
     .     Tn**2d0 - 0.00000015d0 * Tn**3d0 + 0.00000000073d0 * Tn**4d0 ! time of mean phases of the Moon
      E0 = 1d0 - 0.002516d0 * Tn - 0.0000074d0 * Tn**2d0

!------MMR = Mean Anomaly the Sun at time JDE.
      MMR = (2.5534d0+29.1053567d0*K-0.0000014d0*Tn**2d0-
     .      0.00000011d0 * Tn**3d0) * DD2R                            ! radians
!------MLR = Mean Anomaly of the Moon.
      MLR = (201.5643d0 + 385.81693528d0 * K + 0.0107582d0 * Tn**2d0
     .  + 0.00001238d0 * Tn**3d0 - 0.000000058d0 * Tn**4d0) * DD2R    ! radians
!------FFR = Argument Lunar Latitude. 
      FFR = (160.7108d0 + 390.67050284d0 * K - 0.0016118d0 * Tn**2d0
     .  - 0.00000227d0 * Tn**3d0 + 0.000000011d0 * Tn**4d0) * DD2R      ! radians
!------OMR = Longit.of Lunar ascending node .
      OMR =(124.7746d0 - 1.56375588d0 * K + 0.0020672d0 * Tn**2d0 +
     .   0.00000215d0 * Tn**3d0) * DD2R                              ! radians                    
     
!------ Planetary argument.  
      A1 = sin((299.77d0 + 0.107408d0 * K-0.009173d0 * Tn**2d0)*DD2R)
      A2 = sin((251.88d0 + 0.016321d0 * K) * DD2R)
      A3 = sin((251.83d0 + 26.651886d0 * K) * DD2R)
      A4 = sin((349.42d0 + 36.412478d0 * K) * DD2R)
      A5 = sin((84.66d0 + 18.206239d0 * K) * DD2R)
      A6 = sin((141.74d0 + 53.303771d0 * K) * DD2R)
      A7 = sin((207.14d0 + 2.453732d0 * K) * DD2R)
      A8 = sin((154.84d0 + 7.30686d0 * K) * DD2R)
      A9 = sin((34.52d0 + 27.261239d0 * K) * DD2R)
      A10 = sin((207.19d0 + 0.121824d0 * K) * DD2R)
      A11 = sin((291.34d0 + 1.844379d0 * K) * DD2R)
      A12 = sin((161.72d0 + 24.198154d0 * K) * DD2R)
      A13 = sin((239.56d0 + 25.513099d0 * K) * DD2R)
      A14 = sin((331.55d0 + 3.592518d0 * K) * DD2R)

!  'Subroutine per il calcolo della Luna Nuova
      JE0= -0.4072d0 * Sin(MLR)
      JE0= JE0 + 0.17241d0 * E0 * Sin(MMR)
      JE0= JE0 + 0.01608d0 * Sin(2d0 * MLR)
      JE0= JE0 + 0.01039d0 * Sin(2d0 * FFR)
      JE0= JE0 + 0.00739d0 * E0 * Sin(MLR - MMR)
      JE0= JE0 - 0.00514d0 * E0 * Sin(MLR + MMR)
      JE0= JE0 + 0.00208d0 * Sin(2d0 * MMR) * E0**2d0
      JE0= JE0 - 0.00111d0 * Sin(MLR - 2d0 * FFR)
      JE0= JE0 - 0.00057d0 * Sin(MLR + 2d0 * FFR)
      JE0= JE0 + 0.00056d0 * E0 * Sin(2d0 * MLR + MMR)
      JE0= JE0 - 0.00042d0 * Sin(3d0 * MLR)
      JE0= JE0 + 0.00042d0 * E0 * Sin(MMR + 2d0 * FFR)
      JE0= JE0 + 0.00038d0 * E0 * Sin(MMR - 2d0 * FFR)
      JE0= JE0 - 0.00024d0 * E0 * Sin(2d0 * MLR - MMR)
      JE0= JE0 - 0.00017d0 * Sin(OMR)
      JE0= JE0 - 0.00007d0 * Sin(MLR + 2d0 * MMR) ! * E0**2d0
      JE0= JE0 + 0.00004d0 * Sin(2d0 * MLR - 2d0 * FFR)
      JE0= JE0 + 0.00004d0 * Sin(3d0 * MMR)
      JE0= JE0 + 0.00003d0 * Sin(MLR + MMR - 2d0 * FFR)
      JE0= JE0 + 0.00003d0 * Sin(2d0 * MLR + 2d0 * FFR)
      JE0= JE0 - 0.00003d0 * Sin(MLR + MMR + 2d0 * FFR)
      JE0= JE0 + 0.00003d0 * Sin(MLR - MMR + 2d0 * FFR)
      JE0= JE0 - 0.00002d0 * Sin(MLR - MMR - 2d0 * FFR)
      JE0= JE0 - 0.00002d0 * Sin(3d0 * MLR + MMR)
      JE0= JE0 + 0.00002d0 * Sin(4d0 * MLR)
      

! 'Subroutine per il calcolo del 1° quarto di Luna .

      JE1= -0.62801d0 * Sin(MLR)
      JE1= JE1 + 0.17172d0 * E0 * Sin(MMR)
      JE1= JE1 - 0.01183d0 * E0 * Sin(MMR + MLR)
      JE1= JE1 + 0.00862d0 * Sin(2d0 * MLR)
      JE1= JE1 + 0.00804d0 * Sin(2D0 * FFR)
      JE1= JE1 + 0.00454d0 * E0 * Sin(MLR - MMR)
      JE1= JE1 + 0.00204d0 * Sin(2d0 * MMR) * E0**2d0
      JE1= JE1 - 0.00180d0 * Sin(MLR - 2d0 * FFR)     ! <<<<<<<<<<<<<<<<<<<
      JE1= JE1 - 0.00070d0 * Sin(MMR + 2d0 * FFR)
      JE1= JE1 - 0.00040d0 * Sin(3d0 * MLR )
      JE1= JE1 - 0.00034d0 * E0 * Sin(2d0 * MLR - MMR)
      JE1= JE1 + 0.00032d0 * E0 * Sin(MMR + 2d0 * FFR)
      JE1= JE1 + 0.00032d0 * E0 * Sin(MMR - 2d0 * FFR)
      JE1= JE1 - 0.00028d0 * Sin(MLR + 2D0 * MMR) * E0**2D0 
      JE1= JE1 + 0.00027d0 * E0 * Sin(2D0 * MLR + MMR)  
      JE1= JE1 - 0.00017d0 * Sin(OMR)      
      JE1= JE1 - 0.00005d0 * Sin(MLR - MMR - 2D0 * FFR)
      JE1= JE1 + 0.00004d0 * Sin(2d0 * MLR + 2D0 * FFR)
      JE1= JE1 - 0.00004d0 * Sin(MLR + MMR + 2d0 * FFR)
      JE1= JE1 + 0.00004d0 * Sin(MLR - 2d0 * MMR)
      JE1= JE1 + 0.00003d0 * Sin(MLR + MMR - 2d0 * FFR)
      JE1= JE1 + 0.00003d0 * Sin(3D0 * MMR)
      JE1= JE1 + 0.00002d0 * Sin(2D0 * MLR - 2d0 * FFR)
      JE1= JE1 + 0.00002d0 * Sin(MLR - MMR + 2d0 * FFR)
      JE1= JE1 - 0.00002d0 * Sin(3d0 * MLR + MMR)
!      print*,"JE!",JE1
!------- Correction for First and Last Quarters

      W = 0.00306 - 0.00038*E0*cos(MMR) + 0.00026*cos(MLR) 
     .   -0.00002*cos(MLR-MMR) + 0.00002*cos(MLR+MMR)
     .   + 0.00002*cos(2d0*FFR)
!      print*,"W ",W  
      JE1 = JE1 + W    ! +W for First Quarter

      JEU = JE1 - W    ! -W for Last Quarter

!  'Subroutine per il calcolo della Luna Piena
      JEP= -0.40614d0 * Sin(MLR)
      JEP= JEP + 0.17302d0 * E0 * Sin(MMR)
      JEP= JEP + 0.01614d0 * Sin(2d0 * MLR)
      JEP= JEP + 0.01043d0 * Sin(2d0 * FFR)
      JEP= JEP + 0.00734d0 * E0 * Sin(MLR - MMR)
      JEP= JEP - 0.00515d0 * E0 * Sin(MLR + MMR)
      JEP= JEP + 0.00209d0 * Sin(2d0 * MMR) * E0**2d0
      JEP= JEP - 0.00111d0 * Sin(MLR - 2d0 * FFR)
      JEP= JEP - 0.00057d0 * Sin(MLR + 2d0 * FFR)
      JEP= JEP + 0.00056d0 * E0 * Sin(2d0 * MLR + MMR)
      JEP= JEP - 0.00042d0 * Sin(3d0 * MLR)
      JEP= JEP + 0.00042d0 * E0 * Sin(MMR + 2d0 * FFR)
      JEP= JEP + 0.00038d0 * E0 * Sin(MMR - 2d0 * FFR)
      JEP= JEP - 0.00024d0 * E0 * Sin(2d0 * MLR - MMR)
      JEP= JEP - 0.00017d0 * Sin(OMR)
      JEP= JEP - 0.00007d0 * Sin(MLR + 2d0 * MMR) ! * E0**2d0
      JEP= JEP + 0.00004d0 * Sin(2d0 * MLR - 2d0 * FFR)
      JEP= JEP + 0.00004d0 * Sin(3d0 * MMR)
      JEP= JEP + 0.00003d0 * Sin(MLR + MMR - 2d0 * FFR)
      JEP= JEP + 0.00003d0 * Sin(2d0 * MLR + 2d0 * FFR)
      JEP= JEP - 0.00003d0 * Sin(MLR + MMR + 2d0 * FFR)
      JEP= JEP + 0.00003d0 * Sin(MLR - MMR + 2d0 * FFR)
      JEP= JEP - 0.00002d0 * Sin(MLR - MMR - 2d0 * FFR)
      JEP= JEP - 0.00002d0 * Sin(3d0 * MLR + MMR)
      JEP= JEP + 0.00002d0 * Sin(4d0 * MLR)

!------- Additional correction for all phases
      JH= 0.000325d0 * A1
      JH= JH + 0.000165d0 * A2
      JH= JH + 0.000164d0 * A3
      JH= JH + 0.000126d0 * A4
      JH= JH + 0.00011d0 * A5
      JH= JH + 0.000062d0 * A6
      JH= JH + 0.00006d0 * A7
      JH= JH + 0.000056d0 * A8
      JH= JH + 0.000047d0 * A9
      JH= JH + 0.000042d0 * A10
      JH= JH + 0.00004d0 * A11
      JH= JH + 0.000037d0 * A12
      JH= JH + 0.000035d0 * A13
      JH= JH + 0.000023d0 * A14
!      print*,"IH ",JH
      IF (N == 0) THEN
        JH00 = JH 
        JDE0 = JDE       
      ELSE IF (N == 1) THEN
        JH10 = JH
        JDE1 = JDE 
      ELSE IF (N == 2) THEN
        JH20 = JH
        JDEP = JDE 
      ELSE IF (N == 3) THEN 
        JH30 = JH
        JDEU = JDE  
      END IF     
!       print*,"JDE0 + JE0 + JH00",JDE0,JE0,JH00
      JH0 = JDE0 + JE0 + JH00                       !Add a JDE the corrections, for julian day of New Moon
      JH1 = JDE1 + JE1 + JH10                       !Add a JDE the corrections, for julian day of First Quarter
      JHP = JDEP + JEP + JH20                       !Add a JDE the corrections, for julian day of Full Moon  
      JHU = JDEU + JEU + JH30                       !Add a JDE the corrections, for julian day of Last Quarter

      IF (N == 0) THEN
        JEX = JH0        
      ELSE IF (N == 1) THEN
        JEX = JH1
      ELSE IF (N == 2) THEN
        JEX = JHP
      ELSE IF (N == 3) THEN 
        JEX = JHU 
      END IF     

      END

!--------------------------------------------------------------------------------------------------------


      SUBROUTINE MEANLON(E,L1R,M1R,MR,FR,DR,OMEGAR,A01,A02,A03,
     .           OME0,SL,BL,SR)

      IMPLICIT NONE

!   INPUT
!     E     =  Additional Coeff. due at eccentric terrestrial orbit. 
!     M1R   = MEAN ANOMALY OF THE MOON IN RADIANS.
!     MR    = MEAN ANOMALY OF THE SUN IN RADIANS.
!     FR    = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
!             OF THE MOON'S ASCENDING NODE IN RADIANS.
!     DR    = MEAN ELONGATION OF THE MOON FROM THE SUN IN RADIANS.
!     OMEGAR = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE IN RADIANS.
!     A01,A02,A03 = ADDITIONAL PLANETARY ARGOMENT. 
!   OUTPUT 
!     OME0  = LONGITUDE OF THE MOON'S ASCENDING NODE CORRECTED.(degrees)
!     SL    = Diff. MOON'S LONGITUDE. (0.000001 degs)
!     BL    = Diff. MOON'S LATITUDE.  (0.000001 degs) 
!     SR    = Diff. MOON'S DISTANCE (meters)
! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

      DOUBLE PRECISION E,M1R,MR,FR,DR,OMEGAR,A01,A02,A03,OME0,SL
      DOUBLE PRECISION BL,SR,L1R
!  Compute periodic term for the lunar longitude (coeff. is 0.000001 degree)

      SL = 6288774d0 * Sin(M1R)
      SL = SL + 1274027d0 * Sin(2d0 * DR - M1R)
      SL = SL + 658314d0 * Sin(2d0 * DR)
      SL = SL + 213618d0 * Sin(2d0 * M1R)
      SL = SL - 185116d0 * Sin(MR) * E
      SL = SL - 114332d0 * Sin(2d0 * FR)
      SL = SL + 58793d0 * Sin(2d0 * DR - 2d0 * M1R)
      SL = SL + 57066d0 * Sin(2d0 * DR - MR - M1R) * E
      SL = SL + 53322d0 * Sin(2d0 * DR + M1R)
      SL = SL + 45758d0 * Sin(2d0 * DR - MR) * E
      SL = SL - 40923d0 * Sin(MR - M1R) * E
      SL = SL - 34720d0 * Sin(DR)
      SL = SL - 30383d0 * Sin(MR + M1R) * E
      SL = SL + 15327d0 * Sin(2d0 * DR - 2d0 * FR)
      SL = SL - 12528d0 * Sin(M1R + 2d0 * FR)
      SL = SL + 10980d0 * Sin(M1R - 2d0 * FR)
      SL = SL + 10675d0 * Sin(4d0 * DR - M1R)
      SL = SL + 10034d0 * Sin(3d0 * M1R)
      SL = SL + 8548d0 * Sin(4d0 * DR - 2d0 * M1R)
      SL = SL - 7888d0 * Sin(2d0 * DR + MR - M1R) * E
      SL = SL - 6766d0 * Sin(2d0 * DR + MR) * E
      SL = SL - 5163d0 * Sin(DR - M1R)
      SL = SL + 4987d0 * Sin(DR + MR) * E
      SL = SL + 4036d0 * Sin(2d0 * DR - MR + M1R) * E
      SL = SL + 3994d0 * Sin(2d0 * DR + 2d0 * M1R)
      SL = SL + 3861d0 * Sin(4d0 * DR)
      SL = SL + 3665d0 * Sin(2d0 * DR - 3d0 * M1R)
      SL = SL - 2689d0 * Sin(MR - 2d0 * M1R) * E
      SL = SL - 2602d0 * Sin(2d0 * DR - M1R + 2d0 * FR)
      SL = SL + 2390d0 * Sin(2d0 * DR - MR - 2d0 * M1R) * E
      SL = SL - 2348d0 * Sin(DR + M1R)
      SL = SL + 2236d0 * Sin(2d0 * DR - 2d0 * MR) * E**2d0
      SL = SL - 2120d0 * Sin(MR + 2d0 * M1R) * E
      SL = SL - 2069d0 * Sin(2d0 * MR) * E**2d0
      SL = SL + 2048d0 * Sin(2d0 * DR - 2d0 * MR - M1R) * E**2d0
      SL = SL - 1773d0 * Sin(2d0 * DR + M1R - 2d0 * FR)
      SL = SL - 1595d0 * Sin(2d0 * DR + 2d0 * FR)
      SL = SL + 1215d0 * Sin(4d0 * DR - MR - M1R) * E
      SL = SL - 1110d0 * Sin(2d0 * M1R + 2d0 * FR)
      SL = SL - 892d0 * Sin(3d0 * DR - M1R)
      SL = SL - 810d0 * Sin(2d0 * DR + MR + M1R) * E
      SL = SL + 759d0 * Sin(4d0 * DR - MR - 2d0 * M1R) * E
      SL = SL - 713d0 * Sin(2d0 * MR - M1R) * E**2d0
      SL = SL - 700d0 * Sin(2d0 * DR + 2d0 * MR - M1R) * E**2d0
      SL = SL + 691d0 * Sin(2d0 * DR + MR - 2d0 * M1R) * E
      SL = SL + 596d0 * Sin(2d0 * DR - MR - 2d0 * FR) * E
      SL = SL + 549d0 * Sin(4d0 * DR + M1R)
      SL = SL + 537d0 * Sin(4d0 * M1R)
      SL = SL + 520d0 * Sin(4d0 * DR - MR) * E
      SL = SL - 487d0 * Sin(DR - 2d0 * M1R)
      SL = SL - 399d0 * Sin(2d0 * DR + MR - 2d0 * FR) * E
      SL = SL - 381d0 * Sin(2d0 * M1R - 2d0 * FR)
      SL = SL + 351d0 * Sin(DR + MR + M1R) * E
      SL = SL - 340d0 * Sin(3d0 * DR - 2d0 * M1R)
      SL = SL + 330d0 * Sin(4d0 * DR - 3d0 * M1R)
      SL = SL + 327d0 * Sin(2d0 * DR - MR + 2d0 * M1R) * E
      SL = SL - 323d0 * Sin(2d0 * MR + M1R) * E**2d0
      SL = SL + 299d0 * Sin(DR + MR - M1R) * E
      SL = SL + 294d0 * Sin(2d0 * DR + 3d0 * M1R)



!  Compute the periodic term for lunar distance DL (coeff.in 0.001 Km)
 
      SR = -20905355d0 * Cos(M1R)
      SR = SR - 3699111d0 * Cos(2d0 * DR - M1R)
      SR = SR - 2955968d0 * Cos(2d0 * DR)
      SR = SR - 569925d0 * Cos(2d0 * M1R)
      SR = SR + 48888d0 * Cos(MR) * E
      SR = SR - 3149d0 * Cos(2d0 * FR)
      SR = SR + 246158d0 * Cos(2d0 * DR - 2d0 * M1R)
      SR = SR - 152138d0 * Cos(2d0 * DR - MR - M1R) * E
      SR = SR - 170733d0 * Cos(2d0 * DR + M1R)
      SR = SR - 204586d0 * Cos(2d0 * DR - MR) * E
      SR = SR - 129620d0 * Cos(MR - M1R) * E
      SR = SR + 108743d0 * Cos(DR)
      SR = SR + 104755d0 * Cos(MR + M1R) * E
      SR = SR + 10321d0 * Cos(2d0 * DR - 2d0 * FR)
      SR = SR + 79661d0 * Cos(M1R - 2d0 * FR)
      SR = SR - 34782d0 * Cos(4d0 * DR - M1R)
      SR = SR - 23210d0 * Cos(3d0 * M1R)
      SR = SR - 21636d0 * Cos(4d0 * DR - 2d0 * M1R)
      SR = SR + 24208d0 * Cos(2d0 * DR + MR - M1R) * E
      SR = SR + 30824d0 * Cos(2d0 * DR + MR) * E
      SR = SR - 8379d0 * Cos(DR - M1R)
      SR = SR - 16675d0 * Cos(DR + MR) * E
      SR = SR - 12831d0 * Cos(2d0 * DR - MR + M1R) * E
      SR = SR - 10445d0 * Cos(2d0 * DR + 2d0 * M1R)
      SR = SR - 11650d0 * Cos(4d0 * DR)
      SR = SR + 14403d0 * Cos(2d0 * DR - 3d0 * M1R)
      SR = SR - 7003d0 * Cos(MR - 2d0 * M1R) * E
      SR = SR + 10056d0 * Cos(2d0 * DR - MR - 2d0 * M1R) * E
      SR = SR + 6322d0 * Cos(DR + M1R)
      SR = SR - 9884d0 * Cos(2d0 * DR - 2d0 * MR) * E**2d0
      SR = SR + 5751d0 * Cos(MR + 2d0 * M1R) * E
      SR = SR - 4950d0 * Cos(2d0 * DR - 2d0 * MR - M1R) * E**2d0
      SR = SR + 4130d0 * Cos(2d0 * DR + M1R - 2d0 * FR)
      SR = SR - 3958d0 * Cos(4d0 * DR - MR - M1R) * E
      SR = SR + 3258d0 * Cos(3d0 * DR - M1R)
      SR = SR + 2616d0 * Cos(2d0 * DR + MR + M1R) * E
      SR = SR - 1897d0 * Cos(4d0 * DR - MR - 2d0 * M1R) * E
      SR = SR - 2117d0 * Cos(2d0 * MR - M1R) * E**2d0
      SR = SR + 2354d0 * Cos(2d0 * DR + 2d0 * MR - M1R) * E**2d0
      SR = SR - 1423d0 * Cos(4d0 * DR + M1R)
      SR = SR - 1117d0 * Cos(4d0 * M1R)
      SR = SR - 1571d0 * Cos(4d0 * DR - MR) * E
      SR = SR - 1739d0 * Cos(DR - 2d0 * M1R)
      SR = SR - 4421d0 * Cos(2d0 * M1R - 2d0 * FR)
      SR = SR + 1165d0 * Cos(2d0 * MR + M1R) * E**2d0
      SR = SR + 8752d0 * Cos(2d0 * DR - M1R - 2d0 * FR)



!   Compute the periodic term for lunar latitude BL (coeff.in 0.000001 Gradi)

      BL = 5128122d0 * Sin(FR)
      BL = BL + 280602d0 * Sin(M1R + FR)
      BL = BL + 277693d0 * Sin(M1R - FR)
      BL = BL + 173237d0 * Sin(2d0 * DR - FR)
      BL = BL + 55413d0 * Sin(2d0 * DR - M1R + FR)
      BL = BL + 46271d0 * Sin(2d0 * DR - M1R - FR)
      BL = BL + 32573d0 * Sin(2d0 * DR + FR)
      BL = BL + 17198d0 * Sin(2d0 * M1R + FR)
      BL = BL + 9266d0 * Sin(2d0 * DR + M1R - FR)
      BL = BL + 8822d0 * Sin(2d0 * M1R - FR)
      BL = BL + 8216d0 * Sin(2d0 * DR - MR - FR) * E
      BL = BL + 4324d0 * Sin(2d0 * DR - 2d0 * M1R - FR)
      BL = BL + 4200d0 * Sin(2d0 * DR + M1R + FR)
      BL = BL - 3359d0 * Sin(2d0 * DR + MR - FR) * E
      BL = BL + 2463d0 * Sin(2d0 * DR - MR - M1R + FR) * E
      BL = BL + 2211d0 * Sin(2d0 * DR - MR + FR) * E
      BL = BL + 2065d0 * Sin(2d0 * DR - MR - M1R - FR) * E
      BL = BL - 1870d0 * Sin(MR - M1R - FR) * E
      BL = BL + 1828d0 * Sin(4d0 * DR - M1R - FR)
      BL = BL - 1794d0 * Sin(MR + FR) * E
      BL = BL - 1749d0 * Sin(3d0 * FR)
      BL = BL - 1565d0 * Sin(MR - M1R + FR) * E
      BL = BL - 1491d0 * Sin(DR + FR)
      BL = BL - 1475d0 * Sin(MR + M1R + FR) * E
      BL = BL - 1410d0 * Sin(MR + M1R - FR) * E
      BL = BL - 1344d0 * Sin(MR - FR) * E
      BL = BL - 1335d0 * Sin(DR - FR)
      BL = BL + 1107d0 * Sin(3d0 * M1R + FR)
      BL = BL + 1021d0 * Sin(4d0 * DR - FR)
      BL = BL + 833d0 * Sin(4d0 * DR - M1R + FR)
      BL = BL + 777d0 * Sin(M1R - 3d0 * FR)
      BL = BL + 671d0 * Sin(4d0 * DR - 2d0 * M1R + FR)
      BL = BL + 607d0 * Sin(2d0 * DR - 3d0 * FR)
      BL = BL + 596d0 * Sin(2d0 * DR + 2d0 * M1R - FR)
      BL = BL + 491d0 * Sin(2d0 * DR - MR + M1R - FR) * E
      BL = BL - 451d0 * Sin(2d0 * DR - 2d0 * M1R + FR)
      BL = BL + 439d0 * Sin(3d0 * M1R - FR)
      BL = BL + 422d0 * Sin(2d0 * DR + 2d0 * M1R + FR)
      BL = BL + 421d0 * Sin(2d0 * DR - 3d0 * M1R - FR)
      BL = BL - 366d0 * Sin(2d0 * DR + MR - M1R + FR) * E
      BL = BL - 351d0 * Sin(2d0 * DR + MR + FR) * E
      BL = BL + 331d0 * Sin(4d0 * DR + FR)
      BL = BL + 315d0 * Sin(2d0 * DR - MR + M1R + FR) * E
      BL = BL + 302d0 * Sin(2d0 * DR - 2d0 * MR - FR) * E**2d0
      BL = BL - 283d0 * Sin(M1R + 3d0 * FR)
      BL = BL - 229d0 * Sin(2d0 * DR + MR + M1R - FR) * E
      BL = BL + 223d0 * Sin(DR + MR - FR) * E
      BL = BL + 223d0 * Sin(DR + MR + FR) * E
      BL = BL - 220d0 * Sin(MR - 2d0 * M1R - FR) * E
      BL = BL - 220d0 * Sin(2d0 * DR + MR - M1R - FR) * E
      BL = BL - 185d0 * Sin(DR + M1R + FR)
      BL = BL + 181d0 * Sin(2d0 * DR - MR - 2d0 * M1R - FR) * E
      BL = BL - 177d0 * Sin(MR + 2d0 * M1R + FR) * E
      BL = BL + 176d0 * Sin(4d0 * DR - 2d0 * M1R - FR)
      BL = BL + 166d0 * Sin(4d0 * DR - MR - M1R - FR) * E
      BL = BL - 164d0 * Sin(DR + M1R - FR)
      BL = BL + 132d0 * Sin(4d0 * DR + M1R - FR)
      BL = BL - 119d0 * Sin(DR - M1R - FR)
      BL = BL + 115d0 * Sin(4d0 * DR - MR - FR) * E
      BL = BL + 107d0 * Sin(2d0 * DR - 2d0 * MR + FR) * E**2d0

!  Compute the longitude term of the Moon's ascending node.
      OME0 =  OMEGAR * DR2D - 1.4979d0 * Sin(2d0 * (DR - FR))
      OME0 =  OME0 - 0.15d0 * Sin(MR)
      OME0 =  OME0 - 0.1226d0 * Sin(2d0 * DR)
      OME0 =  OME0 + 0.1176d0 * Sin(2d0 * FR)
      OME0 =  OME0 - 0.0801d0 * Sin(2d0 * (M1R - FR))


! Add the terms osculating A01,A02,A03
! for  SL
        SL = SL + 3958d0 * Sin(A01)
        SL = SL + 1962d0 * Sin(L1R - FR)
        SL = SL + 318d0 * Sin(A02)
  
! for BL
        BL = BL - 2235d0 * Sin(L1R)
        BL = BL + 382d0 * Sin(A03)
        BL = BL + 175d0 * Sin(A01 - FR)
        BL = BL + 175d0 * Sin(A01 + FR)
        BL = BL + 127d0 * Sin(L1R - M1R)
        BL = BL - 115d0 * Sin(L1R + M1R)
  


      END
!--------------------------------------------------------------------------------------------------

      SUBROUTINE PER_APOG(TAP,DAP,MPA,FAP,PE,AP,APS)

      implicit none

      double precision TAP,DAP,MPA,FAP,PE,AP,APS
!------- PE = coefficient of PERIGEE  (DAY)

      PE = -1.6769d0 * Sin(2d0 * DAP)
      PE = PE + 0.4589d0 * Sin(4d0 * DAP)
      PE = PE - 0.1856d0 * Sin(6d0 * DAP)
      PE = PE + 0.0883d0 * Sin(8d0 * DAP)
      PE = PE - (0.0773 + 0.00019d0 * TAP) * Sin(2d0 * DAP - MPA)
      PE = PE + (0.0502 - 0.00013d0 * TAP) * Sin(MPA)
      PE = PE - 0.046d0 * Sin(10 * DAP)
      PE = PE + (0.0422 - 0.00011d0 * TAP) * Sin(4d0 * DAP - MPA)
      PE = PE - 0.0256d0 * Sin(6d0 * DAP - MPA)
      PE = PE + 0.0253d0 * Sin(12d0 * DAP)
      PE = PE + 0.0237d0 * Sin(DAP)
      PE = PE + 0.0162d0 * Sin(8d0 * DAP - MPA)
      PE = PE - 0.0145d0 * Sin(14d0 * DAP)
      PE = PE + 0.0129d0 * Sin(2d0 * FAP)
      PE = PE - 0.0112d0 * Sin(3d0 * DAP)
      PE = PE - 0.0104d0 * Sin(10 * DAP - MPA)
      PE = PE + 0.0086d0 * Sin(16d0 * DAP)
      PE = PE + 0.0069d0 * Sin(12d0 * DAP - MPA)
      PE = PE + 0.0066d0 * Sin(5d0 * DAP)
      PE = PE - 0.0053d0 * Sin(2d0 * DAP + 2d0 * FAP)
      PE = PE - 0.0052d0 * Sin(18d0 * DAP)
      PE = PE - 0.0046d0 * Sin(14d0 * DAP - MPA)
      PE = PE - 0.0041d0 * Sin(7d0 * DAP)
      PE = PE + 0.004d0 * Sin(2d0 * DAP + MPA)
      PE = PE + 0.0032d0 * Sin(20 * DAP)
      PE = PE - 0.0032d0 * Sin(DAP + MPA)
      PE = PE + 0.0031d0 * Sin(16d0 * DAP - MPA)
      PE = PE - 0.0029d0 * Sin(4d0 * DAP + MPA)
      PE = PE + 0.0027d0 * Sin(9d0 * DAP)
      PE = PE + 0.0027d0 * Sin(4d0 * DAP + 2d0 * FAP)
      PE = PE - 0.0027d0 * Sin(2d0 * DAP - 2d0 * MPA)
      PE = PE + 0.0024d0 * Sin(4d0 * DAP - 2d0 * MPA)
      PE = PE - 0.0021d0 * Sin(6d0 * DAP - 2d0 * MPA)
      PE = PE - 0.0021d0 * Sin(22d0 * DAP)
      PE = PE - 0.0021d0 * Sin(18d0 * DAP - MPA)
      PE = PE + 0.0019d0 * Sin(6d0 * DAP + MPA)
      PE = PE - 0.0018d0 * Sin(11d0 * DAP)
      PE = PE - 0.0014d0 * Sin(8d0 * DAP + MPA)
      PE = PE - 0.0014d0 * Sin(4d0 * DAP - 2d0 * FAP)
      PE = PE - 0.0014d0 * Sin(6d0 * DAP + 2d0 * FAP)
      PE = PE + 0.0014d0 * Sin(3d0 * DAP + MPA)
      PE = PE - 0.0014d0 * Sin(5d0 * DAP + MPA)
      PE = PE + 0.0013d0 * Sin(13d0 * DAP)
      PE = PE + 0.0013d0 * Sin(20 * DAP - MPA)
      PE = PE + 0.0011d0 * Sin(3d0 * DAP + 2d0 * MPA)
      PE = PE - 0.0011d0 * Sin(4d0 * DAP + 2d0 * FAP - 2d0 * MPA)
      PE = PE - 0.001d0 * Sin(DAP + 2d0 * MPA)
      PE = PE - 0.0009d0 * Sin(22d0 * DAP - MPA)
      PE = PE - 0.0008d0 * Sin(4d0 * FAP)
      PE = PE + 0.0008d0 * Sin(6d0 * DAP - 2d0 * FAP)
      PE = PE + 0.0008d0 * Sin(2d0 * DAP - 2d0 * FAP + MPA)
      PE = PE + 0.0007d0 * Sin(2d0 * MPA)
      PE = PE + 0.0007d0 * Sin(2d0 * FAP - MPA)
      PE = PE + 0.0007d0 * Sin(2d0 * DAP + 4d0 * FAP)
      PE = PE - 0.0006d0 * Sin(2d0 * FAP - 2d0 * MPA)
      PE = PE - 0.0006d0 * Sin(2d0 * DAP - 2d0 * FAP + 2d0 * MPA)
      PE = PE + 0.0006d0 * Sin(24d0 * DAP)
      PE = PE + 0.0005d0 * Sin(4d0 * DAP - 4d0 * FAP)
      PE = PE + 0.0005d0 * Sin(2d0 * DAP + 2d0 * MPA)
      PE = PE - 0.0004d0 * Sin(DAP - MPA)

!------ AP = = coefficient of APOGEE  (DAY).  
      AP = 0.4392d0 * Sin(2d0 * DAP)
      AP = AP + 0.0684d0 * Sin(4d0 * DAP)
      AP = AP + 0.0456d0 * Sin(MPA) - 0.00011d0 * TAP
      AP = AP + 0.0426d0 * Sin(2d0 * DAP - MPA) - 0.00011d0 * TAP
      AP = AP + 0.0212d0 * Sin(2d0 * FAP)
      AP = AP - 0.0189d0 * Sin(DAP)
      AP = AP + 0.0144d0 * Sin(6d0 * DAP)
      AP = AP + 0.0113d0 * Sin(4d0 * DAP - MPA)
      AP = AP + 0.0047d0 * Sin(2d0 * DAP + 2d0 * FAP)
      AP = AP + 0.0036d0 * Sin(DAP + MPA)
      AP = AP + 0.0035d0 * Sin(8d0 * DAP)
      AP = AP + 0.0034d0 * Sin(6d0 * DAP - MPA)
      AP = AP - 0.0034d0 * Sin(2d0 * DAP - 2d0 * FAP)
      AP = AP + 0.0022d0 * Sin(2d0 * DAP - 2d0 * MPA)
      AP = AP - 0.0017d0 * Sin(3d0 * DAP)
      AP = AP + 0.0013d0 * Sin(4d0 * DAP + 2d0 * FAP)
      AP = AP + 0.0011d0 * Sin(8d0 * DAP - MPA)
      AP = AP + 0.001d0 * Sin(4d0 * DAP - 2d0 * MPA)
      AP = AP + 0.0009d0 * Sin(10d0 * DAP)
      AP = AP + 0.0007d0 * Sin(3d0 * DAP + MPA)
      AP = AP + 0.0006d0 * Sin(2d0 * MPA)
      AP = AP + 0.0005d0 * Sin(2d0 * DAP + MPA)
      AP = AP + 0.0005d0 * Sin(2d0 * DAP + 2d0 * MPA)
      AP = AP + 0.0004d0 * Sin(6d0 * DAP + 2d0 * FAP)
      AP = AP + 0.0004d0 * Sin(6d0 * DAP - 2d0 * MPA)
      AP = AP + 0.0004d0 * Sin(10d0 * DAP - MPA)
      AP = AP - 0.0004d0 * Sin(5d0 * DAP)
      AP = AP - 0.0004d0 * Sin(4d0 * DAP - 2d0 * FAP)
      AP = AP + 0.0003d0 * Sin(2d0 * FAP + MPA)
      AP = AP + 0.0003d0 * Sin(12d0 * DAP)
      AP = AP + 0.0003d0 * Sin(2d0 * DAP + 2d0 * FAP - MPA)
      AP = AP - 0.0003d0 * Sin(DAP - MPA)

!------ AP = = coefficient of APOGEE  (arcseconds). 

      APS =3245.251d0
      APS = APS - 9.147d0 * cos(2d0*DAP)
      APS = APS - 0.841d0 * cos(DAP)
      APS = APS + 0.697d0 * cos(2d0*FAP)
      APS = APS -(0.656d0 + 0.0016d0*TAP) * cos(MPA)
      APS = APS + 0.355d0 * cos(4d0*DAP)
      APS = APS + 0.159d0 * cos(2d0*DAP-MPA)
      APS = APS + 0.127d0 * cos(DAP+MPA)
      APS = APS + 0.065d0 * cos(4d0*DAP-MPA)
      APS = APS + 0.052d0 * cos(6d0*DAP)
      APS = APS + 0.043d0 * cos(2d0*DAP+MPA)
      APS = APS + 0.031d0 * cos(2d0*DAP+2d0*FAP) 
      APS = APS + 0.023d0 * cos(2d0*DAP-2d0*FAP)
      APS = APS + 0.022d0 * cos(2d0*DAP-2d0*MPA)
      APS = APS + 0.019d0 * cos(2d0*DAP+2d0*MPA)
      APS = APS - 0.016d0 * cos(2d0*MPA)
      APS = APS + 0.014d0 * cos(6d0*DAP-4d0*MPA)
      APS = APS - 0.010d0 * cos(8d0*DAP)

      END
!-----------------------------------------------------------------------------------------------

      subroutine LIBRATION(E,ELP,EL,F,D,K1,K2,OMEGAR,RO,SIGMA,TAU) 

      implicit none 

      double precision E,ELP,EL,F,D,K1,K2,OMEGAR,RO,SIGMA,TAU 

      RO = -0.02752d0 * Cos(EL )
      RO = RO - 0.02245d0 *  Sin(F )
      RO = RO + 0.00684d0 *  Cos(EL  - 2d0 * F )
      RO = RO - 0.00293d0 *  Cos(2d0 * F )
      RO = RO - 0.00085d0 *  Cos(2d0 * F  - 2d0 * D )
      RO = RO - 0.00054d0 *  Cos(EL  - 2d0 * D )
      RO = RO - 0.0002d0 *  Sin(EL  + F )
      RO = RO - 0.0002d0 *  Cos(EL  + 2d0 * F )
      RO = RO - 0.0002d0 *  Cos(EL  - F )
      RO = RO + 0.00014d0 *  Cos(EL  + 2d0 * F  - 2d0 * D )
     
      SIGMA = -0.02816d0 *  Sin(EL )
      SIGMA = SIGMA + 0.02244d0 *  Cos(F )
      SIGMA = SIGMA - 0.00682d0 *  Sin(EL  - 2d0 * F )
      SIGMA = SIGMA - 0.00279d0 *  Sin(2d0 * F )
      SIGMA = SIGMA - 0.00083d0 *  Sin(2d0 * F  - 2d0 * D )
      SIGMA = SIGMA + 0.00069d0 *  Sin(EL  - 2d0 * D )
      SIGMA = SIGMA + 0.0004d0 *  Cos(EL  + F )
      SIGMA = SIGMA - 0.00025d0 *  Sin(2d0 * EL )
      SIGMA = SIGMA - 0.00023d0 *  Sin(EL  + 2d0 * F )
      SIGMA = SIGMA + 0.0002d0 *  Cos(EL  - F )
      SIGMA = SIGMA + 0.00019d0 *  Sin(EL  - F )
      SIGMA = SIGMA + 0.00013d0 *  Sin(EL  + 2d0 * F  - 2d0 * D )
      SIGMA = SIGMA - 0.0001d0 *  Cos(EL  - 3d0 *  F )
     
     
      TAU = 0.0252d0 * E * Sin(ELP )
      TAU = TAU + 0.00473d0 *  Sin(2d0 * EL  - 2d0 * F )
      TAU = TAU - 0.00467d0 *  Sin(EL )
      TAU = TAU + 0.00396d0 *  Sin(K1)
      TAU = TAU + 0.00276d0 *  Sin(2d0 * EL  - 2d0 * D )
      TAU = TAU + 0.00196d0 *  Sin(OMEGAR)
      TAU = TAU - 0.00183d0 *  Cos(EL  - F )
      TAU = TAU + 0.00115d0 *  Sin(EL  - 2d0 * D )
      TAU = TAU - 0.00096d0 *  Sin(EL  - D )
      TAU = TAU + 0.00046d0 *  Sin(2d0 * F  - 2d0 * D )
      TAU = TAU - 0.00039d0 *  Sin(EL  - F )
      TAU = TAU - 0.00032d0 *  Sin(EL  - ELP  - D )
      TAU = TAU + 0.00027d0 *  Sin(2d0 * EL  - ELP  - 2d0 * D )
      TAU = TAU + 0.00023d0 *  Sin(K2)
      TAU = TAU - 0.00014d0 *  Sin(2d0 * D )
      TAU = TAU + 0.00014d0 *  Cos(2d0 * EL  - 2d0 * F )
      TAU = TAU - 0.00012d0 *  Sin(EL  - 2d0 * F )
      TAU = TAU - 0.00012d0 *  Sin(2d0 * EL )
      TAU = TAU + 0.00011d0 *  Sin(2d0 * EL  - 2d0 * ELP  - 2d0 * D )
     
      END
!-----------------------------------------------------------------------------------------------------
      
      subroutine GEPV (GMB,OBL,J,LON,LAT)

!     This routine rotate an r-matrix about the x-axis
!     and discarge velocity vector , tranform the position
!     vector found  to longitude and latitude .      

!     GMB(3,2)     Helioc. state vector.(input)
!     OBM          Ecliptic obliquity. (input)
!     LON          Heliocentric ecliptic longitude.(output) 
!     LAT          Heliocentric ecliptic latitude. (output) 
!     GEP(3)       Geocentric ecliptic position.
!     iau_C2S      unit vector to spherical coordinates.
!     iau_ANP      normalize angle into range 0 to 2pi.


      implicit none

! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )
! 2Pi
      DOUBLE PRECISION D2PI
      PARAMETER (D2PI = 6.283185307179586476925287D0 )

      double precision GMB(3,2),LON,LAT,OBL,GEP(3),DT,iau_ANP
      integer J

      DT = 0.0d0

      if(J == 1)then
        call iau_RX(OBL,GMB)
      else if (J == 2) then
        call iau_RY(OBL,GMB)
      else if (J == 3) then
        call iau_RZ(OBL,GMB) 
      end if
      call iau_PVUP(DT,GMB,GEP)
      call iau_C2S(GEP,LON,LAT)
      LON = iau_ANP(LON)

      end
!---------------------------------------------------------------------------

      subroutine kepler( Mean,ecce,E)

C ------------solve kepler equation -------- --------------------------
C Equation: E1 = E0 + (M + e*sin(E0) - E0) / (1 - e*cos(E0)) ,
C being  ecce = Orbit Eccentricity E= Eccentric Anomaly, M= Mean Anomaly
C------------------------------------------------------------------
      implicit none

! Degrees to radians
      DOUBLE PRECISION DD2R
      PARAMETER ( DD2R = 1.745329251994329576923691D-2 )
!
      double precision Mean,tol,Xe,ff,gg,Ex,diff,E,ecce,M

      M = Mean*DD2R 
      tol = 1 * 10E-12
      XE = M       ! first approach
      do 
       ff = M + ecce * sin(XE) - XE
       gg = 1d0 - ecce * cos(XE)
       Ex = XE +( ff/gg)
       diff = (Ex-XE)  
       if ( abs(diff) < tol) then
        XE = Ex 
        goto 300
       else
        XE = XE + diff
       end if
      end do
300   continue
      E = XE  
   
      end
C----------------------------------------------------------------------------

!---------------------------------------------------------------------------

      SUBROUTINE AN_HMS ( ANGLE, IH, IM, SEC) 

      ! Transform angle radians in hours, minuts seconds and fraction  

      IMPLICIT NONE
!----  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )
!----- Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )
      DOUBLE PRECISION ANGLE,SEC,IHD,IMD,RD,W,DEGS
      INTEGER IH,IM

      W = MOD(ANGLE,D2PI)
      IF ( W .LT. 0D0 ) W = W + D2PI


      DEGS = W * DR2D         ! radians to degrees
      IHD = DEGS /15d0         ! degrees to hour   
      IH = int(IHD)
      IMD = (IHD - IH) * 60d0  ! fraction hour to minutes
      IM = int(IMD)            
      SEC = (IMD - IM) * 60d0  ! fraction minutes to seconds 
     
      END 
!----------------------------------------------------------------------------

      SUBROUTINE AN_DMS ( ANGLE, sign, DG, DM, DSEC) 

!    Transform radians in degrees, minuts seconds of arc and fraction  

      IMPLICIT NONE
!----  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )
!----- Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )
      DOUBLE PRECISION ANGLE,SEC,IMD,W,DSEC,DEGS
      INTEGER DG,DM
      CHARACTER*1  sign 

      W = MOD(ANGLE,D2PI)
!      IF ( W .LT. 0D0 ) W = W + D2PI


      DEGS = abs(W * DR2D)      ! radians to degrees
      DG = int(DEGS)
      IMD = (DEGS - DG) * 60d0  ! fraction degrees to minutes of arc
      DM = int(IMD)            
      DSEC = (IMD - DM) * 60d0  ! fraction minutes to seconds of arc

      IF ( ANGLE .LT. 0D0  ) THEN
        sign = "-"
       ELSE IF ( ANGLE .GE. 0d0) THEN
        sign = "+"
      END IF
     
      END 

C--------------------------------------------------------------------------------------------

      subroutine TWILIGHT(NTARG,TDB0,GST0,AUDAY,DELTAT,LAT,LONG,
     . DX00,DY00,Ha,tw1,tw2,tw3,tw4,tw5,tw6,n1,n2,n3)

!     This program compute the times of transit,rise and set of body
!     of the solar sistem , this algorithm takes in account the atmospheric
!     refraction and all the parameters relate.     

!     INPUT:
!           NTARG:    is Solar System Body number.              
!           TDB0 :    Barycentric Dinamical Time to 0.0h UT1.
!           GST0 :    Greenwich Apparent Sideral Time at 0.0h UT1 
!           AUDAY:    Speed of light. (AU per day  "TDB")
!           DELTAT:   DeltaT in fraction day.
!           LAT,LONG: Latitude and Longitude of observer. (Rad.)
!           DX00,DY00:Celestial pole offsets (dX, dY; Rad.) 
!           
!     OUTPUT
!           tw1:      Begin civil twilight.
!           tw2:      End civil twilight.
!           tw3:      Begin nutical twilight.
!           tw4:      End nautical twilight.
!           tw5:      Begin astronomical twilight.
!           tw6:      end astronomical twilight.
 
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


      INTEGER NTARG,C,I,n1,n2,n3,flag,L,M,N
      DOUBLE PRECISION TDB0,GST0,AUDAY,LAT,LONG,H(3),Z(3)
      DOUBLE PRECISION RA3,R2,R3,R6,R8,R9,R10,DE3,D2,D3,D6,D8,D9,D10
      DOUBLE PRECISION H00,H0,delta,Ha,m12,m13,m22,m23,tw6
      DOUBLE PRECISION COH00,DELTAT,h_,m1,m2,newm,oldm,m0
      DOUBLE PRECISION DX00,DY00,H02,H03,azm,HN,HAS,COH2,COH3,COH1
      DOUBLE PRECISION HC,H01,m11,m21,tw1,tw2,tw3,tw4,tw5
            
      if(NTARG .EQ. 11) then
       H0 = (90d0 ) * DD2R  !+ 0.833333333d0
       HC = (96d0 )* DD2R   
       HN = 102d0 * DD2R
       HAS = 108d0 * DD2R
      endif
      
      H(1) = HC
      H(2) = HN
      H(3) = HAS

      flag = -1
      if(LAT .GT. 1.4 .OR. LAT .LT. -1.4) goto 300

      call COORD(NTARG,TDB0,AUDAY,DX00,DY00,RA3,R2,R3,R6,R8,R9,       
     .  R10,DE3,D2,D3,D6,D8,D9,D10)
       
      COH00 = (cos(H0) - sin(LAT) * sin(DE3)) / (cos(LAT) * cos(DE3))       !Horizon observer
      m0 = (RA3 - LONG - GST0) / D2PI
      m0 = modulo (m0,1d0)
      H00 = acos(COH00)                                                     !Horizon observer
      H00 = modulo(H00,PI)
      m1 = m0 - H00 / D2PI
      m1 = modulo(m1,1d0)
      m2 = m0 + H00 / D2PI
      m2 = modulo(m2,1d0)  
                                                 
!------Compute Civil Twilight -------------------------------------------------------------------

      COH1 = (cos(HC) - sin(LAT) * sin(DE3)) / (cos(LAT) * cos(DE3))      !Civil twilight
      
       IF(COH1 .LE. -1.0d0) THEN
        tw1 = 0; tw2 = 0; n1 = 1 
        GOTO 300
       ELSE
        tw1 = tw1; tw2 = tw2; n1 = 0
       END IF

      H01 = acos(COH1)                     
      H01 = modulo(H01,PI)
      m11 = m0 - H01 / D2PI                                                 ! Begin civil twilight
      m11 = modulo(m11,1d0)
      m11 = m11
      C = 1
      call Interpol(GST0,m11,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! begin civil twilight = tw1
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HC,tw1,Ha,delta,azm)
      tw1 = tw1 
      
      m21 = m0 + H01 / D2PI                                                 ! End civil twilight
      m21 = modulo(m21,1d0)
      m21 = m21
      C = 1
       call Interpol(GST0,m21,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! end civil twilight = tw2
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HC,tw2,Ha,delta,azm)
      tw2 = tw2 

!------Compute Nautical Twilight -------------------------------------------------------------------


      COH2 = (cos(HN) - sin(LAT) * sin(DE3)) / (cos(LAT) * cos(DE3))      !Nautical twilight
      
       IF(COH2 .LE. -1.0d0) THEN
        tw3 = 0; tw4 = 0;  n2 = 1
        GOTO 300
       ELSE
        tw3 = tw3; tw4 = tw4; n2 = 0
       END IF

      H02 = acos(COH2)                     !
      H02 = modulo(H02,PI)
      m12 = m0 - H02 / D2PI                                                 ! Begin nutical twilight
      m12 = modulo(m12,1d0)
      m12 = m12
      C = 1 
       call Interpol(GST0,m12,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! begin nautical twilight = tw3
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HN,tw3,Ha,delta,azm)
      tw3 = tw3 
       
      m22 = m0 + H02 / D2PI                                                 ! End nautical twilight
      m22 = modulo(m22,1d0)
      m22 = m22
      C = 1
       call Interpol(GST0,m22,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! end nautical twilight = tw4
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HN,tw4,Ha,delta,azm)
      tw4 = tw4 

!------Compute Astronomical Twilight -------------------------------------------------------------------

      COH3 = (cos(HAS) - sin(LAT) * sin(DE3)) / (cos(LAT) * cos(DE3))     !Astronomical twilight
      
       IF(COH3 .LE. -1.0d0) THEN
        tw5 = 0; tw6 = 0; n3 = 1
        GOTO 300
       ELSE
        tw5 = tw5; tw6 = tw6; n3 = 0
       END IF

      H03 = acos(COH3) 
      H03 = modulo(H03,PI)
      m13 = m0 - H03 / D2PI                ! Begin astronomical twilight
      m13 = modulo(m13,1d0)      
      m13 = m13
      C = 1
       call Interpol(GST0,m13,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! begin astronomical twilight = tw5
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HAS,tw5,Ha,delta,azm)
      tw5 = tw5       

      m23 = m0 + H03 / D2PI                ! End astronomical twilight
      m23 = modulo(m23,1d0)
      m23 = m23
      C = 1
       call Interpol(GST0,m23,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,    ! end astronomical twilight = tw6
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,HAS,tw6,Ha,delta,azm)
      tw6 = tw6 
      
300   end
     
             
C----------------------------------------------------------------------------------------------------------------------------------------


      subroutine TRANSIT(NTARG,TDB0,GST0,AUDAY,DELTAT,PL,SD,LAT,LONG,
     . DX00,DY00,RIS,SET,TR,azm1,azm2,Ha,aS,NR,MS,P,LX,LM,LN)
     

!     Reference: Meeus. " Astronomical formulae for calculators"
!               " 2th edition.   pag.101"            

!     This program compute the times of transit,rise and set of body
!     of the solar sistem , this algorithm takes in account the atmospheric
!     refraction and all the parameters relate.     

!     INPUT:
!           NTARG:    is Solar System Body number.              
!           TDB0 :    Barycentric Dinamical Time to 0.0h UT1.
!           GST0 :    Greenwich Apparent Sideral Time at 0.0h UT1 
!           AUDAY:    Speed of light. (AU per day  "TDB")
!           DELTAT:   DeltaT in fraction day.
!           PL,SD:    SD = Moon's semidiameter  PL = Moon's horiz.parallax. (degs)
!           LAT,LONG: Latitude and Longitude of observer. (Rad.)
!           
!     OUTPUT
!           RIS,SET,TR:  Times of rising,set and transit. (UT)
!           azm1,azm2:   Angle of azimuth .(Deg.)
!           Ha:          Altitude of sea level.
!           labelR:      Label of rise time. 
!           labelS:      Label of setting time. 
!           labelT:      Label of transit time. 
!           labelNR,MS,P,LX,LM,LN:  Visual references.
 
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


      INTEGER NTARG,I ,J,II,C,flag,MS,NR,P,LX,LM,LN
      DOUBLE PRECISION TDB0,GST0,AUDAY,RNPB0(3,3),PL,SD,LAT,LONG
      DOUBLE PRECISION RA3,R2,R3,R6,R8,R9,R10,DE3,D2,D3,D6,D8,D9,D10
      DOUBLE PRECISION trn,T,H00,H0,azm,azm1,azm2,delta,Ha,alt,M01
      DOUBLE PRECISION COH00,DELTAT,h_,m,m1,m2,m_,m_m,newm,oldm,m0
      DOUBLE PRECISION rise,s_,set_,settime,trntime,DX00,DY00
      DOUBLE PRECISION T1,T2,T3,RIS,SET,TRL,TR,MMM ,Z(3)       
      CHARACTER * 45 labelR,labelS,labelT,SIDE,labelM,labelN,labelX
      CHARACTER * 45 labelY
      CHARACTER * 1 aS,labelO
      
      if(NTARG .EQ. 11) then
      H0 = (90d0 + 0.833333333d0) * DD2R
       elseif(NTARG .EQ. 10) then
       H0 = (90d0 + 0.5666666667d0 -PL + SD ) * DD2R ! SD= Moon's Semidiameter "deg",   
       elseif(NTARG .LT. 10) then                   ! PL= Moon's parallax "deg"
       H0 = (90d0 + 0.5666666667d0) * DD2R
      endif
      
      flag = -1

      if(LAT .GT. 1.4 .OR. LAT .LT. -1.4) goto 300

      call COORD(NTARG,TDB0,AUDAY,DX00,DY00,RA3,R2,R3,R6,R8,R9,       
     .  R10,DE3,D2,D3,D6,D8,D9,D10)

      COH00 = (cos(H0) - sin(LAT) * sin(DE3)) / (cos(LAT) * cos(DE3))

      if (COH00 .LT. -1d0) then
       LX = 1
       labelX = '   ---- Body circumpolar ----   '
       if ((H0 * DR2D) .GE. 96d0) then
        labelM = ' .... BRIGHT .... '
        LM = 1
       else
        labelM = ' .... NO SET - body always up the horizon '
        LM = 2
       endif
       flag = 0
      elseif ( COH00 .GT. 1d0) then
       if ((H0 * DR2D) .GE. 96d0) then
        labelN = ' .... DARK .... '
        LN = 1
       else
        labelN = '... NO RISE - body always below the horizon ' 
        LN = 2  
       endif
        labelY =' -------------- '
        flag = 0
      endif
C---------------
      m0 = (RA3 - LONG - GST0) / D2PI
      m0 = modulo (m0,1d0)
!      print*,"m0---",m0
      if(flag .EQ. flag) then
       H00 = acos(COH00)
       H00 = modulo(H00,PI)

       if(COH00 .LT. -1d0 .OR. COH00 .GT. 1d0)then
        m0 = m0
        RIS = 0d0
        SET = 0d0
        goto 100
       endif
       m1 = m0 - H00 / D2PI
       m1 = modulo(m1,1d0)
       m2 = m0 + H00 / D2PI
       m2 = modulo(m2,1d0)

C-------- Rising
       oldm = m1
       C = 1

       call Interpol(GST0,oldm,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,H0,newm,Ha,delta,azm)

       if(azm .LT. 0) then
        azm = azm + PI
       elseif (azm .GT. PI) then
        azm = azm - PI
       else
        azm = azm
       endif
       azm1 = azm * DR2D

       m = newm 
       rise = 24d0 * m
       if(rise .GT. 24d0) then
        rise = rise - 24d0
C -----Event occurs the following day
        RIS = rise
        labelR = '-- Rise - event occurs the next day' 
        NR = 1
       elseif(rise .LT. 0d0) then
        rise = rise + 24d0
C------Event occurs the previus day
        RIS = rise
        labelR = '-- Rise - event occurs the previus day' 
        NR = 2
       else
        RIS = rise
        labelR = ' '
        NR = 0
       endif

C------Setting
      
       oldm = m2
       C = 1

       call Interpol(GST0,oldm,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,H0,newm,Ha,delta,azm)

       if(azm .LT. 0) then
        azm = azm + 2d0 * PI
       elseif (azm .LT. PI .and. azm .GT. 0d0)then
        azm = azm + PI
       else
        azm = azm
       endif

       azm2 = azm * DR2D
       m = newm
 
       settime = 24d0 * m
       
       if(settime .GT. 24d0)then
        settime = settime - 24d0
        SET = settime     
        labelS = '--Set - event occurs the next day'
        MS = 1
       elseif (settime .LT. 0d0)then
        settime = settime + 24d0
        SET = settime
        labelS = '--Set- event occurs the previus day'
        MS = 2 
       else
        SET = settime
        labelS = ' ' 
        MS = 0     
       endif
      endif 

C----- transiting
100   oldm = m0
      C = 0d0

       call Interpol(GST0,oldm,DELTAT,RA3,R2,R3,R6,R8,R9,R10,DE3,
     .     D2,D3,D6,D8,D9,D10,LONG,LAT,C,H0,newm,Ha,delta,azm)

      if(delta .GT. LAT)then
       aS ='N'
      elseif(delta .LT. LAT)then
       aS ='S'
      endif

      Ha = Ha * DR2D
      m = newm

      trntime = 24d0 * m
      if(trntime .GT. 24d0)then
       trntime = trntime - 24d0
       TR = trntime
       labelT = '--Transit- event occurs the next day'
       P = 1
      elseif(trntime .LT. 0d0)then
       trntime = trntime + 24d0
       TR = trntime
       labelT = '--Transit- event occurs the previus day'
       P = 2
      else
       TR = trntime
       labelT = ' '
       P = 0
      endif

      goto 500

300   print*,"****  Algorithm used in program is not optimized   ****"
      print*,"****   for locations very close to the poles.      ****"
      print*," Specificy a latitude between 80 degs N and 80 degs S."  

500   end
C----------------------------------------------------------------------------------------------------------------------------------------


      subroutine JPLEPH_0(FINAL,XD,TIME,ET,TT,UT1,NTARG,longit,latit,
     .HM,TC,GP,GV,XPO,YPO,DELTA,RA4,DE4,SGP,MGP,DIST,DX00,DY00,ELEV,AZ)
    

! Radians to degrees
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )

      INTEGER NTARG,ISTEP,NSTEP,NVS,J
      DOUBLE PRECISION TT,RA4,DE4,R(6),SS(3),VALS(400),PVH(3,2),GAST
      DOUBLE PRECISION SGE(3,2),PVE(3,2),EB(3),GPB(3),DT,DIST,PVB(3,2)
      DOUBLE PRECISION BVE(3),SGP(3),HMOD,GMB(3,2),iau_ANP,RA3,ET
      DOUBLE PRECISION iau_GMST00,DJMJD0,DP00, DE00,RA1,DE1,DE3,UT1
      DOUBLE PRECISION iau_EE00,EPSA, RB(3,3), RP(3,3), RPB(3,3),MO2
      DOUBLE PRECISION BODY_G(3),BODY_H(3),GD,HD,XET,PV(3,2),RA5,DE5
      DOUBLE PRECISION BEV(3,2),BEP(3),HPE(3),u1(3),u2,Mo1,u0(3)
      DOUBLE PRECISION u4(3),THD,EPV(3,2),HE0(3,2),u3(3),RN(3,3)
      DOUBLE PRECISION RNPB0(3,3),CVE(3),BCV,BM1,C,AU,DS,K,AUDAY,VC
      DOUBLE PRECISION RA6,DE6,TIME,BGP(3),iau_S06,DX00,DY00,X,Y,S
      DOUBLE PRECISION RC2I(3,3),TC,WL,PHPA,DET,UJD,ELEV
      DOUBLE PRECISION XPO,YPO,longit,latit,HM,RA7H,DE7D,ZD,AZ,GP(3)
      DOUBLE PRECISION RA9,DE9,RA7,DE7,u5(3),u6(3),OB(3),POS_T(3),XD
      DOUBLE PRECISION POS_H(3),RH,TSL,iau_DTDB,DIST2,u01(3),GV(3)
      DOUBLE PRECISION HOB, DOB, ROB,u03(3),OV(3),p(3),DELTA
      CHARACTER*6 NAMS(400) 

!------- Use subroutine  PLEPH ( DE432) 
      
      ISTEP=1
      NSTEP=1
       
      CALL  Const (NAMS, VALS, SS, NVS)
      IF (XD == 0) then
      END IF     
      CALL PLEPH ( ET, 3, 12, PVE)              ! PVE = Earth barycentric state vector
      
      CALL PLEPH ( ET, 3, 11, PVH)              ! PVH = Earth heliocentric state vector

      CALL PLEPH ( ET, NTARG, 12, PVB)          ! PVB = Body  barycentric state vector

      CALL PLEPH ( ET, NTARG, 11, BEV)          ! BEV = Body's heliocentric state vector 

      CALL PLEPH ( ET, 11, 3, SGE )             ! SGE = Sun geocentric state vector  

      CALL PLEPH ( ET, 10, 3, BGV)              ! BGV = Body geocentric state vector  

!------Update a PVH-vector, discarding the velocity component.

      DT = 0.0D0 
      CALL iau_PVUP (DT,PVE,EB)                  ! EB Barycentric  position Earth vector

      CALL iau_PPP (EB,GP,OB)                    ! OB = EB + GP (barycentric position of observer)
      
      DO  I=1,3
        BVE(I) = PVE(I,1) * DT + PVE(I,2)
      END DO                                     ! BVE Barycentric velocity Earth vector
  
      CALL iau_PVUP (DT,BGV,MGP)                 ! MGP Moon geocentric position vector
 
      CALL iau_PVUP (DT,SGE,SGP)                 ! SGP Sun Geocentric position vector 

      CALL iau_PPP (BVE,GV,OV)                   ! OV = BVE + GV (barycentric velocity of observer)
      
      CALL iau_PVUP (DT,PVH,HPE)                 ! HPE Heliocentric position Earth vector 

      CALL iau_PVUP (DT,BEV,BEP)                 ! BEP Body's Heliocentric position
       
      CALL iau_PM(HPE,HMOD)                      ! HMOD  Earth Heliocentric position modulus

      CALL iau_PVMPV(PVB,PVE,GMB)                ! GMB  Body's geocentric mean state vector
      
      CALL iau_PVUP (DT,GMB,GPB)                 ! GPB  Body's geocentric mean J2000 position    
      
      CALL iau_C2S ( GPB,RA1,DE1 )
              
      RA1 = iau_ANP(RA1)                         ! range angle (0,2!pi)
      
      DIST = DSQRT(GPB(1)*GPB(1)+GPB(2)*GPB(2)+GPB(3)*GPB(3))     ! Geocentric Distance 

!---------- Astrometric coordinate      
      
      CALL L_TIME(ET,NTARG,EB,BVE,HMOD,BODY_G,BODY_H,GD,HD,XET,PV,
     . HE0,THD,EPV) 
      CALL  iau_C2S ( BODY_G, RA3, DE3 )  
      
      RA3 = iau_ANP(RA3)                         ! range angle (0,2!pi)
!      print*, "RA3,DE3",RA3 * DR2D,DE3 *DR2D
  

! =========================================================
!            IAU 2000A, equinox based
! =========================================================
      DJMJD0 = 2400000.5d0
!      TT = JD - DJMJD0       ! TT in MJD
      
!---------- Nutation, IAU 2000A. 
      CALL iau_NUT00A ( DJMJD0, TT, DP00, DE00 )
      
!---------- Precession-nutation quantities, IAU 2000. 
      CALL iau_PN00 ( DJMJD0, TT, DP00, DE00, 
     .      EPSA, RB, RP, RPB, RN, RNPB0 )
!===================================================================================
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
!       ERA  =  iau_ERA00  (   UT11,  UT12    ) 

!====================================================================================
 
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun ( for apparent geocentric coordinate)
      
      CALL DEFLIGHT_0(BODY_G,BEP,HPE,u1)
      
      u1 = u1
      
      CALL iau_PN(u1,u2,p)          ! p = unit vector in order MU/C^2, u2 = modulus

      CALL  iau_C2S   ( p, RA5, DE5 )       
      RA5 = iau_ANP(RA5)            ! range angle (0,2!pi)
      
!---------------------------------------------------------------------------       
!  RELATIVISTIC ANNUAL ABERRATION   ( for apparent geocentric coordinate)
      C  = 299792458D0                     ! Speed of light  m/s 
      AU = 149597870700D0                  ! Astronomical Unit (meter)(TDB)
      DS = 86400D0                         ! n° seconds  Day
      K = 0.01720209895D0                  ! K  heliocentric gravitational costant (UA/Day)
      AUDAY = (C * DS / AU)                ! Speed of light (AU per day  "TDB")
      VC = AU / DS / C                     ! Speed of body in unit of C.

         
      
      CALL iau_SXP ( VC,BVE,CVE)    ! CVE body velocity vector in unit of C 
      
      BCV = SQRT(CVE(1)**2d0+CVE(2)**2d0+CVE(3)**2d0) !form module of vector CVE

      BM1 = SQRT(1d0-abs(BCV)**2d0)            ! for use with routine SOFA iau_AB()  

      CALL ABERR (AUDAY,u1,BVE,u0)

      CALL iau_PN(u0,MO1,u3)         !Form unit vector of u0 = u3 ) 
 
           CALL  iau_C2S ( u3, RA6, DE6 )  
     
      RA6 = iau_ANP(RA6)                         ! range angle (0,2!pi)


!---------- Form apparent geocentric coordinate (RA4,DE4)

       CALL   iau_RXP  ( RNPB0, u3, u4 )
!       "Apparent geocentric position true equinox and ecliptic of the date.(IAU 2000 EQUINOX BASED)" 

      CALL  iau_C2S ( u4, RA4, DE4 )       
      RA4 = iau_ANP(RA4)                         ! range angle (0,2!pi)
      
      CALL L_TIME(ET,NTARG,OB,BVE,HMOD,POS_T,POS_H,DIST2,THD,XET,PV,
     .  HE0,THD,EPV) ! RECALL for Topocentric data
      
!====================================================================================
!     RELATIVISTIC DEFLECTION OF LIGHT due the Sun (for Topocentric coordinate)
      
      CALL DEFLIGHT_0(POS_T,BEP,HPE,u01)

!--------- Add in deflection due to Earth

      CALL DEFLIGHT_T(u01,BEP,HPE,u01)
!      u01(1)=0.49725432787932250D0
!      u01(2)=0.80893654831628226D0
!      u01(3)=0.35069964375626600D0
      
      CALL iau_PN(u01,u02,p)          ! p = unit vector in order MU/C^2, u02 = modulus

!---------------------------------------------------------------------------       
!         RELATIVISTIC ANNUAL ABERRATION  (for Topocentric coordinate)
     
      CALL ABERR (AUDAY,u01,OV,u03)
      
      CALL iau_PN(u03,MO2,u5)         !Form unit vector of u03 = u5 ) 
      
!=========================================================================================

       CALL   iau_RXP  ( RNPB0, u5, u6 )
! ......... Call for direction cosines to spherical coordinates
       CALL   iau_C2S   ( u6, RA7, DE7 )
         RA7= iau_ANP(RA7)                       ! range angle (0,2!pi)
!        print*,"u6,RA7,DE7",u6,RA7,DE7  
         RA7H = RA7 * DR2D / 15D0                ! RA7 in hour
         DE7D = DE7 * DR2D
      
!---------- "Topocentric observed position with refraction.(IAU 2000 EQUINOX BASED)"

      UJD = DJMJD0 + UT1              ! UJD = UT1  Julian day       
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

      RH = 0.5D0                                      ! relative humidity at the observer (range 0-1)
      WL = 0.55D0                                     ! wavelength (micrometers)
      TSL = 20.0D0                                    ! Is the approximate sea-level air temperature in K"
      PHPA = 1013.25D0 * exp ( -HM / ( 29.3D0 * TSL ) )! pressure at the observer (hPa = mB

      DET = DELTA
      CALL REFDAT (PHPA,TSL,0D0,WL)
      CALL SETDT ( DET) 
!      DUT = (UT1+2400000.5d0 - TIME)*86400d0  
      IF(TT > FINAL) THEN
       XPO = 0d0; YPO = 0d0
      END IF

      CALL ZDAZ (UJD, XPO, YPO, longit, latit, HM, RA7H, DE7D, 1, 
     .                    ZD, AZ, RA9, DE9 )
       ELEV = 90d0 - ZD 
       AZ = AZ
!       print*,"ZD,AZ",ZD,AZ

       

      END

*----------------------------------------------------------------------------------------------------

      subroutine JPLEPH2 (JD,RA1,DE1)

      INTEGER NTARG,ISTEP,NSTEP
      DOUBLE PRECISION JD,RA1,DE1,R6,SS3,VALS(400),PVH(3,2),PVB(3,2)
      DOUBLE PRECISION SGE(3,2),LEP(3),EB(3),GPB(3),DT,DIST,LES(3,2)
      DOUBLE PRECISION BVE(3),SGP(3),HPE(3),HMOD(3),GMB(3,2),iau_ANP
      DOUBLE PRECISION PVE(3,2)
      CHARACTER*6 NAMS(400)  

!------- Use subroutine  PLEPH ( DE432) 
    

      ISTEP=1
      NSTEP=1
         
      CALL  Const (NAMS, VALS, SS, NVS)

      CALL PLEPH ( JD, 3, 12, PVE)              ! PVE = Earth barycentric state vector   
 
      CALL PLEPH ( JD, 10, 11, LES)              ! LES = Lunar Ecliptic state vector

      
      CALL iau_PVUP (DT,LES,LEP)                 ! LEP  Lunar Elioc. mean J2000 position 

!      CALL iau_PPP (BVE,GV,OV)                   ! OV = BVE + GV (barycentric velocity of observer)

      CALL  iau_C2S   ( LEP, RA1, DE1 )       
      RA1 = iau_ANP(RA1)                         ! range angle (0,2!pi)
      
!      DIST = DSQRT(GPB(1)*GPB(1)+GPB(2)*GPB(2)+GPB(3)*GPB(3))     ! Geocentric Distance 

      END
!-------------------------------------------------------------------------------------------

      subroutine city_calc(num,longit,latit,height,zone,site)

      implicit none

      DOUBLE PRECISION  longit,latit,lo,la        
      INTEGER :: num,I,N,height,he ,ios 
      CHARACTER (len = 19) zone,zo                     
      CHARACTER (len = 70) si,site
C     Open database  "city_coord.txt" 
      open (unit=11, file="city_coord.txt", status="old",         
     .  action="read",form="formatted",position="rewind")
      
      DO I = 1,430 

      read (unit=11,fmt="(i4,f11.4,f11.4,i7,2x,a19,3x,a71)"  
     . ,iostat=ios)N,lo,la,he,zo,si
      
       IF(I == num) THEN 
         longit = lo
         latit  = la
         height = he
         zone   = zo
         site   = si
       END IF
      END DO

      rewind (unit=11)

      end

!----------------------------------------------------------------------------

