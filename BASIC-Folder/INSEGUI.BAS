' Inseguimento di un SATELLITE ARTIFICIALE - Insegui.bas
' =======================================  (Gen.05)
'
  COLOR 14, 1: CLS
  DEFDBL A-Z
  pi = 4 * ATN(1): rad = pi / 180
DEF FNatan (y, x) = (1 + (x = 0)) * ATN(y / (x + (x = 0))) - pi * (x <= 0) - 2 * pi * (x > 0 AND y < 0) + (x = 0) * SGN(y) * pi / 2
DEF FNridu (x) = (x / (2 * pi) - (INT(x / (2 * pi)) - (x < 0))) * 2 * pi
DEF FNasin (x) = (2 + (x < 1 AND x > -1)) * ATN(x / SQR(1 + x * x * (x < 1 AND x > -1)))
'----------------------------------
' Costanti astronomiche
  ae = 6378.14: 'Raggio Equatoriale in Km
  ec = .081819489#
  mu = 398600.5: ' Param.gravit.Terra in km^3/sec^2
'------------------------------------------------------------------------------------------------------------------------------------
COLOR 14, 4
  LOCATE 2, 17: PRINT "ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿"
  LOCATE 3, 17: PRINT "³                                             ³"
  LOCATE 4, 17: PRINT "ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ"
COLOR 14, 5
  LOCATE 3, 18: PRINT "   Moto dei SATELLITI Artificiali TERRESTRI  "
COLOR 14, 1
  LOCATE 6, 17: PRINT "ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿"
  LOCATE 7, 17: PRINT "³                                             ³"
  LOCATE 8, 17: PRINT "³                                             ³"
  LOCATE 9, 17: PRINT "³                                             ³"
 LOCATE 10, 17: PRINT "³                                             ³"
 LOCATE 11, 17: PRINT "ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ"

 LOCATE 8, 20: COLOR 13: PRINT "  (Formato NASA - 3 righe)": COLOR 11
 LOCATE 9, 20: PRINT "URL: http://datameteo.com/sat/track.htm": COLOR 12
LOCATE 10, 20: PRINT "Coordinate Stazione nel file: coo.dat": COLOR 14
 LOCATE 7, 20: INPUT ; "Nome DataBase dei Satelliti: ", z$
CLS
'
' -------- Inizio LETTURA DATI da file esterno (dummy.dat) -----------
OPEN z$ FOR INPUT AS #1
1946 :
' Legge la 1^ riga di "dummy.dat"
INPUT #1, nomesat$: ' Nome del satellite
ON ERROR GOTO 2000
' Legge la 2^ riga di "dummy.dat"
INPUT #1, a$
 ann0 = VAL(MID$(a$, 19, 2)) + 2000: annoepo$ = STR$(ann0) + "-"
' Controlla se l'anno è bisestile e stabilisce il numero dei giorni all'anno
 IF (ann0 / 4 - INT(ann0 / 4)) <> 0 THEN
 feblim = 59
 ELSE
 feblim = 60
 END IF
'------------------------------------
 gio0 = VAL(MID$(a$, 21, 12))
 giornoepo$ = STR$(gio0)
'
' Calcolo del mese
  IF INT(gio0) <= 31 THEN mes0 = 1
  IF INT(gio0) > 31 AND INT(gio0) <= feblim THEN mes0 = 2
  IF INT(gio0) > feblim AND INT(gio0) <= feblim + 31 THEN mes0 = 3
  IF INT(gio0) > feblim + 31 AND INT(gio0) <= feblim + 61 THEN mes0 = 4
  IF INT(gio0) > feblim + 61 AND INT(gio0) <= feblim + 92 THEN mes0 = 5
  IF INT(gio0) > feblim + 92 AND INT(gio0) <= feblim + 122 THEN mes0 = 6
  IF INT(gio0) > feblim + 122 AND INT(gio0) <= feblim + 153 THEN mes0 = 7
  IF INT(gio0) > feblim + 153 AND INT(gio0) <= feblim + 184 THEN mes0 = 8
  IF INT(gio0) > feblim + 184 AND INT(gio0) <= feblim + 214 THEN mes0 = 9
  IF INT(gio0) > feblim + 214 AND INT(gio0) <= feblim + 245 THEN mes0 = 10
  IF INT(gio0) > feblim + 245 AND INT(gio0) <= feblim + 275 THEN mes0 = 11
  IF INT(gio0) > feblim + 275 AND INT(gio0) <= feblim + 306 THEN mes0 = 12

  meseepo$ = STR$(mes0) + "-"
' JD0 dell'Epoca
      u7 = 2 - INT(ann0 / 100) + INT(INT(ann0 / 100) / 4)
      IF mes0 < 3 THEN
          ann0 = ann0 - 1: mes0 = mes0 + 12
          END IF
      jd0 = INT(365.25 * ann0) + INT(30.6001 * (mes0 + 1)) + gio0 + 1720994.5#
      jd0 = jd0 + u7: ' tempo universale TU
'
'ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ

' Legge la 3^ riga di "dummy.dat"
INPUT #1, a$
' Elementi orbitali
 nqq = VAL(MID$(a$, 53, 17)): 'rivoluzioni/giorno
 eqq = VAL(MID$(a$, 27, 7)) / 10000000#: ' eccentricità
 iqq = VAL(MID$(a$, 9, 8)): ' Inclinazione orbitale (sull'equatore terrestre) [gradi]
 Oqq = VAL(MID$(a$, 18, 8)): ' Long.Nodo Asc. [gradi]
 wqq = VAL(MID$(a$, 35, 8)): ' Argom. del Perigeo [gradi]
 M0qq = VAL(MID$(a$, 44, 8)): ' M0 = Amomalia media all'Epoca
'''CLOSE #1
'ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
' Calcolo del semiasse maggiore (a) dell'orbita e del parametro p
' ---------------------------------------------------------------
' Periodo orbitale in [sec]
  Tsec = 86400 / nqq
  aux2 = (Tsec / (2 * pi)) ^ 2 * mu
  aqq = aux2 ^ (1 / 3): 'Semiasse maggiore
  pqq = aqq * (1 - eqq ^ 2): ' Parametro dell'ellisse: p=a*(1-e^2)
'--------------------------------------
OPEN "coo.dat" FOR INPUT AS #3
 INPUT #3, qq$
 lonqq = VAL(qq$)
 INPUT #3, qq$
 latqq = VAL(qq$)
 INPUT #3, qq$
 Hqq = VAL(qq$)
 INPUT #3, qq$
 fusoqq = VAL(qq$)
CLOSE #3
'--------------------------------------
'ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
COLOR 10
LOCATE 2, 8: PRINT "Premere: (c) per cambiare satellite;   (u) per uscire"
COLOR 15
LOCATE 3, 20: PRINT "Nome Sat -> "; : COLOR 12: PRINT nomesat$
COLOR 10
LOCATE 18, 20: PRINT "Elementi Orbitali del Satellite"
COLOR 14
LOCATE 19, 6: PRINT "Epoca:" + annoepo$ + meseepo$ + giornoepo$;
PRINT USING "    M= ###.####ø    n= ##.########ø/g"; M0qq; nqq
LOCATE 20, 6: PRINT USING "e= #.#######   i= ###.####ø    O= ###.####ø    w= ###.####ø"; eqq; iqq; Oqq; wqq
LOCATE 21, 6: PRINT USING "a= ##### km"; aqq
COLOR 15
LOCATE 22, 6: PRINT "Coord. geografiche Stazione:";
COLOR 11: PRINT USING " Lon(est)=+###.####ø  Lat= +##.####ø"; lonqq; latqq
COLOR 14
'ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
' -------- Inizio CICLO -----------
100
DO
COLOR 11
LOCATE 4, 16: PRINT MID$(DATE$, 4, 2) + "-" + MID$(DATE$, 1, 2) + "-" + MID$(DATE$, 7, 4)
'COLOR 11
LOCATE 4, 42: PRINT MID$(TIME$, 1, 2) + ":" + MID$(TIME$, 4, 2) + ":" + MID$(TIME$, 7, 2); : COLOR 14: PRINT " TMEC"

' Riassegnazione elementi orbitali e coo.sito osserv.
  n = nqq: e = eqq: i = iqq: a = aqq: p = pqq
  O = Oqq: w = wqq: M0 = M0qq
  lon = lonqq: lat = latqq: H = Hqq
  fuso = fusoqq
'--------------------------------

' Data dell'effemeride
' ---------------------
oreqq = VAL(MID$(TIME$, 1, 2)) - fuso
minqq = VAL(MID$(TIME$, 4, 2))
secqq = INT(VAL(MID$(TIME$, 7, 2)))
'
gioqq = VAL(MID$(DATE$, 4, 2))
mesqq = VAL(MID$(DATE$, 1, 2))
annqq = VAL(MID$(DATE$, 7, 4))

' Riassegnazione data e tempo
 ore = oreqq: min = minqq: sec = secqq
 gio = gioqq: mes = mesqq: ann = annqq
' ---------------------------------

' ---- Calcoli ----

 'Calcolo di JD
      u7 = 2 - INT(ann / 100) + INT(INT(ann / 100) / 4)
      IF mes < 3 THEN
          ann = ann - 1: mes = mes + 12
          END IF
      jd = INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5#
      IF jd >= 2299160.5# THEN jd = jd + u7

 'Calcolo del T-sid a Greenwich all'ora fissata -
      wn1 = ore + min / 60 + sec / 3600: ' tempo universale TU
     mat1 = .279057273# + 100.002139027# * t + 1.0776E-06 * t ^ 2
     mat2 = mat1 - INT(mat1)
   tsidd = 1.00273790935# * wn1 + mat2 * 24
   tsidd = tsidd - INT(tsidd / 24) * 24: 'in ORE
    ts = tsidd * 15: 'in GRADI
    jd1 = jd + wn1 / 24

' Tempo Delta-t tra le due date: Epoca e Effemeride
  dt = jd1 - jd0: ' in giorni

 'Tempo Siderale Locale (LST) = (Tsid-Green. + Long_est) [gradi]
  LST = ts + lon:
  LST = LST - INT(LST / 360) * 360: ' normalizzazione angolo

' Latitudine e LST in radianti
   Fi = lat * rad: te = LST * rad

' Vettore posizione inerziale (R) -> Centro Terra-Stazione osservativa
  aux1 = SQR(1 - (ec * SIN(Fi)) ^ 2)
  xe = ae / aux1 + H / 1000
  ze = (1 - ec ^ 2) * ae / aux1 + H / 1000

  Ri = xe * COS(Fi) * COS(te)
  Rj = xe * COS(Fi) * SIN(te)
  Rk = ze * SIN(Fi)
 Rg = SQR(Ri ^ 2 + Rj ^ 2 + Rk ^ 2)

' Anomalia Media: M = M0+ngùdt   (ng= moto medio in gradi/giorno)
 ng = 360 * n
 m = M0 + ng * dt
 m = m - INT(m / 360) * 360: ' normalizzazione anomalia media [gradi]
 mgra = m

' Equazione di Keplero
'----------------------
 m = m * rad
anomecc = m    'Valore iniziale di E assunto pari all'anom.media M
 DO
 variaz = (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc))
  anomecc = anomecc + variaz
LOOP UNTIL ABS(variaz) < .0000000001#


 'Anomalia VERA e RAGGIO geocentrico (r)
 '--------------------------------------
 tgv = TAN(anomecc / 2) * SQR((1 + e) / (1 - e))
 v = FNridu(2 * FNatan(tgv, 1))
 r = a * (1 - e * COS(anomecc))

' Versori geocentrici
 ugra = w + v / rad: ' argomento della latitudine (u = w + v) in [gradi]
 ugra = ugra - INT(ugra / 360) * 360: ' normalizzazione (u) in [gradi]

 u = ugra * rad: ' (u) in [RAD]
 nodo = O * rad: ' (O) in [RAD]
 i = i * rad: ' (i) in [RAD]

'...........................................................
 Rx = COS(u) * COS(nodo) - SIN(u) * SIN(nodo) * COS(i)
 Ry = COS(u) * SIN(nodo) + SIN(u) * COS(nodo) * COS(i)
 Rz = SIN(u) * SIN(i)

 xpri = -SIN(u) * COS(nodo) - COS(u) * SIN(nodo) * COS(i)
 ypri = -SIN(u) * SIN(nodo) + COS(u) * COS(nodo) * COS(i)
 zpri = COS(u) * SIN(i)
'...........................................................


' VETTORE POSIZIONE - geocentrico: r=r(x,y,z)
' -----------------
  x = r * Rx: y = r * Ry: z = r * Rz

' VETTORE VELOCITA' - geocentrico: V=V(x',y',z')
' -----------------
 f = v: 'anomalia vera ricopiata su f, in (rad)

 Vx = SQR(mu / p) * (Rx * e * SIN(f) + (1 + e * COS(f)) * xpri)
 Vy = SQR(mu / p) * (Ry * e * SIN(f) + (1 + e * COS(f)) * ypri)
 Vz = SQR(mu / p) * (Rz * e * SIN(f) + (1 + e * COS(f)) * zpri)
 vel = SQR(Vx ^ 2 + Vy ^ 2 + Vz ^ 2)

GOTO 89
  PRINT
  PRINT
  PRINT USING "   Componenti Vett. V [Km/s]:  +#####.###    +#####.###    +#####.###"; Vx; Vy; Vz
  PRINT USING "                         V->    #####.###"; vel
'  PRINT
  PRINT "   -------------------------------------------------------------------"
  PRINT USING "   Componenti Vett.   r [Km]:  +#####.###    +#####.###    +#####.###"; x; y; z
  PRINT USING "                         r->    #####.###"; r
  PRINT
  PRINT USING "   Componenti Vett.   R [Km]:   +####.###     +####.###     +####.###"; Ri; Rj; Rk
  PRINT USING "                         R->    #####.###"; Rg
  PRINT "   -------------------------------------------------------------------"
89 :

' VETTORE TOPOCENTRICO: rho=rho(x-Ri,y-Rj,z-Rk) = rho(rhoi,rhoj,rhok)
' --------------------
  rhoi = x - Ri
  rhoj = y - Rj
  rhok = z - Rk
 rho = SQR(rhoi ^ 2 + rhoj ^ 2 + rhok ^ 2)

GOTO 99
  PRINT USING "   Componenti Vett. rho [Km]:  +#####.###    +#####.###    +#####.###"; rhoi; rhoj; rhok
  PRINT USING "     (rho= r-R)        rho->    #####.###"; rho
  PRINT
99 :


' AR e DECL. Topocentriche
' ------------------------
  AR = FNatan(rhoj, rhoi)
  DE = ATN(rhok / SQR(rhoi ^ 2 + rhoj ^ 2))
GOTO 59
  delgr = AR / (rad * 15)
  GOSUB GradAsess
  PRINT USING "         AR ---> ###h ##m ##.##s"; deltagra; deltamin; deltasec

  delgr = DE / rad
  GOSUB GradAsess
  PRINT USING "         DE ---> +##ø ##' ##.##''"; deltagra; deltamin; deltasec
  PRINT
59 :

' Azimut e Altezza sull'orizzonte (frame SEZ = Sud-Est-Zenit)
' -------------------------------
  rhos = SIN(Fi) * COS(te) * rhoi + SIN(Fi) * SIN(te) * rhoj - COS(Fi) * rhok
  rhoe = -SIN(te) * rhoi + COS(te) * rhoj
  rhoz = COS(Fi) * COS(te) * rhoi + COS(Fi) * SIN(te) * rhoj + SIN(Fi) * rhok


  El = ATN(rhoz / SQR(rhos ^ 2 + rhoe ^ 2))
  aux3 = ATN(ABS(rhoe) / ABS(rhos))

  IF rhos < 0 AND rhoe > 0 THEN Az = aux3
  IF rhos > 0 AND rhoe > 0 THEN Az = pi - aux3
  IF rhos > 0 AND rhoe < 0 THEN Az = pi + aux3
  IF rhos < 0 AND rhoe < 0 THEN Az = 2 * pi - aux3
  Azgra = Az / rad: ' Azimut in gradi
  Elgra = El / rad: ' Altezza in gradi
GOTO 39
  delgr = Az / rad
  PRINT USING "         Azimut ---> ###.####ø"; delgr

  delgr = El / rad
  PRINT USING "         Alt. h ---> +##.####ø"; delgr
  PRINT
39 :

'
' Versori di X = (Rx,Ry,Rz), uguali a quelli del vettore posiz.(r)

'-------------- Procedura (Calcolo ATAN2) ------------------
  aux4 = ABS(ATN(Ry / Rx))

  IF Ry > 0 AND Rx > 0 THEN tet = aux4:        ' 1ø quadr.
  IF Ry > 0 AND Rx < 0 THEN tet = pi - aux4:   ' 2ø quadr.
  IF Ry < 0 AND Rx < 0 THEN tet = pi + aux4:   ' 3ø quadr.
  IF Ry < 0 AND Rx > 0 THEN tet = 2 * pi - aux4: ' 4ø quadr.

  Tetagra = tet / rad
  aux5 = Tetagra - ts: ' in gradi

  LonX = aux5
  IF aux5 < -360 THEN LonX = 360 + aux5
  IF aux5 > 360 THEN LonX = 720 - aux5
'-------------------------------------------
  IF aux5 < -180 THEN LonX = 360 + aux5
  IF aux5 > 180 THEN LonX = aux5 - 360

  LatX = FNasin(Rz) / rad

 ' Calcolo di (Long, Lat) della localit… dove il satellite
 ' si trova allo ZENIT all'istante (t) considerato e la sua
 ' Quota (q) in quel punto
 ' ----------------------------------------------------------
' Quota del Sat. in X
  Rt = ae * SQR((1 - ec ^ 2)) / (1 - ec ^ 2 * SIN(LatX * rad) ^ 2): 'Raggio sfera locale in X
  q = r - Rt
'
' Flag di visibilit…
 IF Elgra < 0 THEN
 wo$ = "NON Visib."
 ELSE
 wo$ = " VISIBILE "
 END IF

  COLOR 13
        PRINT "                                                  Posiz.Geogr.Satellite": COLOR 14
  PRINT USING "       Azimut (N-E):   ###.####ø          Longitudine (est)=  +###.####ø"; Azgra; LonX
  PRINT USING "     Elevazione (h):   +##.####ø"; Elgra; : COLOR 11
  PRINT wo$; : COLOR 14
  PRINT USING "Latitudine       =   +##.####ø"; LatX
  PRINT USING "    Distanza Topoc.: #####.# km           Quota            = #####.# km"; rho; q
        PRINT "   -----------------------------------------------------------------------"
  PRINT USING "   Componenti Vettore V [Km/s]:    +##.###       +##.###       +##.###"; Vx; Vy; Vz
  PRINT USING "      (Veloc. geocentr.)  V->       ##.###"; vel
  PRINT USING "   Componenti Vettore r [Km]:   +#####.###    +#####.###    +#####.###"; x; y; z
  PRINT USING "     (Geocentrico)       r->     #####.###"; r
        PRINT "   -----------------------------------------------------------------------"
  PRINT USING "   Componenti Vettore rho[Km]:  +#####.###    +#####.###    +#####.###"; rhoi; rhoj; rhok
  PRINT USING "     (Topocentrico)     rho->    #####.###"; rho
        PRINT "   -----------------------------------------------------------------------"

b$ = INKEY$
LOOP WHILE b$ = ""
IF UCASE$(b$) = "U" GOTO 2000
IF UCASE$(b$) = "C" GOTO 1946 ELSE 100
'
2000 END
'------------- SUBROUTINES ------------------------------------------------
GradAsess:
' Conversione di GRADI in GRADI,PRIMI,SEC e GRADI,PRIMI,Fraz.PRIMI
'
 deltagra = INT(ABS(delgr))
 deltafraz = (ABS(delgr) - deltagra) * 60
 deltamin = INT(deltafraz)
 deltasec = (deltafraz - deltamin) * 60
 deltagra = deltagra * SGN(delgr)
'
 RETURN


' ******** Fine Programma: Insegui.BAS ************ (Gen.05) ****************

