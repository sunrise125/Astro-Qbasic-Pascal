' Effemeridi di un SATELLITE ARTIFICIALE - Art-Sats.bas
' =======================================  (Gen.05)
'
' Provvisorio ---> Dati inseriti nel LISTATO (non da lettura files esterni)
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
'
' -------- Inizio Inserimento DATI -----------
' nota per dopo..... val(left$(time$,2))-1
'                    val(mid$(date$,4,2))
' Data dell'effemeride
  gio = 11: mes = 1: ann = 2005
  ore = 14: min = 0: sec = 0
  ' ---
' Coordinate geografiche della Stazione di osservazione
  lon = 30.3274: 'Long.est in gradi
  lat = 59.7721: 'Lat.nord in gradi
  H = 65: ' Metri sul livello del mare
' ---
' Elementi orbitali del satellite (LAGEOS 1) -> Geodetic.txt da NASA Dbase
  n = 6.38664605#: 'rivoluzioni/giorno
  e = .0044295#: ' eccentricità
  i = 109.8532: ' Inclinazione orbitale (sull'equatore terrestre) [gradi]
  O = 17.7748: ' Long.Nodo Asc. [gradi]
  w = 179.483: ' Argom. del Perigeo [gradi]
  M0 = 180.5876: ' M0 = Amomalia media all'Epoca
' ------ Epoca: 2005 01 09.92370312 -> seguono i dati dell'epoca
  gio0 = 9: mes0 = 1: ann0 = 2005
  ore0 = 22: min0 = 10: sec0 = 8
' --------   Fine Inserimento DATI -----------
'
' ---- Calcoli ----

 'Calcolo di JD
      u7 = 2 - INT(ann / 100) + INT(INT(ann / 100) / 4)
      IF mes < 3 THEN
          ann = ann - 1: mes = mes + 12
          END IF
      jd = INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5#
      IF jd >= 2299160.5# THEN jd = jd + u7
      jq = jd - 2451545: t = jq / 36525

 'Calcolo del T-sid a Greenwich all'ora fissata -
      wn1 = ore + min / 60 + sec / 3600: ' tempo universale TU
     mat1 = .279057273# + 100.002139027# * t + 1.0776E-06 * t ^ 2
     mat2 = mat1 - INT(mat1)
   tsidd = 1.00273790935# * wn1 + mat2 * 24
   tsidd = tsidd - INT(tsidd / 24) * 24: 'in ORE
    ts = tsidd * 15: 'in GRADI
    jd1 = jd + wn1 / 24

 'Tempo Siderale Locale (LST) = (Tsid-Green. + Long_est) [gradi]
  LST = ts + lon:
  LST = LST - INT(LST / 360) * 360: ' normalizzazione angolo

'  PRINT
'  PRINT
'  PRINT USING "      Tempo Siderale Locale=  ###.#####ø"; LST

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

' Calcolo del semiasse maggiore (a) dell'orbita
' ---------------------------------------------
' Periodo orbitale in [sec]
  Tsec = 86400 / n
  aux2 = (Tsec / (2 * pi)) ^ 2 * mu
  a = aux2 ^ (1 / 3): 'Semiasse maggiore

' JD0 dell'Epoca
 'Calcolo di JD0
      u7 = 2 - INT(ann0 / 100) + INT(INT(ann0 / 100) / 4)
      IF mes0 < 3 THEN
          ann0 = ann0 - 1: mes0 = mes0 + 12
          END IF
      jd0 = INT(365.25 * ann0) + INT(30.6001 * (mes0 + 1)) + gio0 + 1720994.5#
      IF jd0 >= 2299160.5# THEN jd0 = jd0 + u7
      wn1 = ore0 + min0 / 60 + sec0 / 3600: ' tempo universale TU
    jd0 = jd0 + wn1 / 24

' Tempo Delta-t tra le due date: Epoca e Effemeride
  dt = jd1 - jd0: ' in giorni

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
 p = a * (1 - e ^ 2): ' Parametro dell'ellisse p

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
 v = SQR(Vx ^ 2 + Vy ^ 2 + Vz ^ 2)

  PRINT
  PRINT
  PRINT USING "   Componenti Vett. V [Km/s]:  +#####.###    +#####.###    +#####.###"; Vx; Vy; Vz
  PRINT USING "                         V->    #####.###"; v
'  PRINT
  PRINT "   -------------------------------------------------------------------"
  PRINT USING "   Componenti Vett.   r [Km]:  +#####.###    +#####.###    +#####.###"; x; y; z
  PRINT USING "                         r->    #####.###"; r
  PRINT
  PRINT USING "   Componenti Vett.   R [Km]:   +####.###     +####.###     +####.###"; Ri; Rj; Rk
  PRINT USING "                         R->    #####.###"; Rg
  PRINT "   -------------------------------------------------------------------"


' VETTORE TOPOCENTRICO: rho=rho(x-Ri,y-Rj,z-Rk) = rho(rhoi,rhoj,rhok)
' --------------------
  rhoi = x - Ri
  rhoj = y - Rj
  rhok = z - Rk
 rho = SQR(rhoi ^ 2 + rhoj ^ 2 + rhok ^ 2)

  PRINT USING "   Componenti Vett. rho [Km]:  +#####.###    +#####.###    +#####.###"; rhoi; rhoj; rhok
  PRINT USING "     (rho= r-R)        rho->    #####.###"; rho
  PRINT

' AR e DECL. Topocentriche
' ------------------------
  AR = FNatan(rhoj, rhoi)
  DE = ATN(rhok / SQR(rhoi ^ 2 + rhoj ^ 2))

  delgr = AR / (rad * 15)
  GOSUB GradAsess
  PRINT USING "         AR ---> ###h ##m ##.##s"; deltagra; deltamin; deltasec

  delgr = DE / rad
  GOSUB GradAsess
  PRINT USING "         DE ---> +##ø ##' ##.##''"; deltagra; deltamin; deltasec
  PRINT

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

  delgr = Az / rad
  PRINT USING "         Azimut ---> ###.#ø"; delgr

  delgr = El / rad
  PRINT USING "         Alt. h ---> +##.#ø"; delgr
  PRINT

 ' Calcolo di (Long, Lat) della localit… dove il satellite
 ' si trova allo ZENIT all'istante (t) considerato e la sua
 ' Quota (q) in quel punto
 ' ----------------------------------------------------------
' Quota del Sat. in X
  Rt = ae * SQR((1 - ec ^ 2)) / (1 - ec ^ 2 * SIN(LatX * rad) ^ 2): 'Raggio sfera locale in X
  q = r - Rt
'
' Versori di X = (Rx,Ry,Rz), uguali a quelli del vettore posiz.(r)

'-------------- Procedura Toni (Calcolo ATAN2) ------------------
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

'----------------------------------------------------------------
'  IF Rx < 0 THEN tet = -tet
'  IF ts > 180 THEN ts = -360 + ts
' LonX = tet / rad - ts
' IF LonX > 360 THEN LonX = LonX - 360


  LatX = FNasin(Rz)

  PRINT USING "                    Long.X (est) ---> +###.##ø"; LonX
  PRINT USING "                    Lati.X       --->  +##.##ø"; LatX / rad
  PRINT USING "                    Quota in X   --->#####.# km"; q

 END
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


' ******** Fine Programma: ART-SATS.BAS ************ (Gen.05) ****************

