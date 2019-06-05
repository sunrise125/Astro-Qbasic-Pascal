REM ---------  EFFEMERIDI Sole, Luna e Pianeti  ------------
'               (in tempo REALE per Canicattini Bagni)
'
     SCREEN 9: PALETTE 0, 1: COLOR 14
   DEFDBL A-H, J, L-Z
                  pi# = 4 * ATN(1#): rad# = pi# / 180:
    lim = 10: ps = 2
    WINDOW (0, lim)-(lim, 0):
    VIEW SCREEN (15, 15)-(625, 333), , 12: 'Viewport
'
DIM x(125) AS DOUBLE, k(125) AS DOUBLE, c(20) AS DOUBLE, a(7, 3) AS DOUBLE
DIM z(7, 6, 3) AS DOUBLE
DIM l(7) AS DOUBLE, amagg(7) AS DOUBLE, ecc(7) AS DOUBLE
DIM incl(7) AS DOUBLE, periel(7) AS DOUBLE, lnodo(7) AS DOUBLE
DIM manom(7) AS DOUBLE, ekepl(7) AS DOUBLE, vanom(7) AS DOUBLE
DIM arglat(7) AS DOUBLE, lonelio(7) AS DOUBLE, latelio(7) AS DOUBLE
DIM rvett(7) AS DOUBLE, longeo(7) AS DOUBLE, latgeo(7) AS DOUBLE
DIM dpian(7) AS DOUBLE, alfa(7) AS DOUBLE, delta(7) AS DOUBLE
DIM altez(7) AS DOUBLE, azim(7) AS DOUBLE

'
   ' Lettura dei Parametri
   FOR i = 1 TO 20: READ c(i): NEXT
   FOR k = 1 TO 7: FOR i = 0 TO 3: READ a(k, i): NEXT: NEXT
   FOR k = 1 TO 125: READ k(k): NEXT
   FOR k = 1 TO 7: FOR i = 1 TO 6: FOR i1 = 0 TO 3: READ z(k, i, i1): NEXT: NEXT: NEXT
    '
     LINE (1, 6.25)-(9, 9), 11, BF
      LOCATE 4, 21: PRINT "GRAFICO EFFEMERIDI  -   Sole,Luna,Pianeti"
      LOCATE 6, 21: PRINT " (In tempo REALE per CANICATTINI Bagni ) "
      LOCATE 8, 21: PRINT "             Premere un tasto            "
      WHILE INKEY$ = "": WEND
      '
 ab$ = RIGHT$(DATE$, 4): aa = VAL(ab$): m$ = LEFT$(DATE$, 2): m = VAL(m$): g$ = MID$(DATE$, 4, 2): g = VAL(g$)
      '
      gdat = g: mdat = m: adat = aa
      a8 = aa: u7 = 2 - INT(aa / 100) + INT(INT(aa / 100) / 4)
      IF m >= 3 THEN 1100
      a8 = a8 - 1: m = m + 12
1100  jd = INT(365.25 * a8) + INT(30.6001 * (m + 1)) + g + 1720994.5#
      IF jd >= 2299160.5# THEN jd = jd + u7
      n0 = jd - 2415020#
      '
      je = n0: u = 360
      '           '==================================
                        lt = 37: lg = 15: fuso = 1:   'lg = Longit. Est
      '           '==================================
      '
 m$ = LEFT$(TIME$, 2): ha = VAL(m$): g$ = MID$(TIME$, 4, 2): mn = VAL(g$)
      '
      CLS
1200  jq = je + (ha - fuso) / 24 + mn / 1440:     t = jq / 36525
      ore1 = ha: min1 = mn
      gio1 = ore1 / 24 + min1 / 1440
      ore2 = INT(gio1 * 24): min2 = (gio1 * 24 - ore2) * 60
      IF INT(min2 + .5) = 60 THEN ore2 = ore2 + 1: min2 = 0
      IF ore2 >= 24 THEN ore2 = ore2 - 24
      '
    'Calcolo parametri SOLE

    'ja=Long.media Sole(L) ; jb=Anom.media Sole(M) : jc=Eccentr.Orbita Terra
 ja = (c(3) * t + c(2)) * t + c(1): xxx = ja: GOSUB mo: ja = xxx
 jb = ((c(7) * t + c(6)) * t + c(5)) * t + c(4): xxx = jb: GOSUB mo: jb = xxx
 jc = (c(10) * t + c(9)) * t + c(8)
 '
    'sd = Equaz.Centro Sole
 sd = ((c(13) * t + c(12)) * t + c(11)) * SIN(jb * rad#)
 sd = sd + (c(15) * t + c(14)) * SIN(2 * jb * rad#)
 sd = sd + c(16) * SIN(3 * jb * rad#)
 '
    'sf=Long.Vera Sole (teta) ; sg=Anom.Vera Sole (v)
 sf = ja + sd: xxx = sf: GOSUB mo: sf = xxx
 sg = jb + sd
 sh = (1 - jc ^ 2) / (1 + jc * COS(sg * rad#)): ' Raggio Vettore Sole (U.A.)
 si = ((c(20) * t + c(19)) * t + c(18)) * t + c(17): ' Obliquit… Eclitt.(î)
     sl = SIN(si * rad#)
     sm = COS(si * rad#)
 so = 1 / rad# * ATN(SIN(sf * rad#) * sm / COS(sf * rad#))
     IF COS(sf * rad#) < 0 THEN so = so + u / 2
     IF so < 0 THEN so = so + u
 s9 = sl * SIN(sf * rad#): c9 = SQR(1 - s9 ^ 2):   'so =   AR Sole(à-sole)
 sq = 1 / rad# * ATN(s9 / c9):                     'sq = Decl.Sole(ë-sole)
 '
 wn = ha * 15 + mn / 4:  'Trasf. ore+min. in GRADI e fraz. di grado
 '
   'Caricamenti dati orbitali LUNA
       FOR k = 1 TO 7: y(k) = ((a(k, 3) * t + a(k, 2)) * t + a(k, 1)) * t + a(k, 0)
       xxx = y(k): GOSUB mo: y(k) = xxx
       NEXT
 '
 s1 = 51.2 + 20.2 * t: s2 = 346.56 + 132.87 * t - .0091731 * t * t
 s3 = y(6) + 275.05 - 2.3 * t: t1 = SIN(s1 * rad#)
 t2 = SIN(s2 * rad#) * .003964: t3 = SIN(s3 * rad#): t4 = SIN(y(6) * rad#)
 '
    y(1) = y(1) + .000233 * t1 + t2 + .001964 * t4
    y(2) = y(2) - .001778 * t1
    y(3) = y(3) + .000817 * t1 + t2 + .002541 * t4
    y(4) = y(4) + .002011 * t1 + t2 + .001964 * t4
    y(5) = y(5) + t2 - .024691 * t4 - .004328 * t3
 '
       l1 = y(1): m = y(2): m1 = y(3): d = y(4): f = y(5)
 '
    d2 = 2 * d: m2 = 2 * m: n2 = 2 * m1: f2 = 2 * f
    d3 = 3 * d: m3 = 3 * m: n3 = 3 * m1: f3 = 3 * f
    d4 = 4 * d: m4 = 4 * m: n4 = 4 * m1: f4 = 4 * f
 '-------------------------------------------------------------------
x(1) = m1:   x(2) = d2 - m1:  x(3) = d2:     x(4) = n2:     x(5) = f2
x(6) = d2 - n2: x(7) = d2 + m1: x(8) = d:    x(9) = d2 - f2: x(10) = f2 + m1
x(11) = f2 - m1: x(12) = d4 - m1: x(13) = n3: x(14) = d4 - n2: x(15) = m1 - d
x(16) = n2 + d2: x(17) = d4:  x(18) = d2 - n3: x(19) = m1 - f2 - d2: x(20) = m1 + d
x(21) = m1 + d2 - f2: x(22) = f2 + d2: x(23) = n2 + f2: x(24) = m1 - d3: x(25) = m1 + d4
x(26) = n4:  x(27) = n2 - d:  x(28) = m:     x(29) = d2 - m - m1: x(30) = d2 - m
x(31) = m1 - m: x(32) = m + m1: x(33) = m - m1 + d2: x(34) = d2 + m: x(35) = m + d
x(36) = m1 - m + d2: x(37) = n2 - m: x(38) = d2 - m - n2: x(39) = n2 + m: x(40) = d4 - m - m1
x(41) = m + m1 + d2: x(42) = d4 - m - n2: x(43) = m - n2 + d2: x(44) = d2 - m - f2: x(45) = d4 - m
x(46) = d2 - m2: x(47) = m2: x(48) = d2 - m1 - m2: x(49) = m1 - m2: x(50) = m1 - m2 - d2
x(51) = f:   x(52) = m1 + f: x(53) = m1 - f: x(54) = d2 - f: x(55) = d2 + f - m1
x(56) = d2 - f - m1: x(57) = d2 + f: x(58) = n2 + f: x(59) = d2 + m1 - f: x(60) = n2 - f
x(61) = d2 - f - n2: x(62) = d2 + f + m1: x(63) = d4 - f - m1: x(64) = f3: x(65) = f + d
x(66) = f - d: x(67) = f + n3: x(68) = d4 - f: x(69) = f + d4: x(70) = m1 - f3
x(71) = f + d4 - n2: x(72) = d2 - f3: x(73) = d2 + n2 - f: x(74) = n2 - f - d2: x(75) = n3 - f
x(76) = f + d2 + n2: x(77) = d2 - f - n3: x(78) = f + d4: x(79) = m1 + f3: x(80) = d2 - m - f
x(81) = f - m - d2: x(82) = d2 + f - m - m1: x(83) = d2 + f - m: x(84) = d2 - f - m - m1: x(85) = f - m + m1
x(86) = f + m: x(87) = m1 - m - f: x(88) = f + m + m1: x(89) = f - m - m1: x(90) = f - m
x(91) = d2 + m1 - m - f: x(92) = m + f + d2 - m1: x(93) = m + f + d2: x(94) = d2 + f - m + m1: x(95) = d2 - m2 - f
x(96) = m1:    x(97) = d2 - m1: x(98) = d2:  x(99) = n2:    x(100) = d2 + m1
x(101) = d:    x(102) = f2 - m1: x(103) = n3: x(104) = d4 - m1: x(105) = d4 - n2
x(106) = n2 - d2: x(107) = d2 + n2: x(108) = d4: x(109) = n3 - d2: x(110) = m1 + d
x(111) = f2 - d2: x(112) = f2 - d2 + m1: x(113) = d2 - m: x(114) = d2 - m - m1: x(115) = m1 - m
x(116) = m + m1: x(117) = m: x(118) = d2 + m: x(119) = d2 - m + m1: x(120) = d2 + m - m1
x(121) = m + d: x(122) = n2 - m: x(123) = n2 + m: x(124) = d4 - m - m1: x(125) = d2 - m2
'-----------------------------------------------------------------------
b1 = 0: p1 = .950724
   'Correzioni LONGITUDINE Luna
   FOR k = 1 TO 27: l1 = l1 + k(k) * SIN(x(k) * rad#):    NEXT
   FOR k = 28 TO 45: l1 = l1 + k(k) * SIN(x(k) * rad#) * y(7): NEXT
   FOR k = 46 TO 50: l1 = l1 + k(k) * SIN(x(k) * rad#) * y(7) * y(7): NEXT
   'Correzioni LATITUDINE Luna
   FOR k = 51 TO 79: b1 = b1 + k(k) * SIN(x(k) * rad#):   NEXT
   FOR k = 80 TO 94: b1 = b1 + k(k) * SIN(x(k) * rad#) * y(7): NEXT
   FOR k = 95 TO 95: b1 = b1 + k(k) * SIN(x(k) * rad#) * y(7) * y(7): NEXT
   'Correzioni PARALLASSE Luna
   FOR k = 96 TO 112: p1 = p1 + k(k) * COS(x(k) * rad#):  NEXT
   FOR k = 113 TO 124: p1 = p1 + k(k) * COS(x(k) * rad#) * y(7): NEXT
   FOR k = 125 TO 125: p1 = p1 + k(k) * COS(x(k) * rad#) * y(7) * y(7): NEXT
'-----------------------------------------------------------------------
 w1 = .0004664 * COS(y(6) * rad#): w2 = .0000754 * COS(s3 * rad#)
     B = b1 * (1 - w1 - w2): l = l1: p = p1
 '--------- l = LONG.  Geocentrica LUNA ----------- Accuratezza: 10"
 '--------- b = LATIT. Geocentrica LUNA ----------- Accuratezza:  3"
 '--------- p = PARALLASSE Orizz.Equatoriale LUNA - Accuratezza:0.2"
   eps = si
ar = ATN((SIN(l * rad#) * COS(eps * rad#) - TAN(B * rad#) * SIN(eps * rad#)) / COS(l * rad#))
   IF COS(l * rad#) < 0 THEN ar = ar + pi#
   IF ar < 0 THEN ar = ar + 2 * pi#
ar = ar / rad#:                               ' ar = AR Luna (à-luna)
 ' ------------------------------------------------------------------
 x = SIN(B * rad#) * COS(eps * rad#) + COS(B * rad#) * SIN(eps * rad#) * SIN(l * rad#)
   de = ATN(x / SQR(-x * x + 1)): de = de / rad#: sl = p * 3600 * .272476

      ' de = Decl.Luna (ë-luna) ; sl = Raggio geocentrico Luna (in " d'arco)

  'Calcolo del T-sid a Greenwich all'ora fissata
   wn1 = (ha - fuso) + mn / 60: torig = wn1 / 876600
   mat1 = .276919398# + 100.0021359# * (t - torig) + 1.075E-06 * (t - torig) ^ 2
   mat2 = mat1 - INT(mat1)
   tsid = 1.002737908# * wn1 + mat2 * 24
   tsidgr = tsid * 15
   '
   CLS 2
   COLOR 15

   LOCATE 2, 5: PRINT USING "##/##/####   Diff.Fuso Orario:## h"; gdat; mdat; adat; fuso
   LOCATE 2, 44: PRINT USING "Latit.= ###.###   Long.Est=####.###"; lt; lg
   LOCATE 23, 68: PRINT USING "ORE ##:##"; ore2; min2
   COLOR 13
   aogreen = tsidgr - so:   ' ANGOLO ORARIO del SOLE a Greenwich
       decl = sq
   GOSUB altaz: hs1 = hs: azs1 = azs
  '
   aogreen = tsidgr - ar:   ' ANGOLO ORARIO della  LUNA a Greenwich
       decl = de
   GOSUB altaz: hs2 = hs: azs2 = azs
  '
 LOCATE 22, 6: PRINT USING "SOLE => h=###.##   Az-SW=####.##"; hs1; azs1
 LOCATE 22, 42: PRINT USING "LUNA => h=###.##   Az-SW=####.##"; hs2; azs2
  '
  'WHILE INKEY$ = "": WEND
  '*********************************************************************
  'Calcolo Parametri orbitali medi
               FOR k = 1 TO 7
     l(k) = ((z(k, 1, 3) * t + z(k, 1, 2)) * t + z(k, 1, 1)) * t + z(k, 1, 0)
                           xxx = l(k): GOSUB mo: l(k) = xxx
 amagg(k) = ((z(k, 2, 3) * t + z(k, 2, 2)) * t + z(k, 2, 1)) * t + z(k, 2, 0)
   ecc(k) = ((z(k, 3, 3) * t + z(k, 3, 2)) * t + z(k, 3, 1)) * t + z(k, 3, 0)
  incl(k) = ((z(k, 4, 3) * t + z(k, 4, 2)) * t + z(k, 4, 1)) * t + z(k, 4, 0)
periel(k) = ((z(k, 5, 3) * t + z(k, 5, 2)) * t + z(k, 5, 1)) * t + z(k, 5, 0)
                           xxx = periel(k): GOSUB mo: periel(k) = xxx
 lnodo(k) = ((z(k, 6, 3) * t + z(k, 6, 2)) * t + z(k, 6, 1)) * t + z(k, 6, 0)
                           xxx = lnodo(k): GOSUB mo: lnodo(k) = xxx
   '------------------------------------------------------------------
 manom(k) = l(k) - periel(k) - lnodo(k) + 360
                           xxx = manom(k): GOSUB mo: manom(k) = xxx
 ekepl(k) = manom(k)
99  ekepl(k) = ekepl(k) + (manom(k) + ecc(k) / rad# * SIN(ekepl(k) * rad#) - ekepl(k)) / (1 - ecc(k) * COS(ekepl(k) * rad#))

    IF ABS(ekepl(k) - flag) > .000001 THEN flag = ekepl(k): GOTO 99
  '
  vanom(k) = 2 * ATN(TAN(ekepl(k) / 2 * rad#) * SQR((1 + ecc(k)) / (1 - ecc(k)))) / rad#
  IF vanom(k) < 0 THEN vanom(k) = vanom(k) + 360
  '
  arglat(k) = l(k) + vanom(k) - manom(k) - lnodo(k) + 360
                           xxx = arglat(k): GOSUB mo: arglat(k) = xxx
  '
  nv = ATN(COS(incl(k) * rad#) * TAN(arglat(k) * rad#)) / rad#
    IF COS(arglat(k) * rad#) < 0 THEN nv = nv + 180
  lonelio(k) = nv + lnodo(k) + 360
                           xxx = lonelio(k): GOSUB mo: lonelio(k) = xxx
  sb = SIN(arglat(k) * rad#) * SIN(incl(k) * rad#)
  latelio(k) = ATN(sb / SQR(1 - sb ^ 2)) / rad#
  '
  rvett(k) = amagg(k) * (1 - ecc(k) * COS(ekepl(k) * rad#))
  '
  numer = rvett(k) * COS(latelio(k) * rad#) * SIN((lonelio(k) - sf) * rad#)
  denom = rvett(k) * COS(latelio(k) * rad#) * COS((lonelio(k) - sf) * rad#) + sh
  tango = 1 / rad# * ATN(numer / denom)
  IF denom < 0 THEN tango = tango + 180
  longeo(k) = sf + tango
  ' IF longeo(k) > 360 THEN longeo(k) = longeo(k) - 360
  '
  dpian(k) = SQR(numer ^ 2 + denom ^ 2 + (rvett(k) * SIN(latelio(k) * rad#)) ^ 2)
  '
  sbet = rvett(k) / dpian(k) * SIN(latelio(k) * rad#)
  latgeo(k) = ATN(sbet / SQR(1 - sbet ^ 2)) / rad#
 ' ------------------------------------------------------------------
  aret = ATN((SIN(longeo(k) * rad#) * COS(eps * rad#) - TAN(latgeo(k) * rad#) * SIN(eps * rad#)) / COS(longeo(k) * rad#))
  IF COS(longeo(k) * rad#) < 0 THEN aret = aret + pi#
   IF aret < 0 THEN aret = aret + 2 * pi#
  alfa(k) = aret / rad#:                   ' AR Pianeta
 '                                        ---------------------------
 x = SIN(latgeo(k) * rad#) * COS(eps * rad#) + COS(latgeo(k) * rad#) * SIN(eps * rad#) * SIN(longeo(k) * rad#)
 decln = ATN(x / SQR(-x * x + 1))
  delta(k) = decln / rad#:                 ' Decl. Pianeta
 '
   aogreen = tsidgr - alfa(k):   ' ANGOLO ORARIO del Pianeta a Greenwich
       decl = delta(k)
   GOSUB altaz: altez(k) = hs: azim(k) = azs
 '
   NEXT
  '*********************************************************************
   'CLS :
   lim = 21: delta = lim / 10
     SCREEN 9: PALETTE 0, 0:
     COLOR 14:
    REM -----------------------------------------------------------------
        vid = .735:             'Fattore deformaz. video
        sta = .735 * 1.13:      'Fattore deformaz. stampa
    REM -----------------------------------------------------------------
     qw = vid:                   'Opzione VIDEO/STAMPA
     '
     IF qw = sta THEN kk = sta ELSE kk = vid
     IF qw = sta THEN qw = vid / 1.13: 'Costante Aspect Ratio del CERCHIO
            fs = 100 / 3.52: 'Fattore di scala video-stampante
            '
            st = 150:       '[1 cm = 100/3.52 = 28.409 pixels]
         '
          r1 = 100: 'Raggio equatore
           sc = fs / st:               'Scala disegno
            r = r1 * sc: xorig = lim / 2: yorig = lim / 2.6
          REM -----------------------------------------------------------

1    IF lim <= 0 THEN lim = 10
    WINDOW (0, lim)-(lim, 0):
3   VIEW SCREEN (15, 15)-(625, 333), , 12: 'Viewport

          REM -----------------------------------------------------------
                                              'Disegno Box-Mercator
           fi = -72 * rad#: lam = -180 * rad#: GOSUB subr: x1 = xx: y1 = y
           fi = 72 * rad#: lam = 180 * rad#: GOSUB subr: x2 = xx: y2 = y
           LINE (x1, y1)-(x2, y2), 11, B
                                               'Disegno Paralleli
          FOR j = -70 TO 70 STEP 10
          fi = j * rad#: GOSUB subr: y3 = y
          IF j = 0 THEN colo = 9:  ELSE colo = 2
          LINE (x1, y3)-(x2, y3), colo
          NEXT
            dx = (x2 - x1) / 24
                                                'Disegno Meridiani
          kon = 0
          FOR xxx = x1 + dx TO x2 - dx STEP dx
          kon = kon + 1
          LINE (xxx, y1)-(xxx, y2), 5
          IF kon = 18 OR kon = 6 THEN LINE (xxx, y1)-(xxx, y2), 12
          IF kon = 12 THEN LINE (xxx, y1)-(xxx, y2), 9
          NEXT
          COLOR 11
              LOCATE 3, 41: PRINT "S"
              LOCATE 3, 5: PRINT "N"
              LOCATE 3, 76: PRINT "N"
         COLOR 12
              LOCATE 3, 23: PRINT "E"
              LOCATE 3, 58: PRINT "W"
           COLOR 14
                                                ' PLOT Posiz. Sole
          fi = hs1 * rad#: lam = -azs1 * rad#
          GOSUB subr: x1 = xx: y1 = y
          CIRCLE (x1, y1), .3, 12, , , kk
          PAINT (x1, y1), 12
          PSET (x1, y1), 14
                                                ' PLOT Posiz. Luna
          fi = hs2 * rad#: lam = -azs2 * rad#
          GOSUB subr: x1 = xx: y1 = y
          CIRCLE (x1, y1), .25, 14, , , kk
          PAINT (x1, y1), 14
          PSET (x1, y1), 0
 '===================================================================
      FOR k = 1 TO 7: pittur = k + 8
                                                ' PLOT Posiz. Pianeti
          fi = altez(k) * rad#: lam = -azim(k) * rad#
         ' PRINT TAB(10); altez(k), azim(k)
          GOSUB subr: x1 = xx: y1 = y
          CIRCLE (x1, y1), .12, pittur, , , kk
          PAINT (x1, y1), pittur
          PSET (x1, y1), 15
          NEXT
 '===================================================================
   LOCATE 23, 6: PRINT "Me   Ve    Ma    Gi   Sa    Ur    Ne     Sole     Luna        "
   '
   cl = 9: x21 = lim / 14: y21 = lim / 14
   CIRCLE (x21, y21), .12, cl, , , kk: PAINT (x21, y21), cl
   cl = 10: x22 = x21 + lim / 14:
   CIRCLE (x22, y21), .12, cl, , , kk: PAINT (x22, y21), cl
   cl = 11: x23 = x22 + lim / 14:
   CIRCLE (x23, y21), .12, cl, , , kk: PAINT (x23, y21), cl
   cl = 12: x24 = x23 + lim / 14:
   CIRCLE (x24, y21), .12, cl, , , kk: PAINT (x24, y21), cl
   cl = 13: x25 = x24 + lim / 14:
   CIRCLE (x25, y21), .12, cl, , , kk: PAINT (x25, y21), cl
   cl = 14: x26 = x25 + lim / 14:
   CIRCLE (x26, y21), .12, cl, , , kk: PAINT (x26, y21), cl
   cl = 15: x27 = x26 + lim / 14:
   CIRCLE (x27, y21), .12, cl, , , kk: PAINT (x27, y21), cl
   '
   cl = 12: x28 = x27 + 1.5 * lim / 14:
   CIRCLE (x28, y21), .3, cl, , , kk: PAINT (x28, y21), cl
   cl = 14: x29 = x28 + 1.5 * lim / 14:
   CIRCLE (x29, y21), .25, cl, , , kk: PAINT (x29, y21), cl
 '===================================================================
                                   'Plot ECLITTICA per punti
     FOR i = 0 TO 360 STEP 5
 '
  alfaec = ATN(COS(eps * rad#) * TAN(i * rad#))
  IF COS(i * rad#) < 0 THEN alfaec = alfaec + pi#
   IF alfaec < 0 THEN alfaec = alfaec + 2 * pi#
  alfaec = alfaec / rad#:                   ' AR di un punto Eclittica
 '                                        ---------------------------
 sendec = SIN(eps * rad#) * SIN(i * rad#)
 cosdec = SQR(1 - sendec ^ 2):
 deltec = ATN(sendec / cosdec) / rad#:     '  Decl. di un punto Eclittica
 '
   aogreen = tsidgr - alfaec:   ' ANGOLO ORARIO 1 punto Eclitt. a Greenwich
       decl = deltec
   GOSUB altaz: hs3 = hs: azs3 = azs
 '
   fi = hs3 * rad#: lam = -azs3 * rad#
   GOSUB subr: x1 = xx: y1 = y
   PSET (x1, y1), 15
   IF i MOD 90 = 0 THEN CIRCLE (x1, y1), .05, 14, , , kk
   IF i = 0 THEN CIRCLE (x1, y1), .1, 13, , , kk '
   NEXT
 '===================================================================
                                               'Plot Equatore Celeste
   decl = 0
   FOR i = 0 TO 360 STEP 4
   aogreen = tsidgr - i
   GOSUB altaz: hs3 = hs: azs3 = azs
 '
   fi = hs3 * rad#: lam = -azs3 * rad#
   GOSUB subr: x1 = xx: y1 = y
   PSET (x1, y1), 11
   NEXT
'===========================================================================
  'COLOR 14
15 LOCATE 1, 1: INPUT "         (q)=+60';(w)=-60';(n)=+10';(m)=-10';(i)=-1';(j)=+1';(f)=fine"; q$
 COLOR 15
 IF q$ <> "f" AND q$ <> "n" AND q$ <> "i" AND q$ <> "m" AND q$ <> "j" AND q$ <> "q" AND q$ <> "w" THEN 15
 '
 IF q$ = "n" THEN mn = mn + 10: GOTO 1200
 IF q$ = "m" THEN mn = mn - 10: GOTO 1200
 '
 IF q$ = "i" THEN mn = mn - 1: GOTO 1200
 IF q$ = "j" THEN mn = mn + 1: GOTO 1200
 '
 IF q$ = "w" THEN mn = mn - 60: GOTO 1200
 IF q$ = "q" THEN mn = mn + 60: GOTO 1200

 '

  END
subr:
     p = LOG(TAN(pi# / 4 - fi / 2))
     x9 = -r * lam / 2 / pi#
     y9 = -p * r / (2 * pi# * 1)
     kx = -tx: ky = -ty
     xx = COS(te) * x9 - SIN(te) * y9 + kx + xorig
     y = (SIN(te) * x9 + COS(te) * y9 + ky + yorig) / qw
     RETURN
     '
altaz:
 aoloc = aogreen + lg
 s6 = SIN(lt * rad#) * SIN(decl * rad#) + COS(lt * rad#) * COS(decl * rad#) * COS(aoloc * rad#)
 c6 = SQR(1 - s6 ^ 2)
       hs = 1 / rad# * ATN(s6 / c6): 'hs = Altezza
 '
 s9 = SIN(decl * rad#): c9 = COS(decl * rad#)
 '
 tg1 = SIN(aoloc * rad#)
 tg2 = COS(aoloc * rad#) * SIN(lt * rad#) - s9 / c9 * COS(lt * rad#)
   tgaz = tg1 / tg2: azs = ATN(tgaz)
 '
   IF tg2 < 0 THEN azs = azs + pi#
   IF azs < 0 THEN azs = azs + 2 * pi#
   IF azs > pi# THEN azs = azs - 2 * pi#
   azs = azs / rad#:                  ' azs = Azimut Sud-Ovest
   RETURN
'
mo: xxx = xxx - INT(xxx / u) * u: RETURN
  '
    REM ============== Vettore C(i) ===============================
DATA 279.69668,   36000.76892,   .0003025:             ' L = long.media Sole
DATA 358.475833,  35999.04975,  -.00015,    -3.3e-6:   ' M = anom.media Sole
DATA .01675104,   -.0000418,    -1.26e-7:              ' e = eccentr.orb.Ter.
DATA 1.91946,-.004789,-.000014,.020094,-.0001,.000293: ' Eq.centro
DATA 23.452294,   -.0130125,   -16.4e-7,  5.03e-7:     ' î = Incl.asse terr.
    REM ============== Vettore A(k,i) =============================
DATA 270.434164, 481267.8831,  -1133e-6,  19e-7:       ' L'= long.media Luna
DATA 358.475833,  35999.0498,  -15e-5,   -33e-7:       ' M = anom.media Sole
DATA 296.104608, 477198.8491,   9192e-6,  144e-7:      ' M'= anom.media Luna
DATA 350.737486, 445267.1142,  -1436e-6,  19e-7:       ' D = elon.media Luna
DATA  11.250889, 483202.0251,  -3211e-6,  -3e-7:       ' F = Dist.Luna dal Nodo Ascend.
DATA 259.183275,  -1934.142,    2078e-6,  22e-7:       ' ê = long.Nodo Ascend.
DATA   1,        -2495e-6,     -752e-8,   0:           ' Coeff.moltipl.
    REM ============== Vettore K(k) ===============================
DATA 6.28875,   1.274018,   .658309,   .213616,  -.114336: ' Termini correttivi Long.LUNA
DATA  .058793,   .05332,   -.034718,   .015326,  -.012528: ' LON.10
DATA -.01098,    .010674,   .010034,   .008548,   .005162
DATA  .003996,   .003862,   .003665,   .002602,  -.002349: ' LON.20
DATA -.001773,  -.001595,  -.00111,    .000892,   .00055
DATA  .000538,   .000486,  -.185596,   .057212,   .045874: ' LON.30
DATA  .041024,  -.030465,  -.00791,   -.006783,   .005
DATA  .004049,   .002695,   .002396,  -.002125,   .00122:  ' LON.40
DATA -.000811,   .000761,   .000693,   .000598,   .000521
DATA  .002249,  -.002079,   .002059,   .000717,   .000704: ' LON.50
    REM ===========================================================
DATA 5.128189,   .280606,   .277693,   .173238,   .055413: ' Termini correttivi Lat. LUNA
DATA  .046272,   .032573,   .017198,   .009267,   .008823: ' Lat.10
DATA  .004323,   .0042,     .001828,  -.00175,   -.001487
DATA  .00133,    .001106,   .00102,    .000833,   .000781: ' Lat.20
DATA  .00067,    .000606,   .000597,   .00045,    .000439
DATA  .000423,   .000422,   .000331,  -.000283,   .008247: ' Lat.30
DATA  .003372,   .002472,   .002222,   .002072,   .001877
DATA -.001803,   .00157,   -.001481,   .001417,   .00135:  ' Lat.40
DATA  .000492,  -.000367,  -.000353,   .000317,   .000306
    REM ===========================================================
DATA  .051818,   .009531,   .007843,   .002824,   .000857: ' Termini correttivi Parall. Luna
DATA -.000271,  -.000198,   .000173,   .000167,   .000103: ' Par.10
DATA -.000084,   .000079,   .000072,  -.000033,  -.00003
DATA -.000029,  -.000023,   .000533,   .000401,   .00032:  ' Par.20
DATA -.000264,  -.000111,  -.000083,   .000064,  -.000063
DATA  .000041,   .000035,  -.000029,   .000019,   .000026: ' Par.30
    REM ===========================================================
' Dati orbitali Pianeti
DATA 178.179078,  149474.07078,   .0003011,   0:     'Mercurio
DATA   0.3870986,      0,         0,          0
DATA   0.20561421,  2046e-8,     -3e-8,       0
DATA   7.002881,   .0018608,   -183e-7,       0
DATA  28.753753,   .3702806,   1208e-7,       0
DATA  47.145944,  1.1852083,   1739e-7,       0
'
DATA 342.767053,   58519.21191,   .0003097,   0:     'Venere
DATA   0.7233316,      0,         0,          0
DATA   0.00682069, -4774e-8,     91e-9,       0
DATA   3.393631,   .0010058,     -1e-6,       0
DATA  54.384186,   .5081861,   -.0013864,     0
DATA  75.779647,   .89985,       41e-5,       0
'
DATA 293.737334,   19141.69551,   .0003107,   0:     'Marte
DATA   1.5236883,      0,         0,          0
DATA   0.09331290, 92064e-9,    -77e-9,       0
DATA   1.850333,  -.0006750,    126e-7,       0
DATA 285.431761,  1.0697667,    .0001313,     414e-8
DATA  48.786442,   .7709917,    -14e-7,      -533e-8
'
DATA 238.049257,    3036.301986,  .0003347,  -165e-8:    'Giove
DATA   5.202561,       0,         0,          0
DATA   0.04833475,  .00016418, -4676e-10,    -17e-10
DATA   1.308736,  -.0056961,     39e-7,       0
DATA 273.277558,   .5994317,    .00070405,    508e-8
DATA  99.443414,  1.0105300,    .00035222,   -851e-8
'
DATA 266.564377,    1223.509884,  .0003245,  -58e-8:     'Saturno
DATA   9.554747,       0,         0,          0
DATA   0.05589232, -.00034550,  -728e-9,      74e-11
DATA   2.492519,  -.0039189,   -1549e-8,      4e-8
DATA 338.307800,  1.0852207,    .00097854,    992e-8
DATA 112.790414,   .8731951,   -.00015218,   -531e-8
'
DATA 244.197470,     429.863546,  .0003160,  -60e-8:     'Urano
DATA  19.218140,       0,         0,          0
DATA   0.04634440, -.00002658,    77e-9,      0
DATA    .772464,   .0006253,     395e-7,      0
DATA  98.071581,   .9857650,   -.00107445,   -61e-8
DATA  73.477111,   .4986678,    .00131170,    0
'
DATA  84.457994,     219.885914,  .0003205,  -60e-8:     'Nettuno
DATA  30.109570,       0,         0,          0
DATA   0.00899704,  .00000633,   -2e-9,       0
DATA   1.779242,  -.0095436,    -91e-7,       0
DATA 276.045975,   .3256394,    .00014095,  4113e-9
DATA 130.681389,  1.0989350,    .00024987, -4718e-9
'
   REM ============== FINE PROGRAMMA: PIA_CB.bas ==================

