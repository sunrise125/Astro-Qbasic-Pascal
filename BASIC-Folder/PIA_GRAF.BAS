REM ---------  EFFEMERIDI  TABULATE dei Pianeti  ------------
'
     SCREEN 9: PALETTE 0, 1: COLOR 14
   DEFDBL A-H, J, L-Z
                  pi# = 4 * ATN(1#): rad# = pi# / 180:
    lim = 10: ps = 5: lim0 = lim
    WINDOW (0, lim)-(lim, 0):
    VIEW SCREEN (15, 15)-(625, 333), , 12: 'Viewport
'
DIM c(20) AS DOUBLE
DIM z(8, 6, 3) AS DOUBLE
DIM l(8) AS DOUBLE, amagg(8) AS DOUBLE, ecc(8) AS DOUBLE
DIM incl(8) AS DOUBLE, periel(8) AS DOUBLE, lnodo(8) AS DOUBLE
DIM manom(8) AS DOUBLE, ekepl(8) AS DOUBLE, vanom(8) AS DOUBLE
DIM arglat(8) AS DOUBLE, lonelio(8) AS DOUBLE, latelio(8) AS DOUBLE
DIM rvett(8) AS DOUBLE, longeo(8) AS DOUBLE, latgeo(8) AS DOUBLE
DIM dpian(8) AS DOUBLE, alfa(8) AS DOUBLE, delta(8) AS DOUBLE
DIM x(8), y(8)
'
   ' Lettura dei Parametri
   FOR i = 1 TO 20: READ c(i): NEXT
   FOR i = 1 TO 8: READ pia$(i): NEXT
   FOR k = 1 TO 8: FOR i = 1 TO 6: FOR i1 = 0 TO 3: READ z(k, i, i1): NEXT: NEXT: NEXT
    '
     LINE (1, 8.25)-(9, 9), 11, BF
      LOCATE 4, 18: PRINT " EFFEMERIDI Planetarie & GRAFICO  [ 1 ö 8 ]   "
      LINE (1, 5)-(5, 8), 2, BF: COLOR 12
      '
34    LOCATE 7, 15: INPUT "GIORNO (gg)"; g: IF g < 1 OR g > 31 THEN 34
35    LOCATE 9, 15: INPUT "MESE   (mm)"; m: IF m < 1 OR m > 12 THEN 35
     LOCATE 11, 15: INPUT "ANNO (aaaa)"; aa
     gdat = g: mdat = m: adat = aa
      a8 = aa: u7 = 2 - INT(aa / 100) + INT(INT(aa / 100) / 4)
      IF m >= 3 THEN 1100
      a8 = a8 - 1: m = m + 12
1100  jd = INT(365.25 * a8) + INT(30.6001 * (m + 1)) + g + 1720994.5#
      IF jd >= 2299160.5# THEN jd = jd + u7
      n0 = jd - 2415020#
      '
      je = n0: u = 360
      '
      LINE (5, 5)-(9, 8), 5, BF: COLOR 14
44    LOCATE 7, 45: INPUT "ORE    [T.U.]"; ha: IF ha < 0 OR ha > 23 THEN 44
45    LOCATE 9, 45: INPUT "MINUTI [T.U.]"; mn: IF mn < 0 OR mn > 59 THEN 45
      'LOCATE 10, 45: INPUT "Ore FUSO da GREENW."; fuso
      CLS
1200  jq = je + (ha - fuso) / 24 + mn / 1440:     t = jq / 36525
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
 'anecc=Anomalia eccentrica (E)
  anecc = 2 * ATN(TAN(sg / 2 * rad#) * SQR((1 - jc) / (1 + jc))) / rad#
  IF anecc < 0 THEN anecc = anecc + 360
  '
  si = ((c(20) * t + c(19)) * t + c(18)) * t + c(17): ' Obliquit… Eclitt.(î)
     sl = SIN(si * rad#)
     sm = COS(si * rad#)
 so = 1 / rad# * ATN(SIN(sf * rad#) * sm / COS(sf * rad#))
     IF COS(sf * rad#) < 0 THEN so = so + u / 2
     IF so < 0 THEN so = so + u
 s9 = sl * SIN(sf * rad#): c9 = SQR(1 - s9 ^ 2):   'so =   AR Sole(à-sole)
 sq = 1 / rad# * ATN(s9 / c9):                     'sq = Decl.Sole(ë-sole)
  '*********************************************************************
  CLS 2
  GOSUB logo
  '*********************************************************************
  'Calcolo Parametri orbitali Pianeti
               FOR k = 1 TO 8
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
  longeo(k) = sf + tango: xxx = longeo(k): GOSUB mo: longeo(k) = xxx
  dpian(k) = SQR(numer ^ 2 + denom ^ 2 + (rvett(k) * SIN(latelio(k) * rad#)) ^ 2)
  '
  sbet = rvett(k) / dpian(k) * SIN(latelio(k) * rad#)
  latgeo(k) = ATN(sbet / SQR(1 - sbet ^ 2)) / rad#
 ' ------------------------------------------------------------------
  eps = si
  aret = ATN((SIN(longeo(k) * rad#) * COS(eps * rad#) - TAN(latgeo(k) * rad#) * SIN(eps * rad#)) / COS(longeo(k) * rad#))
  IF COS(longeo(k) * rad#) < 0 THEN aret = aret + pi#
   IF aret < 0 THEN aret = aret + 2 * pi#
  alfa(k) = aret / rad#:                   ' AR Pianeta
 '                                        ---------------------------
 x = SIN(latgeo(k) * rad#) * COS(eps * rad#) + COS(latgeo(k) * rad#) * SIN(eps * rad#) * SIN(longeo(k) * rad#)
 decln = ATN(x / SQR(-x * x + 1))
  delta(k) = decln / rad#:                 ' Decl. Pianeta
  '
  IF k = 3 THEN longeo(k) = 0
  LOCATE k + 4, 16
  PRINT USING "³ \  \³###.#####  ###.#####³ ###.#####  ###.#####³"; pia$(k); lonelio(k); latelio(k); longeo(k); latgeo(k)
  '
   NEXT
 GOSUB logo1
 alfa(3) = 0: delta(3) = 0
  FOR k = 1 TO 8
  LOCATE k + 15, 16
  PRINT USING "³ \  \³###.#####  ###.#####³ ###.#####  ###.#####³"; pia$(k); rvett(k); dpian(k); alfa(k); delta(k)
  NEXT

  '*********************************************************************
 WHILE INKEY$ = "": WEND
'
CLS : delta = lim / 10
     SCREEN 9: COLOR 14
    
    REM -----------------------------------------------------------------
        vid = .735:             'Fattore deformaz. video
        sta = .735 * 1.13:      'Fattore deformaz. stampa
    REM -----------------------------------------------------------------
     qw = vid:                   'Opzione VIDEO/STAMPA
     '
     IF qw = sta THEN kk = sta ELSE kk = vid
     IF qw = sta THEN qw = vid / 1.13: 'Costante Aspect Ratio del CERCHIO
        REM ==========================================================
              xorig = lim / 2: yorig = lim / 2 * qw
              'yorig = lim / 3
        REM -----------------------------------------------------------
         '
 '
1   IF lim <= 0 THEN lim = lim0
    WINDOW (0, lim)-(lim, 0):
    xorig = lim / 2: yorig = lim / 2 * qw
3   VIEW SCREEN (15, 15)-(625, 333), , 12: 'Viewport
'
 '===================================================================
   COLOR 9
   LOCATE 23, 6: PRINT "Me    Ve   Te    Ma    Gi   Sa    Ur    Ne"
   sr = lim / 14: dm = .06
   '
   cl = 9: x21 = sr: y21 = sr
   CIRCLE (x21, y21), dm, cl, , , kk: PAINT (x21, y21), cl
   cl = 10: x22 = x21 + sr
   CIRCLE (x22, y21), dm, cl, , , kk: PAINT (x22, y21), cl
   cl = 6: x23 = x22 + sr:
   CIRCLE (x23, y21), dm, cl, , , kk: PAINT (x23, y21), cl
   cl = 11: x24 = x23 + sr
   CIRCLE (x24, y21), dm, cl, , , kk: PAINT (x24, y21), cl
   cl = 12: x25 = x24 + sr
   CIRCLE (x25, y21), dm, cl, , , kk: PAINT (x25, y21), cl
   cl = 13: x26 = x25 + sr
   CIRCLE (x26, y21), dm, cl, , , kk: PAINT (x26, y21), cl
   cl = 14: x27 = x26 + sr
   CIRCLE (x27, y21), dm, cl, , , kk: PAINT (x27, y21), cl
   cl = 15: x28 = x27 + sr
   CIRCLE (x28, y21), dm, cl, , , kk: PAINT (x28, y21), cl
   '===================================================================
   COLOR 15
   LOCATE 23, 52: PRINT USING "##/##/#### - ##h ##m T.U."; gdat; mdat; adat; ha; mn
  '
    'ESECUZIONE Disegno ......... per punti
            x = 0: y = 0: GOSUB subr
            CIRCLE (x, y), .2, 12, , , kk
            PAINT (x, y), 12
         '--------------------------------
             kcolo = 1
             '
            FOR k = 1 TO 8
           x(k) = rvett(k) * COS(lonelio(k) * rad#)
           y(k) = rvett(k) * SIN(lonelio(k) * rad#)
           x = x(k): y = y(k)
           GOSUB subr
             x2 = x: y2 = y
             x = 0: y = 0: GOSUB subr: x1 = x: y1 = y
          LINE (x2, y2)-(x1, y1), 2, , &HAA00
           CIRCLE (x1, y1), rvett(k), 13, , , kk
             pittur = k + 7 + kcolo
             IF k = 3 THEN pittur = 6: kcolo = 0
             x = x(k): y = y(k)
          GOSUB subr
          CIRCLE (x, y), .1, pittur, , , kk
          PAINT (x, y), pittur
          PSET (x, y), 15
           NEXT

        COLOR 14
        LOCATE 1, 73: PRINT "Pas="; ps: LOCATE 1, 64: PRINT "Quad="; lim
15 LOCATE 1, 1: INPUT "Opz:(+);(-);(d)=dx;(s)=sx;(a)=alto;(b)=bas;(r)=rot;(f)=fine"; q$
 COLOR 15
 IF q$ <> "+" AND q$ <> "-" AND q$ <> "f" AND q$ <> "d" AND q$ <> "s" AND q$ <> "a" AND q$ <> "b" AND q$ <> "r" THEN 15
 IF q$ = "-" THEN lim = lim + ps: CLS 2: GOTO 1
 IF q$ = "+" THEN lim = lim - ps: CLS 2: GOTO 1
 '
 IF q$ = "s" THEN dd = dd + delta: tx = dd: CLS 2: GOTO 3
 IF q$ = "d" THEN dd = dd - delta: tx = dd: CLS 2: GOTO 3
 '
 IF q$ = "b" THEN dd1 = dd1 + delta: ty = dd1: CLS 2: GOTO 3
 IF q$ = "a" THEN dd1 = dd1 - delta: ty = dd1: CLS 2: GOTO 3
 '
 IF q$ = "r" THEN te = te + pi# / 8: CLS 2: GOTO 3

      END
'
subr:
     x9 = x
     y9 = y
     kx = -tx: ky = -ty
     x = COS(te) * x9 - SIN(te) * y9 + kx + xorig
     y = (SIN(te) * x9 + COS(te) * y9 + ky + yorig) / qw
     RETURN
'
mo: xxx = xxx - INT(xxx / u) * u: RETURN
  '
logo: SCREEN 0
      COLOR 10, 4: PALETTE 0, 1
      LOCATE 1, 23: PRINT USING "Data : ##/##/####   ore T.U.-> ##:##"; gdat; mdat; adat; ha; mn
      LOCATE 3, 16: PRINT "³Pian.³ LG_elio    LT_elio ³  LG_geo      LT_geo ³"
      LOCATE 4, 16: PRINT "³     ³       [gradi]      ³       [gradi]       ³"
      COLOR 14, 1
      RETURN
'
logo1: SCREEN 0
      COLOR 10, 4: PALETTE 0, 1
     LOCATE 14, 16: PRINT "³Pian.³ d_sole     d_terra ³   A.R.       DECL.  ³"
     LOCATE 15, 16: PRINT "³     ³      [ U.A.]       ³       [gradi]       ³"
      COLOR 14, 1
      RETURN
'


    REM ============== Vettore C(i) ===============================
DATA 279.69668,   36000.76892,   .0003025:             ' L = long.media Sole
DATA 358.475833,  35999.04975,  -.00015,    -3.3e-6:   ' M = anom.media Sole
DATA .01675104,   -.0000418,    -1.26e-7:              ' e = eccentr.orb.Ter.
DATA 1.91946,-.004789,-.000014,.020094,-.0001,.000293: ' Eq.centro
DATA 23.452294,   -.0130125,   -16.4e-7,  5.03e-7:     ' î = Incl.asse terr.
    REM ===========================================================
   'Nomi Pianeti
DATA "MER","VEN","TER","MAR","GIO","SAT","URA","NET"
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
DATA  99.69668,   36000.76892,   .0003025,    0:     'Terra
DATA   1.0000002,      0,         0,          0
DATA   0.01673104, -418e-7,    126e-9,        0
DATA   0,              0,         0,          0
DATA 101.22085,   1.71917,    4525e-7,        0
DATA   0,              0,         0,          0
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
   REM ============== FINE PROGRAMMA: PIA_GRAF.bas ==================
     '

