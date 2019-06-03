'Soluzione dell'EQUAZIONE di KEPLERO: E = M + e * SIN(E)
' Con due metodi tradizionali e con le Variabili Universali
'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
COLOR 14, 6: CLS
DEFDBL A-Z
DEF FNatan (y, x) = (1 + (x = 0)) * ATN(y / (x + (x = 0))) - pi * (x <= 0) - 2 * pi * (x > 0 AND y < 0) + (x = 0) * SGN(y) * pi / 2
DEF FNridu (x) = (x / (2 * pi) - (INT(x / (2 * pi)) - (x < 0))) * 2 * pi
DEF FNacos (x) = pi / 2 - (2 + (x < 1 AND x > -1)) * ATN(x / SQR(1 + x * x * (x < 1 AND x > -1)))

   pi = ATN(1) * 4: rad = pi / 180: kgauss = .01720209895#

'     mgra = 29.568230006#: e = .9: a = 1: 'q=0.1 t=30 days
'     mgra = 18.47902673#: e = .987: a = 7.692307692#: 'q=0.1 t=400 days
'     mgra = 149.7463873#: e = .3978305: a = .756656
'     mgra = 19.952457#:   e = .00684458: a = .7233322169
' ******************* Epoca: 1.5/1/2002 (Hansen) **************************
     mgra = 286.53759248#: e = .20560396#: a = .387093595#: 'Mercury
     mgra = 141.58474327#: e = .00684458#: a = .723331617#: 'Venus
     mgra = 358.49294412#: e = .01677126#: a = 1.000000219#: 'Earth
     mgra = 43.11829743#: e = .09326685#: a = 1.52368838#: 'Mars
     mgra = 81.09098333999999#: e = .04825382#: a = 5.202755032#: 'Jupiter
     mgra = 342.33429146#: e = .05606075#: a = 9.54285512#: 'Saturn
     mgra = 151.48239878#: e = .0469055#: a = 19.1909792#: 'Uranus
     mgra = 263.83036278#: e = .0085082#: a = 30.070947#:  'Neptune
' *************************************************************************

''     mgra = 204.79768#: e = .0044295: a = 12271: '(a) in km; Sat Lageos-1


m = mgra * rad
anomecc = m    'Valore iniziale di E assunto pari all'anom.media M
kount = 0
'GOTO 12    'SELEZIONA Laguerre(goto12-OUT) - Mikkola(goto12-IN)
 n = 8 ' Con n=8 -> Laguerre ð Newton
 DO
  F0 = anomecc - e * SIN(anomecc) - m
  F1 = 1 - e * COS(anomecc)
  F2 = e * SIN(anomecc)
  RAX = (n - 1) ^ 2 * F1 ^ 2 - n * (n - 1) * F0 * F2
  IF RAX < 0 THEN RAX = -RAX
    RAD1 = F1 - SQR(RAX): RAD2 = F1 + SQR(RAX)
  IF ABS(RAD1) > ABS(RAD2) THEN RADI = RAD1 ELSE RADI = RAD2
  variaz = -n * F0 / RADI
   kount = kount + 1
   anomecc = anomecc + variaz
   PRINT kount, anomecc, variaz
LOOP UNTIL ABS(variaz) < .0000000001#
GOTO 15
12 :  ' MIKKOLA
  alfa = (1 - e) / (4 * e + .5): beta = m / (8 * e + 1)
  z = SGN(m) * ABS(beta + SQR(beta ^ 2 + alfa ^ 3)) ^ (1 / 3)
s0 = z - alfa / 2: s = s0 - .078 * s0 * s0 ^ 2 * s0 ^ 2 / (1 + e)
anomecc = m + e * SGN(m) * (3 * s - 4 * s * s * s)'Valore iniziale di E
 DO
  F0 = anomecc - e * SIN(anomecc) - m
  F1 = 1 - e * COS(anomecc)
  variaz = -F0 / F1
   kount = kount + 1
   anomecc = anomecc + variaz
   PRINT kount, anomecc, variaz
LOOP UNTIL ABS(variaz) < .0000000001#
15 :
 tgv = TAN(anomecc / 2) * SQR((1 + e) / (1 - e))
 v = FNridu(2 * FNatan(tgv, 1))
 r = a * (1 - e * COS(anomecc))
'
   PRINT "                                DATI"
   PRINT USING " Anom. media M = ###.######  Eccentr.(e) = #.#######   Sem.mag (a) = ##.######"; mgra#; e#; a#
   PRINT "                              RISULTATI"
   PRINT " Anomalia Eccentrica E = ", anomecc / rad; "ø = "; anomecc; " rad"
   PRINT " Anomalia Vera       é = ", v / rad
   PRINT " Raggio Vettore      r = ", r
 PRINT USING " X-orb=##.########    Y-orb=##.########"; r * COS(v); r * SIN(v)
PRINT
PRINT "******************************** V.U. **********************************"
'
22 :
   deltaTp = m * a ^ 1.5 / kgauss: 'Giorni dal Perielio

  q = a * (1 - e): tau = kgauss * deltaTp: r0 = q
' Starting VALUE (x0) derived by this cubic equation -> [ x^3+6ùr0ùx-6ç=0 ]
  c1 = 3 * tau + SQR(9 * tau ^ 2 + 8 * r0 ^ 3): c2 = 3 * tau - SQR(9 * tau ^ 2 + 8 * r0 ^ 3)
  cc1 = SGN(c1) * ABS(c1) ^ (1 / 3): cc2 = SGN(c2) * ABS(c2) ^ (1 / 3)
  x0 = cc1 + cc2
'Generalized Conic Functions: S(z)=1/3!-z/5!+zý/7!...; C(z)=1/2!-z/4!+zý/6!....
' z=xý/a; x=unknown=Generalized Eccentric Anomaly
 DEF FNS (z) = 1 / 6 * ((((((((z / 342 - 1) * z / 272 + 1) * z / 210 - 1) * z / 156 + 1) * z / 110 - 1) * z / 72 + 1) * z / 42 - 1) * z / 20 + 1)
 DEF FNC (z) = 1 / 2 * ((((((((z / 306 - 1) * z / 240 + 1) * z / 182 - 1) * z / 132 + 1) * z / 90 - 1) * z / 56 + 1) * z / 30 - 1) * z / 12 + 1)
kount = 0
10 :  kount = kount + 1
      z0 = x0 ^ 2 / a
      F = x0 ^ 3 * FNS(z0) + r0 * x0 * (1 - z0 * FNS(z0)) - tau 'Kepler Eq.
      r = x0 ^ 2 * FNC(z0) + r0 * (1 - z0 * FNC(z0)) 'Radius Vector
      
     x = x0 - F / r: ' Newton-Raphson iteration
     
      IF ABS(F / r) < .000000001# THEN 20
      x0 = x
    PRINT USING " It-No=##   r= ##.########  x= ##.######### (+#.##^^^^)"; kount; r; x; ABS(F / r)
 GOTO 10
'
20 : teta = FNacos((q * (1 + e) / r - 1) / e) / rad
     IF tau < 0 THEN teta = 360 - teta
     COLOR 10: PRINT
     PRINT USING "     True Anom. v= ###.#####ø     r=####.######## UA  (ëx=+#.##^^^^)"; teta; r; ABS(F / r): COLOR 14
END


