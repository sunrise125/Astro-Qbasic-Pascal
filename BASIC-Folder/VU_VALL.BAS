' ---------------- Adattamento a moto sat geocentrico ----- (Vallado) -- 14.4.2013
COLOR 14, 1: CLS
CLEAR , , 2000
PRINT "# gauss # preliminary orbital elements"
PRINT "# by the method of gauss"
PRINT
DEFDBL A-Z
DEFINT I-K, N
'
DEF fnf (x) = c + x * x * x * (b + x * x * x * (a + x * x)): ' f(x)=equaz. di lagrange
DEF fndf (x) = x * x * (3 * b + x * x * x * (6 * a + 8 * x * x)): ' f'(x)
DEF fnp (x) = aa + m * bb / (x * x * x): ' rho = a + b/(r2*r2*r2)
kgauss# = .0012394466678# ' moto geoc.    '= .01720209895#
RT = 6378.14#: uv0 = 7.905364437#
m = 1:  ' æ = 1 -> unit… canoniche
    pi = 4 * ATN(1#): rad = pi / 180
ab = 0#   ' coeff. per velocità luce -> ab = 1 / 173.1446#  c=173.1446 UA/g
nr = 0: p1 = 0: p2 = 0: p3 = 0
g$ = "####.#######"
'
DIM XXX(30)
READ n$, e$
  READ t1, Alfa1, Delt1, XS1, YS1, ZS1
  READ t2, Alfa2, Delt2, XS2, YS2, ZS2
  READ t3, Alfa3, Delt3, XS3, YS3, ZS3
IF kgauss# > .07 THEN ab = 0
'
PRINT n$ + " - "; e$: PRINT
'
  a1 = COS(Delt1 * rad) * COS(Alfa1 * rad)
  b1 = COS(Delt1 * rad) * SIN(Alfa1 * rad)
  c1 = SIN(Delt1 * rad)
'
  a2 = COS(Delt2 * rad) * COS(Alfa2 * rad)
  b2 = COS(Delt2 * rad) * SIN(Alfa2 * rad)
  c2 = SIN(Delt2 * rad)
'
  a3 = COS(Delt3 * rad) * COS(Alfa3 * rad)
  b3 = COS(Delt3 * rad) * SIN(Alfa3 * rad)
  c3 = SIN(Delt3 * rad)

PRINT USING " t1 = ####.#####  à1 = ##.######    ë1 = +##.######"; t1; Alfa1; Delt1
PRINT USING "                  XS1= +#.#######   YS1=  +#.#######   ZS1= +#.#######"; XS1; YS1; ZS1
PRINT
PRINT USING " t2 = ####.#####  à2 = ##.######    ë2 = +##.######"; t2; Alfa2; Delt2
PRINT USING "                  XS2= +#.#######   YS2=  +#.#######   ZS2= +#.#######"; XS2; YS2; ZS2
PRINT
PRINT USING " t3 = ####.#####  à3 = ##.######    ë3 = +##.######"; t3; Alfa3; Delt3
PRINT USING "                  XS3= +#.#######   YS3=  +#.#######   ZS3= +#.#######"; XS3; YS3; ZS3
'LINE INPUT ""; l$
PRINT
PRINT "                    *** numerical coefficients ***"
tt1 = kgauss# * (t1 - t2)
tt3 = kgauss# * (t3 - t2)
tt = tt3 - tt1
'------------------------------------------------------------------
e1 = a1 * (b2 * c3 - b3 * c2)
e2 = b1 * (a2 * c3 - a3 * c2)
e3 = c1 * (a2 * b3 - a3 * b2):       d0 = e1 - e2 + e3
'------------------------------------------------------------------
  e1 = XS1 * (b2 * c3 - b3 * c2)
  e2 = YS1 * (a2 * c3 - a3 * c2)
  e3 = ZS1 * (a2 * b3 - a3 * b2):   d11 = e1 - e2 + e3
  e1 = XS2 * (b2 * c3 - b3 * c2)
  e2 = YS2 * (a2 * c3 - a3 * c2)
  e3 = ZS2 * (a2 * b3 - a3 * b2):   d12 = e1 - e2 + e3
  e1 = XS3 * (b2 * c3 - b3 * c2)
  e2 = YS3 * (a2 * c3 - a3 * c2)
  e3 = ZS3 * (a2 * b3 - a3 * b2):   d13 = e1 - e2 + e3
'------------------------------------------------------------------
  e1 = a1 * (YS1 * c3 - b3 * ZS1)
  e2 = b1 * (XS1 * c3 - a3 * ZS1)
  e3 = c1 * (XS1 * b3 - a3 * YS1):  d21 = e1 - e2 + e3
  e1 = a1 * (YS2 * c3 - b3 * ZS2)
  e2 = b1 * (XS2 * c3 - a3 * ZS2)
  e3 = c1 * (XS2 * b3 - a3 * YS2):  d22 = e1 - e2 + e3
  e1 = a1 * (YS3 * c3 - b3 * ZS3)
  e2 = b1 * (XS3 * c3 - a3 * ZS3)
  e3 = c1 * (XS3 * b3 - a3 * YS3):  d23 = e1 - e2 + e3
'------------------------------------------------------------------
  e1 = a1 * (b2 * ZS1 - YS1 * c2)
  e2 = b1 * (a2 * ZS1 - XS1 * c2)
  e3 = c1 * (a2 * YS1 - XS1 * b2):  d31 = e1 - e2 + e3
  e1 = a1 * (b2 * ZS2 - YS2 * c2)
  e2 = b1 * (a2 * ZS2 - XS2 * c2)
  e3 = c1 * (a2 * YS2 - XS2 * b2):  d32 = e1 - e2 + e3
  e1 = a1 * (b2 * ZS3 - YS3 * c2)
  e2 = b1 * (a2 * ZS3 - XS3 * c2)
  e3 = c1 * (a2 * YS3 - XS3 * b2):  d33 = e1 - e2 + e3
'------------------------------------------------------------------
'
ee = (a2 * XS2 + b2 * YS2 + c2 * ZS2) * (-2)
ff = XS2 ^ 2 + YS2 ^ 2 + ZS2 ^ 2
'
aq1 = tt3 / tt
bq1 = aq1 * (tt * tt - tt3 * tt3) / 6
aq3 = -tt1 / tt
bq3 = aq3 * (tt * tt - tt1 * tt1) / 6
'
aa = (aq1 * d21 - d22 + aq3 * d23) / (-d0)
bb = (bq1 * d21 + bq3 * d23) / (-d0)
'
a = -(aa * aa + aa * ee + ff):     ' a,b,c = coeff. dell'equaz. di lagrange
b = -(m) * (2 * aa * bb + bb * ee)
c = -(m * m) * (bb * bb)
'
PRINT "d0"; TAB(9); : PRINT USING g$; d0

PRINT "aa"; TAB(9); : PRINT USING g$; aa
PRINT "bb"; TAB(9); : PRINT USING g$; bb

PRINT "ee"; TAB(9); : PRINT USING g$; ee
PRINT "ff"; TAB(9); : PRINT USING g$; ff

PRINT "         K"; TAB(29); : PRINT USING g$; a
PRINT "         L"; TAB(29); : PRINT USING g$; b
PRINT "         M"; TAB(29); : PRINT USING g$; c
'PRINT
LINE INPUT ""; l$
CLS
iflag = 0
 FOR qw = .1 TO 4.1 STEP .2
iflag = iflag + 1
x = qw
'
2340 IF ABS(fnf(x)) < .0000001# THEN 2400
  x0 = x
  x = x - fnf(x) / fndf(x): res = ABS(x0 - x)
GOTO 2340
2400 XXX(iflag) = x
PRINT USING "           Radice r2(##) = #.########    Residuo = #.##^^^^"; iflag; x; res
NEXT
PRINT : COLOR 12
INPUT "                         Scelta indice = ", ind: COLOR 14
x = XXX(ind)
'
 u2 = m / (x * x * x):
'
fw1 = 1 - u2 * tt1 * tt1 / 2
fw3 = 1 - u2 * tt3 * tt3 / 2
gw1 = tt1 * (1 - u2 * tt1 * tt1 / 6)
gw3 = tt3 * (1 - u2 * tt3 * tt3 / 6)
'
fg = fw1 * gw3 - fw3 * gw1
'
cp1 = gw3 / fg
cp2 = -1
cp3 = -gw1 / fg
'
dn1 = -fw3 / fg
dn3 = fw1 / fg
'
2620 CLS
PRINT "*** preliminary orbital elements ***"
PRINT
'-----------------------------------------------------------------------------
 pp1 = (cp1 * d11 + cp2 * d12 + cp3 * d13) / (cp1 * d0)
  x1 = pp1 * a1 - XS1: y1 = pp1 * b1 - YS1: z1 = pp1 * c1 - ZS1
  r1 = SQR(x1 ^ 2 + y1 ^ 2 + z1 ^ 2)
'
 pp2 = (cp1 * d21 + cp2 * d22 + cp3 * d23) / (cp2 * d0)
  x2 = pp2 * a2 - XS2: y2 = pp2 * b2 - YS2: z2 = pp2 * c2 - ZS2
  r2 = SQR(x2 ^ 2 + y2 ^ 2 + z2 ^ 2)
'
 pp3 = (cp1 * d31 + cp2 * d32 + cp3 * d33) / (cp3 * d0)
  x3 = pp3 * a3 - XS3: y3 = pp3 * b3 - YS3: z3 = pp3 * c3 - ZS3
  r3 = SQR(x3 ^ 2 + y3 ^ 2 + z3 ^ 2)
'-----------------------------------------------------------------------------
  V2x = dn1 * x1 + dn3 * x3
  V2y = dn1 * y1 + dn3 * y3
  V2z = dn1 * z1 + dn3 * z3: V2 = SQR(V2x ^ 2 + V2y ^ 2 + V2z ^ 2)
'-----------------------------------------------------------------------------
dp1 = pp1 - p1
dp2 = pp2 - p2
dp3 = pp3 - p3
'
  tabe1 = t1 - ab * pp1
  tabe2 = t2 - ab * pp2
  tabe3 = t3 - ab * pp3
'
tt1 = kgauss# * (tabe1 - tabe2)
tt3 = kgauss# * (tabe3 - tabe2)
tt = tt3 - tt1
'
PRINT USING " Rho-1 = ##.#######    Rho-2 = ##.#######    Rho-3 = ##.#######"; pp1; pp2; pp3
PRINT USING " Rho-1 = #####.### km  Rho-2 = #####.### km  Rho-3 = #####.### km"; pp1 * RT; pp2 * RT; pp3 * RT
PRINT
PRINT USING "    t2 = ####.#####       r2 = ##.#######       V2 = ##.#######"; t2; r2; V2
PRINT USING "                          r2 = #####.### km     V2 = #.#### km/s"; r2 * RT; V2 * uv0
PRINT
PRINT USING "  x2= +##.#######       y2= +##.#######        z2= +##.#######"; x2; y2; z2
PRINT USING "  x2= +#####.### km     y2= +#####.### km      z2= +#####.### km"; x2 * RT; y2 * RT; z2 * RT
PRINT
PRINT USING " V2x= +##.#######      V2y= +##.#######       V2z= +##.#######"; V2x; V2y; V2z
PRINT USING " V2x= +#.#### km/s     V2y= +#.#### km/s      V2z= +#.#### km/s"; V2x * uv0; V2y * uv0; V2z * uv0

PRINT : COLOR 10
PRINT USING "                 ##"; nr: COLOR 14
PRINT USING "       ë(Rho-1)= ##.########"; dp1
PRINT USING "       ë(Rho-2)= ##.########"; dp2
PRINT USING "       ë(Rho-3)= ##.########"; dp3
'
PRINT : COLOR 12
INPUT "    Risultati SODDISFACENTI (s/ )"; a$: COLOR 14
IF a$ = "s" OR a$ = "S" THEN 3520
   nr = nr + 1
p1 = pp1
p2 = pp2
p3 = pp3
'
f1 = fw1
f3 = fw3
g1 = gw1
g3 = gw3
'
   h = tt1: GOSUB 16010: ' fg/sub
   fw1 = f: gw1 = g
'
   h = tt3: GOSUB 16010
   fw3 = f: gw3 = g
'
fw1 = (fw1 + f1) / 2
fw3 = (fw3 + f3) / 2
gw1 = (gw1 + g1) / 2
gw3 = (gw3 + g3) / 2
'
fg = fw1 * gw3 - fw3 * gw1
'
cp1 = gw3 / fg
cp3 = -gw1 / fg
'
dn1 = -fw3 / fg
dn3 = fw1 / fg
'
GOTO 2620
'
3520 PRINT n$
PRINT "method of gauss"
END
 STOP
16010 ' # ufg/sub # universal f&g expressions
 '     unit… di misura
 '     secondo herrick -> [æ]=[l^3ùt^(-2)]   parametro gravitazionale ð goodyear
 '                                           = 1 in unit… canoniche
 '     [ç]=[ûæùt] = [l^(3/2)ùt^(-1)ùt] ==> [l^(3/2)]  unit… di tempo
 '     [r]                             ==> [l]        unit… di lunghezza
 '     [r/ç] = [l/l^(3/2)]             ==> [l^(-1/2)] unit… di velocit…
 '

 sigma2 = (x2 * V2x + y2 * V2y + z2 * V2z) / SQR(m)
         ' å2 = (r2ùv2)/ûæ = (x2ùV2x+y2ùV2y+z2ùV2z)/ûæ = prodotto scalare
 ai = 2 / r2 - V2 ^ 2 / m: ' ai = 1/a = 2/r2-V2ý/æ
 co = 1 - r2 * ai: ' c0=1-r2/a
 ww = h * SQR(m):  ' ç ð ëtùûæ
 xx = ww / r2:  ' xðxx = anomalia eccentrica generalizzata=  x=ç/r0 (1^approx.)
                ' convenzione goodyear: psi=anom.ecc.gener.= x/ûæ

16090    xq2 = xx * xx: xq3 = xq2 * xx:  '  xý e x^3
         xa = -xq2 * ai: ' = -xý/a   --> corrisponde a á=àùpsiý di goodyear

   'funzioni coniche generalizzate ->  uu=s3ùæ^(3/2), cc=s2ùæ, ss=s1ùûæ
 uu = xq3 / 6 * ((((((((xa / 342 + 1) * xa / 272 + 1) * xa / 210 + 1) * xa / 156 + 1) * xa / 110 + 1) * xa / 72 + 1) * xa / 42 + 1) * xa / 20 + 1)
 cc = xq2 / 2 * ((((((((xa / 306 + 1) * xa / 240 + 1) * xa / 182 + 1) * xa / 132 + 1) * xa / 90 + 1) * xa / 56 + 1) * xa / 30 + 1) * xa / 12 + 1)
 ss = xx - uu * ai
    fx = r2 * xx + co * uu + sigma2 * cc - ww: 'equaz. keplero generalizzata
                                         ' f(x)= r0ùx + c0ùuu + d0ùcc - ç
    IF ABS(fx) < .0000001# THEN 16200
    df = r2 + co * cc + sigma2 * ss:  ' r=f'(x)=df(x)/dx = r0 + c0ùcc + d0ùss
                                   '    (cc=duu/dx; ss=dcc/dx)
    xx = xx - fx / df: ' iterazione di newton-raphson
    GOTO 16090
16200 f = 1 - cc / r2:                    ' f = 1 - cc/r0
      g = (r2 * ss + sigma2 * cc) / SQR(m):  ' g = (ç-uu)/ûæ = (r0ùss+d0ùcc)/ûæ
'           ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿
'           ³  ->    ->     -> ³
'           ³  r = fùr0 + gùv0 ³        raggio vettore al tempo t (=t0+ç)
'           ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ
 RETURN

'
' ----------------------------- RAAN & Dec in gradi ----------------
 DATA "sat447", "Vallado"
 '
 DATA  0.0#
 DATA  0.939913#, 18.667717#
 DATA  -0.63369117#,-0.432729918#,-0.639570627#
 '
 DATA  480.0#
 DATA  45.025748#, 35.664741#
 DATA  -0.618159682#,-0.454641008#,-0.639570627#
 '
 DATA  720.0#
 DATA  67.886635#, 36.996583#
 DATA  -0.610108746#,-0.465389282#,-0.639570627#

