'Soluzione di un'EQUAZIONE trascendente F(x)=0 per Succ.approx. NEWTON
'applicata al calcolo dell'anomalia di un settore ellittico nota l'area
'                       Settorel.bas (Feb.05)
'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
COLOR 14, 1: CLS
DEFDBL A-Z
pi = 4 * ATN(1): rad = pi / 180
DEF fnc (x) = 2 / (1 - e ^ 2) ^ (1.5) * ATN(SQR((1 - e) / (1 + e)) * TAN(x / 2)) - e * SIN(x) / (1 - e ^ 2) / (1 + e * COS(x)) - K0: 'Funzione
DEF fnc1 (x) = 1 / (1 + e * COS(x)) ^ 2: 'Derivata prima
'-------------------
 e = .2:  'eccentricita'
 eta = .35: ' frazione di mezza ellisse (0 < eta < 1)
 K0 = pi * eta * (1 - e ^ 2) ^ (-1.5)
'-------------------
' Ciclo iterativo
x = 1:       'Valore iniziale in radianti
'             (da modificare se F(x) diverge; varia da 0 a pi)
 DO
variaz = -fnc(x) / fnc1(x)
   kount = kount + 1
   x = x + variaz
   PRINT kount, x, variaz
LOOP UNTIL ABS(variaz) < .00000001#
 PRINT : PRINT USING "  x= ##.####### rad =  ###.####ø"; x; x / rad

END
'**************** Programma: Settorel.bas (Feb.05) ************************

