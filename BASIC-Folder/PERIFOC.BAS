' ----------- Versori Perifocali (P,Q) orbita ellittica ------ Gen.05 ------
'                   (Sistema di riferimento eclittico)      Perifoc.bas
' ----------------------------------------------------------------------------
 COLOR 14, 1, 6: CLS
 DEFDBL A-Z
 pi = 4 * ATN(1): rad = pi / 180
'
' Px = COSwùCOSê - SINwùSINêùCOSi =[Rx]; Qx = -SINwùCOSê - COSwùSINêùCOSi =[x']
' Py = COSwùSINê + SINwùCOSêùCOSi =[Ry]; Qy = -SINwùSINê + COSwùCOSêùCOSi =[y']
' Pz = SINwùSINi                  =[Rz]; Qz =  COSwùSINi                  =[z']
'
'  Nota IMPORTANTE: Se al posto di (w) si sostituisce l'argomento della
'       ÜÜÜÜÜÜÜÜÜÜ  latitudine (u=w+f) con f=anomalia vera del corpo,
'                   il versore P diventa il versore posizione (Rx,Ry,Rz),
'                   mentre entrambi (P,Q) concorrono a determinare
'                   il versore velocità (x',y',z').
'
'                   Quindi:  x= rùRx
'    Vettore POSIZIONE --->  y= rùRy  con r= p/(1+eùCOSf)=p=aù(1-eý)/(1+eùCOSf)
'                            z= rùRz
'
'                            Vx= k/ûpù[RxùeùSINf+(1+eùCOSf)ùx']
'    Vettore VELOCITA' --->  Vy= k/ûpù[RyùeùSINf+(1+eùCOSf)ùy']
'                            Vz= k/ûpù[RzùeùSINf+(1+eùCOSf)ùz']
'
'           con  k = ûæ= COSTANTE DINAMICA DEL MOTO = 1  (unit… canoniche)
'           oppure:
'           1) moto eliocentrico: k=k_gauss= 0.017 2020 9895  [V]=[UA/giorno]
'           2) moto  geocentrico: k= 0.00 123 944 6678
'
' ------------------- COSTANTI del MOTO ELIOCENTRICO -------------------------
' k = .01720209895#: 'costante dinamica (di Gauss)
' uv0 = 29.78469169#: ' costante per trasform. velocit… da Un.Canon.in [km/s]
' uag0 = uv0 / k:    '(= 1731.456828) per trasform. velocit… da UA/g in [km/s]
' wo$ = "UA"
' ----------------------------------------------------------------------------
'
' ------------------- COSTANTI del MOTO GEOCENTRICO -------------------------
   k = .001239446678#: 'costante dinamica
   uv0 = 7.905364437#: ' costante per trasform. velocit… da Un.Canon.in [km/s]
   uag0 = uv0 / k: '(= 6378.14 = DU) costante per trasform. velocit… in [km/s]
   wo$ = "DU"
' ----------------------------------------------------------------------------
'
' ---------- Dati del problema -------- Sistema GEOcentrico ------
  a = 1.923915122
  e = .0044295
  i = 109.8532 * rad
  peri = 179.4830 * rad
  nodo = 17.7748# * rad
  f = 204.58586# * rad
  u = peri + f: ' argomento della latitudine (u = w + f)
' ---------- FINE .... Dati del problema -------------------------------
'
' ---------- Dati del problema: (1) -------- Sistema Eliocentrico ------
'  a = 1.000001017#
'  e = .01671166845#
'  i = 0 * rad
'  peri = 102.81262678# * rad
'  nodo = 0 * rad
'  f = 264.4766308# * rad
'  u = peri + f: ' argomento della latitudine (u = w + f)
' ---------- FINE .... Dati del problema -------------------------------
'
' ---------- Dati del problema: (2) -------- Sistema Eliocentrico ------
' a = 1.523679341#
' e = .093394934206#
' i = 1.8497638141# * rad
' peri = 286.43497485# * rad
' nodo = 49.50957884# * rad
' f = 246.63325923# * rad
' u = peri + f: ' argomento della latitudine (u = w + f)
' ---------- FINE .... Dati del problema -------------------------------
'
' ......... Calcoli .........
'
  p = a * (1 - e ^ 2)
  r = p / (1 + e * COS(f)): ' raggio vettore

 Rx = COS(u) * COS(nodo) - SIN(u) * SIN(nodo) * COS(i)
 Ry = COS(u) * SIN(nodo) + SIN(u) * COS(nodo) * COS(i)
 Rz = SIN(u) * SIN(i)

 xpri = -SIN(u) * COS(nodo) - COS(u) * SIN(nodo) * COS(i)
 ypri = -SIN(u) * SIN(nodo) + COS(u) * COS(nodo) * COS(i)
 zpri = COS(u) * SIN(i)

'PRINT
'PRINT USING "  r= +#.######### UA"; r
 PRINT
 PRINT USING "  Px= +#.#########   Py= +#.#########   Pz= +#.#########"; Rx; Ry; Rz
 PRINT USING "  Qx= +#.#########   Qy= +#.#########   Qz= +#.#########"; xpri; ypri; zpri
 PRINT
' Vettore posizione
 x = r * Rx: y = r * Ry: z = r * Rz
 V = SQR(Vx ^ 2 + Vy ^ 2 + Vz ^ 2)

' Vettore velocita'
 Vx = k / SQR(p) * (Rx * e * SIN(f) + (1 + e * COS(f)) * xpri)
 Vy = k / SQR(p) * (Ry * e * SIN(f) + (1 + e * COS(f)) * ypri)
 Vz = k / SQR(p) * (Rz * e * SIN(f) + (1 + e * COS(f)) * zpri)
 V = SQR(Vx ^ 2 + Vy ^ 2 + Vz ^ 2)

 PRINT USING "  Vett.Posizione (" + wo$ + "):  x= +#.#########   y= +#.#########   z= +#.#########"; x; y; z
 PRINT USING "                        r= +#.#########"; r
PRINT
 PRINT USING "  Vett.Velocita' (km/s):  Vx= +##.#####   Vy= +##.#####   Vz= +##.#####"; Vx * uag0; Vy * uag0; Vz * uag0
 PRINT USING "                           V= ###.#####"; V * uag0
  PRINT
 r = SQR(x ^ 2 + y ^ 2 + z ^ 2)
 PRINT "   r="; r
  PRINT "                            F I N E"

