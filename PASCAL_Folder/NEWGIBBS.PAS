  Program NEWGIBBS; { Determinazione della traiettoria (iperbolica) di
                     una METEORA, noti i dati fotografici (à,ë) di 3
                     suoi punti da 2 stazioni terrestri di rilevamento di
                     date coordinate geografiche applicando il Metodo di
                     GIBBS: risoluzione di una conica noti 3 raggi vettori
                     (in questo caso geocentrici)
   ------------------------------------------------------------------------
 }
{$N+}
  Uses  Crt;

  Const
         rad=pi/180;
         mu=1; {unit… di misura CANONICHE }
         DU=6378.14;  {Unit… canonica di lunghezza per campo gravitaz.Terra}
         ameri=6378.388; emeri=0.081991889; {Ellisse di Hayford, Colombo:p.91
                                             f=1/297 -> eý-1=(1-f)ý -> e=emeri}
         Kgauss=0.001239446678; {Cost.gaussiana per centro di gravit… terrestre}
         UV0=DU*Kgauss;  {Per trasf. DU/DT in km/sec ->attraz. Terra}

         {Kgauss=0.01720209895; Cost.gaussiana per centro di gravit… SOLARE}
         {UV0=Kgauss*149597870/86400; Per traf. UA/UT0 in Km/sec gravit.solare}
  Var
      fl,ft                 : Text;
         count,
         ann,mes,gio,ore,min,sec: Integer;
         x1,y1,z1,x2,y2,z2, V1x,V1y,V1z, V2x,V2y,V2z,
         x3,y3,z3, V3x,V3y,V3z, Om,Incl,Nodo, p1,p2,p3, q1,q2,q3,
         r1,r2,r3, Det,
         fi1,lo1,fi2,lo2,
         al1A,de1A,al1B,de1B, teta0, AA,CC,DD,
         al2A,de2A,al2B,de2B,
         al3A,de3A,al3B,de3B,
         ro1,ro2,
         x11, y11, z11, Xg1, Yg1, Zg1,   x21, y21, z21, Xg2, Yg2, Zg2,
         x12, y12, z12, x22, y22, z22,
         x13, y13, z13, x23, y23, z23,
         r1x,r1y,r1z,  r1xver,r1yver,r1zver, r1xdef,r1ydef,r1zdef,
         r2x,r2y,r2z,  r2xver,r2yver,r2zver, r2xdef,r2ydef,r2zdef,
         r3x,r3y,r3z,  r3xver,r3yver,r3zver, r3xdef,r3ydef,r3zdef,
         RT, alfarad,deltrad, alfavero,deltvero,
         a,p,e,q, V1,V2,V3
                         : Double;

   Function ArcoSIN(valore:Double): Double;
     Begin
     ArcoSIN:=arctan(valore/SQRT(1-valore*valore));
     End;

   Function ArcoCOS(valore:Double): Double;
     Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(1-valore*valore));
     End;

   Function ArcoTANG(num,den:Double) : Double;
     Var atan : Double;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
   End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
   End;

  Procedure LogoDatiGIBBS;
  Var
     u7,jd,jq,t,wn1,mat1,mat2,tsidd :Double;
   Begin
   ASSIGN (fl,'NEWGIBBS.DAT');
   Reset(fl);
   ASSIGN (ft,'NEWGIBBS.RIS');
   Rewrite(ft);
    readln(fl,ann,mes,gio,ore,min,sec);
    readln(fl,fi1,lo1,fi2,lo2);
    readln(fl,al1A,de1A,al1B,de1B);
    readln(fl,al2A,de2A,al2B,de2B);
    readln(fl,al3A,de3A,al3B,de3B);
   Close(fl);

     textcolor(10);
     gotoxy(5,2); writeln('Elementi Orbitali di una METEORA fotografata da 2 stazioni terrestri,');
     gotoxy(5,3); writeln('note le coordinate (à,ë) di 3 punti omologhi  -  Metodo astrodinamico');
     gotoxy(5,4); writeln('---------------------------------------------------------------------');
     writeln(ft,'     Elementi Orbitali di una METEORA fotografata da 2 stazioni terrestri,');
     writeln(ft,'     note le coordinate (à,ë) di 3 punti omologhi  -  Metodo astrodinamico');
     writeln(ft,'     ---------------------------------------------------------------------');


 { Raggio della sfera terrestre locale -> RT, dall'ellissoide di Bessel}
    RT:=ameri*sqrt((1-emeri*emeri)/(1-emeri*emeri*SQR(SIN((fi1+fi2)/2*rad))));
 writeln(ft,'   RT= ',RT:8:3);
   {Calcolo di JD}
      u7:= 2-INT(ann/100)+INT(INT(ann/100)/4);
      IF mes < 3 THEN begin
          ann:= ann - 1; mes:= mes + 12;
          end;
      jd:= INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5;
      IF jd >= 2299160.5 THEN jd:= jd + u7;
      jq:= jd - 2415020; t:= jq / 36525;

  {Calcolo del T-sid a Greenwich all'ora fissata}
      wn1:= ore + min/60+sec/3600;
     mat1:= 0.276919398 + 100.0021359 * (t) + 1.075E-06 * SQR(t);
     mat2:= mat1 - INT(mat1);
   tsidd:= 1.002737908 * wn1 + mat2 * 24;
   tsidd:= tsidd - INT(tsidd / 24) * 24;  {in ORE}
   teta0:= tsidd  * 15;  {in GRADI}

 {Calcolo componenti 1ø raggio vettore r1 }
         x11:=COS(de1A*rad)*COS((al1A-teta0)*rad);
         y11:=COS(de1A*rad)*SIN((al1A-teta0)*rad);
         z11:=SIN(de1A*rad);
         Xg1:=RT*COS(fi1*rad)*COS(lo1*rad);
         Yg1:=RT*COS(fi1*rad)*SIN(lo1*rad);
         Zg1:=RT*SIN(fi1*rad);

         x21:=COS(de1B*rad)*COS((al1B-teta0)*rad);
         y21:=COS(de1B*rad)*SIN((al1B-teta0)*rad);
         z21:=SIN(de1B*rad);
         Xg2:=RT*COS(fi2*rad)*COS(lo2*rad);
         Yg2:=RT*COS(fi2*rad)*SIN(lo2*rad);
         Zg2:=RT*SIN(fi2*rad);
       {--------------------------------------------------}
 {Calcolo componenti 2ø raggio vettore r2 }
         x12:=COS(de2A*rad)*COS((al2A-teta0)*rad);
         y12:=COS(de2A*rad)*SIN((al2A-teta0)*rad);
         z12:=SIN(de2A*rad);

         x22:=COS(de2B*rad)*COS((al2B-teta0)*rad);
         y22:=COS(de2B*rad)*SIN((al2B-teta0)*rad);
         z22:=SIN(de2B*rad);
       {--------------------------------------------------}
 {Calcolo componenti 3ø raggio vettore r3 }
         x13:=COS(de3A*rad)*COS((al3A-teta0)*rad);
         y13:=COS(de3A*rad)*SIN((al3A-teta0)*rad);
         z13:=SIN(de3A*rad);

         x23:=COS(de3B*rad)*COS((al3B-teta0)*rad);
         y23:=COS(de3B*rad)*SIN((al3B-teta0)*rad);
         z23:=SIN(de3B*rad);

{------------------- 1ø PUNTO ------------------------------------------------}
    {Parametri AA,CC,DD per il calcolo di (ro1,ro2), che sono le
     distanze topoc. del Punto1 della meteora dalle staz. A e B }

   AA:= -(x11*x21+y11*y21+z11*z21);
   CC:= (Xg2-Xg1)*x11+(Yg2-Yg1)*y11+(Zg2-Zg1)*z11;
   DD:= (Xg2-Xg1)*(-x21)+(Yg2-Yg1)*(-y21)+(Zg2-Zg1)*(-z21);
      ro1:=(CC-DD*AA)/(1-AA*AA);
      ro2:=(DD-CC*AA)/(1-AA*AA);
   r1x:=ro1*x11+Xg1;  r1xver:=ro2*x21+Xg2;
   r1y:=ro1*y11+Yg1;  r1yver:=ro2*y21+Yg2;
   r1z:=ro1*z11+Zg1;  r1zver:=ro2*z21+Zg2;
       r1xdef:=(r1x/ro1+r1xver/ro2)/(1/ro1+1/ro2);
       r1ydef:=(r1y/ro1+r1yver/ro2)/(1/ro1+1/ro2);
       r1zdef:=(r1z/ro1+r1zver/ro2)/(1/ro1+1/ro2);
{------------------- 2ø PUNTO ------------------------------------------------}
   AA:= -(x12*x22+y12*y22+z12*z22);
   CC:= (Xg2-Xg1)*x12+(Yg2-Yg1)*y12+(Zg2-Zg1)*z12;
   DD:= (Xg2-Xg1)*(-x22)+(Yg2-Yg1)*(-y22)+(Zg2-Zg1)*(-z22);
      ro1:=(CC-DD*AA)/(1-AA*AA);
      ro2:=(DD-CC*AA)/(1-AA*AA);
   r2x:=ro1*x12+Xg1;  r2xver:=ro2*x22+Xg2;
   r2y:=ro1*y12+Yg1;  r2yver:=ro2*y22+Yg2;
   r2z:=ro1*z12+Zg1;  r2zver:=ro2*z22+Zg2;
       r2xdef:=(r2x/ro1+r2xver/ro2)/(1/ro1+1/ro2);
       r2ydef:=(r2y/ro1+r2yver/ro2)/(1/ro1+1/ro2);
       r2zdef:=(r2z/ro1+r2zver/ro2)/(1/ro1+1/ro2);
{------------------- 3ø PUNTO ------------------------------------------------}
   AA:= -(x13*x23+y13*y23+z13*z23);
   CC:= (Xg2-Xg1)*x13+(Yg2-Yg1)*y13+(Zg2-Zg1)*z13;
   DD:= (Xg2-Xg1)*(-x23)+(Yg2-Yg1)*(-y23)+(Zg2-Zg1)*(-z23);
      ro1:=(CC-DD*AA)/(1-AA*AA);
      ro2:=(DD-CC*AA)/(1-AA*AA);
   r3x:=ro1*x13+Xg1;  r3xver:=ro2*x23+Xg2;
   r3y:=ro1*y13+Yg1;  r3yver:=ro2*y23+Yg2;
   r3z:=ro1*z13+Zg1;  r3zver:=ro2*z23+Zg2;
       r3xdef:=(r3x/ro1+r3xver/ro2)/(1/ro1+1/ro2);
       r3ydef:=(r3y/ro1+r3yver/ro2)/(1/ro1+1/ro2);
       r3zdef:=(r3z/ro1+r3zver/ro2)/(1/ro1+1/ro2);

 { Raggi geocentrici (r1,r2,r3) in UNITA' CANONICHE: 1ùDU= 6378.14 km }
     x1:=r1xdef/DU;   y1:=r1ydef/DU;   z1:=r1zdef/DU;
     x2:=r2xdef/DU;   y2:=r2ydef/DU;   z2:=r2zdef/DU;
     x3:=r3xdef/DU;   y3:=r3ydef/DU;   z3:=r3zdef/DU;
   End;

  Procedure GIBBS;
  Label  Skip10;
  Var
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       d,d1,d2,d3, nn,nn1,nn2,nn3, s,s1,s2,s3,
       co1,se1,co2,se2,co3,se3,
       sigma2,
       cos12,sin12,cos13,sin13,
       R11,R12,R21,R22,R31,R32,
       pvers,qvers,xvers,yvers,zvers : Double;

   Begin

      r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));  r2:=SQRT(SQR(x2)+SQR(y2)+SQR(z2));
      r3:=SQRT(SQR(x3)+SQR(y3)+SQR(z3));

writeln(ft,'   h1= ',r1*DU-RT:4:1,'   h2= ',r2*DU-RT:4:1,'   h3= ',r3*DU-RT:4:1);

{ Calcolo del RADIANTE APPARENTE àApp,ëApp, cioŠ direzione vettore (P1-P3)}
    alfarad:=ArcoTANG(y1-y3,x1-x3);
    deltrad:=ArcoSIN((z1-z3)/sqrt(SQR(x1-x3)+SQR(y1-y3)+SQR(z1-z3)));

{ Verifica COMPLANARITA' dei Vettori r2,r2,r3 -> Det=0 }
    Det:=x1*(y2*z3-y3*z2)+y1*(x3*z2-z3*x2)+z1*(x2*y3-x3*y2);

 {----------- Inizio Algoritmo di GIBBS --------------------------------------}
      { Calcolo del Vettore D }
  d1:=(y1*z2-z1*y2)+(y2*z3-z2*y3)+(y3*z1-z3*y1);
  d2:=(z1*x2-x1*z2)+(z2*x3-x2*z3)+(z3*x1-x3*z1);
  d3:=(x1*y2-y1*x2)+(x2*y3-y2*x3)+(x3*y1-y3*x1);
       d:= SQRT(d1*d1+d2*d2+d3*d3);

      { Calcolo del Vettore N }
nn1:= r3*(y1*z2-z1*y2)+r1*(y2*z3-z2*y3)+r2*(y3*z1-z3*y1);
nn2:= r3*(z1*x2-x1*z2)+r1*(z2*x3-x2*z3)+r2*(z3*x1-x3*z1);
nn3:= r3*(x1*y2-y1*x2)+r1*(x2*y3-y2*x3)+r2*(x3*y1-y3*x1);
       nn:= SQRT(nn1*nn1+nn2*nn2+nn3*nn3);
  { Parametro p }
     p:=nn/d;
      { Calcolo del Vettore S }
s1:=(r2-r3)*x1+(r3-r1)*x2+(r1-r2)*x3;
s2:=(r2-r3)*y1+(r3-r1)*y2+(r1-r2)*y3;
s3:=(r2-r3)*z1+(r3-r1)*z2+(r1-r2)*z3;
       s:= SQRT(s1*s1+s2*s2+s3*s3);

  { Eccentricit…,Semiasse magg.,Dist.Perielio -> e,a,q }
         e:=s/d; a:=p/(1-e*e); q:=a*(1-e);

{ Calcolo del VERSORE P}
 p1:= (s2 * nn3 - s3 * nn2) / (s * nn);
 p2:= (s3 * nn1 - s1 * nn3) / (s * nn);
 p3:= (s1 * nn2 - s2 * nn1) / (s * nn);

{ Calcolo del VERSORE Q}
 q1:= s1/s; q2:= s2/s; q3:= s3/s;
{----------------------------------------------------------------------------
       ->                ->          ->
       V = û(æ/p)ù[-SINéùP +(e+COSé)ùQ ]
---------------------------------------------------------------------------- }
 { Componenti della VELOCITA' V1 }
 co1:= (x1 * p1 + y1 * p2 + z1 * p3) / r1;  { COS(é1) }
 se1:= SQRT(1 - SQR(co1)); { SIN(é1) }
  V1x:= SQRT(1 / p) * (-se1 * p1 + (e + co1) * q1);
  V1y:= SQRT(1 / p) * (-se1 * p2 + (e + co1) * q2);
  V1z:= SQRT(1 / p) * (-se1 * p3 + (e + co1) * q3);
     V1:= SQRT(V1x*V1x+V1y*V1y+V1z*V1z);

 { Componenti della VELOCITA' V2 }
 co2:= (x2 * p1 + y2 * p2 + z2 * p3) / r2;  { COS(é2) }
 se2:= SQRT(1 - SQR(co2)); { SIN(é2) }
  V2x:= SQRT(1 / p) * (-se2 * p1 + (e + co2) * q1);
  V2y:= SQRT(1 / p) * (-se2 * p2 + (e + co2) * q2);
  V2z:= SQRT(1 / p) * (-se2 * p3 + (e + co2) * q3);
     V2:= SQRT(V2x*V2x+V2y*V2y+V2z*V2z);

 { Componenti della VELOCITA' V3 }
 co3:= (x3 * p1 + y3 * p2 + z3 * p3) / r3;  { COS(é3) }
 se3:= SQRT(1 - SQR(co3)); { SIN(é3) }
  V3x:= SQRT(1 / p) * (-se3 * p1 + (e + co3) * q1);
  V3y:= SQRT(1 / p) * (-se3 * p2 + (e + co3) * q2);
  V3z:= SQRT(1 / p) * (-se3 * p3 + (e + co3) * q3);
     V3:= SQRT(V3x*V3x+V3y*V3y+V3z*V3z);
writeln(ft,'   é1= ',ArcoCOS(co1)/rad:8:4,'   é2= ',ArcoCOS(co2)/rad:8:4,'   é3= ',ArcoCOS(co3)/rad:8:4);
 { Componenti Vettori h,n,e }
  hi:=y2*V2z-z2*V2y;   hj:=z2*V2x-x2*V2z;   hk:=x2*V2y-y2*V2x;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);
  ei:=V2y*hk-V2z*hj-x2/r2; ej:=V2z*hi-V2x*hk-y2/r2; ek:=V2x*hj-V2y*hi-z2/r2;

      { Parametri Orbitali i,ê,ì }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

     {  Matrice di Trasformazione }
 R11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 R12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 R21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 R22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 R31:= SIN(incl) * SIN(om);
 R32:= SIN(incl) * COS(om);


    if e<=1 then goto Skip10;

 { Versori Perifocali dell'asintoto della traiettoria iperbolica METEORA}
  pvers:=-1/e;
  qvers:=-sqrt(1-SQR(pvers));

 { Versori Equatoriali dell'asintoto}
  xvers:=R11*pvers+R12*qvers;
  yvers:=R21*pvers+R22*qvers;
  zvers:=R31*pvers+R32*qvers;

{ RADIANTE VERO àV,ëV, cioŠ direzione dell'asintoto dell'iperbole }
    alfavero:=ArcoTANG(yvers,xvers);
    deltvero:=ArcoSIN(zvers);
 Skip10:

 {----------- Fine Procedura GIBBS --------------------------------------}
    End;

  Procedure StampaRisultati;
  Label  Skip11,Skip20;
   Begin
   {-------------- STAMPA RISULTATI -------------------------}
     textcolor(15); writeln('                DATI (Vett_POS) e RISULTATI (Elem.Conica+Vett_VEL)'); textcolor(14);
     writeln('                      Semiasse magg.  a = ',a:11:8,' DU  = ',a*DU:5:0,' Km');
     writeln('                      Eccentricit…    e = ',e:11:8);
     writeln('          Verifica Complanarit… Vettori Pos. (Det=0) = ',Det);
       textcolor(10);
       if e<=1 then
     writeln('                              ORBITA ELLITTICA ')
       else
     writeln('                              ORBITA IPERBOLICA');
       textcolor(14);
     writeln('                      i= ',Incl/rad:8:5,'ø   w= ',Om/rad:9:5,'ø  ê= ',Nodo/rad:9:5,'ø');
  textcolor(12);
     writeln('                      Sistema EQUATORIALE Geocentrico');
  textcolor(15); 
     writeln('POS/VEL_1    r1=',r1:13:10,' DU = ',r1*DU:5:0,' km  V1 = ',V1:6:5,' DU/TU = ',V1*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x1= ',x1:13:10,' DU = ',x1*DU:6:0,' km      V1x= ',V1x:8:5,' DU/TU = ',V1x*UV0:7:3,' Km/s');
     writeln('     y1= ',y1:13:10,' DU = ',y1*DU:6:0,' km      V1y= ',V1y:8:5,' DU/TU = ',V1y*UV0:7:3,' Km/s');
     writeln('     z1= ',z1:13:10,' DU = ',z1*DU:6:0,' km      V1z= ',V1z:8:5,' DU/TU = ',V1z*UV0:7:3,' Km/s');
  textcolor(15); 
     writeln('POS/VEL_2    r2=',r2:13:10,' DU = ',r2*DU:5:0,' km  V2 = ',V2:6:5,' DU/TU = ',V2*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x2= ',x2:13:10,' DU = ',x2*DU:6:0,' km      V2x= ',V2x:8:5,' DU/TU = ',V2x*UV0:7:3,' Km/s');
     writeln('     y2= ',y2:13:10,' DU = ',y2*DU:6:0,' km      V2y= ',V2y:8:5,' DU/TU = ',V2y*UV0:7:3,' Km/s');
     writeln('     z2= ',z2:13:10,' DU = ',z2*DU:6:0,' km      V2z= ',V2z:8:5,' DU/TU = ',V2z*UV0:7:3,' Km/s');
  textcolor(15); 
     writeln('POS/VEL_3    r3=',r3:13:10,' DU = ',r3*DU:5:0,' km  V3 = ',V3:6:5,' DU/TU = ',V3*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x3= ',x3:13:10,' DU = ',x3*DU:6:0,' km      V3x= ',V3x:8:5,' DU/TU = ',V3x*UV0:7:3,' Km/s');
     writeln('     y3= ',y3:13:10,' DU = ',y3*DU:6:0,' km      V3y= ',V3y:8:5,' DU/TU = ',V3y*UV0:7:3,' Km/s');
     writeln('     z3= ',z3:13:10,' DU = ',z3*DU:6:0,' km      V3z= ',V3z:8:5,' DU/TU = ',V3z*UV0:7:3,' Km/s');
  textcolor(10);
    write('     RADIANTE Apparente: à=',alfarad/rad:5:1,'ø   ë=',deltrad/rad:5:1,'ø');

    if e<=1 then goto Skip11;
     writeln('     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
  Skip11:
  textcolor(14);

     writeln(ft,'           Verifica Complanarit… Vettori Pos. (Det=0) = ',Det);
     writeln(ft,'           DATI (Vett_POS) e RISULTATI (Elem.Conica+Vett_VEL)');
     writeln(ft,'                      Semiasse magg.  a = ',a:11:8,' DU  = ',a*DU:5:0,' Km');
     writeln(ft,'                      Eccentricit…    e = ',e:11:8);
       if e<=1 then
     writeln(ft,'                              ORBITA ELLITTICA ')
       else
     writeln(ft,'                              ORBITA IPERBOLICA');
     writeln(ft,'                      i= ',Incl/rad:8:5,'ø   w= ',Om/rad:9:5,'ø  ê= ',Nodo/rad:9:5,'ø');
     writeln(ft,'                      Sistema EQUATORIALE Geocentrico');
     writeln(ft,' POS/VEL_1   r1=',r1:13:10,' DU = ',r1*DU:5:0,' km  V1 = ',V1:6:5,' DU/TU = ',V1*UV0:6:3,' Km/s');
     writeln(ft,'     x1= ',x1:13:10,' DU = ',x1*DU:6:0,' km      V1x= ',V1x:8:5,' DU/TU = ',V1x*UV0:7:3,' Km/s');
     writeln(ft,'     y1= ',y1:13:10,' DU = ',y1*DU:6:0,' km      V1y= ',V1y:8:5,' DU/TU = ',V1y*UV0:7:3,' Km/s');
     writeln(ft,'     z1= ',z1:13:10,' DU = ',z1*DU:6:0,' km      V1z= ',V1z:8:5,' DU/TU = ',V1z*UV0:7:3,' Km/s');
     writeln(ft);
     writeln(ft,' POS/VEL_2   r2=',r2:13:10,' DU = ',r2*DU:5:0,' km  V2 = ',V2:6:5,' DU/TU = ',V2*UV0:6:3,' Km/s');
     writeln(ft,'     x2= ',x2:13:10,' DU = ',x2*DU:6:0,' km      V2x= ',V2x:8:5,' DU/TU = ',V2x*UV0:7:3,' Km/s');
     writeln(ft,'     y2= ',y2:13:10,' DU = ',y2*DU:6:0,' km      V2y= ',V2y:8:5,' DU/TU = ',V2y*UV0:7:3,' Km/s');
     writeln(ft,'     z2= ',z2:13:10,' DU = ',z2*DU:6:0,' km      V2z= ',V2z:8:5,' DU/TU = ',V2z*UV0:7:3,' Km/s');
     writeln(ft);
     writeln(ft,' POS/VEL_3   r3=',r3:13:10,' DU = ',r3*DU:5:0,' km  V3 = ',V3:6:5,' DU/TU = ',V3*UV0:6:3,' Km/s');
     writeln(ft,'     x3= ',x3:13:10,' DU = ',x3*DU:6:0,' km      V3x= ',V3x:8:5,' DU/TU = ',V3x*UV0:7:3,' Km/s');
     writeln(ft,'     y3= ',y3:13:10,' DU = ',y3*DU:6:0,' km      V3y= ',V3y:8:5,' DU/TU = ',V3y*UV0:7:3,' Km/s');
     writeln(ft,'     z3= ',z3:13:10,' DU = ',z3*DU:6:0,' km      V3z= ',V3z:8:5,' DU/TU = ',V3z*UV0:7:3,' Km/s');
     writeln(ft);
  write(ft,'     RADIANTE Apparente: à=',alfarad/rad:5:1,'ø   ë=',deltrad/rad:5:1,'ø');
    if e<=1 then goto Skip20;
    writeln(ft,'     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
  Skip20:
  writeln(ft);
  writeln(ft,'      [Vettore (P1-P3)]');
   Close(ft);
    End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiGIBBS;
       GIBBS;
       StampaRisultati;
       End.
   {FINE Programma}
