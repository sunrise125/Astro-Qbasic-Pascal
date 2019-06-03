  Program METE2ORB; {               Þ The Best Þ
                      Metodo vettoriale per la determinazione delle ORBITE
                      GEOCENTRICA ed ELIOCENTRICA di un METEOROIDE tramite
                      il teorema di LAMBERT (Metodo Battin)
   ------------------------------------------------------------------------
                 Sono NOTE:  le QUOTE d'inizio  e fine traccia,
                 le loro coordinate terrestri e la durata ç del
                 fenomeno -> Lambert Š completamente DEFINITO
   ------------------------------------------------------------------------
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
     METE2ORB.pas   (Nov.00)
              (Input Manuali)
              Vers.KILLER  up to 2001.Nov.01
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ

    Gen.01 ->  Inserite LATITUDINI Geocentriche í', pi— precise rispetto a
                quelle Geografiche í
               Inserito Raggio sfera terrestre secondo ellissoide
                 usato da prof. Ceplecha (esperto mondiale in Meteore)
 }

{$N+}
  Uses Crt,Dos;

  Const
         rad=pi/180;
         mu=1; {unit… di misura CANONICHE }
         DU=6378.14;  {Unit… canonica di lunghezza per campo gravitaz.Terra}
         ameri=6378.14; emeri=0.08181922; {a,e -> ellisse meridiano terrestre}
         Kgaut=0.001239446678; {Cost.gaussiana per centro di gravit… terrestre}
         UV0=DU*Kgaut;
         Kgauss=0.01720209895; uv1=29.784691694;
  Var
          fr                  : file of BYTE;
          NVOLTE              : Byte;
          nome                : String[12];
          yea,mon,day,dow,yea1,mon1,day1     : Word;
          nrif,nrif2,io       : LongInt;
        {...... Fine Variabili THRILLING }

         ddata      : String;
         ch: String[1];
         count,ann,mes,gio,ore,min,sec: Integer;

         x1,y1,z1,x2,y2,z2, Vx1,Vy1,Vz1,Vx2,Vy2,Vz2, costeta,aux,
         r1,r2,teta,tau, d1,ha1,Az1, d2,ha2,Az2, hh1,hh2, t1,t2,
         fi1,fi2,lo1,lo2, Q1,Q2, num1,num2,den1,den2, d1ang,d2ang,
         c,s,lam,Tadim, x1km,y1km,z1km, x2km,y2km,z2km, fi,lo,tsid1,tsid2,
         sinf,cosf,tgf2,L,rop,m, x1eq,y1eq,z1eq, x2eq,y2eq,z2eq, RT,
         x,y,eta,csi, alfarad,deltrad, alfavero,deltvero, Om,Incl,Nodo,
         gs4,gs3, gs2, gs1, gs0,
         fs4,fs3, fs2, fs1, fs0,ku,  jd, Epsdata, T, xvers,yvers,zvers,
         hden,h1,h2,B,u,yit,
         a,p,e, V1,V2,
         f,g,fpunto,gpunto,  teta1,teta2,  qx,qy,qz, lop,fip, tsidd,
    Xdata,Ydata,Zdata, z, re,
    VXdata,VYdata,VZdata, vx,vy,vz, ve, alf, motog, Perio, Tp,
    sigma, ex,ey,ez, hx,hy,hz, nx,ny,nz, h, lnodo, nn, corda,
    ecc,amag,q, sinE,cosE,E2,anoM
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

  Function MODULO(x,y:Double) : Double;
   Var modu : Double;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
   End;

  Procedure LogoDatiLAMBERT;
  Label  Salta;
  Var
     u7,jq,t,torig,wn1,mat1,mat2,tsloc1,tsloc2  :Double;
   Begin
     textcolor(10);
     gotoxy(5,2); writeln('PROBLEMA di LAMBERT per traiettorie di METEORE Terrestri (Metodo BATTIN)');
     gotoxy(5,3); writeln('------------------------------------------------------------------------');
{
gotoxy(22,4); write('ÚÄÄÄÄÄÄ Data di Rilevamento ÄÄÄÄÄ¿');
gotoxy(22,5); write('³                                ³');
gotoxy(22,6); write('³                                ³');
gotoxy(22,7); write('³                                ³');
gotoxy(22,8); write('³                                ³');
gotoxy(22,9); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(23,5);write('  Giorno = ');readln(gio);
      until (gio<31) and (g>=0);
repeat gotoxy(23,6);write('    Mese = ');readln(mes);
      until mes in [1..12];
       gotoxy(23,7);write('    Anno = ');readln(ann);
repeat gotoxy(40,5);write('  Ore TU = ');readln(ore);
      until ore in [0..23];
repeat gotoxy(40,6);write('  Min    = ');readln(min);
      until min in [0..59];
repeat gotoxy(40,7);write('  Sec    = ');readln(sec);
      until sec in [0..59];
  t1:=0;
      gotoxy(23,8);write('  Durata traccia ç [Sec]  = ');readln(t2);

gotoxy(20,11);write('ÚÄÄÄ Coordinate Punto INIZIO scia ÄÄÄ¿');
gotoxy(20,12);write('³                                    ³');
gotoxy(20,13);write('³                                    ³');
gotoxy(20,14);write('³                                    ³');
gotoxy(20,15);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
 gotoxy(25,12);write('      Quota in KM = ');readln(Q1);
 gotoxy(25,13);write('Latitud. in GRADI = ');readln(fi1);
 gotoxy(25,14);write('Long-Est in GRADI = ');readln(lo1);

gotoxy(20,17);write('ÚÄÄÄ Coordinate Punto  FINE  scia ÄÄÄ¿');
gotoxy(20,18);write('³                                    ³');
gotoxy(20,19);write('³                                    ³');
gotoxy(20,20);write('³                                    ³');
gotoxy(20,21);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
 gotoxy(25,18);write('      Quota in KM = ');readln(Q2);
 gotoxy(25,19);write('Latitud. in GRADI = ');readln(fi2);
 gotoxy(25,20);write('Long-Est in GRADI = ');readln(lo2);

     textcolor(15);
 gotoxy(32,22);writeln('Premere un tasto'); read (ch);

 clrscr;
}
{   gio:=10; mes:=8; ann:=1972; ore:=20; min:=28; sec:=29;
   t1:=0; t2:=101.2;
   Q1:=75.98; fi1:=39.289; lo1:=247.425;
   Q2:=101.55;fi2:=52.331; lo2:=246.242;

   gio:=14; mes:=3; ann:=1978; ore:=2; min:=46; sec:=38;
   t1:=0; t2:=1.392;
   Q1:=100;  fi1:=48.6569; lo1:=13.2290;
   Q2:=77.84;fi2:=48.8587; lo2:=13.6245;   }

      {Fireball Borovicka}
   gio:=13; mes:=10; ann:=1990; ore:=3; min:=27; sec:=16;
   t1:=0; t2:=9.8;
   Q1:=103.7;  fi1:=49.2166; lo1:=17.65;
   Q2:=100.4;  fi2:=52.6833; lo2:=17.0666;


 {------------------ INIZIO MODIFICA -----------------}
  {Latitudini geocentriche, usate stesse variab. di latit. geograf.}
  fi1:= fi1-0.1924240867*SIN(2*fi1*rad)+
            0.0003231220*SIN(4*fi1*rad)-
            0.0000007235*SIN(6*fi1*rad);

  fi2:= fi2-0.1924240867*SIN(2*fi2*rad)+
            0.0003231220*SIN(4*fi2*rad)-
            0.0000007235*SIN(6*fi2*rad);

 { Raggio della sfera terrestre locale -> RT, dall'ellissoide di Bessel}
    fi:=(fi1+fi2)/2;  {latit. MEDIA tra le 2 localit…}
{    RT:=ameri*sqrt((1-emeri*emeri)/(1-emeri*emeri*SQR(SIN(fi*rad)))); }

  { ... raggio terrestre secondo Ceplecha, BAIC 30,350 (1979) }
   RT:=sqrt((1-0.0133439554*SQR(SIN(fi*rad)))/(1-0.006694385096*SQR(SIN(fi*rad))));
   RT:=RT*sqrt(40680669.86);

   {Tempo di volo}
   tau:=(t2-t1)*Kgaut;

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

  {Calcolo dei 2 T-sid LOCALI}
 tsloc1:= (tsidd+lo1/15)-INT((tsidd+lo1/15)/24)*24;  {Ore e fraz.NORMALIZZATE}
 tsloc2:= (tsidd+lo2/15)-INT((tsidd+lo2/15)/24)*24;

   teta1:=15*(tsloc1);              {GRADI}
   teta2:=15*tsloc2+(t2-t1)/3600;

   x1eq:=(RT+Q1)/RT*COS(fi1*rad)*COS(teta1*rad);
   y1eq:=(RT+Q1)/RT*COS(fi1*rad)*SIN(teta1*rad);
   z1eq:=(RT+Q1)/RT*SIN(fi1*rad);

      x2eq:=(RT+Q2)/RT*COS(fi2*rad)*COS(teta2*rad);
      y2eq:=(RT+Q2)/RT*COS(fi2*rad)*SIN(teta2*rad);
      z2eq:=(RT+Q2)/RT*SIN(fi2*rad);

 {------------------   FINE MODIFICA -----------------}

   { Coordinate EQUATORIALI rinominate come quelle orizzontali per
     non cambiare la subroutine LAMBERT}
     x1:=x1eq;  y1:=y1eq;  z1:=z1eq;
     x2:=x2eq;  y2:=y2eq;  z2:=z2eq;
   End;

  Procedure LAMBERT;
  Label  Skip10,Skip11,Skip20;
  Var
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       ye1,ze1,ye2,ze2,ye3,ze3,
       sigma2,
       R11,R12,R21,R22,R31,R32,  xs,ys,zs, an,cosa,sina, ra, va, tmenoTp,
       aux1,aux2,csi1,sinh,
       pvers,qvers : Double;

   Begin
      r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));  r2:=SQRT(SQR(x2)+SQR(y2)+SQR(z2));
      costeta:=(x1*x2+y1*y2+z1*z2)/(r1*r2);
      aux:=(1-costeta)/2;  aux:=sqrt(aux);
      teta:=2*ArcoSIN(aux);
    {  Calcolo PARAMETRI Intermedi }
 c:= SQRT(SQR(r1) + SQR(r2) - 2 * r1 * r2 * COS(teta));
 corda:=c;
     s:= (r1 + r2 + c) / 2; { SEMI-PERIMETRO }
   lam:= SQRT(r1 * r2) / s * COS(teta / 2);
 Tadim:= SQRT(8 * mu / s) * tau / s;
{-----------------------------------------------------------------------}
 cosf:= lam * s * 2 / (r1 + r2); sinf:= SQRT(1 - SQR(cosf));
 tgf2:= (1 - cosf) / sinf;
 {Parametri Elle(l), rop, m }
   L:= SQR(tgf2);
   rop:= (r1 + r2 + 2 * lam * s) / 4;
   m:= mu * SQR(tau) / (8 * rop * rop * rop);
{-----------------------------------------------------------------------}
 { Valore iniziale di x}
  x:= L; GOTO Skip11;
Skip10:
    x:= SQRT(SQR((1 - L) / 2) + m / SQR(y)) - (1 + L) / 2; {1^ equaz. di GAUSS}
Skip11:
count:= count + 1;
 { Param. ausiliario eta ->     -1 < eta < +1     formula (52)}
 eta:= x / SQR(1 + SQRT(1 + x));

{Calcolo delle Funzione Ipergeometrica  Csi = f(x,eta)}
 gs4:=0;
 gs3:= 169 / 25 / 27 * eta / (1 + 196 / 27 / 29 * eta / (1 + 225 / 29 / 31 * eta / (1 + gs4)));
 gs2:= 100 / 19 / 21 * eta / (1 + 121 / 21 / 23 * eta / (1 + 144 / 23 / 25 * eta / (1 + gs3)));
 gs1:= 49 / 13 / 15 * eta / (1 + 64 / 15 / 17 * eta / (1 + 81 / 17 / 19 * eta / (1 + gs2)));
 gs0:= 16 / 7 / 9 * eta / (1 + 25 / 9 / 11 * eta / (1 + 36 / 11 / 13 * eta / (1 + gs1)));

 csi:= 8 * (1 + SQRT(1 + x)) / (3 + 1 / (5 + eta + 9 / 7 * eta / (1 + gs0))); {Formula (53)}

{ Parametri h1,h2 }
  hden:= (1 + 2 * x + L) * (4 * x + csi * (3 + x));
  h1:= SQR(L + x) * (1 + 3 * x + csi) / hden; h2:= m * (x - L + csi) / hden;
{ Parametri ausiliari }
 B:= 27 * h2 / (4 * (1 + h1)*SQR(1 + h1));
 u:= -B / 2 / (1 + SQRT(1 + B));

{ Calcolo della Funzione Ipergeometrica K(u) }
  fs3:= 928 / 3591 * u / (1 - 1054 / 4347 * u / (1 - 266 / 1035 * u / (1 - fs4)));
  fs2:= 418 / 1755 * u / (1 - 598 / 2295 * u / (1 - 700 / 2907 * u / (1 - fs3)));
  fs1:= 22 / 81 * u / (1 - 208 / 891 * u / (1 - 340 / 1287 * u / (1 - fs2)));
  fs0:= 4 / 27 * u / (1 - 8 / 27 * u / (1 - 2 / 9 * u / (1 - fs1)));

  ku:= 1 / 3 / (1 - fs0);        {Formula (58)}

{ 2^ Equaz. di Gauss -> Radice positiva y }
 
   yit:= (1 + h1) / 3 * (2 + SQRT(1 + B) / (1 - 2 * u * ku*ku));

  textcolor(11);
{  gotoxy(13,3); writeln('  ',count:2,')  yit= ',yit:11:8,'       Residuo = ',yit-y); }
   IF ABS(yit - y) < 1e-10 THEN goto Skip20
       ELSE y:= yit; GOTO Skip10;
Skip20:
textcolor(14);

{ Calcolo dei PARAMETRI ORBITALI a,p,e }
  a:= m * s * SQR(1 + lam) / (8 * x * SQR(yit));
  p:= 2 * r1 * r2 * SQR((y * (1 + x) * SIN(teta / 2))) / SQR(1 + lam) / m / s;
  e:= SQRT(1 - p / a);

{ Calcolo delle VELOCITA' V1, V2 tramite i parametri lagrangiani f,g,f',g'}

  f:= 1 - r2 / p * 2 * SQR(SIN(teta / 2));
  g:= r1 * r2 * SIN(teta) / SQRT(mu * p);
  fpunto:= (2 / p * SQR(SIN(teta / 2)) - 1 / r1 - 1 / r2) * SIN(teta/2)/COS(teta/2) * SQRT(mu / p);
  gpunto:= 1 - r1 / p * 2 * SQR(SIN(teta / 2));

    { Componenti delle VELOCITA' ai 2 estremi dell'arco }
  Vx1:=(x2-f*x1)/g;   Vy1:=(y2-f*y1)/g;   Vz1:=(z2-f*z1)/g;
 Vx2:=fpunto*x1+gpunto*Vx1; Vy2:=fpunto*y1+gpunto*Vy1; Vz2:=fpunto*z1+gpunto*Vz1;
    { Moduli delle VELOCITA' }
  V1:=SQRT(SQR(Vx1)+SQR(Vy1)+SQR(Vz1));
  V2:=SQRT(SQR(Vx2)+SQR(Vy2)+SQR(Vz2));

{ Calcolo del RADIANTE APPARENTE àApp,ëApp, cioŠ direzione vettore (P1-P2) }
    alfarad:=ArcoTANG(y1-y2,x1-x2);
    deltrad:=ArcoSIN((z1-z2)/sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1)));

     { Prodotto Scalare r2ùV2 = å2 }
      sigma2:= x2*Vx2+y2*Vy2+z2*Vz2;

 { Componenti Vettori h,n,e }
  hi:=y2*Vz2-z2*Vy2;   hj:=z2*Vx2-x2*Vz2;   hk:=x2*Vy2-y2*Vx2;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);
  ei:=Vy2*hk-Vz2*hj-x2/r2; ej:=Vz2*hi-Vx2*hk-y2/r2; ek:=Vx2*hj-Vy2*hi-z2/r2;
            e:=SQRT(ei*ei+ej*ej+ek*ek);

      { Parametri Orbitali i,ê,w }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

     {  Matrice di Trasformazione }
 R11:=  COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 R12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 R21:=  SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 R22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 R31:=  SIN(incl) * SIN(om);
 R32:=  SIN(incl) * COS(om);

 { Versori Perifocali dell'asintoto della traiettoria iperbolica METEORA}
  pvers:=-1/e; qvers:=-sqrt(1-SQR(pvers));

 { Versori Equatoriali dell'asintoto}
  xvers:=R11*pvers+R12*qvers;
  yvers:=R21*pvers+R22*qvers;
  zvers:=R31*pvers+R32*qvers;

{ RADIANTE VERO àV,ëV, cioŠ direzione dell'asintoto dell'iperbole }
    alfavero:=ArcoTANG(yvers,xvers);
    deltvero:=ArcoSIN(zvers);

  { Coordinate geografiche PERIGEO e distanza }
    fip:=ArcoSIN(R31)/rad;
    lop:=MODULO(ArcoTANG(R21,R11)/rad - tsidd*15,360);
      q:=a*(1-e);

{  Versori generici al variare dell'anomalia à e coord. geografiche
   write('  Anomalia in GRADI moto geocentrico = '); readln(an);

   cosa:=COS(an*rad);  sina:=SIN(an*rad);

   xs:= R11*cosa + R12*sina;
   ys:= R21*cosa + R22*sina;
   zs:= R31*cosa + R32*sina;
     ra:=a*(1-e*e)/(1+e*cosa);        Raggio vettore
     va:=sqrt(v1*v1+2*(1/ra-1/r1));   Velocit…
    aux1:=(ra-a)/(-a*e);
    csi1:=LN(aux1+sqrt(aux1*aux1-1));
    aux2:=sqrt(-a*a*a*RT*RT*RT/398600);
    sinh:=(EXP(csi1)-EXP(-csi1))/2;
    tmenoTp:=(e*sinh-csi1)*aux2;       Tempo di percorrenza arco da Tp

    fip:=ArcoSIN(zs)/rad;
    lop:=MODULO(ArcoTANG(ys,xs)/rad - tsidd*15,360);

  writeln('    à=',an:6:2,'   Quota_à=',(ra-1)*RT:6:1,' Km');
  writeln('    í=',fip:8:2,'ø   lon=',lop:8:2,'ø');
  writeln('    V=',va*uv0:7:3,' Km/s');
  writeln('    û(a^3/æ)=',aux2:7:3,' sec       Csi=',csi1:10:8);
  writeln('    t-Tp=',tmenoTp:5:1,' sec');
  writeln('    á=',ArcoCOS(-1/e)/rad:8:2,'ø');
  halt;
}
 {----------- Fine Procedura Lambert --------------------------------------}
    End;

  Procedure StampaRisultati;
   Begin
   {-------------- STAMPA RISULTATI -------------------------}
   textcolor(11); gotoxy(19,1);writeln('O R B I T A   G E O C E N T R I C A'); textcolor(14);
{
     writeln(' Lat_pgeo=',fip:6:2,'ø    Semiasse magg.  a = ',a:11:8,' RT  = ',a*RT:5:0,' Km');
     writeln(' Lon_pgeo=',lop:6:2,'ø    Eccentricit…    e = ',e:11:8,'   RT=',RT:6:1,' Km');
}
     writeln('    Semiasse magg. a=',a:11:8,' RT = ',a*RT:5:0,' Km    q=',q:11:8,' RT =',q*RT:5:0,' Km');
     writeln('    Eccentricit…   e = ',e:11:8,'   RT=',RT:6:1,' Km    Quota Perigeo =',(q-1)*RT:5:0,' Km');
     writeln('                      i= ',Incl/rad:7:4,'ø    w= ',Om/rad:8:4,'ø   ê= ',Nodo/rad:8:4,'ø');
     writeln('                      Angolo di trasf.= ',ArcoCOS(costeta)/rad:7:4,'ø     L_scia=',corda*DU:5:1,' Km');
  textcolor(12);
     writeln('                      Sistema EQUATORIALE Geocentrico');
  textcolor(15); 
     writeln('POS/VEL_1       r1=',r1:9:6,' RT = ',r1*RT:5:0,' km  V1 = ',V1:6:5,' DU/TU = ',V1*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x1= ',x1:9:6,' RT = ',x1*RT:6:0,' km      V1x= ',Vx1:8:5,' DU/TU = ',Vx1*UV0:7:3,' Km/s');
     writeln('     y1= ',y1:9:6,' RT = ',y1*RT:6:0,' km      V1y= ',Vy1:8:5,' DU/TU = ',Vy1*UV0:7:3,' Km/s');
     writeln('     z1= ',z1:9:6,' RT = ',z1*RT:6:0,' km      V1z= ',Vz1:8:5,' DU/TU = ',Vz1*UV0:7:3,' Km/s');
  textcolor(15); 
     writeln('POS/VEL_2       r2=',r2:9:6,' RT = ',r2*RT:5:0,' km  V2 = ',V2:6:5,' DU/TU = ',V2*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x2= ',x2:9:6,' RT = ',x2*RT:6:0,' km      V2x= ',Vx2:8:5,' DU/TU = ',Vx2*UV0:7:3,' Km/s');
     writeln('     y2= ',y2:9:6,' RT = ',y2*RT:6:0,' km      V2y= ',Vy2:8:5,' DU/TU = ',Vy2*UV0:7:3,' Km/s');
     writeln('     z2= ',z2:9:6,' RT = ',z2*RT:6:0,' km      V2z= ',Vz2:8:5,' DU/TU = ',Vz2*UV0:7:3,' Km/s');
  textcolor(10);
     write('     RADIANTE Apparente: à=',alfarad/rad:5:1,'ø   ë=',deltrad/rad:5:1,'ø');
     write('     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
   writeln; writeln;
   textcolor(11); writeln('                  O R B I T A   E L I O C E N T R I C A'); textcolor(14);

    End;

   Procedure Kepler(m,e:Double; VAR v:Double);
   Var anomecc,variaz,tgv: Double;
   Begin
   anomecc:= m;  {Valore iniziale di E assunto pari all'anom.media M}
 repeat
    variaz:= (m + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 UNTIL ABS(variaz) < 1e-10;
   tgv:= SIN(anomecc / 2)/COS(anomecc / 2) * SQRT((1 + e) / (1 - e));
   v:= MODULO(2 * ArcoTANG(tgv, 1), 2 * pi);
   End;

  Procedure PosVelTerra;

 Const  c1=102.937348*rad; c2=1.7195269*rad; c3=4.5962e-4*rad; c33=499e-9*rad;
        c4=357.529101*rad; c5=35999.0502962*rad; c6=-1.5594e-4*rad; c66=-289e-9*rad;
        c7=0.01670862; c8=-4.2037e-5; c9=-1.236e-7; c99=4e-11;
        c16=1.000001018;   {Semiasse maggiore}

 Var Lp,M,e,R,te,p,pp   :Double;
   Begin
         { Vettore Posizione }
    Lp:=MODULO(c1+c2*T+c3*T*T+c33*T*T*T,2*pi);  {Longitudine Perielio}
     M:=MODULO(c4+c5*T+c6*T*T+c66*T*T*T,2*pi);  {Anomalia Media}
     e:=c7+c8*T+c9*T*T+c99*T*T*T;               {Obliquit… Eclittica}

    Kepler(M,e,te);
     p:=c16*(1-e*e); pp:=sqrt(p);
     R:=p/(1+e*COS(te));
   Xdata:=COS(Lp)*R*COS(te)-SIN(Lp)*R*SIN(te);
   Ydata:=SIN(Lp)*R*COS(te)+COS(Lp)*R*SIN(te);

         { Vettore Velocit… }
   VXdata:=COS(Lp)/pp*(-SIN(te))-SIN(Lp)/pp*(e+COS(te)); {in unit… Canoniche}
   VYdata:=SIN(Lp)/pp*(-SIN(te))+COS(Lp)/pp*(e+COS(te));
   VZdata:=0;
    End;

   Procedure ElementiOrbitali;
  Label Skip;
  Var yb,bb,sx,F1,dd: Double;
   Begin
  re:=SQRT(x*x+y*y+z*z);   ve:=SQRT(vx*vx+vy*vy+vz*vz);  {Moduli r,V elioc.}
  v2:=ve*ve;
  sigma:=x*vx+y*vy+z*vz; { Prodotto scalare rùV }

  {Componenti vettore eccentricit…, e}
  ex:= (v2-1/re)*x-sigma*vx;  ey:= (v2-1/re)*y-sigma*vy;
  ez:= (v2-1/re)*z-sigma*vz;

  {Componenti vettore momento angolare, h}
  hx:= y*vz - z*vy;  hy:= z*vx - x*vz;
  hz:= x*vy - y*vx;  h:= sqrt(hx*hx+hy*hy+hz*hz);

  {Componenti vettore longitudine nodo, n}
  nx:= -hy;   ny:= hx;   nz:= 0;     nn:= sqrt(nx*nx+ny*ny);

   ecc:=sqrt(ex*ex+ey*ey+ez*ez);  {eccentricit… orbitale eliocentrica}
   alf:= 2/re - v2;  {costante energetica}
   amag:= 1/alf;     { semiasse maggiore }
      q:= h*h/(1+ecc);  { dist_perielia }
   incl:= ArcoCOS(hz/h);
  lnodo:= ArcoCOS(nx/nn);
  if ny<0 then lnodo:=2*pi-lnodo;
    om:= ArcoCOS((nx*ex+ny*ey)/(nn*ecc));
    if ez<0 then om:=2*pi-om;

   textcolor(12);
   writeln('                      Sistema ECLITTICO Eliocentrico'); textcolor(14);

   if abs(1-ecc)<0.001 then begin            {Orbita PARABOLICA}
   dd:=sigma;
   anoM:=q * dd + dd * dd * dd / 6;
   motog:=Kgauss;
   Tp:= -anoM/motog;
 writeln('          e=',ecc:10:8,'    q=',q:10:8,' UA      V=',ve*uv1:7:3,' km/s');
 writeln('          i=',incl/rad:7:4,'ø      ê=',lnodo/rad:7:4,'ø     w=',om/rad:7:4,'ø');
 writeln('          M=',anoM/rad:7:4,'ø      Tp=',Tp:8:3,' giorni');
   goto skip;
   end;

   if ecc<1 then begin               {Orbita ELLITTICA}
   sinE:=sigma/(ecc*sqrt(amag));
   cosE:=1/ecc*(1-re/amag);
   E2:=ArcoTANG(sinE,cosE);
   anoM:=(E2-ecc*sinE);
   motog:=Kgauss/rad*alf*sqrt(alf);  {moto angolare medio= n}
   Perio:=2*pi/kgauss*amag*sqrt(amag); {P= periodo orbitale in giorni}
 writeln('          e=',ecc:10:8,'    a=',amag:10:8,' UA   n=',motog:10:8,' ø/giorno');
 writeln('          i=',incl/rad:7:4,'ø      ê=',lnodo/rad:7:4,'ø       w=',om/rad:7:4,'ø');
 writeln('          M=',anoM/rad:7:4,'ø      V=',ve*uv1:7:3,' Km/s');
   end;

   if ecc>1 then begin               {Orbita IPERBOLICA}
   yb:= sigma * h/ecc;
   bb:= (-1 / alf) * sqrt(ecc * ecc - 1);
   sx:= yb/bb;
   F1:=LN(sx + sqrt(sx * sx + 1));     {anomalia eccentr. iperbolica}
   anoM:=ecc*sx - F1;
   motog:= Kgauss/rad*(-alf)*sqrt(-alf);
 writeln('          e=',ecc:10:8,'    a=',amag:10:8,' UA    M=',anoM/rad:7:4,'ø      V=',ve*uv1:7:3,' km/s');
 writeln('          i=',incl/rad:8:5,'ø      ê=',lnodo/rad:8:5,'ø         w=',om/rad:8:5,'ø');
   end;
 skip:
   End;

  Procedure OrbitaEliocentrica;
  Var U,Ux,Uy,Uz, Uxecl,Uyecl,Uzecl
                   :Double;
   Begin
    jd:=jd+ore/24+min/1440+sec/86400;  {Giorno Giuliano}
    T:=(jd-2451545)/36525;             {Secolo Giuliano}
   {Calcolo dell'obliquit… eclittica alla data}
   Epsdata:=23.4392911-0.013004166*T-1.64e-7*T*T;

      {Componenti del vettore geocentrico EQUATORIALE Vì ð U in km/s}
            U:=sqrt(V1*V1-2/r1)*UV0;   {Vì in km/s}

      {Sotto: cambiato di segno dei versori perchŠ il verso dell'asintoto,
              che va dal centro dell'iperbole a Pì, Š SEMPRE OPPOSTO al
              verso della velocit… U }

                     Ux:= -U* xvers;
                     Uy:= -U* yvers;
                     Uz:= -U* zvers;
      {Componenti del vettore geocentrico ECLITTICO Vì ð U in km/s}
                     Uxecl:= Ux;
                     Uyecl:= COS(Epsdata*rad)*Uy+SIN(Epsdata*rad)*Uz;
                     Uzecl:=-SIN(Epsdata*rad)*Uy+COS(Epsdata*rad)*Uz;

     PosVelTerra;

 {Componenti del vettore VELOCITA' eliocentrica ECLITTICA Meteora in UA/UT0}
      Vx:= Vxdata + Uxecl/uv1;
      Vy:= Vydata + Uyecl/uv1;
      Vz:= Vzdata + Uzecl/uv1;

 {Componenti del vettore POSIZIONE eliocentrica ECLITTICA Meteora in UA/UT0
                 coincidente con quello dela Terra}
      x:= xdata;   y:= ydata;     z:= zdata;
{
  writeln('         Vettore POSIZIONE e VELOCITA'' TERRA [Unit… Canoniche]');
  writeln(xdata:13:8,ydata:13:8,zdata:13:8,'³',vxdata:12:9,vydata:13:9,vzdata:13:9);
}
  writeln('         Vettore POSIZIONE e VELOCITA'' Eclittica METEORA [Unit… Canoniche]');
  writeln(x:13:8,y:13:8,z:13:8,'³',vx:12:9,vy:13:9,vz:13:9);

    ElementiOrbitali;
    End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
 {...... Inizio ZONA THRILLING }
   GetDate(yea,mon,day,dow);
   yea1:=yea; mon1:=mon; day1:=day;
   nrif:= trunc(365.25*yea)+trunc(mon*30.6)+day;

   SetDate(2001,11,1);       {DATA Stabilita oltre la quale c'Š il CRASH}
   GetDate(yea,mon,day,dow);
   nrif2:= trunc(365.25*yea)+trunc(mon*30.6)+day;

   SetDate(yea1,mon1,day1);  {Ritorno alla Data di Sistema}

      nome:='METE2ORB.exe';   {Nome di QUESTO Programma ... da distruggere}
     NVOLTE:=32;              {Numero byte a caso; 32=carattere spazio}

     ASSIGN(fr,nome);
          Reset(fr);

      if (nrif2-nrif) <=0 then begin
           for io:=1 to 16 do begin    {Primi 16 Bytes ...corrotti}
             seek(fr,io-1);
             write(fr,NVOLTE);
             end;
          halt;
          end;
    Close(fr);
 {...... fine ZONA THRILLING }
       ColoriSchermo;
       LogoDatiLAMBERT;
       LAMBERT;
       StampaRisultati;
       OrbitaEliocentrica;
 END.
   {FINE Programma}
