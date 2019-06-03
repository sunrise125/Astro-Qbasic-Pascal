  Program LAMB4MET; { PROBLEMA di LAMBERT applicato alla traiettoria di
                     una METEORA: Metodo di BATTIN -> Bettbest.bas: Ago_95 )
   ------------------------------------------------------------------------
   Lamb1met.pas rispetto a Lamb_met.pas ->
                   INPUTati (dist.in km, alt.in km, Azim.in ø, tempo in sec)
                              cioŠ (d1,h1,Az1,t1) per raggio vett. r1
                                   (d2,h2,Az2,t2) per raggio vett. r2
   Conviene adottare: t1=0, cosicchŠ t2= tempo di volo= ç
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
   Lamb4met.pas
 Copia Personalizzata per .......... (Input Manuali)
              Vers.KILLER  up to 2001.Ott.01
   Lamb4met.pas risp. a Lam3met -> Altezza in gradi invece che in Km
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ
   ------------------------------------------------------------------------
   Lamb2met.pas rispetto a Lamb1met.pas ->
         Usato il sistema di coordinate inerziali EQUATORIALE (I,J,K),
       trasformando con l'apposita matrice le coord. orizzontali (S,E,Z)
       di centro Cðcentro Terra in coord. equatoriali terrestri
   ------------------------------------------------------------------------
 .....  Dati in coordinate cartesiane
             Determinazione degli elementi ORBITALI di una CONICA, NOTI:
               1   Vettori posizione r1(x1,y1,z1) e r2(x2,y2,z2)
               2 - Il tempo di Percorrenza ç=t2-t1
                L'angolo di trasferimento é Š calcolato AUTOMATICAMENTE
                tramite le coordinate dei 2 punti }
  Uses Crt,Dos;

  Const
         rad=pi/180;
         mu=1; {unit… di misura CANONICHE }
         DU=6378.14;  {Unit… canonica di lunghezza per campo gravitaz.Terra}
         ameri=6378.14; emeri=0.08181922; {a,e -> ellisse meridiano terrestre}
         Kgauss=0.001239446678; {Cost.gaussiana per centro di gravit… terrestre}
         UV0=DU*Kgauss;
  Var
          fr                  : file of BYTE;
          NVOLTE              : Byte;
          nome                : String[12];
          yea,mon,day,dow,yea1,mon1,day1     : Word;
          nrif,nrif2,io       : LongInt;
        {...... Fine Variabili THRILLING }

         ch: String[1];
         count,ann,mes,gio,ore,min,sec: Integer;

         x1,y1,z1,x2,y2,z2, Vx1,Vy1,Vz1,Vx2,Vy2,Vz2, costeta,aux,
         r1,r2,teta,tau, d1,ha1,Az1, d2,ha2,Az2, hh1,hh2, t1,t2,
         c,s,lam,Tadim, x1km,y1km,z1km, x2km,y2km,z2km, fi,lo,tsid1,tsid2,
         sinf,cosf,tgf2,L,rop,m, x1eq,y1eq,z1eq, x2eq,y2eq,z2eq, RT,
         x,y,eta,csi, alfarad,deltrad, alfavero,deltvero, Om,Incl,Nodo,
         gs4,gs3, gs2, gs1, gs0,
         fs4,fs3, fs2, fs1, fs0,ku,
         hden,h1,h2,B,u,yit,
         a,p,e, V1,V2,
         f,g,fpunto,gpunto        : Real;

   Function ArcoSIN(valore:Real): Real;
     Begin
     ArcoSIN:=arctan(valore/SQRT(1-valore*valore));
     End;

   Function ArcoCOS(valore:Real): Real;
     Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(1-valore*valore));
     End;

   Function ArcoTANG(num,den:Real) : Real;
     Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
   End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
   End;

  Procedure MatriceRotaz(tsid,x,y,z:Real; VAR xeq,yeq,zeq:Real);
   Var
       d11,d12,d13,d21,d22,d23,d31,d32,d33, sf,cf, sts,cts: Real;
   Begin
      { Matrice di Rotazione per passaggio da Coord. (S,E,Z centro C_terra)
        a Coord. Equatoriali (I,J,K stesso centro) -> (xeq,yeq,zeq) }
                  sf:=SIN(fi*rad);     cf:=COS(fi*rad);
                 sts:=SIN(tsid*rad);  cts:=COS(tsid*rad);
   d11:=sf*cts;   d12:=-sts;   d13:=cf*cts;
   d21:=sf*sts;   d22:=cts;    d23:=cf*sts;
   d31:=-cf;      d32:=0;      d33:=sf;
   xeq:=d11*x+d12*y+d13*z;  yeq:=d21*x+d22*y+d23*z;  zeq:=d31*x+d32*y+d33*z;
   End;

  Procedure LogoDatiLAMBERT;
  Var
     u7,jd,jq,t,torig,wn1,mat1,mat2,tsidd,tsloc  :Real;
   Begin
     textcolor(10);
     gotoxy(5,2); writeln('PROBLEMA di LAMBERT per traiettorie di METEORE Terrestri (Metodo BATTIN)');
     gotoxy(5,3); writeln('------------------------------------------------------------------------');

gotoxy(5,4); write('ÚÄÄÄÄÄÄ Data di Rilevamento ÄÄÄÄÄ¿');
gotoxy(5,5); write('³                                ³');
gotoxy(5,6); write('³                                ³');
gotoxy(5,7); write('³                                ³');
gotoxy(5,8); write('³                                ³');
gotoxy(5,9); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(6,5);write('  Giorno = ');readln(gio);
      until (gio<31) and (g>=0);
repeat gotoxy(6,6);write('    Mese = ');readln(mes);
      until mes in [1..12];
       gotoxy(6,7);write('    Anno = ');readln(ann);
repeat gotoxy(23,5);write('  Ore TU = ');readln(ore);
      until ore in [0..23];
repeat gotoxy(23,6);write('  Min    = ');readln(min);
      until min in [0..59];
repeat gotoxy(23,7);write('  Sec    = ');readln(sec);
      until sec in [0..59];
  t1:=0;
      gotoxy(6,8);write('  Durata traccia ç [Sec] = ');readln(t2);

gotoxy(41,4); write('ÚÄÄÄÄÄ Coord.Geogr. POSTAZIONE ÄÄÄÄ¿');
gotoxy(41,5); write('³                                  ³');
gotoxy(41,6); write('³                                  ³');
gotoxy(41,7); write('³                                  ³');
gotoxy(41,8); write('³                                  ³');
gotoxy(41,9); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');

gotoxy(43,6);write('Longitudine Est [Gradi]= ');readln(lo);
gotoxy(43,7);write('Latitudine Nord [Gradi]= ');readln(fi);

gotoxy(20,11);write('ÚÄÄÄ Coord. altazimutali Punto PI ÄÄÄ¿');
gotoxy(20,12);write('³                                    ³');
gotoxy(20,13);write('³                                    ³');
gotoxy(20,14);write('³                                    ³');
gotoxy(20,15);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
 gotoxy(25,12);write('   Distanza in KM = ');readln(d1);
 gotoxy(25,13);write('Altezza  in GRADI = ');readln(ha1);
 gotoxy(25,14);write('  Azimut in GRADI = ');readln(Az1);

gotoxy(20,17);write('ÚÄÄÄ Coord. altazimutali Punto PF ÄÄÄ¿');
gotoxy(20,18);write('³                                    ³');
gotoxy(20,19);write('³                                    ³');
gotoxy(20,20);write('³                                    ³');
gotoxy(20,21);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
 gotoxy(25,18);write('   Distanza in KM = ');readln(d2);
 gotoxy(25,19);write('Altezza  in GRADI = ');readln(ha2);
 gotoxy(25,20);write('  Azimut in GRADI = ');readln(Az2);
     textcolor(15);
 gotoxy(32,22);writeln('Premere un tasto'); read (ch);

 clrscr;
     textcolor(10);
     gotoxy(5,1); writeln('PROBLEMA di LAMBERT per traiettorie di METEORE Terrestri (Metodo BATTIN)');
     gotoxy(5,2); writeln('------------------------------------------------------------------------');

 { Raggio della sfera terrestre locale -> RT, dall'ellissoide di Bessel}
    RT:=ameri*sqrt((1-emeri*emeri)/(1-emeri*emeri*SQR(SIN(fi*rad))));

 {---------------------------------}
   hh1:=ha1*rad;  hh2:=ha2*rad;
 {---------------------------------}
   x1km:=-d1*COS(hh1)*COS(Az1*rad);        x1:=x1km/DU;
   y1km:= d1*COS(hh1)*SIN(Az1*rad);        y1:=y1km/DU;
   z1km:= (d1*SIN(hh1)+RT);                z1:=z1km/DU;
       x2km:=-d2*COS(hh2)*COS(Az2*rad);      x2:=x2km/DU;
       y2km:= d2*COS(hh2)*SIN(Az2*rad);      y2:=y2km/DU;
       z2km:=(d2*SIN(hh2)+RT);               z2:=z2km/DU;
   {Tempo di volo}
   tau:=(t2-t1)*Kgauss;

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
  {Calcolo del T-sid LOCALE all'ora fissata}
   tsloc:= tsidd + lo / 15; {Tsid locale in ORE e fraz. di ora}
   tsloc:= tsloc - INT(tsloc / 24) * 24;  {in ORE}

   MatriceRotaz(15*tsloc,x1,y1,z1,x1eq,y1eq,z1eq);
   MatriceRotaz(15*(tsloc+(t2-t1)/3600),x2,y2,z2,x2eq,y2eq,z2eq);

   { Coordinate EQUATORIALI rinominate come quelle orizzontali per
     non cambiare la subroutine LAMBERT}
     x1:=x1eq;  y1:=y1eq;  z1:=z1eq;
     x2:=x2eq;  y2:=y2eq;  z2:=z2eq;
 {-------------------------------------------------------------------------}
     textcolor(14);
     gotoxy(05,4); write('Vett_Pos  r1  S1 [Km] = '); write(x1km:5:0);
     gotoxy(05,5); write('              E1 [Km] = '); write(y1km:5:0);
     gotoxy(05,6); write('              Z1 [Km] = '); write(z1km:5:0);

     gotoxy(43,4); write('Vett_Pos  r2  S2 [Km] = '); write(x2km:5:0);
     gotoxy(43,5); write('              E2 [Km] = '); write(y2km:5:0);
     gotoxy(43,6); write('              Z2 [Km] = '); write(z2km:5:0);

     gotoxy(25,7); write(' Tempo di volo   ç [sec] = '); write(tau/Kgauss:6:3);
{--------------------------------------------------------------------------------}

   End;

  Procedure LAMBERT;
  Label  Skip10,Skip11,Skip20;
  Var
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       ye1,ze1,ye2,ze2,ye3,ze3,
       sigma2,
       R11,R12,R21,R22,R31,R32,
       pvers,qvers,xvers,yvers,zvers : Real;

   Begin
      r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));  r2:=SQRT(SQR(x2)+SQR(y2)+SQR(z2));
      costeta:=(x1*x2+y1*y2+z1*z2)/(r1*r2);
      aux:=(1-costeta)/2;  aux:=sqrt(aux);
      teta:=2*ArcoSIN(aux);
     gotoxy(25,8); write(' Angolo Trasfer. é [ gradi] = ',teta/rad:9:5);
    gotoxy(25,9); write(' r1= ',r1:11:8,' DU  r2= ',r2:11:8,' DU');

    {  Calcolo PARAMETRI Intermedi }
 c:= SQRT(SQR(r1) + SQR(r2) - 2 * r1 * r2 * COS(teta));
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
  gotoxy(13,3); writeln('  ',count:2,')  yit= ',yit:11:8,'       Residuo = ',yit-y);
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
 R11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 R12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 R21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 R22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 R31:= SIN(incl) * SIN(om);
 R32:= SIN(incl) * COS(om);

 { Versori Perifocali dell'asintoto della traiettoria iperbolica METEORA}
  pvers:=-1/e; qvers:=-sqrt(1-SQR(pvers));

 { Versori Equatoriali dell'asintoto}
  xvers:=R11*pvers+R12*qvers;
  yvers:=R21*pvers+R22*qvers;
  zvers:=R31*pvers+R32*qvers;

{ RADIANTE VERO àV,ëV, cioŠ direzione dell'asintoto dell'iperbole }
    alfavero:=ArcoTANG(yvers,xvers);
    deltvero:=ArcoSIN(zvers);

 {----------- Fine Procedura Lambert --------------------------------------}
    End;

  Procedure StampaRisultati;
   Begin
   {-------------- STAMPA RISULTATI -------------------------}
     textcolor(1); gotoxy(32,10);writeln('RISULTATI FINALI'); textcolor(14);
     writeln('                      Semiasse magg.  a = ',a:11:8,' DU  = ',a*DU:5:0,' Km');
     writeln('                      Eccentricit…    e = ',e:11:8);
     writeln('                      i= ',Incl/rad:7:4,'ø    w= ',Om/rad:8:4,'ø   ê= ',Nodo/rad:8:4,'ø');
  textcolor(12);
     writeln('                      Sistema EQUATORIALE Geocentrico');
  textcolor(15); 
     writeln('POS/VEL_1       r1=',r1:9:6,' DU = ',r1*DU:5:0,' km  V1 = ',V1:6:5,' DU/TU = ',V1*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x1= ',x1:9:6,' DU = ',x1*DU:6:0,' km      V1x= ',Vx1:8:5,' DU/TU = ',Vx1*UV0:7:3,' Km/s');
     writeln('     y1= ',y1:9:6,' DU = ',y1*DU:6:0,' km      V1y= ',Vy1:8:5,' DU/TU = ',Vy1*UV0:7:3,' Km/s');
     writeln('     z1= ',z1:9:6,' DU = ',z1*DU:6:0,' km      V1z= ',Vz1:8:5,' DU/TU = ',Vz1*UV0:7:3,' Km/s');
  textcolor(15); 
     writeln('POS/VEL_2       r2=',r2:9:6,' DU = ',r2*DU:5:0,' km  V2 = ',V2:6:5,' DU/TU = ',V2*UV0:6:3,' Km/s');
  textcolor(14);
     writeln('     x2= ',x2:9:6,' DU = ',x2*DU:6:0,' km      V2x= ',Vx2:8:5,' DU/TU = ',Vx2*UV0:7:3,' Km/s');
     writeln('     y2= ',y2:9:6,' DU = ',y2*DU:6:0,' km      V2y= ',Vy2:8:5,' DU/TU = ',Vy2*UV0:7:3,' Km/s');
     writeln('     z2= ',z2:9:6,' DU = ',z2*DU:6:0,' km      V2z= ',Vz2:8:5,' DU/TU = ',Vz2*UV0:7:3,' Km/s');
  textcolor(10);
     write('     RADIANTE Apparente: à=',alfarad/rad:5:1,'ø   ë=',deltrad/rad:5:1,'ø');
     write('     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
  textcolor(14);

    End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
 {...... Inizio ZONA THRILLING }
   GetDate(yea,mon,day,dow);
   yea1:=yea; mon1:=mon; day1:=day;
   nrif:= trunc(365.25*yea)+trunc(mon*30.6)+day;

   SetDate(2021,10,1);       {DATA Stabilita oltre la quale c'Š il CRASH}
   GetDate(yea,mon,day,dow);
   nrif2:= trunc(365.25*yea)+trunc(mon*30.6)+day;

   SetDate(yea1,mon1,day1);  {Ritorno alla Data di Sistema}

     nome:='LAMB4MET.exe';   {Nome di QUESTO Programma ... da distruggere}
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
       End.
   {FINE Programma}
