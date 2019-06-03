  Program MY_DUBY;  { Determinazione traiettoria geocentrica di una METEORA,
                       note, tramite lastre fotografiche, le coordinate
                       (à11,ë11), (à21,ë21), (à12,ë12), (à22,ë22) cioŠ:
            [Staz.1,Punto1, Staz.2,Punto1, Staz.1,Punto2, Staz.2,Punto2]
                       di 4 punti da 2 stazioni terrestri di rilevamento
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ    di date coordinate geografiche (Metodo DUBYAGO).
 Cambiati i COS_DIR    --------------------------------------------------
   a11,a12, etc.          NON Š necessario che i punti (A1,B1) e (A2,B2)
 ERRATI nel Dubyago       siano omologhi - E' qui il VANTAGGIO del metodo
 ÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜÜ    --------------------------------------------------  }

{$N+}
  Uses  Crt;

  Const
         rad=pi/180;
         mu=1; {unit… di misura CANONICHE }
         DU=6378.14;  {Unit… canonica di lunghezza per campo gravitaz.Terra}
         ameri=6378.14; emeri=0.08181922; {a,e -> ellisse meridiano terrestre}
         Kgauss=0.001239446678; {Cost.gaussiana per centro di gravit… terrestre}
         UV0=DU*Kgauss;

  Var
      fl,ft                 : Text;
         w1: String;
         ch: String[1];
         ann,mes,gio,ore,min,sec,count: Integer;
         x21,y21,z21, x22,y22,z22, x11,y11,z11, x12,y12,z12,
         Xg1,Yg1,Zg1, Xg2,Yg2,Zg2, hA1,hA2, hB1,hB2, DistA,DistB,
         D0,D1,D2,D3, DX,DY,DZ,  d11,d21,s21, d22,s22, d12,s12, Contr,
         fi1,lo1,fi2,lo2,
         al11,de11,al12,de12, teta0,
         al21,de21,al22,de22,
         a11,b11,c11, a12,b12,c12,
         a21,b21,c21, a22,b22,c22,
         l1,m1,n1, l2,m2,n2, apri,bpri, denom, a,b,c,
         RT
                         : Double;

         x1,y1,z1,x2,y2,z2, Vx1,Vy1,Vz1,Vx2,Vy2,Vz2, costeta,aux,
         r1,r2,teta,tau,t1,t2,
         s,lam,
         sinf,cosf,tgf2,L,rop,m,
         x,y,eta,csi, alfarad,deltrad, alfavero,deltvero,
         gs4,gs3, gs2, gs1, gs0,
         fs4,fs3, fs2, fs1, fs0,ku,
         hden,h1,h2,u,yit,
         p,e, V1,V2,
         f,g,fpunto,gpunto        : Double;

   Function ArcoTANG(num,den:Double) : Double;
     Var atan : Double;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
   End;

      Function ArcoSIN(valore:Double): Double;
   Begin
     ArcoSIN:=arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

      Function ArcoCOS(valore:Double): Double;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
   End;

  Procedure LogoDati;
  Label Skip77;
  Var
     u7,jd,jq,t,wn1,mat1,mat2,tsidd,
     je,ore1,gio1,ore2,torig, min1,min2
      :Double;
   Begin
   ASSIGN (fl,'DUBY.DAT');
   Reset(fl);
   ASSIGN (ft,'DUBY.RIS');
   Rewrite(ft);
    readln(fl,ann,mes,gio,ore,min,sec);
    readln(fl,fi1,lo1,fi2,lo2);
    readln(fl,t1,t2);
    readln(fl,al11,de11,al21,de21); { Staz.1,Punto1, Staz.2,Punto1 }
    readln(fl,al12,de12,al22,de22); { Staz.1,Punto2, Staz.2,Punto2 }
   Close(fl);

     textcolor(10);
     gotoxy(3,2); writeln('Vettore geocentrico di una METEORA fotografata da 2 stazioni terrestri, note');
     gotoxy(3,3); writeln('le coord. (à,ë) di 2+2 punti NON necessariamente OMOLOGHI  -  Metodo DUBYAGO');
     gotoxy(3,4); writeln('----------------------------------------------------------------------------');
           writeln(ft,'  Vettore geocentrico di una METEORA fotografata da 2 stazioni terrestri, note');
           writeln(ft,'  le coord. (à,ë) di 2+2 punti NON necessariamente OMOLOGHI  -  Metodo DUBYAGO');
           writeln(ft,'  ----------------------------------------------------------------------------');

 { Raggio della sfera terrestre locale -> RT, dall'ellissoide di Bessel}
    RT:=ameri*sqrt((1-emeri*emeri)/(1-emeri*emeri*SQR(SIN((fi1+fi2)/2*rad))));
 writeln(ft,' RT= ',RT:8:3);

  textcolor(12);
 gotoxy(60,4);writeln(' RT=',RT:8:3,' km ');
  textcolor(14);

   {Calcolo di JD}
      u7:= 2-INT(ann/100)+INT(INT(ann/100)/4);
      IF mes < 3 THEN begin
          ann:= ann - 1; mes:= mes + 12;
          end;
      jd:= INT(365.25 * ann) + INT(30.6001 * (mes + 1)) + gio + 1720994.5;
      IF jd >= 2299160.5 THEN jd:= jd + u7;

      je:= jd - 2415020;
      jq:= je + ore/24 + (min + sec/60)/1440; t:= jq / 36525;
      ore1:= ore; min1:= min;
      gio1:= ore1 / 24 + min1 / 1440;
      ore2:= INT(gio1 * 24); min2:= (gio1 * 24 - ore2) * 60;

      IF INT(min2 + 0.5) = 60 THEN begin
              ore2:= ore2 + 1; min2:= 0;
                end;
      IF ore2 < 24 THEN goto Skip77;
      IF( ore2 >= 24) AND (ore2 < 48) THEN begin
         ore2:= ore2 - 24; end
           ELSE ore2:= ore2 - 48;
Skip77:
  { Calcolo del T-sid a Greenwich all'ora fissata }
    wn1:= ore + min/60 + sec/3600; torig:= wn1/876600;
   mat1:= 0.276919398 + 100.0021359 * (t-torig) + 1.075E-06 * SQR(t-torig);
   mat2:= mat1 - INT(mat1);
  tsidd:= 1.002737908 * wn1 + mat2 * 24;
  tsidd:= tsidd - INT(tsidd / 24) * 24;  {in ORE}
   teta0:= tsidd*15;     {Tsid-Greenw. in GRADI}

 {Calcolo coseni direttori aik,bik,cik dei 4 punti
  ------------------------------------------------}
   a11:=COS(de11*rad)*COS(al11*rad);
   b11:=COS(de11*rad)*SIN(al11*rad);
   c11:=SIN(de11*rad);
     a12:=COS(de12*rad)*COS(al12*rad);
     b12:=COS(de12*rad)*SIN(al12*rad);
     c12:=SIN(de12*rad);

   a21:=COS(de21*rad)*COS(al21*rad);
   b21:=COS(de21*rad)*SIN(al21*rad);
   c21:=SIN(de21*rad);
     a22:=COS(de22*rad)*COS(al22*rad);
     b22:=COS(de22*rad)*SIN(al22*rad);
     c22:=SIN(de22*rad);

  { Parametri ausiliari }
    l1:=b11*c12-c11*b12;   m1:=c11*a12-a11*c12;   n1:=a11*b12-b11*a12;
    l2:=b21*c22-c21*b22;   m2:=c21*a22-a21*c22;   n2:=a21*b22-b21*a22;

   apri:=(m1*n2-m2*n1)/(l1*m2-m1*l2);
   bpri:=(n1*l2-n2*l1)/(l1*m2-m1*l2);  denom:=sqrt(SQR(apri)+SQR(bpri)+1);
   a:=apri/denom;    b:=bpri/denom;    c:=   1/denom;
 {-------------------------------------------------------------------------}
{  writeln(' Tsid-Greenw.= ',teta0:9:6,'ø'); }

   { Coordinate geocentriche delle Stazioni 1,2 }
   Xg1:= RT*Cos(fi1*rad)*Cos((lo1+teta0)*rad);     {in Km}
   Yg1:= RT*Cos(fi1*rad)*Sin((lo1+teta0)*rad);     {in Km}
   Zg1:= RT*Sin(fi1*rad);                          {in Km}

     Xg2:= RT*Cos(fi2*rad)*Cos((lo2+teta0)*rad);   {in Km}
     Yg2:= RT*Cos(fi2*rad)*Sin((lo2+teta0)*rad);   {in Km}
     Zg2:= RT*Sin(fi2*rad);                        {in Km}

   DX:=Xg2-Xg1;   DY:=Yg2-Yg1;   DZ:=Zg2-Zg1;

   { Calcolo delle coordinate geocentriche di 2 punti della seconda traccia
     ----------------------------------------------------------------------}
   { 1ø sistema lineare 3x3 }
   D0:= a11*(b21*c-b*c21)+a21*(b*c11-c*b11)+a*(b11*c21-b21*c11);
   D1:= DX*(b21*c-b*c21)+a21*(b*DZ-c*DY)+a*(DY*c21-b21*DZ);
   D2:= a11*(-DY*c+b*DZ)+DX*(-b*c11+c*b11)-a*(b11*DZ-c11*DY);
   D3:= a11*(-DZ*b21+c21*DY)-a21*(DY*c11-DZ*b11)+DX*(-b11*c21+b21*c11);

   d11:=D1/D0;    d21:=D2/D0;     s21:=D3/D0;
 {---------------------------------------------------------------------}
   { 2ø sistema lineare 3x3 }
   D0:= a11*(b22*c-b*c22)+a22*(b*c11-c*b11)+a*(b11*c22-b22*c11);
   D1:= DX*(b22*c-b*c22)+a22*(b*DZ-c*DY)+a*(DY*c22-b22*DZ);
   D2:= a11*(-DY*c+b*DZ)+DX*(-b*c11+c*b11)-a*(b11*DZ-c11*DY);
   D3:= a11*(-DZ*b22+c22*DY)-a22*(DY*c11-DZ*b11)+DX*(-b11*c22+b22*c11);

   d11:=D1/D0;    d22:=D2/D0;     s22:=D3/D0;
 {---------------------------------------------------------------------}
   { Coordinate geocentriche di A1, Staz_1 punto 1 }
   x11:=Xg1+a11*d11;    y11:=Yg1+b11*d11;    z11:=Zg1+c11*d11;

   { Coordinate geocentriche di B1 e B2, Staz_2 punti 1,2 }
   x21:=x11-a*s21;   y21:=y11-b*s21;   z21:=z11-c*s21;
   x22:=x11-a*s22;   y22:=y11-b*s22;   z22:=z11-c*s22;

   { Coordinate geocentriche di A2, Staz_1 punto 2 }
   { 3ø sistema lineare 2x2}
     d12:=d11*(a11*b-a*b11)/(a12*b-a*b12);
     s12:=d11*(a12*b11-a11*b12)/(a12*b-a*b12);
   {la 3^ equaz. come Controllo }
           Contr:=c12*d12-c11*d11+c*s12;
     writeln('         Controllo= ',Contr);
   writeln(ft,'         Controllo= ',Contr);
      x12:=x11-a*s12;    { x_A2 }
      y12:=y11-b*s12;    { y_A2 }
      z12:=z11-c*s12;    { z_A2 }
   writeln;
   writeln(ft);

       { Quote dei Punti }
   hA1:=sqrt(SQR(x11)+SQR(y11)+SQR(z11)) - RT;
   hA2:=sqrt(SQR(x12)+SQR(y12)+SQR(z12)) - RT;
     hB1:=sqrt(SQR(x21)+SQR(y21)+SQR(z21)) - RT;
     hB2:=sqrt(SQR(x22)+SQR(y22)+SQR(z22)) - RT;

       { LUNGHEZZE scie luminose }
   DistA:=sqrt(SQR(x12-x11)+SQR(y12-y11)+SQR(z12-z11));
   DistB:=sqrt(SQR(x22-x21)+SQR(y22-y21)+SQR(z22-z21));

   writeln(' QUOTE met. ->   A1= ',hA1:4:1,' km     A2= ',hA2:4:1,' km');
   writeln('                 B1= ',hB1:4:1,' km     B2= ',hB2:4:1,' km');
   writeln(' LUNGH.scie ->   lA= ',DistA:4:1,' km     lB= ',DistB:4:1,' km');
   writeln;
   writeln(' ... Per punti OMOLOGHI Q_A1=Q_B1, Q_A2=Q_B2 + L_scie Uguali');
   writeln;

   writeln(ft,' QUOTE met. ->   A1= ',hA1:4:1,' km     A2= ',hA2:4:1,' km');
   writeln(ft,'                 B1= ',hB1:4:1,' km     B2= ',hB2:4:1,' km');
   writeln(ft,' LUNGH.scie ->   lA= ',DistA:4:1,' km     lB= ',DistB:4:1,' km');
   writeln(ft);
   writeln(ft,' ... Per punti OMOLOGHI Q_A1=Q_B1, Q_A2=Q_B2 + L_scie Uguali');
   writeln(ft);

   { Coordinate in unit… CANONICHE dei 4 Punti}
   x11:=x11/DU;   y11:=y11/DU;   z11:=z11/DU;  { A1 }
   x12:=x12/DU;   y12:=y12/DU;   z12:=z12/DU;  { A2 }
   x21:=x21/DU;   y21:=y21/DU;   z21:=z21/DU;  { B1 }
   x22:=x22/DU;   y22:=y22/DU;   z22:=z22/DU;  { B2 }

   textcolor(11);
   writeln('        Coordinate geocentriche di A1 e A2 (punti visti dalla Staz_1)');
   textcolor(14);
   writeln(ft,'        Coordinate geocentriche di A1 e A2 (punti visti dalla Staz_1)');
   writeln('           x11=',x11:11:8,'    y11=',y11:11:8,'     z11=',z11:11:8);
   writeln('           x12=',x12:11:8,'    y12=',y12:11:8,'     z12=',z12:11:8);
   writeln(ft,'           x11=',x11:11:8,'    y11=',y11:11:8,'     z11=',z11:11:8);
   writeln(ft,'           x12=',x12:11:8,'    y12=',y12:11:8,'     z12=',z12:11:8);
                        {------------------------------------}
   textcolor(11);
   writeln('        Coordinate geocentriche di B1 e B2 (punti visti dalla Staz_2)');
   textcolor(14);
   writeln(ft,'        Coordinate geocentriche di B1 e B2 (punti visti dalla Staz_2)');
   writeln('           x21=',x21:11:8,'    y21=',y21:11:8,'     z21=',z21:11:8);
   writeln('           x22=',x22:11:8,'    y22=',y22:11:8,'     z22=',z22:11:8);
   writeln(ft,'           x21=',x21:11:8,'    y21=',y21:11:8,'     z21=',z21:11:8);
   writeln(ft,'           x22=',x22:11:8,'    y22=',y22:11:8,'     z22=',z22:11:8);

       {Punti A1 e A2 della 1^ traccia scelti per procedura Lambert }
   w1:='                      Per Lambert: punti 1^ traccia';
   x1:=x11;   y1:=y11;   z1:=z11;
   x2:=x12;   y2:=y12;   z2:=z12;

       {Punti B1 e B2 della 2^ traccia scelti per procedura Lambert
   w1:='                      Per Lambert: punti 2^ traccia';
   x1:=x21;   y1:=y21;   z1:=z21;
   x2:=x22;   y2:=y22;   z2:=z22; }

   textcolor(12);
   writeln;
   writeln(w1);
   writeln(ft);
   writeln(ft,w1);

     textcolor(15);
 writeln;
 writeln('                              Premere un tasto'); read (ch);
 clrscr;
     textcolor(14);

   End;

  Procedure LAMBERT;
  Label  Skip10,Skip11,Skip20;
  Var
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       ye1,ze1,ye2,ze2,ye3,ze3,
       Om,Incl,Nodo,sigma2,
       R11,R12,R21,R22,R31,R32, msambo,lsambo, mlimite, maux,laux,
       pvers,qvers,xvers,yvers,zvers : Double;

   Begin
      r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));  r2:=SQRT(SQR(x2)+SQR(y2)+SQR(z2));
      costeta:=(x1*x2+y1*y2+z1*z2)/(r1*r2);
      aux:=(1-costeta)/2;  aux:=sqrt(aux);
      teta:=2*ArcoSIN(aux);
   tau:=(t2-t1)*Kgauss;

     textcolor(10);
     gotoxy(5,2); writeln('PROBLEMA di LAMBERT per traiettorie di METEORE Terrestri (Metodo BATTIN)');
     gotoxy(5,3); writeln('------------------------------------------------------------------------');
     writeln(ft);
     writeln(ft,'  PROBLEMA di LAMBERT per traiettorie di METEORE Terrestri (Metodo BATTIN)');
     writeln(ft,'  ------------------------------------------------------------------------');
     gotoxy(22,6); write(' Tempo di volo   ç [sec] = '); write(tau/Kgauss:6:3);
     gotoxy(22,7); write(' Angolo Trasfer. é [ gradi] = ',teta/rad:9:5);
     gotoxy(22,8); write(' r1= ',r1:11:8,' DU  r2= ',r2:11:8,' DU');
     writeln(ft,'                        Tempo di volo   ç [sec] = ',tau/Kgauss:6:3);
                   writeln(ft,'                        Angolo Trasfer. é [ gradi] = ',teta/rad:9:5);
                   writeln(ft,'                        r1= ',r1:11:8,' DU  r2= ',r2:11:8,' DU');

  { Controllo se l'orbita Š IPEBOLICA, Sambo: m<4/3ùl^(1.5)+ûl
    ---------------------------------------------------------- }
    maux:=sqrt(2)*sqrt(sqrt(r1*r2))*sqrt(COS(teta/2));
    msambo:=tau/(maux*maux*maux);
      laux:= 4*sqrt(r1*r2)*COS(teta/2);
      lsambo:= (r1+r2)/laux-1/2;
    mlimite:= 4/3*lsambo*sqrt(lsambo)+sqrt(lsambo);

    if msambo >= mlimite then begin
       textcolor(15);textbackground(4);
       writeln;
       writeln('          Attenzione: con i dati immessi la traiettoria NON Š iperbolica       ');
       writeln('          Diminuire il tempo t2, cioŠ la durata in Sec della scia luminosa      ');
       textcolor(14);textbackground(6);
       writeln(ft,'          Attenzione: con i dati immessi la traiettoria NON Š iperbolica');
       writeln(ft,'          Diminuire il tempo t2, cioŠ la durata in Sec della scia luminosa      ');
       close(ft);
    halt;
    end;

    {  Calcolo PARAMETRI Intermedi }
 c:= SQRT(SQR(r1) + SQR(r2) - 2 * r1 * r2 * COS(teta));
     s:= (r1 + r2 + c) / 2; { SEMI-PERIMETRO }
   lam:= SQRT(r1 * r2) / s * COS(teta / 2);
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
  gotoxy(13,4); writeln('  ',count:2,')  yit= ',yit:11:8,'       Residuo = ',yit-y);
                writeln(ft,'        N_iter. ',count:2,')  yit= ',yit:11:8,'       Residuo = ',yit-y);
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
 writeln(ft,'       i= ',Incl/rad:7:4,'   w= ',Om/rad:8:4,'  ê= ',Nodo/rad:8:4);

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
     textcolor(11); gotoxy(32,11);writeln('RISULTATI FINALI'); textcolor(14);
     writeln('                      Semiasse magg.  a = ',a:11:8,' DU  = ',a*DU:5:0,' Km');
     writeln('                      Eccentricit…    e = ',e:11:8);
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
     writeln('     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
  textcolor(14);

     writeln(ft,'                               RISULTATI FINALI');
     writeln(ft,'                      Semiasse magg.  a = ',a:11:8,' DU  = ',a*DU:5:0,' Km');
     writeln(ft,'                      Eccentricit…    e = ',e:11:8);
     writeln(ft,'                      Sistema EQUATORIALE Geocentrico');
     writeln(ft,' POS/VEL_1       r1=',r1:9:6,' DU = ',r1*DU:5:0,' km  V1 = ',V1:6:5,' DU/TU = ',V1*UV0:6:3,' Km/s');
     writeln(ft,'     x1= ',x1:9:6,' DU = ',x1*DU:6:0,' km      V1x= ',Vx1:8:5,' DU/TU = ',Vx1*UV0:7:3,' Km/s');
     writeln(ft,'     y1= ',y1:9:6,' DU = ',y1*DU:6:0,' km      V1y= ',Vy1:8:5,' DU/TU = ',Vy1*UV0:7:3,' Km/s');
     writeln(ft,'     z1= ',z1:9:6,' DU = ',z1*DU:6:0,' km      V1z= ',Vz1:8:5,' DU/TU = ',Vz1*UV0:7:3,' Km/s');
     writeln(ft);
     writeln(ft,' POS/VEL_2       r2=',r2:9:6,' DU = ',r2*DU:5:0,' km  V2 = ',V2:6:5,' DU/TU = ',V2*UV0:6:3,' Km/s');
     writeln(ft,'     x2= ',x2:9:6,' DU = ',x2*DU:6:0,' km      V2x= ',Vx2:8:5,' DU/TU = ',Vx2*UV0:7:3,' Km/s');
     writeln(ft,'     y2= ',y2:9:6,' DU = ',y2*DU:6:0,' km      V2y= ',Vy2:8:5,' DU/TU = ',Vy2*UV0:7:3,' Km/s');
     writeln(ft,'     z2= ',z2:9:6,' DU = ',z2*DU:6:0,' km      V2z= ',Vz2:8:5,' DU/TU = ',Vz2*UV0:7:3,' Km/s');
     writeln(ft);
     write(ft,'     RADIANTE Apparente: à=',alfarad/rad:5:1,'ø   ë=',deltrad/rad:5:1,'ø');
     writeln(ft,'     Vero: à=',alfavero/rad:5:1,'ø   ë=',deltvero/rad:5:1,'ø');
     writeln(ft);
   Close(ft);
  End;

{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDati;
       LAMBERT;
       StampaRisultati;

   End.
   {FINE Programma}
