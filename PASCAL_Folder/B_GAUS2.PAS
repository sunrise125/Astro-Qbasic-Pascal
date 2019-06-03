Program B_GAUS2;      { Determinazione dei Parametri Orbitali di
                        un Corpo Celeste, note 3 Osservazioni,
                        con il metodo di GAUSS, con correzione dell'
                        orbita preliminare; si differenzia dal Gauss-Merton
                        nel 'refining' degli elementi orbitali }

Uses  Graph,Crt;

Const
         np=3; npp=3;
kgauss=0.01720209895;
rad=pi/180;

Var
           ASAT: Array[1..npp,1..npp+1] of Real;
w,tr,ch  : String;
nme  : String[16];
eps,si,co,
TG1,TG2,TG3,tau1,tau2,tau3,tau13,taux,ff,gg,alfa,
y2t,p,
Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,Vx2,Vy2,Vz2,
XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
r2,V2,u2,p2,q2,sigma2,
x2,y2,z2,
B11,B12,B21,B22,B31,B32,
a,e,Triv,q,Tper,Om,
l1,l2,l3,b1,b2,b3,n1,n2,n3,m1,m2,m3,
      ll1,ll2,ll3,nn1,nn2,nn3,mm1,mm2,mm3,
x1,y1,z1,x3,y3,z3,r1,r3,V1,
R2S,R2CosPsi2,KT,LT,MT,REL2,RQ,RC,
RO2,Xel2,Yel2,Zel2,RO1,RO3,yyy: Real;

Function ArcoTANG(num,den:Real): Real;
Var atan: Real;
Begin
atan:=arctan(num/den);
if den<0 then atan:=atan+pi;
if atan<0 then atan:=atan+2*pi;
ArcoTANG:=atan;
End;

Function ArcoSIN(valore:Real): Real;
Begin
ArcoSIN:=arctan(valore/SQRT(1-valore*valore));
End;

Function ArcoCOS(valore:Real): Real;
Begin
ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
End;

      Function MODULO(x,y:Real) : Real;
   Var modu : Real;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;


Procedure ColoriSchermo;
Begin
textcolor(14);textbackground(1);clrscr;
End;

Procedure Titolo;
Begin;
textcolor(10);
gotoxy(8,2); writeln('CALCOLO DI UN''ORBITA PRELIMINARE ELLITTICA CON IL METODO DI GAUSS');
gotoxy(8,3); writeln('-----------------------------------------------------------------');
textcolor(14);
End;

Procedure InputDatiTastiera;
Begin
{
textcolor(12);
gotoxy(32,4); writeln('Introduzione Dati'); textcolor(14);
gotoxy( 8,5); write('1^ Osservazione: Giorno = '); readln(TG1);
gotoxy(25,6); write('Ascensione Retta [Gradi] = ');
readln(ALFA1);
gotoxy(25,7); write('Declinazione     [Gradi] = ');
readln(DELT1);

gotoxy( 8,8);  write('2^ Osservazione: Giorno = '); readln(TG2);
gotoxy(25,9);  write('Ascensione Retta [Gradi] = ');
readln(ALFA2);
gotoxy(25,10); write('Declinazione     [Gradi] = ');
readln(DELT2);

gotoxy( 8,11);  write('3^ Osservazione: Giorno = ');
readln(TG3);
gotoxy(25,12);  write('Ascensione Retta [Gradi] = ');
readln(ALFA3);
gotoxy(25,13);  write('Declinazione     [Gradi] = ');
readln(DELT3);

gotoxy(8,15);  writeln('Coordinate Equatoriali del SOLE alle Date delle 3 Osservazioni');
gotoxy(13,16); writeln('1^ Data                 2^ Data             3^ Data');
gotoxy(8 ,17); write('XS1 = '); readln(XS1);
gotoxy(8 ,18); write('YS1 = '); readln(YS1);
gotoxy(8 ,19); write('ZS1 = '); readln(ZS1);
gotoxy(32,17); write('XS2 = '); readln(XS2);
gotoxy(32,18); write('YS2 = '); readln(YS2);
gotoxy(32,19); write('ZS2 = '); readln(ZS2);
gotoxy(54,17); write('XS3 = '); readln(XS3);
gotoxy(54,18); write('YS3 = '); readln(YS3);
gotoxy(54,19); write('ZS3 = '); readln(ZS3);
gotoxy(28,20); write('Obliquit… Eclittica î = '); readln(EPS);
                 BAADE.dat
}
   TG1:= 23.309398; TG2:= 28.477754; TG3:= 34.313831;
   ALFA1:= 316.298333; ALFA2:= 322.613458; ALFA3:= 329.633791;
   DELT1:= 15.476055; DELT2:= 12.075694; DELT3:= 8.408166;
   XS1:=-0.860286; YS1:=-0.457649; ZS1:=-0.198547;
   XS2:=-0.810811; YS2:=-0.526063; ZS2:=-0.228209;
   XS3:=-0.747072; YS3:=-0.598100; ZS3:=-0.259468;
   EPS:= 23.4497171;
End;

  Procedure TempiCanoniciGAUSS;
   Begin
      tau1:=(TG1-TG2)*kgauss; tau3:=(TG3-TG2)*kgauss;
      tau13:=tau3-tau1;
   End;

  Procedure InversioneMATRICE3x3;
   Label Skip1,Skip2;
   Var
      i,k,j        : Integer;
      SOMMAT       : Real;
     Begin
          FOR K:= 1 TO NPP DO
             begin
          FOR J:= 1 TO NPP DO
             begin
             IF J <> K THEN ASAT[K, J] := ASAT[K, J] / ASAT[K,K];
             end;
          FOR I:= 1 TO NPP DO
             begin
             IF I= K THEN goto Skip2;
          FOR J:= 1 TO NPP DO
             begin
             IF J = K THEN goto Skip1;
                    ASAT[I, J]:= ASAT[I, J] - ASAT[K, J] * ASAT[I, K];
Skip1: end;
Skip2: end;
          FOR I:= 1 TO NPP DO
             begin
            IF I <> K THEN ASAT[I, K]:= -ASAT[I, K] / ASAT[K, K];
             end;
          ASAT[K, K]:= 1 / ASAT[K, K];
             end;

   End;

  Procedure VettoriL123GAUSS;
   Begin
l1:= COS(Delt1*rad)*COS(Alfa1*rad); m1:= COS(Delt1*rad)*SIN(Alfa1*rad); n1:= SIN(Delt1*rad);
l2:= COS(Delt2*rad)*COS(Alfa2*rad); m2:= COS(Delt2*rad)*SIN(Alfa2*rad); n2:= SIN(Delt2*rad);
l3:= COS(Delt3*rad)*COS(Alfa3*rad); m3:= COS(Delt3*rad)*SIN(Alfa3*rad); n3:= SIN(Delt3*rad);

  {Trasferimento Coefficienti nella Matrice ASAT[3,3] per la sua INVERSIONE }
     ASAT[1, 1]:= l1;  ASAT[1, 2]:= l2;   ASAT[1, 3]:= l3;
     ASAT[2, 1]:= m1;  ASAT[2, 2]:= m2;   ASAT[2, 3]:= m3;
     ASAT[3, 1]:= n1;  ASAT[3, 2]:= n2;   ASAT[3, 3]:= n3;

    InversioneMATRICE3x3;

{ Coefficienti nella Matrice INVERTITA }
       ll1:=ASAT[1, 1];    ll2:=ASAT[1, 2];       ll3:=ASAT[1, 3];
       mm1:=ASAT[2, 1];    mm2:=ASAT[2, 2];       mm3:=ASAT[2, 3];
       nn1:=ASAT[3, 1];    nn2:=ASAT[3, 2];       nn3:=ASAT[3, 3];
   End;


  Procedure LAGUERRE;
   Label Salta1;
   Const niter=5;
   Var
       XX      : Array [1..niter] of Real;
      n,fla,ind,kount     : Integer;
      z,zit,INCR,RAD1,RAD2,RADI,xrad1,xrad2,
      zax,RAX,F0,F1,F2, GW0,GW1,GW2, HW0,HW1,HW2      : Real;

  Function FX(x:Real): Real;
   Begin
    FX:=SQR(SQR(x)*SQR(x))+KT*SQR(x)*SQR(x)*SQR(x)+LT*x*SQR(x)+MT;
  End;

    Function FPriX(x:Real): Real;
   Begin
    FPriX:=8*x*SQR(x)*SQR(x)*SQR(x)+6*KT*x*SQR(x)*SQR(x)+3*LT*SQR(x);
  End;

      Function FSecX(x:Real): Real;
   Begin
    FSecX:=56*SQR(x)*SQR(x)*SQR(x)+30*KT*SQR(x)*SQR(x)+6*LT*x;
  End;

   Begin
  ClrScr;textcolor(10);
  writeln('                      RADICI della Distanza Eliocentrica r2');textcolor(14);
  fla:=0;incr:=0;
  n:= 8;  {Costante iterativa di LAGUERRE; Se n=1 Metodo Laguerre ð N-R}
    {Valore iniziale della 1^RADICE}
           z:= 10; fla:= 1;
       { Calcolo 1^ Radice }
  repeat
    F0:= FX(z); F1:= FPriX(z); F2:= FSecX(z); { f, f' e f" }

  RAX:= SQR(n - 1) * SQR(F1) - n * (n - 1) * F0 * F2;
  if (RAX < 0) then  RAX:=0;
RAD1:= F1 - SQRT(RAX); RAD2:= F1 + SQRT(RAX);
  if (ABS(RAD1) > ABS(RAD2)) THEN
   begin
    RADI:= RAD1;
   end
     ELSE
    begin
    RADI:= RAD2;
   end;
         zit:= z - n * F0 / RADI;     { Iterazione di LAGUERRE }

     zax:= ABS(zit - z);
     z:= zit;
   until zax < 0.000001;
         if (zit<0) then goto Salta1;
 writeln('             r2(',fla:1,')=',zit:16:12,'   Residuo = ',FX(zit));
 XX[fla]:= zit;

            { Calcolo 2^ Radice }
         kount:=0;
          z:=zit+0.5;
          fla:=fla+1;
            xrad1:= zit; { Prima radice }
    repeat
    F0:= FX(z); F1:= FPriX(z); F2:= FSecX(z); { f, f' e f" }
      GW0:= F0/(z-xrad1);
      GW1:=(F1*(z-xrad1)-F0)/SQR(z-xrad1);
      GW2:=(F2*SQR(z-xrad1)-2*F1*(z-xrad1)+2*F0)/((z-xrad1)*SQR(z-xrad1));

  RAX:= SQR(n - 1) * SQR(GW1) - n * (n - 1) * GW0 * GW2;
       if (RAX < 0) then begin
              RAX:= 0; {Se Radicando=NEGATIVO -> Radicando = 0 }
              kount:= kount + 1;
              end;
          IF (kount > 8) then begin
               textcolor(10);
      writeln('                               Radice UNICA ');
          ind:=1;
          goto Salta1;
          end;

RAD1:= GW1 - SQRT(RAX); RAD2:= GW1 + SQRT(RAX);
  if (ABS(RAD1) > ABS(RAD2)) THEN
   begin
    RADI:= RAD1;
   end
     ELSE
    begin
    RADI:= RAD2;
   end;
         zit:= z - n * GW0 / RADI;     { Iterazione di LAGUERRE }

     zax:= ABS(zit - z);
     z:= zit;
   until zax < 0.000001;
         if (zit<0) then goto Salta1;
 writeln('             r2(',fla:1,')=',zit:16:12,'   Residuo = ',FX(zit));
 XX[fla]:= zit;

            { Calcolo 3^ Radice }
          z:=zit+0.5;
          fla:=fla+1;
            xrad2:= zit; { Seconda radice }
    repeat
    F0:= FX(z); F1:= FPriX(z); F2:= FSecX(z); { f, f' e f" }
      GW0:= F0/(z-xrad1);
      GW1:=(F1*(z-xrad1)-F0)/SQR(z-xrad1);
      GW2:=(F2*SQR(z-xrad1)-2*F1*(z-xrad1)+2*F0)/((z-xrad1)*SQR(z-xrad1));

      HW0:= GW0/(z-xrad2);
      HW1:=(GW1*(z-xrad2)-GW0)/SQR(z-xrad2);
      HW2:=(GW2*SQR(z-xrad2)-2*GW1*(z-xrad2)+2*GW0)/((z-xrad2)*SQR(z-xrad2));

  RAX:= SQR(n - 1) * SQR(HW1) - n * (n - 1) * HW0 * HW2;
       if (RAX < 0) then RAX:= 0;
RAD1:= HW1 - SQRT(RAX); RAD2:= HW1 + SQRT(RAX);
  if (ABS(RAD1) > ABS(RAD2)) THEN
   begin
    RADI:= RAD1;
   end
     ELSE
    begin
    RADI:= RAD2;
   end;
         zit:= z - n * HW0 / RADI;     { Iterazione di LAGUERRE }

     zax:= ABS(zit - z);
     z:= zit;
   until zax < 0.000001;
 writeln('             r2(',fla:1,')=',zit:16:12,'   Residuo = ',FX(zit));
 XX[fla]:= zit;

Salta1:
 textcolor(10);
     if fla<=2 then
 writeln('                               RADICE Unica') ELSE
 writeln('                 Scegliere tra 2 RADICI aventi r2 > 1.0 UA ');
 textcolor(14);
 gotoxy(24,wherey);write('Scelta di r2 tramite l''indice: '); read(ind);
  readln;
    REL2:= XX[ind];
    End;

  Procedure ParametriOrbitaCONICA_gauss;
   Var
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       ye1,ze1,ye2,ze2,ye3,ze3,  Vye2,Vze2,
       Incl,Nodo                      : Real;

     Begin
     a:=-1/(V2*V2-2/r2);
 
 { Passaggio dalle Coord. EQUATORIALI a ECLITTICHE
   Le ascisse x2,Vx2 restano INVARIATE }

     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }
  ye2:= y2*co+z2*si; ze2:=-y2*si+z2*co;
   y2:=ye2; z2:=ze2;
  Vye2:= Vy2*co+Vz2*si; Vze2:=-Vy2*si+Vz2*co;
   Vy2:=Vye2; Vz2:=Vze2;

 { Componenti Vettori h,n,e }
  hi:=y2*Vz2-z2*Vy2;   hj:=z2*Vx2-x2*Vz2;   hk:=x2*Vy2-y2*Vx2;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);
  ei:=Vy2*hk-Vz2*hj-x2/r2; ej:=Vz2*hi-Vx2*hk-y2/r2; ek:=Vx2*hj-Vy2*hi-z2/r2;
           e:=SQRT(ei*ei+ej*ej+ek*ek);
      { Parametri Orbitali i,ê,ì }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

     {  Matrice di Trasformazione --> Serve per il GRAFICO dell'Orbita }
 B11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 B12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 B21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 B22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 B31:= SIN(incl) * SIN(om); B32:= SIN(incl) * COS(om);

 textcolor(10);clrscr;
 gotoxy(23,3);writeln('RISULTATI FINALI - Parametri Orbitali'); textcolor(1);
 gotoxy(24,4);writeln(nme+' - ',tr);textcolor(14);
 gotoxy(23,5);writeln('Semiasse magg. a = ',a:10:8,' UA');
 gotoxy(23,6);writeln('Eccentricit…   e = ',e:10:8);
           if a>0 then
             begin
 gotoxy(23,7);writeln('Periodo Rivol. P = ',Exp(1.5*LN(a)):7:5,' Anni Siderali');
             end;
 gotoxy(23,9); writeln('Dist.Perielio  q = ',a*(1-e):10:8,' UA');
 writeln;
 gotoxy(23,11);writeln('Inclinazione   i = ',Incl/rad:12:8,'ø');
 gotoxy(23,12);writeln('Long-Nodo-Asc  ê = ',Nodo/rad:12:8,'ø');
 gotoxy(23,13);writeln('Argom-Perielio w = ',Om/rad:12:8,'ø');
 {repeat until keypressed;}
    End;

  Procedure VARIABILIUNIVERSALI;
  Var
      sg: Integer;
      psi,psia,psir,psi0,trr,b,deltat,
      Sq3,Sq2,Sq1,Sq0,fs,fs1     :Real;

    begin
    deltat:=taux;
     psia:= deltat/r2;  { Valore iniziale di Psi }

   IF (deltat>=0) THEN sg:= 1 ELSE sg:= -1;
                          { sg= Segno di ët per considerare Tempi NEGATIVI
                            Caso di RIVOLUZIONI multiple in Orbite Ellittiche }
 IF (alfa < 0) THEN begin
    psir:= 2*pi/SQRT(-alfa);
    trr:=2*pi/((-alfa)*SQRT(-alfa));
    deltat:= MODULO(deltat, trr);
      psia:= MODULO(psia, psir);
          IF ABS(deltat) >= trr / 2 THEN begin
             deltat:= (ABS(deltat) - trr) * sg;
             psia:= (ABS(psia) - psir) * sg;
          END;
 END;
  repeat
    b:=alfa*psia*psia;
Sq3:= psia*psia*psia/6*((((((((b/342+1)*b/272+1)*b/210+1)*b/156+1)*b/110+1)*b/72+1)*b/42+1)*b/20+1);
Sq2:= psia*psia/2*((((((((b/306+1)*b/240+1)*b/182+1)*b/132+1)*b/90+1)*b/56+1)*b/30+1)*b/12+1);
Sq1:= psia+alfa*Sq3;
Sq0:= 1+alfa*Sq2;

   fs:= r2 * sq1 + sigma2 * sq2 + sq3 - deltat;  { Tempo t }
  fs1:= r2 * sq0 + sigma2 * sq1 + sq2;           { r al tempo t }
      psi:= psia - fs / fs1;   { Iterazione di NEWTON }
   psi0:=ABS(psi-psia);
   psia:=psi;
  until psi0 < 0.00000001;

ff:= 1 - sq2 / r2; gg:= deltat - sq3; { Parametri f,g }

   end;

  Procedure Iterazioner2GAUSS;
   Var
       c1,c3,RO1,RO2,RO3,
       A1,B1,A3,B3,B2ast,A2ast,G1,G2,G3,ff1,gg1,ff3,gg3,aprec,
       c3ro3,c1ro1,R2S,R2CosPsi2: Real;

   Begin
       R2cosPsi2:=-l2*XS2-m2*YS2-n2*ZS2;
       R2S:=SQRT(XS2*XS2+YS2*YS2+ZS2*ZS2);
       A1:=tau3/tau13;
       A3:=-tau1/tau13;
       B1:=tau3/(6*tau13)*(SQR(tau13)-SQR(tau3));
       B3:=-tau1/(6*tau13)*(SQR(tau13)-SQR(tau1));

  A2ast:= mm1*(-XS1*A1+XS2-XS3*A3)+ mm2*(-YS1*A1+YS2-YS3*A3)+ mm3*(-ZS1*A1+ZS2-ZS3*A3);
  B2ast:= mm1*(-XS1*B1    -XS3*B3)+ mm2*(-YS1*B1    -YS3*B3)+ mm3*(-ZS1*B1    -ZS3*B3);

  KT:=-2*A2ast*R2cosPsi2-A2ast*A2ast-R2S*R2S;
  LT:=-2*B2ast*R2cosPsi2-2*A2ast*B2ast;
  MT:=-B2ast*B2ast;
      IF (KT>0) or (LT<0) or (MT>0) THEN
   begin writeln;
  writeln('                 Errore nei Dati INIZIALI ...Ricominciare daccapo !!!');
   writeln('   KT = ',KT:10:8,'    LT = ',LT:10:8,'   MT = ',MT:10:8);
     readln;
     halt;
      end;
       LAGUERRE;
       r2:=REL2;
    C1:= A1+B1/(r2*r2*r2);
    C3:= A3+B3/(r2*r2*r2);
  ff1:=1-tau1*tau1/(2*r2*r2*r2);
  gg1:=tau1-tau1*tau1*tau1/(6*r2*r2*r2);

repeat
   G1:=-XS1*C1+XS2-XS3*C3;
   G2:=-YS1*C1+YS2-YS3*C3;
   G3:=-ZS1*C1+ZS2-ZS3*C3;

    c3ro3:=-nn1*G1-nn2*G2-nn3*G3; RO3:=c3ro3/C3;
    c1ro1:=-ll1*G1-ll2*G2-ll3*G3; RO1:=c1ro1/C1;
    RO2:=mm1*G1+mm2*G2+mm3*G3; {...perchŠ C2=-1 }

     x1:=RO1*l1-XS1;  y1:=RO1*m1-YS1;  z1:=RO1*n1-ZS1;
     x2:=RO2*l2-XS2;  y2:=RO2*m2-YS2;  z2:=RO2*n2-ZS2;
     x3:=RO3*l3-XS3;  y3:=RO3*m3-YS3;  z3:=RO3*n3-ZS3;

      r1:=SQRT(x1*x1+y1*y1+z1*z1);  r3:=SQRT(x3*x3+y3*y3+z3*z3);
      r2:=SQRT(x2*x2+y2*y2+z2*z2);
    Vx2:=(x1-ff1*x2)/gg1;
    Vy2:=(y1-ff1*y2)/gg1;
    Vz2:=(z1-ff1*z2)/gg1;
      V2:=SQRT(Vx2*Vx2+Vy2*Vy2+Vz2*Vz2);
  sigma2:=x2*Vx2+y2*Vy2+z2*Vz2;
       aprec:=a;
       ParametriOrbitaCONICA_gauss;
       alfa:=-1/a;
      taux:=tau3; VARIABILIUNIVERSALI; ff3:=ff; gg3:=gg;
      taux:=tau1; VARIABILIUNIVERSALI; ff1:=ff; gg1:=gg;
    C1:=gg3/(ff1*gg3-ff3*gg1);  C3:=-gg1/(ff1*gg3-ff3*gg1);

{readln;  }
  until ABS(a-aprec)<0.00000001
   End;

 Procedure PassaggioAlPerielio;
  Label Skip1,Skip2;
  Var
      flagg             : Integer;
      psi,psia,psi0,b,tp,rprim,rsec,period,
      S3,S2,S1,S0     :Real;
 Begin
    if e <= 0.9995 then tr:='Traiettoria ELLITTICA'
    ELSE if (e < 1.0005) and (e > 0.9995) then tr:='Traiettoria PARABOLICA'
    ELSE tr:='Traiettoria IPERBOLICA';

      if sigma2<0 then psia:=1
         ELSE psia:=-1;
 Skip1:
  repeat
    b:=alfa*psia*psia;
S3:= psia*psia*psia/6*((((((((b/342+1)*b/272+1)*b/210+1)*b/156+1)*b/110+1)*b/72+1)*b/42+1)*b/20+1);
S2:= psia*psia/2*((((((((b/306+1)*b/240+1)*b/182+1)*b/132+1)*b/90+1)*b/56+1)*b/30+1)*b/12+1);
S1:= psia+alfa*S3;
S0:= 1+alfa*S2;
   rprim:=(1+alfa*r2)*S1+sigma2*S0;      {  r'  }
    rsec:=(1+alfa*r2)*S0+sigma2*alfa*S1; {  r'' }
         psi:= psia - rprim / rsec;   { Iterazione di NEWTON }
   if flagg=1 then goto Skip2;
           if (alfa<0) and (rsec<0) then begin
           flagg:=1;
              if sigma2<0 then psia:=-1 ELSE psia:=1;
           goto Skip1;
           end;
  Skip2:
   psi0:=ABS(psi-psia);
   psia:=psi;
  until psi0 < 0.00000001;

    tp:=r2*S1+sigma2*S2+S3;

   if flagg=1 then begin
     period:=2*pi/(-alfa*SQRT(-alfa));
     tp:=tp-period/2;
   end;

 textcolor(1); gotoxy(24,4);writeln(nme+' - ',tr);textcolor(14);
 gotoxy(23,8);    writeln('Pass. Perielio Tp= ',tp/kgauss:10:6,' giorni da t2');

    End;

Procedure GRAFICOconica;
Label vai;
Var
ch                      : String;
gd,gm,k                 : Integer;
xc,yc,x,y,colore        : Word;
ka,xorbp,yorbp,ro,
xeclp,yeclp,zeclp,rsole : Real;
Begin
gd:=detect;initgraph(gd,gm,'');
if graphresult <> grOk then halt(1);
ka:=30;
vai:
SetFillStyle(1,6);      { Campitura [0..11] e Colore SFONDO }
FloodFill(0,0,12);
SetColor(14);

SetTextStyle(4,0,5);
OutText('            ORBITA');
SetFillStyle(1,12);
xc:=GetmaxX div 2; yc:=GetmaxY div 2;
SetTextStyle(4,0,3); {Font,Direz.,Size}
OutTextXY(xc+250,yc-30,'Punto');OutTextXY(xc+240,yc-8,'Gamma');
SetTextStyle(0,0,1); {Font,Direz.,Size}
OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l''Eclittica; LINEA ROSSA=>Orbita SOTTO l''Eclittica');
OutTextXY(xc+176,yc+200,'SOLE Fuori Scala');
SetLineStyle(0,0,1);
Rectangle(xc-315,yc+160,xc-75,yc+208);
SetColor(10);
OutTextXY(xc-310,yc+165,' + => Ingrandisce     Grafico');
OutTextXY(xc-310,yc+180,' - => Rimpicciolisce  Grafico');
OutTextXY(xc-310,yc+195,'USCITA  con Qualsiasi Freccia');
SetColor(14);
rsole:=0.3*ka;
SetLineStyle(3,0,1);
line(xc,yc,2*xc,yc);
FillEllipse(xc,yc,round(rsole),round(rsole));
{PLOT Orbita}
colore:=14;
for k:=0 to 360 do
Begin
ro:=ka*a*(1-e*e)/(1+e*cos(k*rad));
if (e < 1.0005) and (e > 0.9995) then
begin
ro:=ka*2*q/(1+cos(k*rad));
end;
XORBP:= ro * COS(k*rad); YORBP:= ro * SIN(k*rad);
{Coordinate ORBITALI - Coordinate ECLITTICHE di P }
XECLP:= B11 * XORBP + B12 * YORBP;
YECLP:= B21 * XORBP + B22 * YORBP;
ZECLP:= B31 * XORBP + B32 * YORBP;
x:=xc+round(xeclp);  y:=yc-round(yeclp);
if zeclp>0 then PutPixel(x,y,10) ELSE PutPixel(x,y,12);
end;
repeat
ch:=readkey;
if ch='+' then
begin
ClearDevice;
ka:=ka*2;
goto vai;
end;

if ch='-' then
begin
ClearDevice;
ka:=ka/2;
goto vai;
end;
until keypressed;

closegraph;
restoreCrtMode;
End;

{--------------------------------------------------------------
INIZIO Programma}
Begin
ColoriSchermo;
Titolo;
InputDatiTastiera;
       TempiCanoniciGAUSS;
       VettoriL123GAUSS;
       Iterazioner2GAUSS;
       PassaggioAlPerielio;
       writeln;
       textcolor(12);
       writeln('                      Per il GRAFICO: premere un tasto');
       textcolor(14);
       readln;
       GRAFICOconica;
End.
{FINE Programma}

