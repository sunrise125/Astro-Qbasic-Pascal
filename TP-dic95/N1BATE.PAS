  Program N_BATE;
                   { Determinazione dei Parametri Orbitali di un Corpo Celeste
                     con il metodo di LAPLACE-BATE, note 3 Osservazioni.}
  Uses  Graph,Crt;

    Const
         np=3; npp=6;
         kgauss=0.01720209895;
         rad=pi/180;
   Type
        MatrINTER = Array [0..np] of Integer;
         MatrREAL = Array [0..np] of Real;
   Var
           ASAT: Array[1..npp,1..npp+1] of Real;
           B,XP: Array[1..npp] of Real;
           nme    : String[16];
         w,tr,ch  : String;
                                eps,si,co,
                    TG1,TG2,TG3,tau1,tau3,
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
      L1x,L1y,L1z,L2x,L2y,L2z,L3x,L3y,L3z,
      R1x,R1y,R1z,R2x,R2y,R2z,R3x,R3y,R3z,
                    r2,V2,u2,p2,q2,sigma2,
                          f,g,f1,g1,f3,g3,
                     x2,y2,z2,Vx2,Vy2,Vz2,
                  B11,B12,B21,B22,B31,B32,
              a,e,Triv,q,Tper,Om           : Real;


      Function ArcoTANG(num,den:Real) : Real;
   Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

      Function ArcoCOS(valore:Real): Real;
   Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(ABS(1-valore*valore)));
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure LogoDatiLAPL_BATE;
   Begin
     textcolor(10);
     gotoxy(9,2); writeln('CALCOLO DI UN''ORBITA PRELIMINARE CON IL METODO DI LAPLACE-BATE');
     gotoxy(9,3); writeln('--------------------------------------------------------------');
     textcolor(14);
     gotoxy(28,4); writeln('Nome del File-Dati = ');gotoxy(49,4);
     read (w);
     textcolor(10);
     gotoxy(23,5);writeln('RIEPILOGO Dati delle 3 OSSERVAZIONI');textcolor(14);
  End;

  Procedure SLICE;
   Var
      fl                       : Text;
      XS,YS,ZS,ARs,DEs,giorno  : MatrREAL;
      ARh,ARm,DEg,DEp          : MatrINTER;
      segno,i                  : Integer;

    Begin
     ASSIGN(fl,w+'.dat');
          Reset(fl);
     for
     i:=1 to np do
      Begin
       readln(fl,giorno[i]);
       readln(fl,XS[i],YS[i],ZS[i]);
       readln(fl,ARh[i],ARm[i],ARs[i],DEg[i],DEp[i],DEs[i]);
  
      write(giorno[i]:11:5,XS[i]:14:10,YS[i]:14:10,ZS[i]:14:10,'³');
      write(ARh[i]:2,ARm[i]:3,ARs[i]:6:2,' ',DEg[i]:3,DEp[i]:3,DEs[i]:5:1);
      writeln;
     End;
      readln(fl,nme,eps);
      write('   '+nme,' î=',eps:9:6);

      TG1:=giorno[1];TG2:=giorno[2];TG3:=giorno[3];

       XS1:=XS[1];YS1:=YS[1];ZS1:=ZS[1];
       XS2:=XS[2];YS2:=YS[2];ZS2:=ZS[2];
       XS3:=XS[3];YS3:=YS[3];ZS3:=ZS[3];

    Alfa1:=(ARs[1]/3600+ARm[1]/60+ARh[1])*15;
    Alfa2:=(ARs[2]/3600+ARm[2]/60+ARh[2])*15;
    Alfa3:=(ARs[3]/3600+ARm[3]/60+ARh[3])*15;
       if  DEg[1]  < 0 then segno:=-1 else segno:=1;
    Delt1:=segno*(DEs[1]/3600+DEp[1]/60+ABS(DEg[1]));
       if  DEg[2]  < 0 then segno:=-1 else segno:=1;
    Delt2:=segno*(DEs[2]/3600+DEp[2]/60+ABS(DEg[2]));
       if  DEg[3]  < 0 then segno:=-1 else segno:=1;
    Delt3:=segno*(DEs[3]/3600+DEp[3]/60+ABS(DEg[3]));
      readln;readln;
   End;


  Procedure TempiCanoniciLAPL_BATE;
   Begin
      Tau1:=(TG1-TG2)*kgauss;Tau3:=(TG3-TG2)*kgauss;
   End;

  Procedure VettoriL123;
   Begin
L1x:= COS(Delt1*rad)*COS(Alfa1*rad); L1y:= COS(Delt1*rad)*SIN(Alfa1*rad); L1z:= SIN(Delt1*rad);
L2x:= COS(Delt2*rad)*COS(Alfa2*rad); L2y:= COS(Delt2*rad)*SIN(Alfa2*rad); L2z:= SIN(Delt2*rad);
L3x:= COS(Delt3*rad)*COS(Alfa3*rad); L3y:= COS(Delt3*rad)*SIN(Alfa3*rad); L3z:= SIN(Delt3*rad);
   End;

  Procedure VettoriR123;
   Begin
R1x:= -XS1; R1y:= -YS1; R1z:= -ZS1;
R2x:= -XS2; R2y:= -YS2; R2z:= -ZS2;
R3x:= -XS3; R3y:= -YS3; R3z:= -ZS3;
   End;

  Procedure FG(tau:Real);
Var
    FF0,FF1,FF2,FF3,FF4,FF5,FF6,FF7,FF8:  Real;
    GG0,GG1,GG2,GG3,GG4,GG5,GG6,GG7,GG8:  Real;
   Begin
{    8   n
  f = ä (ç /n!ùF )
     n=0        n  }
 FF0:= 1; FF1:= 0; FF2:= -u2; FF3:= 3*u2*p2; FF4:= u2*(u2-15*p2*p2+3*q2);
 FF5:= 15*u2*p2*(7*p2*p2-u2-3*q2);
 FF6:= 105*u2*p2*p2*(-9*p2*p2+6*q2+2*u2)-u2*(45*q2*q2+24*u2*p2+u2*u2);
 FF7:= 315*u2*p2*p2*p2*(33*p2*p2-30*q2-10*u2)+63*u2*p2*(25*q2*q2+14*u2*p2+u2*u2);
 FF8:= 10395*u2*p2*p2*p2*p2*(-13*p2*SQR(SQR(p2))+15*q2+5*u2);
 FF8:= FF8-315*u2*p2*p2*(15*q2+7*u2)*(9*q2+u2);
 FF8:= FF8+u2*(1575*q2*q2*q2+1107*u2*q2*q2+117*u2*u2*q2+u2*u2*u2);

f:= (((((((FF8*tau/40320+FF7/5040)*tau+FF6/720)*tau+FF5/120)*tau+FF4/24)*tau+FF3/6)*tau+FF2/2)*tau+FF1)*tau+FF0;

{     8   n
  g = ä (ç /n!ùG )
     n=0        n   }
 GG0:= 0; GG1:= 1; GG2:= 0; GG3:= -u2; GG4:= 6*u2*p2;
 GG5:= u2*(u2-45*p2*p2+9*q2); GG6:= 30*u2*p2*(14*p2*p2-6*q2-u2);
 GG7:= 315*u2*p2*p2*(-15*p2*p2+10*q2+2*u2)-u2*(225*q2*q2+54*u2*q2+u2*u2);
 GG8:= 630*u2*p2*p2*p2*(99*p2*p2-90*q2-20*u2)+126*u2*p2*(75*q2*q2+24*u2*q2+u2*u2);

g:= (((((((GG8*tau/40320+GG7/5040)*tau+GG6/720)*tau+GG5/120)*tau+GG4/24)*tau+GG3/6)*tau+GG2/2)*tau+GG1)*tau+GG0;

   End;

  Procedure SistemaNxN;
   Label Skip1,Skip2;
   Var
       i,k,j,p,za: Integer;
       SOMMAT    : Real;
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

           FOR ZA:= 1 TO NPP DO
           begin
           SOMMAT:= 0;
           FOR P:= 1 TO NPP DO
           begin
           SOMMAT:= SOMMAT + ASAT[ZA, P] * B[P];
           end;
           XP[ZA]:= SOMMAT;
           end;

   End;

  Procedure Iterazioner2LAPL_BATE;
   Var
       r2it,V2it,r2ab        : Real;
     Begin
      r2:=1.25;
  repeat
      u2:=1/(r2*r2*r2);

      FG(tau1); f1:=f; g1:=g;
      FG(tau3); f3:=f; g3:=g;

 { MATRICE delle incognite 6x6 }
ASAT[1,1]:=f1*L1z; ASAT[1,3]:=-f1*L1x; ASAT[1,4]:=g1*L1z; ASAT[1,6]:=-g1*L1x;
ASAT[2,2]:=f1*L1z; ASAT[2,3]:=-f1*L1y; ASAT[2,5]:=g1*L1z; ASAT[2,6]:=-g1*L1y;
ASAT[3,1]:=L2z; ASAT[3,3]:=-L2x;
ASAT[4,2]:=L2z; ASAT[4,3]:=-L2y;
ASAT[5,1]:=f3*L3z; ASAT[5,3]:=-f3*L3x; ASAT[5,4]:=g3*L3z; ASAT[5,6]:=-g3*L3x;
ASAT[6,2]:=f3*L3z; ASAT[6,3]:=-f3*L3y; ASAT[6,5]:=g3*L3z; ASAT[6,6]:=-g3*L3y;

ASAT[1,2]:=0; ASAT[1,5]:=0; ASAT[2,1]:=0; ASAT[2,4]:=0;
ASAT[3,2]:=0; ASAT[3,4]:=0; ASAT[3,5]:=0; ASAT[3,6]:=0;
ASAT[4,1]:=0; ASAT[4,4]:=0; ASAT[4,5]:=0; ASAT[4,6]:=0;
ASAT[5,2]:=0; ASAT[5,5]:=0; ASAT[6,1]:=0; ASAT[6,4]:=0;

 { VETTORE Matriciale dei 6 TERMINI NOTI }
B[1]:= R1x * L1z - R1z * L1x;
B[2]:= R1y * L1z - R1z * L1y;
B[3]:= R2x * L2z - R2z * L2x;
B[4]:= R2y * L2z - R2z * L2y;
B[5]:= R3x * L3z - R3z * L3x;
B[6]:= R3y * L3z - R3z * L3y;

   SistemaNxN;

 x2:= XP[1]; y2:= XP[2]; z2:= XP[3];
Vx2:= XP[4]; Vy2:= XP[5]; Vz2:= XP[6];

 r2it:= SQRT(x2*x2+y2*y2+z2*z2); V2it:= SQRT(Vx2*Vx2+Vy2*Vy2+Vz2*Vz2);

{ Calcolo dei paramatri p2,q2 per la definizione completa di f,g
  å2 = Prodotto Scalare dei 2 Vettori r2 e V2 = (x2ùVx2+y2ùVy2+z2ùVz2) }
sigma2:= x2 * Vx2 + y2 * Vy2 + z2 * Vz2; p2:= sigma2 / SQR(r2it);
q2:= SQR(V2it/r2it) - u2;
  r2ab:= ABS(r2it-r2); r2:=r2it;

    writeln('                     r2-iter = ',r2it:11:8,'   î = ',r2ab);

  until r2ab < 0.000001;

    textcolor(10); writeln;
    writeln('                                r2-finale = ',r2it:11:8); textcolor(14);

     r2:= r2it;  { Valore FINALE distanza Eliocentrica (Punto 2) }
     V2:= V2it;  { Valore FINALE VELOCITA'Eliocentrica (Punto 2) }
    repeat until keypressed;
   End;

  Procedure ELLISSE;
   Var
      motoang,M2,E2,SinE2,CosE2:  Real;
     Begin
     tr:='Traiettoria ELLITTICA';
     Triv:=a*SQRT(a);
     motoang:=kgauss/Triv;
     SinE2:=sigma2/(e*SQRT(a));
     CosE2:=1/e*(1-r2/a);
     E2:=ArcoTANG(SinE2,CosE2);
     M2:=E2-e*SinE2;
     Tper:=-M2/motoang;
 End;

  Procedure IPERBOLE;
   Var
      Sinhf2,Coshf2,hf2: Real;
     Begin
     tr:='Traiettoria IPERBOLICA';
     Sinhf2:=sigma2/(e*SQRT(-a));
     Coshf2:=1/e*(1-r2/a);
     hf2:=LN(coshf2+SQRT(coshf2*coshf2-1));
     Tper:=-SQRT(-a*a*a)*(e*Sinhf2-hf2)/kgauss;
 End;

  Procedure PARABOLA;
   Var
      d,ppp : Real;
     Begin
     tr:='Traiettoria PARABOLICA';
     d:=sigma2; ppp:=2*r2-d*d;
     Tper:=-0.5*(ppp*d-d*d*d/3);
 End;

  Procedure ParametriOrbita;
   Var
        x,y,z,Vx,Vy,Vz,V,hi,hj,hk,h,r,
        ni,nj,nk,n,ei,ej,ek,Incl,Nodo         : Real;

    Begin
     si:=Sin(eps*rad); co:=Cos(eps*rad);    { Sin(î) & Cos(î) }

     a:=1/(2/r2-V2*V2);    {Semiasse maggiore a }

 { Vettori Posizione e Velocita' in Coord. Eclittiche }
    x:=x2; y:=y2*co+z2*si; z:=-y2*si+z2*co; r:=r2;
   Vx:=Vx2; Vy:=Vy2*co+Vz2*si; Vz:=-Vy2*si+Vz2*co; V:= V2;

 { Componenti Vettori h (momento angolare) ed n(nodo) }
  hi:=y*Vz-z*Vy; hj:=z*Vx-x*Vz; hk:=x*Vy-y*Vx;
  h:= SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; nk:= 0; n:= SQRT(hi*hi+hj*hj);

 { Componenti Vettore e (eccentricit…, diretto secondo la linea perifocale) }
  ei:= (Vy * hk - Vz * hj) - x/r;
  ej:= (Vz * hi - Vx * hk) - y/r;
  ek:= (Vx * hj - Vy * hi) - z/r;
    e:=SQRT(ei*ei+ej*ej+ek*ek);  q:=a*(1-e);  { Eccentr.,Dist.Perielio -> e,q }

 { Calcolo degli elementi -> i,ê,ì }
  Incl:= ArcoCOS(hk/h);
  Nodo:= ArcoCOS(ni/n); if hi < 0 then Nodo:= 2*pi-Nodo;
    Om:= ArcoCOS((ni*ei+nj*ej)/(n*e)); if ek < 0 then Om:= 2*pi-Om;

             if e <= 0.9995 then ELLISSE
                ELSE if (e < 1.0005) and (e > 0.9995) then PARABOLA
                  ELSE IPERBOLE;

       {Matrice di Trasformazione --> Serve per il GRAFICO dell'Orbita}
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
 gotoxy(23,7);writeln('Periodo Rivol. P = ',Triv:7:5,' Anni Siderali');
             end;
 gotoxy(23,8);writeln('Pass.Perielio Tp = ',Tper:12:6,' giorni da t2');
 gotoxy(23,9); writeln('Dist.Perielio  q = ',a*(1-e):10:8,' UA');
 writeln;
 gotoxy(23,11);writeln('Inclinazione   i = ',Incl/rad:12:8,'ø');
 gotoxy(23,12);writeln('Long-Nodo-Asc  ê = ',Nodo/rad:12:8,'ø');
 gotoxy(23,13);writeln('Argom-Perielio ì = ',Om/rad:12:8,'ø');
 readln;  repeat until keypressed;
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

    SetTextStyle(4,0,5);
    OutText('            ORBITA');
     xc:=GetmaxX div 2; yc:=GetmaxY div 2;
     SetTextStyle(4,0,3); {Font,Direz.,Size}
     OutTextXY(xc+250,yc-30,'Punto');OutTextXY(xc+240,yc-8,'Gamma');
     SetTextStyle(0,0,1); {Font,Direz.,Size}
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l''Eclittica;  LINEA ROSSA=>Orbita SOTTO l''Eclittica');
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
      XORBP:= ro * COS(k*rad); YORBP:= ro * SIN(k*rad);{Coordinate ORBITALI
      Coordinate ECLITTICHE di P }
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

{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiLAPL_BATE;
       SLICE;
       TempiCanoniciLAPL_BATE;
       VettoriL123;
       VettoriR123;
       Iterazioner2LAPL_BATE;
       ParametriOrbita;
       GRAFICOconica;
       End.
   {FINE Programma}
