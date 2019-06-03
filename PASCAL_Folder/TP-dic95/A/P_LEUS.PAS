   Program P_LEUS;   { Determinazione dei Parametri Orbitali di un Corpo Celeste
                       con il metodo di LAPLACE-LEUSCHNER, note 3 Osservazioni}

  Uses  Graph,Crt;
  
   Const
         np=3; npp=4;
         kgauss=0.01720209895;
         srad=206265;
         rad=pi/180;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
   Var
      nme      : String[16];
                                      w,tr,ch,
                          rxa1,rxa3,rxb1,rxb3,
      rxc1,rxc2,rxc3,rxc4,rxc5,rxc6,rxc7,rxc8,
      rxc9,rxc10,rxc11,rxc12,rxc13,rxc14,rxc00,
      rxre1,rxre2,rxre3,rxre4,
      rxs1,rxs2,rxs3,rxs4,rxt1,rxt2,rxt3,rxt4 : String;

           ASAT: Array[1..npp,1..npp+1] of Real;
           B,XP: Array[1..npp] of Real;
                                eps,si,co,
                     TG1,TG2,TG3,T1,T2,T3,
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
      LX1,LY1,LZ1,LX2,LY2,LZ2,LX3,LY3,LZ3,
      AX1,AY1,AZ1,AX2,AY2,AZ2,AX3,AY3,AZ3,
      DX1,DY1,DZ1,DX2,DY2,DZ2,DX3,DY3,DZ3,
                  B11,B12,B21,B22,B31,B32,
                     GG1,GG2,GG3,J1,J2,J3,
                  FX1,FY1,FZ1,FX3,FY3,FZ3,
LpriX2,LpriY2,LpriZ2,LsecX2,LsecY2,LsecZ2,
      Xpri2,Ypri2,Zpri2,Xsec2,Ysec2,Zsec2,
    Det,Delta,Delta1,Delta2,Delta3,Delta4,
     R2,R2Q,R2CosPsi2,KT,LT,MT,REL2,RQ,RC,
        ROPRI2,RO2,Xel2,Yel2,Zel2,RO1,RO3,
       XelPri2,YelPri2,ZelPri2,V2,RelPri2,
                     FF1,FF2,FF3,G1,G3,G4,
              ResAL1,ResDE1,ResAL3,ResDE3,
       a,e,Triv,q,Tper,Om,sigma2          : Real;

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

  Procedure Titolo;
   Begin
     textcolor(10);
     gotoxy(7,2); writeln('CALCOLO DI UN''ORBITA PRELIMINARE ELLITTICA CON IL METODO DI LAPLACE');
     gotoxy(7,3); writeln('-------------------------------------------------------------------');
     textcolor(14);
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure Pausa;
   Begin
     textcolor(12);
     writeln;
     write('                               Premere un tasto');
     repeat until keypressed;
     textcolor(14);
  End;

  Procedure LogoDatiLAPLACE;
   Begin
     Titolo;
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
   End;

  Procedure Video1;
   Begin
     gotoxy(10, 4); write('Vettore L: Lx1 = ',LX1:9:6,'   Ly1 = ',LY1:9:6,'   Lz1 = ',LZ1:9:6);
     gotoxy(10, 5); write('           Lx2 = ',LX2:9:6,'   Ly2 = ',LY2:9:6,'   Lz2 = ',LZ2:9:6);
     gotoxy(10, 6); write('           Lx3 = ',LX3:9:6,'   Ly3 = ',LY3:9:6,'   Lz3 = ',LZ3:9:6);

     gotoxy(10, 8); write('Vettore A: Ax1 = ',AX1:9:6,'   Ay1 = ',AY1:9:6,'   Az1 = ',AZ1:9:6);
     gotoxy(10, 9); write('           Ax2 = ',AX2:9:6,'   Ay2 = ',AY2:9:6,'   Az2 = ',AZ2:9:6);
     gotoxy(10,10); write('           Ax3 = ',AX3:9:6,'   Ay3 = ',AY3:9:6,'   Az3 = ',AZ3:9:6);

     gotoxy(10,12); write('Vettore D: Dx1 = ',DX1:9:6,'   Dy1 = ',DY1:9:6,'   Dz1 = ',DZ1:9:6);
     gotoxy(10,13); write('           Dx2 = ',DX2:9:6,'   Dy2 = ',DY2:9:6,'   Dz2 = ',DZ2:9:6);
     gotoxy(10,14); write('           Dx3 = ',DX3:9:6,'   Dy3 = ',DY3:9:6,'   Dz3 = ',DZ3:9:6);

 gotoxy(6,16); write('DIST.ELIOC.-2: R2ý= ',R2Q:8:6,'  R2= ',R2:8:6,'  R2ùCOS(Psi)= ',R2COSPSI2:8:6);
 gotoxy(5,17); write('Del-T CANONICI:  T1 = ',T1:9:6,'    T2 = ',T2:9:6,'    T3 = ',T3:9:6);
 gotoxy(8,19); write('Costanti Gi:  G1= ',GG1:10:6,'    G2= ',GG2:10:6,'    G3= ',GG3:10:6);
 gotoxy(8,20); write('Costanti Ji:  J1= ',J1:10:6,'    J2= ',J2:10:6,'    J3= ',J3:10:6);

 readln; readln;
   End;

  Procedure Video2;
   Begin

 gotoxy(6, 4); write('Vettore L`2 : Lx`2 = ',LPRIX2:9:6,'    Ly`2 = ',LPRIY2:9:6,'    Lz`2 = ',LPRIZ2:9:6);
 gotoxy(6, 5); write('Vettore L``2: Lx``2= ',LSECX2:9:6,'    Ly``2= ',LSECY2:9:6,'    Lz``2= ',LSECZ2:9:6);
 gotoxy(6, 7); write('Vettore R`2 :  X`2 = ',XPRI2 :9:6,'     Y`2 = ',YPRI2 :9:6,'     Z`2 = ',ZPRI2 :9:6);
 gotoxy(6, 8); write('Vettore R``2:  X``2= ',XSEC2 :9:6,'     Y``2= ',YSEC2 :9:6,'     Z``2= ',ZSEC2 :9:6);

 gotoxy(6,10); write('DETERMINANTI:    D = ',DELTA:9:6,'       D1= ',DELTA1:9:6,'       D2= ',DELTA2:9:6);
 gotoxy(6,11); write('CoefPolinom.:    K = ',KT:9:6,'        L= ',LT:9:6,'        M= ',MT:9:6);

 readln;
   End;

  Procedure Video4;
   Begin
 gotoxy(20,4); write('Radice Scelta  r2 = ',REL2:11:8,' UA');
 gotoxy(6,5); write('DETERMINANTI:    D = ',DELTA:9:6,'       D3= ',DELTA3:9:6,'       D4= ',DELTA4:9:6);
 gotoxy(6,6); write('R/Vel.GEOCEN:  ro2 = ',RO2:11:8,' UA     ro`2 = ',ROPRI2:11:8,' UA/UT0');
 gotoxy(3,8); write('(x,y,z)-2 ELIOC:   x2 = ',XEL2:11:8,'    y2 = ',YEL2:11:8,'    z2 = ',ZEL2:11:8);
 gotoxy(3,9); write('(Vx,Vy,Vz)-2 EL:  x`2 = ',XELPRI2:11:8,'   y`2 = ',YELPRI2:11:8,'   z`2 = ',ZELPRI2:11:8);
 gotoxy(16,10); write('R2-V2 ELIOC:  r2 = ',REL2:11:8,'    V2 = ',V2:11:8);
gotoxy(23,11); write('VEL/RADIALE/2  r`2 = ',RELPRI2:11:8);
 gotoxy(3,12); write('DIST/GEOC-1/2/3:  ro1 = ',RO1:11:8,'   ro2 = ',RO2:11:8,'   ro3 = ',RO3:11:8);

 writeln;writeln;
   End;

  Procedure Video5;
   Begin
 gotoxy(12, 4); write('Parametri  f1 = ',FF1:9:6,'       g1 = ',G1:9:6);
 gotoxy(12, 5); write('           f3 = ',FF3:9:6,'       g3 = ',G3:9:6);
 gotoxy(12, 6); write('Grandezze  a1 = ',rxa1,'     b1 = ',rxb1);
 gotoxy(12, 7); write('           a3 = ',rxa3,'     b3 = ',rxb3);
 gotoxy(12, 8); write('   fx1 = ',FX1:9:6,'    fy1 = ',FY1:9:6,'    fz1 = ',FZ1:9:6);
 gotoxy(12, 9); write('   fx3 = ',FX3:9:6,'    fy3 = ',FY3:9:6,'    fz3 = ',FZ3:9:6);
 gotoxy(16, 11); write('Matrice di LEUSCHNER (4x4) estesa ai termini Noti');
 gotoxy(6, 12); write('³ ',rxc1,'    ',rxc2,'    ',rxc3,'    ',rxc00,' ³ ',rxre1);
 gotoxy(6, 13); write('³ ',rxc4,'    ',rxc5,'    ',rxc6,'    ',rxc00,' ³ ',rxre2);
 gotoxy(6, 14); write('³ ',rxc7,'    ',rxc8,'    ',rxc9,'    ',rxc10,' ³ ',rxre3);
 gotoxy(6, 15); write('³ ',rxc11,'    ',rxc12,'    ',rxc13,'    ',rxc14,' ³ ',rxre4);
 gotoxy(18, 17); write('RISULTATI del Sistema Lineare e loro VERIFICA');
 gotoxy(15, 18); write('ëro2 = ',rxs1,'    Residuo 1^Equaz. = ',rxt1);
 gotoxy(15, 19); write('ëx`2 = ',rxs2,'    Residuo 2^Equaz. = ',rxt2);
 gotoxy(15, 20); write('ëy`2 = ',rxs3,'    Residuo 3^Equaz. = ',rxt3);
 gotoxy(15, 21); write('ëz`2 = ',rxs4,'    Residuo 4^Equaz. = ',rxt4);

 readln;
   End;

  Procedure TempiCanoniciLAPLACE;
   Begin
      T1:=(TG1-TG2)*kgauss; T2:=(TG3-TG1)*kgauss; T3:=(TG3-TG2)*kgauss;
   End;

  Procedure VettoriL123;
   Begin
LX1:= COS(Delt1*rad)*COS(Alfa1*rad); LY1:= COS(Delt1*rad)*SIN(Alfa1*rad); LZ1:= SIN(Delt1*rad);
LX2:= COS(Delt2*rad)*COS(Alfa2*rad); LY2:= COS(Delt2*rad)*SIN(Alfa2*rad); LZ2:= SIN(Delt2*rad);
LX3:= COS(Delt3*rad)*COS(Alfa3*rad); LY3:= COS(Delt3*rad)*SIN(Alfa3*rad); LZ3:= SIN(Delt3*rad);
   End;

  Procedure VettoriA123;
   Begin
AX1:=-SIN(Alfa1*rad); AY1:= COS(Alfa1*rad); AZ1:= 0;
AX2:=-SIN(Alfa2*rad); AY2:= COS(Alfa2*rad); AZ2:= 0;
AX3:=-SIN(Alfa3*rad); AY3:= COS(Alfa3*rad); AZ3:= 0;
   End;

  Procedure VettoriD123;
   Begin
DX1:= -SIN(Delt1*rad)*COS(Alfa1*rad); DY1:=-SIN(Delt1*rad)*SIN(Alfa1*rad); DZ1:= COS(Delt1*rad);
DX2:= -SIN(Delt2*rad)*COS(Alfa2*rad); DY2:=-SIN(Delt2*rad)*SIN(Alfa2*rad); DZ2:= COS(Delt2*rad);
DX3:= -SIN(Delt3*rad)*COS(Alfa3*rad); DY3:=-SIN(Delt3*rad)*SIN(Alfa3*rad); DZ3:= COS(Delt3*rad);
   End;

  Procedure CostantiGJ;
   Begin
    GG1:=-T3/(T1*T2);GG3:=-T1/(T3*T2);GG2:=GG1-GG3;
     J1:=-2/(T1*T2);J3:=2/(T3*T2);J2:=J1+J3;
   End;

  Procedure L2R2derivate;
   Begin
    LpriX2:=-GG1*LX1+GG2*LX2+GG3*LX3;
    LpriY2:=-GG1*LY1+GG2*LY2+GG3*LY3;
    LpriZ2:=-GG1*LZ1+GG2*LZ2+GG3*LZ3;
    LsecX2:=J1*LX1-J2*LX2+J3*LX3;
    LsecY2:=J1*LY1-J2*LY2+J3*LY3;
    LsecZ2:=J1*LZ1-J2*LZ2+J3*LZ3;

    Xpri2:=-GG1*XS1+GG2*XS2+GG3*XS3;
    Ypri2:=-GG1*YS1+GG2*YS2+GG3*YS3;
    Zpri2:=-GG1*ZS1+GG2*ZS2+GG3*ZS3;
    Xsec2:=J1*XS1-J2*XS2+J3*XS3;
    Ysec2:=J1*YS1-J2*YS2+J3*YS3;
    Zsec2:=J1*ZS1-J2*ZS2+J3*ZS3;
   End;

  Procedure Sistema4x4;
   Label Skip1,Skip2;
   Var
      i,k,j,p,za   : Integer;
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

  Procedure R2Sigma2; {R2 e L2ùR2 = Prodotto Scalare}
    Begin
     R2Q:= SQR(XS2)+ SQR(YS2)+SQR(ZS2); R2:= SQRT(R2Q);
     R2CosPsi2:=LX2*XS2+LY2*YS2+LZ2*ZS2;
   End;

  Procedure Determinanti(A11,A12,A13,A21,A22,A23,A31,A32,A33:Real);
   Begin
    Det:= A11*(A22*A33-A32*A23)-A12*(A21*A33-A31*A23)+A13*(A21*A32-A22*A31);
   End;

  Procedure CoefEquazLagrange;
   Begin
    KT:=2*R2CosPsi2*Delta1/Delta-SQR(Delta1/Delta)-SQR(R2);
    LT:=-2*Delta2/Delta*(Delta1/Delta-R2CosPsi2);
    MT:=-SQR(Delta2/Delta);
   End;

  Procedure LAGUERRE;
   Label Salta,Salta1;
   Const niter=25;
   Var
       XX      : Array [1..niter] of Real;
      n,fla,ind,kount     : Integer;
      z,ziniz,zit,INCR,RAD1,RAD2,RADI, xrad1,xrad2,
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
      writeln('                                   Radice UNICA ');
          ind:=1;
          goto Salta;
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
 writeln('                 Scegliere tra 2 RADICI aventi r2 > 1.0 UA ');
 textcolor(14);
 readln;
 gotoxy(24,wherey);write('Scelta di r2 tramite l''indice: '); read(ind);

Salta:
  readln;
  textcolor(12);
    REL2:= XX[ind];
 gotoxy(25,wherey+1); write(' r2 = ',REL2:16:12,' UA'); readln;
 textcolor(14);
    End;

  Procedure DistGeoElio2;
   Begin
   RQ:= SQR(REL2); RC:= REL2 * RQ; { Quadrato e Cubo di r2
 Calcolo di ro2 = distanza geocentrica al tempo t2 e della sua Deriv.Prima}
   RO2:= Delta1 / Delta + Delta2 / (Delta * RC);
   ROPri2:= Delta3 / (2 * Delta) + Delta4 / (2 * Delta * RC);

Xel2:= RO2 * LX2 - XS2;
Yel2:= RO2 * LY2 - YS2;
Zel2:= RO2 * LZ2 - ZS2;

XelPri2:= ROPri2 * LX2 + RO2 * LPriX2 - XPri2;
YelPri2:= ROPri2 * LY2 + RO2 * LPriY2 - YPri2;
ZelPri2:= ROPri2 * LZ2 + RO2 * LPriZ2 - ZPri2;
  End;

  Procedure V2rpri2;
   Begin
     V2:=SQRT(SQR(XelPri2)+SQR(YelPri2)+SQR(ZelPri2));
   RelPri2:=(Xel2*XelPri2+Yel2*YelPri2+Zel2*ZelPri2)/REL2;
    End;

  Procedure FG;
   Begin
    FF2:=-1/(2*RC); FF3:=RelPri2/(2*RQ*RQ); G3:=FF2/3; G4:=FF3/2;

    FF1:=1+FF2*SQR(T1)+FF3*T1*SQR(T1); G1:=T1+G3*T1*SQR(T1)+G4*SQR(SQR(T1));
    FF3:=1+FF2*SQR(T3)+FF3*T3*SQR(T3); G3:=T3+G3*T3*SQR(T3)+G4*SQR(SQR(T3));
    End;

  Procedure AFFIN;
   Var
      ResAL1T,ResDE1T,ResAL3T,ResDE3T,
                          A1,B1,A3,B3,
       C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,
                      C11,C12,C13,C14,
              SCAR1,SCAR2,SCAR3,SCAR4,
   DELTARO2,DELTAXELPRI2,DELTAYELPRI2,
                  DELTAZELPRI2,RO2OLD : Real;
   Begin
   writeln;
 {  Subroutine di AFFINAMENTO Residui in à e ë  }
     RESAL1T:= RO1 * RESAL1 / srad; RESDE1T:= RO1 * RESDE1 / srad;
     RESAL3T:= RO3 * RESAL3 / srad; RESDE3T:= RO3 * RESDE3 / srad;
 FG;
     A1:= (3 * (1 - FF1) * (RO2 - R2COSPSI2)) / RQ; B1:= (3 * (T1 - G1) * (RO2 - R2COSPSI2)) / RQ;
     A3:= (3 * (1 - FF3) * (RO2 - R2COSPSI2)) / RQ; B3:= (3 * (T3 - G3) * (RO2 - R2COSPSI2)) / RQ;
  str(A1,rxa1); delete(rxa1,8,6);  str(A3,rxa3); delete(rxa3,8,6);
  str(B1,rxb1); delete(rxb1,8,6);  str(B3,rxb3); delete(rxb3,8,6);

     FX1:=FF1 * LX2 + A1 * XEL2 + B1 * XELPRI2;
     FY1:=FF1 * LY2 + A1 * YEL2 + B1 * YELPRI2;
     FZ1:=FF1 * LZ2 + A1 * ZEL2 + B1 * ZELPRI2;
     FX3:=FF3 * LX2 + A3 * XEL2 + B3 * XELPRI2;
     FY3:=FF3 * LY2 + A3 * YEL2 + B3 * YELPRI2;
     FZ3:=FF3 * LZ2 + A3 * ZEL2 + B3 * ZELPRI2;
        {14 Coefficienti della Matrice di LEUSCHNER}
      C1:= AX1 * FX1 + AY1 * FY1; C2:= AX1 * G1; C3:= AY1 * G1;
      C4:= AX3 * FX3 + AY3 * FY3; C5:= AX3 * G3; C6:= AY3 * G3;
      C7:= DX1 * FX1 + DY1 * FY1 + DZ1 * FZ1;  C8:= DX1 * G1;  C9:= DY1 * G1;
     C11:= DX3 * FX3 + DY3 * FY3 + DZ3 * FZ3; C12:= DX3 * G3; C13:= DY3 * G3;
     C10:= DZ1 * G1; C14:= DZ3 * G3;
  {Trasferimento Coefficienti nella Matrice ASAT[4,4] e termini noti nel Vettore B[4] }
     ASAT[1, 1]:= C1;  ASAT[1, 2]:= C2;   ASAT[1, 3]:= C3;  ASAT[1, 4]:= 0;
     ASAT[2, 1]:= C4;  ASAT[2, 2]:= C5;   ASAT[2, 3]:= C6;  ASAT[2, 4]:= 0;
     ASAT[3, 1]:= C7;  ASAT[3, 2]:= C8;   ASAT[3, 3]:= C9;  ASAT[3, 4]:= C10;
     ASAT[4, 1]:= C11; ASAT[4, 2]:= C12;  ASAT[4, 3]:= C13; ASAT[4, 4]:= C14;
     B[1]:= RESAL1T; B[2]:= RESAL3T; B[3]:= RESDE1T; B[4]:= RESDE3T;

  str(C1,rxc1); delete(rxc1,8,6);  str(C2,rxc2); delete(rxc2,8,6);
  str(C3,rxc3); delete(rxc3,8,6);  str(C4,rxc4); delete(rxc4,8,6);
  str(C5,rxc5); delete(rxc5,8,6);  str(C6,rxc6); delete(rxc6,8,6);
  str(C7,rxc7); delete(rxc7,8,6);  str(C8,rxc8); delete(rxc8,8,6);
  str(C9,rxc9); delete(rxc9,8,6);  str(C10,rxc10); delete(rxc10,8,6);
  str(C11,rxc11); delete(rxc11,8,6);  str(C12,rxc12); delete(rxc12,8,6);
  str(C13,rxc13); delete(rxc13,8,6);  str(C14,rxc14); delete(rxc14,8,6);
  str(ASAT[1,4],rxc00); delete(rxc00,8,6);

  str(RESAL1T,rxre1); delete(rxre1,8,6); str(RESAL3T,rxre2); delete(rxre2,8,6);
  str(RESDE1T,rxre3); delete(rxre3,8,6); str(RESDE3T,rxre4); delete(rxre4,8,6);

 Sistema4x4;

  { Risultati del sistema lineare -> ëro2, ëx'2, ëy'2, ëz'2 }
     DELTARO2:= XP[1]; DELTAXELPRI2:= XP[2]; DELTAYELPRI2:= XP[3]; DELTAZELPRI2:= XP[4];
     RO2OLD:= RO2;
    { Parametri CORRETTI }
     RO2:= RO2 + DELTARO2;

     XELPRI2:= XELPRI2 + DELTAXELPRI2;
     YELPRI2:= YELPRI2 + DELTAYELPRI2;
     ZELPRI2:= ZELPRI2 + DELTAZELPRI2;

     XEL2:= XEL2 + LX2 * DELTARO2;
     YEL2:= YEL2 + LY2 * DELTARO2;
     ZEL2:= ZEL2 + LZ2 * DELTARO2;
        REL2:= SQRT(XEL2*XEL2 + YEL2*YEL2 + ZEL2*ZEL2);
        RELPRI2:= (XEL2 * XELPRI2 + YEL2 * YELPRI2 + ZEL2 * ZELPRI2) / REL2;
     RQ:= SQR(REL2); RC:= REL2 * RQ; { Quadrato e Cubo di r2 }

   SCAR1:= C1*XP[1]+ C2*XP[2]+ C3*XP[3]+  0*XP[4]-RESAL1T;
   SCAR2:= C4*XP[1]+ C5*XP[2]+ C6*XP[3]+  0*XP[4]-RESAL3T;
   SCAR3:= C7*XP[1]+ C8*XP[2]+ C9*XP[3]+C10*XP[4]-RESDE1T;
   SCAR4:=C11*XP[1]+C12*XP[2]+C13*XP[3]+C14*XP[4]-RESDE3T;

  str(DELTARO2,rxs1); delete(rxs1,8,6);  str(DELTAXELPRI2,rxs2); delete(rxs2,8,6);
  str(DELTAYELPRI2,rxs3); delete(rxs3,8,6);  str(DELTAZELPRI2,rxs4); delete(rxs4,8,6);

  str(SCAR1,rxt1); delete(rxt1,6,8);  str(SCAR2,rxt2); delete(rxt2,6,8);
  str(SCAR3,rxt3); delete(rxt3,6,8);  str(SCAR4,rxt4); delete(rxt4,6,8);

  End;

  Procedure CooElioResidui;
   Label Salta1;
   Var
      Xel1,Yel1,Zel1,Xel3,Yel3,Zel3,
      Csi1,Eta1,Zit1,Csi3,Eta3,Zit3: Real;
   Begin
Xel1:= FF1*Xel2+G1*XelPri2;  Xel3:= FF3*Xel2+G3*XelPri2;
Yel1:= FF1*Yel2+G1*YelPri2;  Yel3:= FF3*Yel2+G3*YelPri2;
Zel1:= FF1*Zel2+G1*ZelPri2;  Zel3:= FF3*Zel2+G3*ZelPri2;

Csi1:=Xel1+XS1; Eta1:=Yel1+YS1; Zit1:=Zel1+ZS1;
Csi3:=Xel3+XS3; Eta3:=Yel3+YS3; Zit3:=Zel3+ZS3;

RO1:= LX1*Csi1+LY1*Eta1+Lz1*Zit1; RO3:= LX3*Csi3+LY3*Eta3+Lz3*Zit3;

  ColoriSchermo; Titolo; Video4;

ResAL1:=(-AX1*Csi1-AY1*Eta1)*srad; ResAL3:=(-AX3*Csi3-AY3*Eta3)*srad;
ResDE1:=(-DX1*Csi1-DY1*Eta1-DZ1*Zit1)*srad; ResDE3:=(-DX3*Csi3-DY3*Eta3-DZ3*Zit3)*srad;
textcolor(10);
writeln('                        RESIDUI in AR e Declinazione');textcolor(14);
writeln('        Residui in 1:  R(à1) = ',ResAL1:8:2,'``    R(ë1) = ',ResDE1:8:2,'``');
writeln('        Residui in 3:  R(à3) = ',ResAL3:8:2,'``    R(ë3) = ',ResDE3:8:2,'``');
textcolor(12);
      write('                        Risultati SODDISFACENTI? (S/ ) ');
textcolor(14);
      ch:=readkey;
        if (ch='s') or (ch='S') then
         begin
          goto Salta1;
         end
         else
         begin
         AFFIN; ColoriSchermo; Titolo; Video5;
         V2rpri2; FG; CooElioResidui;
         end;
  Salta1:
    End;

  Procedure ELLISSE;
   Var
      motoang,M2,E2,SinE2,CosE2:  Real;
    Begin
     tr:='Traiettoria ELLITTICA';
     Triv:=a*SQRT(a);
     motoang:=kgauss/Triv;
     SinE2:=sigma2/(e*SQRT(a));
     CosE2:=1/e*(1-REL2/a);
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
     Coshf2:=1/e*(1-REL2/a);
     hf2:=LN(coshf2+SQRT(coshf2*coshf2-1));
     Tper:=-SQRT(-a*a*a)*(e*Sinhf2-hf2)/kgauss;
 End;

  Procedure PARABOLA;
   Var
      d,ppp : Real;
     Begin
     tr:='Traiettoria PARABOLICA';
     d:=sigma2; ppp:=2*REL2-d*d;
     Tper:=-0.5*(ppp*d-d*d*d/3);
 End;

  Procedure ParametriOrbita;
   Var
        x,y,z,Vx,Vy,Vz,V,hi,hj,hk,h,r,
        ni,nj,nk,n,ei,ej,ek,Incl,Nodo  : Real;

    Begin
     sigma2:=Relpri2*REL2;   { r2'= å2/r2 = Velocit… Radiale }
     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î)
   Semiasse magg.,Eccentricit…,Dist.Perielio -> a,e,q }
     a:=REL2/(2-REL2*SQR(V2));
     e:=SQRT(SQR(1-REL2/a)+SQR(REL2*RelPri2)/a);
     q:=a*(1-e);

 { Vettori Posizione e Velocita' in Coord. Eclittiche }
    x:=xel2; y:=yel2*co+zel2*si; z:=-yel2*si+zel2*co; r:=rel2;
   Vx:=xelpri2; Vy:=yelpri2*co+zelpri2*si; Vz:=-yelpri2*si+zelpri2*co; V:= V2;

 { Componenti Vettori h (momento angolare) ed n(nodo) }
  hi:=y*Vz-z*Vy; hj:=z*Vx-x*Vz; hk:=x*Vy-y*Vx;
  h:= SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; nk:= 0; n:= SQRT(hi*hi+hj*hj);

 { Componenti Vettore e (eccentricit…, diretto secondo la linea perifocale) }
  ei:= (Vy * hk - Vz * hj) - x/r;
  ej:= (Vz * hi - Vx * hk) - y/r;
  ek:= (Vx * hj - Vy * hi) - z/r;

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
   repeat until keypressed;
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

{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
       ColoriSchermo;
       LogoDatiLAPLACE;
       SLICE; Pausa; ColoriSchermo;
       TempiCanoniciLAPLACE;
       VettoriL123;
       VettoriA123;
       VettoriD123;
       CostantiGJ;
       R2Sigma2;
  Titolo; Video1; Pausa; ColoriSchermo;
       L2R2derivate;
       Determinanti(LsecX2,LsecY2,LsecZ2,LpriX2,LpriY2,LpriZ2,LX2,LY2,LZ2); Delta:=-Det;
       Determinanti(Xsec2,Ysec2,Zsec2,LpriX2,LpriY2,LpriZ2,LX2,LY2,LZ2); Delta1:=-Det;
       Determinanti(XS2,YS2,ZS2,LpriX2,LpriY2,LpriZ2,LX2,LY2,LZ2); Delta2:=-Det;
       CoefEquazLagrange;
  Titolo; Video2; Pausa; ColoriSchermo;
       LAGUERRE;
       Determinanti(Xsec2,Ysec2,Zsec2,LX2,LY2,LZ2,LsecX2,LsecY2,LsecZ2); Delta3:=-Det;
       Determinanti(XS2,YS2,ZS2,LX2,LY2,LZ2,LsecX2,LsecY2,LsecZ2); Delta4:=-Det;
       DistGeoElio2;V2rpri2;
       FG; CooElioResidui;
       ParametriOrbita;
       GRAFICOconica;
    End.
   {FINE Programma}

