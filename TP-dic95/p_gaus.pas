Program P_GAUS;       { Determinazione dei Parametri Orbitali di un Corpo Celeste
                        con il metodo di GAUSS, note 3 Osservazioni.}
  Uses  Graph,Crt;

    Const
         np=3;
         kgauss=0.01720209895;
         rad=pi/180;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
   Var
      w,tr,ch  : String;
          nme  : String[16];
      eps,si,co,
               TG1,TG2,TG3,tau1,tau2,tau3,
                                    y2t,p,
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
                    r2,V2,u2,p2,q2,sigma2,
                                 x2,y2,z2,
                  B11,B12,B21,B22,B31,B32,
                       a,e,Triv,q,Tper,Om,
      l1,l2,l3,b1,b2,b3,n1,n2,n3,m1,m2,m3,
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

  Procedure SLICE;
   Var
      fl                       : Text;
      XS,YS,ZS,ARs,DEs,giorno  : MatrREAL;
      ARh,ARm,DEg,DEp          : MatrINTER;
      segno,i                  : Integer;

    Begin
     gotoxy(28,4); writeln('Nome del File-Dati = ');gotoxy(49,4);
     read (w);
     textcolor(10);
     gotoxy(23,5);writeln('RIEPILOGO Dati delle 3 OSSERVAZIONI');textcolor(14);
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

  Procedure TempiCanoniciGAUSS;
   Begin
      tau1:=(TG3-TG2)*kgauss; tau2:=(TG3-TG1)*kgauss; tau3:=(TG2-TG1)*kgauss;
   End;

  Procedure VettoriL123GAUSS;
   Begin
l1:= COS(Delt1*rad)*COS(Alfa1*rad); m1:= COS(Delt1*rad)*SIN(Alfa1*rad); n1:= SIN(Delt1*rad);
l2:= COS(Delt2*rad)*COS(Alfa2*rad); m2:= COS(Delt2*rad)*SIN(Alfa2*rad); n2:= SIN(Delt2*rad);
l3:= COS(Delt3*rad)*COS(Alfa3*rad); m3:= COS(Delt3*rad)*SIN(Alfa3*rad); n3:= SIN(Delt3*rad);
   End;

  Procedure FunzioneIpergeometricaHANSEN( h: Real);
    Var
        gs4,gs3,gs2,gs1,gs0  : Real;
   Begin

 gs3:= 11/9*h/(1+11/9*h/(1+11/9*h/(1+gs4)));
 gs2:= 11/9*h/(1+11/9*h/(1+11/9*h/(1+gs3)));
 gs1:= 11/9*h/(1+11/9*h/(1+11/9*h/(1+gs2)));
 gs0:= 11/9*h/(1+11/9*h/(1+11/9*h/(1+gs1)));

  yyy:= 1+10/9*h/(1+gs0);

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
 gotoxy(24,wherey);write('Scelta di r2 tramite l''indice: '); read(ind);

Salta:
  readln;
  textcolor(12);
    REL2:= XX[ind];
 gotoxy(25,wherey+1); write(' r2 = ',REL2:16:12,' UA'); readln;
 textcolor(14);
    End;

  Procedure Iterazioner2GAUSS;
   Label Salta1,Salta2;
   Var
       c1zero,c3zero,c1,c3,nu1,nu3,R2cosPsi2,
       Ag,Bg,Eg,F1,F2,F3,c3ro3,c1ro1,
       a1,a2,a3,a4,b1,b2,b3,b4,b5,b6,
       cc,
       K1,K2,K3,h1,h2,h3,y1t,y3t,
       c1ab,c3ab,c1it,c3it               : Real;
       ridux1,ridux3:  String;

   Begin
       cc:=2*SQRT(2)/3;
       R2cosPsi2:=l2*XS2+m2*YS2+n2*ZS2;
       R2S:=SQRT(XS2*XS2+YS2*YS2+ZS2*ZS2);
       a1:=m1/l1*l2-m2; a2:=m3-m1/l1*l3; a3:=n1/l1*l2-n2; a4:=n3-n1/l1*l3;
       b1:=YS1-m1/l1*XS1; b2:=m1/l1*XS2-YS2; b3:=YS3-m1/l1*XS3;
       b4:=ZS1-n1/l1*XS1; b5:=n1/l1*XS2-ZS2; b6:=ZS3-n1/l1*XS3;
       Eg:=a3-a4/a2*a1;
       F1:=b4-a4/a2*b1; F2:=b5-a4/a2*b2; F3:=b6-a4/a2*b3;
       c1zero:=tau1/tau2; c3zero:=tau3/tau2;
       nu1:= tau1*tau3/6*(1+c1zero);
       nu3:= tau1*tau3/6*(1+c3zero);
       Ag:=(F1*c1zero+F2+F3*c3zero)/Eg;

Salta1:
       Bg:=(F1*nu1+F3*nu3)/Eg;
  KT:=2*Ag*R2cosPsi2-Ag*Ag-R2S*R2S; LT:=2*Bg*R2cosPsi2-2*Ag*Bg; MT:=-Bg*Bg;
      IF (KT>0) or (LT<0) or (MT>0) THEN
   begin writeln;
  writeln('                 Errore nei Dati INIZIALI ...Ricominciare daccapo !!!');
   writeln('   KT = ',KT:10:8,'    LT = ',LT:10:8,'   MT = ',MT:10:8);
     readln;   repeat until keypressed;
     exit;
      end;
       LAGUERRE;
       r2:=REL2;
    c1:= c1zero+nu1/(r2*r2*r2);
    c3:= c3zero+nu3/(r2*r2*r2);
    RO2:= Ag+Bg/(r2*r2*r2);
    c3ro3:=(b1*c1+b2+b3*c3-a1*RO2)/a2; RO3:=c3ro3/c3;
    c1ro1:=(XS1*c1-XS2+XS3*c3+l2*RO2-l3*c3ro3)/l1; RO1:=c1ro1/c1;

     x1:=RO1*l1-XS1;  y1:=RO1*m1-YS1;  z1:=RO1*n1-ZS1;
     x2:=RO2*l2-XS2;  y2:=RO2*m2-YS2;  z2:=RO2*n2-ZS2;
     x3:=RO3*l3-XS3;  y3:=RO3*m3-YS3;  z3:=RO3*n3-ZS3;
      r1:=SQRT(x1*x1+y1*y1+z1*z1);  r3:=SQRT(x3*x3+y3*y3+z3*z3);

     K1:=SQRT(r2*r3+x2*x3+y2*y3+z2*z3);
     K2:=SQRT(r1*r3+x1*x3+y1*y3+z1*z3);
     K3:=SQRT(r1*r2+x1*x2+y1*y2+z1*z2);

     h1:=tau1*tau1/(K1*K1*(cc*K1+r2+r3));
     h2:=tau2*tau2/(K2*K2*(cc*K2+r1+r3));
     h3:=tau3*tau3/(K3*K3*(cc*K3+r1+r2));

   FunzioneIpergeometricaHANSEN(h1); y1t:=yyy;
   FunzioneIpergeometricaHANSEN(h2); y2t:=yyy;
   FunzioneIpergeometricaHANSEN(h3); y3t:=yyy;


     c1it:=c1zero*y2t/y1t;  c3it:=c3zero*y2t/y3t;

     c1ab:=ABS(c1it-c1);  c3ab:=ABS(c3it-c3);
   str(c1ab,ridux1);  delete(ridux1,6,8);
   str(c3ab,ridux3);  delete(ridux3,6,8);

   writeln('                   nu1 = ', nu1:10:8,'  nu3 = ', nu3:10:8);
   writeln('                    c1 = ',c1it:10:8,'   c3 = ',c3it:10:8);
   writeln('                    r1 = ',r1  :10:8,'   r3 = ',r3  :10:8);
   writeln('         y1t = ',y1t :10:8,'  y2t = ',y2t :10:8,'     y3t = ',y3t:10:8);
   writeln('          KT = ',KT:10:8,'    LT = ',LT:10:8,'   MT = ',MT:10:8);
   writeln('         r2-iter = ',r2:10:8,'   î(c1)=',ridux1,'   î(c3)=',ridux3);

   nu1:=(c1it-c1zero)*r2*r2*r2;
   nu3:=(c3it-c3zero)*r2*r2*r2;

   textcolor(12);
      write('                    Risultati î(c1),î(c3) SODDISFACENTI? (S/ ) ');
   textcolor(14);

      ch:=readkey;
        if (ch='s') or (ch='S') then goto Salta2
         ELSE goto Salta1;

  Salta2:
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

  Procedure ParametriOrbitaCONICA_gauss;
   Var
       d1,d2,d3,d,nn1,nn2,nn3,nn,s1,s2,s3,s,
       Px,Py,Pz,Qx,Qy,Qz, V2x,V2y,V2z,
       bb1,bb2,bb3,lb1,lb2,lb3,lb,
       hi,hj,hk,h,ni,nj,nk,n,ei,ej,ek,
       ye1,ze1,ye2,ze2,ye3,ze3,
       Incl,Nodo                      : Real;

     Begin
 { Passaggio dalle Coord. EQUATORIALI a ECLITTICHE
   Le ascisse x1,x2,x3 restano INVARIATE }

     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }
  ye1:= y1*co+z1*si; ze1:=-y1*si+z1*co;
  ye2:= y2*co+z2*si; ze2:=-y2*si+z2*co;
  ye3:= y3*co+z3*si; ze3:=-y3*si+z3*co;
     y1:=ye1; z1:=ze1; y2:=ye2; z2:=ze2; y3:=ye3; z3:=ze3;

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

     { Calcolo Versori Perifocali P,Q
       Per ora NON SERVONO }
  Px:=(s2*nn3-s3*nn2)/(s*nn); Py:=(s3*nn1-s1*nn3)/(s*nn); Pz:=(s1*nn2-s2*nn1)/(s*nn);
  Qx:=s1/s;  Qy:=s2/s; Qz:=s3/s;

    { Calcolo del Vettore B2 }
bb1:=d2*z2-d3*y2; bb2:=d3*x2-d1*z2; bb3:=d1*y2-d2*x2;
  lb:= SQRT(1/(d*nn));

 { Componenti della VELOCITA' V2 }
V2x:= lb*(bb1/r2+s1); V2y:= lb*(bb2/r2+s2); V2z:= lb*(bb3/r2+s3);
    V2:= SQRT(V2x*V2x+V2y*V2y+V2z*V2z);

     { Prodotto Scalare r2ùV2 = å2 }
      sigma2:= x2*V2x+y2*V2y+z2*V2z;

 { Componenti Vettori h,n,e }
  hi:=y2*V2z-z2*V2y;   hj:=z2*V2x-x2*V2z;   hk:=x2*V2y-y2*V2x;
            h:=SQRT(hi*hi+hj*hj+hk*hk);
  ni:=-hj; nj:=hi; n:=SQRT(ni*ni+nj*nj);
  ei:=V2y*hk-V2z*hj-x2/r2; ej:=V2z*hi-V2x*hk-y2/r2; ek:=V2x*hj-V2y*hi-z2/r2;
            e:=SQRT(ei*ei+ej*ej+ek*ek);

      { Parametri Orbitali i,ê,ì }
   Incl:=ArcoCOS(hk/h);
   Nodo:=ArcoCOS(ni/n);
         if hi<0 then Nodo:=2*pi-Nodo;
   Om:=ArcoCOS((ni*ei+nj*ej)/(n*e));
        if ek<0 then Om:=2*pi-Om;

             if e <= 0.9995 then ELLISSE
                ELSE if (e < 1.0005) and (e > 0.9995) then PARABOLA
                  ELSE IPERBOLE;

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
       Titolo;
       SLICE;
       TempiCanoniciGAUSS;
       VettoriL123GAUSS;
       Iterazioner2GAUSS;
       ParametriOrbitaCONICA_gauss;
       GRAFICOconica;
     End.
   {FINE Programma}

