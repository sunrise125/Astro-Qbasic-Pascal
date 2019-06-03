 Unit MIG;        {UNITA' per Calcolo PARAMETRI ORBITALI da 3 Osservazioni

                          METODI di -> Laplace-Leuschner
                                    -> Laplace-Bate
                                    -> Olbers
                                    -> Gauss
}
 Interface

  Uses  Graph,Crt;

    Const
         np=3;
         kgauss=0.01720209895;
         rad=pi/180;
   Type
        MatrINTER = Array [1..np] of Integer;
         MatrREAL = Array [1..np] of Real;
   Var
           ASAT: Array[1..6,1..7] of Real;
           B,XP: Array[1..6] of Real;
      w,tr,ch  : String;
          nme  : String[16];
      eps,si,co,                          
   TG1,TG2,TG3,dt1,dt2,dt3,tau1,tau2,tau3,
      y2t,anom1,anom3,p,         T1,T2,T3,
      Alfa1,Alfa2,Alfa3,Delt1,Delt2,Delt3,
      XS1,YS1,ZS1,XS2,YS2,ZS2,XS3,YS3,ZS3,
      L1x,L1y,L1z,L2x,L2y,L2z,L3x,L3y,L3z,
      R1x,R1y,R1z,R2x,R2y,R2z,r3x,R3y,R3z,
      AX1,AY1,AZ1,AX2,AY2,AZ2,AX3,AY3,AZ3,
      DX1,DY1,DZ1,DX2,DY2,DZ2,DX3,DY3,DZ3,
                    r2,V2,u2,p2,q2,sigma2,
                          f,g,f1,g1,f3,g3,
                     x2,y2,z2,Vx2,Vy2,Vz2,
                  B11,B12,B21,B22,B31,B32,
   a,e,Triv,q,Tper,Om,UQ2,SinV2,CosV2,VV2,
               lam1,lam2,lam3,mu1,mu2,mu3,
                           nu1,nu2,nu3,O2,
      l1,l2,l3,b1,b2,b3,n1,n2,n3,m1,m2,m3,
               x1,y1,z1,x3,y3,z3,r1,r3,V1,
                     GG1,GG2,GG3,J1,J2,J3,
LpriX2,LpriY2,LpriZ2,LsecX2,LsecY2,LsecZ2,
      Xpri2,Ypri2,Zpri2,Xsec2,Ysec2,Zsec2,
    Det,Delta,Delta1,Delta2,Delta3,Delta4,
    R2S,R2CosPsi2,KT,LT,MT,REL2,RQ,RC,
               RO2,Xel2,Yel2,Zel2,RO1,RO3,
          XelPri2,YelPri2,ZelPri2,RelPri2,
                           FF1,FF2,FF3,G4,
              ResAL1,ResDE1,ResAL3,ResDE3,
      aa1,aa2,aa3,bb1,bb2,bb3,cc1,cc2,cc3 : Real;

  Procedure ColoriSchermo;
  Procedure SLICE;
  Procedure TempiCanoniciOLBERS;
  Procedure TempiCanoniciLAPL_BATE;
  Procedure TempiCanoniciLAPL_LEUS;
  Procedure TempiCanoniciGAUSS;
  Procedure VettoriL123;
  Procedure VettoriL123GAUSS;
  Procedure VettoriA123;
  Procedure VettoriD123;
  Procedure VettoriR123;
  Procedure CosDirettoriRo123;
  Procedure LongLatECLITTICA;
  Procedure IterazioneRo2OLBERS;
  Procedure Iterazioner2LAPL_BATE;
  Procedure Iterazioner2GAUSS;
  Procedure ParametriOrbitaPARABOLICA;
  Procedure ParametriOrbitaCONICA;
  Procedure ParametriOrbitaCONICA_LEUS;
  Procedure ParametriOrbitaCONICA_gauss;
  Procedure GRAFICOparabola;
  Procedure GRAFICOconica;
  Procedure Sistema6x6;
  Procedure Sistema4x4;
  Procedure FG(tau:Real);
  Procedure ELLISSE;
  Procedure IPERBOLE;
  Procedure PARABOLA;
  Procedure CostantiGJ;
  Procedure L2R2derivate;
  Procedure R2Sigma2;
  Procedure Determinanti(A11,A12,A13,A21,A22,A23,A31,A32,A33:Real);
  Procedure CoefEquazLagrange;
  Procedure LAGUERRE;
  Procedure DistGeoElio2;
  Procedure V2rpri2;
  Procedure FG_LEUS;
  Procedure AFFIN;
  Procedure CooElioResidui;

      Function ArcoTANG(num,den:Real): Real;
      Function ArcoSIN(valore:Real): Real;
      Function ArcoCOS(valore:Real): Real;


 Implementation


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
     ArcoCOS:=pi/2-arctan(valore/SQRT(1-valore*valore));
  End;



  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
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


  Procedure TempiCanoniciOLBERS;
   Begin
      dt1:=(TG2-TG1)*kgauss;dt2:=(TG3-TG2)*kgauss; dt3:=(TG3-TG1)*kgauss;
   End;

  Procedure TempiCanoniciLAPL_BATE;
   Begin
      Tau1:=(TG1-TG2)*kgauss;Tau3:=(TG3-TG2)*kgauss;
   End;

  Procedure TempiCanoniciLAPL_LEUS;
   Begin
      T1:=(TG1-TG2)*kgauss; T2:=(TG3-TG1)*kgauss; T3:=(TG3-TG2)*kgauss;
   End;

  Procedure TempiCanoniciGAUSS;
   Begin
      tau1:=(TG3-TG2)*kgauss; tau2:=(TG3-TG1)*kgauss; tau3:=(TG2-TG1)*kgauss;
   End;

  Procedure VettoriL123;
   Begin
L1x:= COS(Delt1*rad)*COS(Alfa1*rad); L1y:= COS(Delt1*rad)*SIN(Alfa1*rad); L1z:= SIN(Delt1*rad);
L2x:= COS(Delt2*rad)*COS(Alfa2*rad); L2y:= COS(Delt2*rad)*SIN(Alfa2*rad); L2z:= SIN(Delt2*rad);
L3x:= COS(Delt3*rad)*COS(Alfa3*rad); L3y:= COS(Delt3*rad)*SIN(Alfa3*rad); L3z:= SIN(Delt3*rad);
   End;

  Procedure VettoriL123GAUSS;
   Begin
l1:= COS(Delt1*rad)*COS(Alfa1*rad); m1:= COS(Delt1*rad)*SIN(Alfa1*rad); n1:= SIN(Delt1*rad);
l2:= COS(Delt2*rad)*COS(Alfa2*rad); m2:= COS(Delt2*rad)*SIN(Alfa2*rad); n2:= SIN(Delt2*rad);
l3:= COS(Delt3*rad)*COS(Alfa3*rad); m3:= COS(Delt3*rad)*SIN(Alfa3*rad); n3:= SIN(Delt3*rad);
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

  Procedure VettoriR123;
   Begin
R1x:= -XS1; R1y:= -YS1; R1z:= -ZS1;
R2x:= -XS2; R2y:= -YS2; R2z:= -ZS2;
R3x:= -XS3; R3y:= -YS3; R3z:= -ZS3;
   End;

  Procedure CosDirettoriRo123;
   Begin
     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }
   lam1:=L1x; lam2:=L2x; lam3:=L3x;
    mu1:=L1z*si+L1y*co; mu2:=L2z*si+L2y*co; mu3:=L3z*si+L3y*co;
    nu1:=L1z*co-L1y*si; nu2:=L2z*co-L2y*si; nu3:=L3z*co-L3y*si;
   O2:=ArcoTANG(YS2,co*XS2);
   End;

  Procedure LongLatECLITTICA;
   Begin
    l1:=ArcoTANG(mu1,lam1);
    l2:=ArcoTANG(mu2,lam2);
    l3:=ArcoTANG(mu3,lam3);
    b1:=ArcoSIN(nu1);
    b2:=ArcoSIN(nu2);
    b3:=ArcoSIN(nu3);
   End;

  Procedure IterazioneRo2OLBERS;
  Label Skip;
  Const srad=206265;
  Var
      aux1,aux2,aux3,QQ,          
      Alfa2c,Delt2c,ResAl2,ResDe2,
      ro1,ro3,passo,ro2,          
      c1,c2,cab,fi,csi,eta,fteta  : Real;
   Begin
    aux1:=Sin(l2-O2)/Sin(b2)*Cos(b2);
    aux2:=Sin(l1-O2)/Sin(b1)*Cos(b1);
    aux3:=Sin(l3-O2)/Sin(b3)*Cos(b3);

  QQ:= dt2/dt1*nu1/nu3*(aux1-aux2)/(aux3-aux1);

  ro1:=8.1; passo:=0.1;
     {Inizio CICLO Iterativo}
 repeat
  ro1:=ro1-passo; ro3:= QQ*ro1;
     if ro1<=0.5 then goto Skip;
  x1:=ro1*L1x-XS1; y1:=ro1*L1y-YS1; z1:=ro1*L1z-ZS1;
  x3:=ro3*L3x-XS3; y3:=ro3*L3y-YS3; z3:=ro3*L3z-ZS3;
     r1:=SQRT(x1*x1+y1*y1+z1*z1);   r3:=SQRT(x3*x3+y3*y3+z3*z3);

  c1:=SQRT(SQR(x3-x1)+SQR(y3-y1)+SQR(z3-z1)); {1ø valore della corda c13}

  fi:=2*dt3/((r1+r3)*SQRT(r1+r3));
  csi:=ArcoSIN(0.75*SQRT(2)*fi)/3;
  eta:=ArcoSIN(SQRT(2)*Sin(csi));
  fteta:=Sin(2*eta)/fi;

  c2:=2*dt3/SQRT(r1+r3)*fteta;              {2ø valore della corda c13}
          cab:= ABS(c2-c1);
   
   textcolor(10);
   gotoxy(18,11);writeln('  ro1 = ',ro1:8:6,'    ëc = ',c2-c1:9:7); textcolor(14);
   gotoxy(58,11);writeln('PASSO=',passo:8:6);
  
     if c2-c1 > 0 then
         begin
         ro1:=ro1+passo;
          passo:=passo/10;
         end;

 until cab < 0.0000001;
     {Fine CICLO Iterativo}

Skip:
 textcolor(12);writeln;
  gotoxy(32,13);writeln('c-finale=',c1:11:8); textcolor(14);

  repeat until keypressed;
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

   Sistema6x6;

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
    writeln('                   r2-finale = ',r2it:11:8); textcolor(14);

     r2:= r2it;  { Valore FINALE distanza Eliocentrica (Punto 2) }
     V2:= V2it;  { Valore FINALE VELOCITA'Eliocentrica (Punto 2) }
    repeat until keypressed;
   End;

  Procedure Iterazioner2GAUSS;
    Var
       c1zero,c3zero,c1,c3,nu1,nu3,R2cosPsi2,
       Ag,Bg,Eg,F1,F2,F3,c3ro3,c1ro1,
       a1,a2,a3,a4,b1,b2,b3,b4,b5,b6,
       mu1,mu2,mu3,cc,
       K1,K2,K3,h1,h2,h3,y1t,y3t,
       c1ab,c3ab,c1it,c3it               : Real;

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

   repeat
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

    y1t:=1+10/11*(11/9*h1)/(1+11/9*h1);
    y2t:=1+10/11*(11/9*h2)/(1+11/9*h2);
    y3t:=1+10/11*(11/9*h3)/(1+11/9*h3);

     c1it:=c1zero*y2t/y1t;  c3it:=c3zero*y2t/y3t;

     c1ab:=ABS(c1it-c1);  c3ab:=ABS(c3it-c3);
   writeln('                   nu1 = ', nu1:10:8,'  nu3 = ', nu3:10:8);
   writeln('                    c1 = ',c1it:10:8,'   c3 = ',c3it:10:8);
   writeln('                    r1 = ',r1  :10:8,'   r3 = ',r3  :10:8);
   writeln('          KT = ',KT:10:8,'    LT = ',LT:10:8,'   MT = ',MT:10:8);
   writeln('     r2-iter = ',r2:10:8,'   î(c1)=',c1ab,'   î(c3)=',c3ab);
   repeat until keypressed;

   nu1:=(c1it-c1zero)*r2*r2*r2;
   nu3:=(c3it-c3zero)*r2*r2*r2;

   until (c1ab<0.000001) and (c3ab<0.000001);

   End;

  Procedure ParametriOrbitaPARABOLICA;
  Var
      aux4,aux5,aux6,aux7,V3om,V1om,
      atilde,btilde,dtilde,etilde,   
      xe1,ye1,ze1,xe3,ye3,ze3,       
      sig1,sig2,sig3,S3,             
      Incl,Nodo,Om,                  
      V3,s1,M1,ss3,M3,Gp,Tp          : Real;
     Begin
  {Passaggio da coord.Equatoriali a Eclittiche}
    xe1:=x1; ye1:=y1*co+z1*si; ze1:=-y1*si+z1*co;
    xe3:=x3; ye3:=y3*co+z3*si; ze3:=-y3*si+z3*co;

  sig1:=(xe1*ye3-xe3*ye1)/2;
  sig2:=(ye1*ze3-ye3*ze1)/2;
  sig3:=(xe1*ze3-xe3*ze1)/2;
      S3:=SQRT(sig1*sig1+sig2*sig2+sig3*sig3);
   { Calcolo i,ê }
  Nodo:=ArcoTANG(sig2,sig3);
  Incl:=ArcoCOS(sig1/S3);
    aux4:=ze3/(r3*Sin(Incl)); aux5:=(xe3*Cos(Nodo)+ye3*Sin(Nodo))/r3;
    aux6:=ze1/(r1*Sin(Incl)); aux7:=(xe1*Cos(Nodo)+ye1*Sin(Nodo))/r1;
  V3om:=ArcoTANG(aux4,aux5);
  V1om:=ArcoTANG(aux6,aux7);
   { Differenza Anomalie V3-V1 }
  atilde:= V3om-V1om;
   { Distanza al perielio q }
  q:=r1*r3*SQR(Sin(atilde/2))/(r1+r3-2*SQRT(r1*r3)*Cos(atilde/2));
  btilde:=2*ArcoSIN(q*(r3-r1)*Cos(atilde/2)/S3);
  V1:=(btilde-atilde)/2; V3:=(atilde+btilde)/2;
   { Argomento Perielio ì }
  Om:=V3om-V3;
   { Giorno del Passaggio al Perielio Gp }
  s1:=Sin(V1/2)/Cos(V1/2); M1:=75*s1+25*s1*SQR(s1);
 ss3:=Sin(V3/2)/Cos(V3/2); M3:=75*ss3+25*ss3*SQR(ss3);
  dtilde:=SQRT(2)/(kgauss*75)*q*SQRT(q)*M1;
  etilde:=SQRT(2)/(kgauss*75)*q*SQRT(q)*M3;
  Gp:= (TG1+TG3-dtilde-etilde)/2;
  { Giorni tra il Passaggio al Perielio e la posiz. 2 }
  Tp:= Gp-TG2;

       {Matrice di Trasformazione --> Serve per il GRAFICO dell'Orbita}
 B11:= COS(nodo) * COS(om) - SIN(nodo) * COS(incl) * SIN(om);
 B12:= -COS(nodo) * SIN(om) - SIN(nodo) * COS(incl) * COS(om);
 B21:= SIN(nodo) * COS(om) + COS(nodo) * COS(incl) * SIN(om);
 B22:= -SIN(nodo) * SIN(om) + COS(nodo) * COS(incl) * COS(om);
 B31:= SIN(incl) * SIN(om); B32:= SIN(incl) * COS(om);
 textcolor(10);clrscr;
 gotoxy(18,3);writeln('RISULTATI FINALI - Parametri Orbita PARABOLICA'); textcolor(1);
 gotoxy(34,4);writeln(nme);textcolor(14);
 gotoxy(23,6);writeln('Pass.Perielio  Tp = ',Tp:8:6,' giorni da t2');
 gotoxy(23,7);writeln('Dist.Perielio   q = ',q:10:8,' UA');
 writeln;
 gotoxy(23,9);writeln('Inclinazione    i = ',Incl/rad:12:8,'ø');
 gotoxy(23,10);writeln('Long-Nodo-Asc   ê = ',Nodo/rad:12:8,'ø');
 gotoxy(23,11);writeln('Argom-Perielio  ì = ',Om/rad:12:8,'ø');
  
   readln; repeat until keypressed;
    End;

  Procedure ParametriOrbitaCONICA;
   Var
       alfa,areolk,r2prim      : Real;
       SinI,CosI,SinL2,CosL2,SinU2,CosU2,Incl,Nodo,L2,
               UX2,UY2,UZ2,VVX2,VVY2,VVZ2,WX2,WY2,WZ2,
                      PP,XX2T,XY2T,XZ2T,YX2T,YY2T,YZ2T : Real;
     Begin
     r2prim:=(x2*Vx2+y2*Vy2+z2*Vz2)/r2; { r2'= å2/r2 = Velocit… Radiale }
     sigma2:=r2prim*r2;
     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }

     alfa:=V2*V2-2/r2; a:= -1/alfa;
     areolk:=SQRT(SQR(V2*r2)-SQR(sigma2));
     e:=SQRT(1+alfa*SQR(areolk));
     q:=a*(1-e);

          { Calcolo Parametri i,ê }
     UX2:=x2/r2; UY2:=y2/r2; UZ2:=z2/r2;
     VVX2:=r2*Vx2-sigma2*UX2;
     VVY2:=r2*Vy2-sigma2*UY2;
     VVZ2:=r2*Vz2-sigma2*UZ2; PP:=SQRT(VVX2*VVX2+VVY2*VVY2+VVZ2*VVZ2);
     WX2:=VVX2/PP; WY2:=VVY2/PP; WZ2:=VVZ2/PP;
     XX2T:=UX2; XY2T:=UY2*co+UZ2*si; XZ2T:=UZ2*co-UY2*si;
     YX2T:=WX2; YY2T:=WY2*co+WZ2*si; YZ2T:=WZ2*co-WY2*si;
          SinI:=SQRT(XZ2T*XZ2T+YZ2T*YZ2T);
          CosI:=SQRT(SQR(XX2T+YY2T)+SQR(XY2T-YX2T))-1;
     Incl:=ArcoTANG(SinI,CosI);
          SinL2:=(XY2T-YX2T)/(1+CosI);
          CosL2:=(XX2T+YY2T)/(1+CosI);
     L2:=ArcoTANG(SinL2,CosL2);
          SinU2:=XZ2T/SinI; CosU2:=YZ2T/SinI;
     UQ2:=ArcoTANG(SinU2,CosU2);
     Nodo:=2*pi+L2-UQ2;

             if e <= 0.9995 then ELLISSE
                ELSE if (e < 1.0005) and (e > 0.9995) then PARABOLA
                  ELSE IPERBOLE;
     Om:=UQ2-VV2;

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
   readln; repeat until keypressed;
    End;

  Procedure ParametriOrbitaCONICA_LEUS;
   Var
       SinI,CosI,SinL2,CosL2,SinU2,CosU2,Incl,Nodo,L2,
               UX2,UY2,UZ2,VVX2,VVY2,VVZ2,WX2,WY2,WZ2,
                      PP,XX2T,XY2T,XZ2T,YX2T,YY2T,YZ2T : Real;
    Begin
     sigma2:=Relpri2*REL2;   { r2'= å2/r2 = Velocit… Radiale }
     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }

     a:=REL2/(2-REL2*SQR(V2));
     e:=SQRT(SQR(1-REL2/a)+SQR(REL2*RelPri2)/a);
     q:=a*(1-e);  r2:=REL2;

          { Calcolo Parametri i,ê }
     UX2:=Xel2/REL2; UY2:=Yel2/REL2; UZ2:=Zel2/REL2;
     VVX2:=REL2*XelPri2-sigma2*UX2;
     VVY2:=REL2*YelPri2-sigma2*UY2;
     VVZ2:=REL2*ZelPri2-sigma2*UZ2; PP:=SQRT(VVX2*VVX2+VVY2*VVY2+VVZ2*VVZ2);
     WX2:=VVX2/PP; WY2:=VVY2/PP; WZ2:=VVZ2/PP;
     XX2T:=UX2; XY2T:=UY2*co+UZ2*si; XZ2T:=UZ2*co-UY2*si;
     YX2T:=WX2; YY2T:=WY2*co+WZ2*si; YZ2T:=WZ2*co-WY2*si;
          SinI:=SQRT(XZ2T*XZ2T+YZ2T*YZ2T);
          CosI:=SQRT(SQR(XX2T+YY2T)+SQR(XY2T-YX2T))-1;
     Incl:=ArcoTANG(SinI,CosI);
          SinL2:=(XY2T-YX2T)/(1+CosI);
          CosL2:=(XX2T+YY2T)/(1+CosI);
     L2:=ArcoTANG(SinL2,CosL2);
          SinU2:=XZ2T/SinI; CosU2:=YZ2T/SinI;
     UQ2:=ArcoTANG(SinU2,CosU2);
     Nodo:=2*pi+L2-UQ2;

             if e <= 0.9995 then ELLISSE
                ELSE if (e < 1.0005) and (e > 0.9995) then PARABOLA
                  ELSE IPERBOLE;
     Om:=UQ2-VV2;

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
  {writeln;
   writeln(' Xel2=',Xel2:10:8,'  Yel2=',Yel2:10:8,'   Zel2=',Zel2:10:8);}
   repeat until keypressed;
    End;

  Procedure ParametriOrbitaCONICA_gauss;
   Var
       sig,xt,yt,zt,rt,q1,q3,V3menoV1,up,down,xe2,ye2,ze2,
       SinV1,CosV1,Px,Py,Pz,Qx,Qy,Qz,
       SinI,CosI,SinL2,CosL2,SinU2,CosU2,Incl,Nodo,L2 : Real;
     Begin

     si:=Sin(eps*rad); co:=Cos(eps*rad);{ Sin(î) & Cos(î) }
     sig:=(x1*x3+y1*y3+z1*z3)/SQR(r1);
     xt:=x3-sig*x1;  yt:=y3-sig*y1;  zt:=z3-sig*z1;
     rt:=SQRT(xt*xt+yt*yt+zt*zt);
      {Parametri p,e,a,q}
     p:=SQR(r1/tau2*rt*y2t);
     q1:=p/r1-1;  q3:=p/r3-1;
     V3menoV1:=ArcoCOS((x1*x3+y1*y3+z1*z3)/r1/r3);
     up:=q1*COS(V3menoV1)-q3; down:=SIN(V3menoV1);

     e:=SQRT(SQR(up/down)+SQR(q1));
     a:=p/(1-e*e);
     q:=a*(1-e);

 { Calcolo Parametri Px,Py,Pz e Qx,Qy,Qz }
     CosV1:=q1/e; SinV1:=SQRT(1-CosV1*CosV1);
   anom1:=ArcoCOS(CosV1); anom3:=V3menoV1+anom1;
   Px:= x1*CosV1/r1-xt*SinV1/rt;
   Py:= y1*CosV1/r1-yt*SinV1/rt;
   Pz:= z1*CosV1/r1-zt*SinV1/rt;

   Qx:= x1*SinV1/r1+xt*CosV1/rt;
   Qy:= y1*SinV1/r1+yt*CosV1/rt;
   Qz:= z1*SinV1/r1+zt*CosV1/rt;

        sigma2:=SQRT((e*e-SQR(1-r2/a))*a);

             if e <= 0.9995 then ELLISSE
                ELSE if (e < 1.0005) and (e > 0.9995) then PARABOLA
                  ELSE IPERBOLE;

       { Parametri orbitali ì,ê,i }
   Om:=arcoTANG(Pz*co-Py*si,Qz*co-Qy*si);
   Nodo:=arcoTANG((Py*COS(om)-Qy*SIN(om))/co,Px*COS(om)-Qx*SIN(om));
   Incl:=arcoCOS(-(Px*SIN(om)+Qx*COS(om))/SIN(Nodo));

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
 readln;  repeat until keypressed;
    End;

  Procedure GRAFICOparabola;
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
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l`Eclittica;  LINEA ROSSA=>Orbita SOTTO l`Eclittica');
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
               ro:=ka*2*q/(1+cos(k*rad-V1));
      XORBP:= ro * COS(k*rad-V1); YORBP:= ro * SIN(k*rad-V1);{Coordinate ORBITALI

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
     OutTextXY(xc-310,yc+220,'LINEA VERDE=>Orbita SOPRA l`Eclittica;  LINEA ROSSA=>Orbita SOTTO l`Eclittica');
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

  Procedure Sistema6x6;
   Label Skip1,Skip2;
   Const npp=6;
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

  Procedure Sistema4x4;
   Label Skip1,Skip2;
   Const npp=4;
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

     SinV2:=a/r2*SQRT(1-e*e)*SinE2; CosV2:=a/r2*(CosE2-e);
     VV2:=ArcoTANG(SinV2,CosV2);
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

     SinV2:=-a/r2*SQRT(e*e-1)*Sinhf2; CosV2:=a/r2*(Coshf2-e);
     VV2:=ArcoTANG(SinV2,CosV2);
 End;


  Procedure PARABOLA;
   Var
      d,ppp : Real;
     Begin
     tr:='Traiettoria PARABOLICA';
     d:=sigma2; ppp:=2*r2-d*d;
     Tper:=-0.5*(ppp*d-d*d*d/3);

     SinV2:=d; CosV2:=SQRT(ppp);
     VV2:=2*ArcoTANG(SinV2,CosV2);
 End;

  Procedure CostantiGJ;
   Begin
    GG1:=-T3/(T1*T2);GG3:=-T1/(T3*T2);GG2:=GG1-GG3;
     J1:=-2/(T1*T2);J3:=2/(T3*T2);J2:=J1+J3;
   End;

  Procedure L2R2derivate;
   Begin
    LpriX2:=-GG1*L1x+GG2*L2x+GG3*L3x;
    LpriY2:=-GG1*L1y+GG2*L2y+GG3*L3y;
    LpriZ2:=-GG1*L1z+GG2*L2z+GG3*L3z;
    LsecX2:=J1*L1x-J2*L2x+J3*L3x;
    LsecY2:=J1*L1y-J2*L2y+J3*L3y;
    LsecZ2:=J1*L1z-J2*L2z+J3*L3z;

    Xpri2:=-GG1*XS1+GG2*XS2+GG3*XS3;
    Ypri2:=-GG1*YS1+GG2*YS2+GG3*YS3;
    Zpri2:=-GG1*ZS1+GG2*ZS2+GG3*ZS3;
    Xsec2:=J1*XS1-J2*XS2+J3*XS3;
    Ysec2:=J1*YS1-J2*YS2+J3*YS3;
    Zsec2:=J1*ZS1-J2*ZS2+J3*ZS3;
   End;

  Procedure R2Sigma2; {R2 e L2ùR2 = Prodotto Scalare}
   Var R2Q   : Real;
    Begin
     R2Q:= SQR(XS2)+ SQR(YS2)+SQR(ZS2); R2S:= SQRT(R2Q);
     R2CosPsi2:=L2x*XS2+L2y*YS2+L2z*ZS2;
   End;

  Procedure Determinanti(A11,A12,A13,A21,A22,A23,A31,A32,A33:Real);
   Begin
    Det:= A11*(A22*A33-A32*A23)-A12*(A21*A33-A31*A23)+A13*(A21*A32-A22*A31);
   End;

  Procedure CoefEquazLagrange;
   Begin
    KT:=2*R2CosPsi2*Delta1/Delta-SQR(Delta1/Delta)-SQR(R2S);
    LT:=-2*Delta2/Delta*(Delta1/Delta-R2CosPsi2);
    MT:=-SQR(Delta2/Delta);
      IF (KT>0) or (LT<0) or (MT>0) THEN
   begin writeln;
  writeln('                 Errore nei Dati INIZIALI ...Ricominciare daccapo !!!');
   writeln('   KT = ',KT:10:8,'    LT = ',LT:10:8,'   MT = ',MT:10:8);
      repeat until keypressed;
         halt;   {Torna al DOS}
      End;
   End;

  Procedure LAGUERRE;
   Label Salta;
   Const niter=25; n=8;
   Var
       XX      : Array [1..niter] of Real;
      fla,ind       : Integer;
      z,ziniz,zit,INCR,RAD1,RAD2,RADI,
      zax,RAX,F0,F1,F2                         : Real;

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

   BEGIN
  ClrScr;textcolor(10);
  writeln('                      RADICI della Distanza Eliocentrica r2');textcolor(14);
  fla:=0;incr:=0; RADI:=0; zit:=0; z:=0;
  RAX:=0; RAD1:=0; RAD2:=0; zax:=0;

     ziniz:= 0.1;    {Valore iniziale della 1^RADICE}
repeat
   z:= ziniz+INCR;

  repeat
    F0:= FX(z); F1:= FPriX(z); F2:= FSecX(z); { f, f' e f" }

  RAX:= SQR(n - 1) * SQR(F1) - n * (n - 1) * F0 * F2;
  if (RAX < 0) then
    begin
   GOTO Salta;  {Radicando NEGATIVO -> Salta}
   end;
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

 fla:= fla + 1; {Contatore}
 XX[fla]:= zit;
 writeln('             r2(',fla:2,')=',zit:16:12,'   Residuo = ',FX(zit));
Salta:
 INCR:= INCR + 0.2; {INCREMENTO DI 0.2 UA}
 until INCR > 4.8;

  textcolor(10);

   IF FPriX(1) < 0 THEN
      begin
      writeln('                   Scegliere l`UNICA Radice r2 > 1.0 UA');
      end
      ELSE
      begin
      writeln('                 Scegliere tra 2 RADICI aventi r2 > 1.0 UA');
     end;
    
  textcolor(14);
      write('                        Scelta di r2 tramite l`indice: '); readln(ind);
  textcolor(0);
    REL2:= XX[ind];
 writeln('                           r2 = ',REL2:16:12,' UA'); textcolor(14);
    End;

  Procedure DistGeoElio2;
   Var
       ROPri2   : Real;
   Begin
   RQ:= SQR(REL2); RC:= REL2 * RQ; { Quadrato e Cubo di r2
 Calcolo di ro2 = distanza geocentrica al tempo t2 e della sua Deriv.Prima}
   RO2:= Delta1 / Delta + Delta2 / (Delta * RC);
   ROPri2:= Delta3 / (2 * Delta) + Delta4 / (2 * Delta * RC);

Xel2:= RO2 * L2x - XS2;
Yel2:= RO2 * L2y - YS2;
Zel2:= RO2 * L2z - ZS2;

XelPri2:= ROPri2 * L2x + RO2 * LPriX2 - XPri2;
YelPri2:= ROPri2 * L2y + RO2 * LPriY2 - YPri2;
ZelPri2:= ROPri2 * L2z + RO2 * LPriZ2 - ZPri2;
  End;

  Procedure V2rpri2;
   Begin
     V2:=SQRT(SQR(XelPri2)+SQR(YelPri2)+SQR(ZelPri2));
   RelPri2:=(Xel2*XelPri2+Yel2*YelPri2+Zel2*ZelPri2)/REL2;
    End;

  Procedure FG_LEUS;
   Begin
    FF2:=-1/(2*RC); FF3:=RelPri2/(2*RQ*RQ); G3:=FF2/3; G4:=FF3/2;

    FF1:=1+FF2*SQR(T1)+FF3*T1*SQR(T1); G1:=T1+G3*T1*SQR(T1)+G4*SQR(SQR(T1));
    FF3:=1+FF2*SQR(T3)+FF3*T3*SQR(T3); G3:=T3+G3*T3*SQR(T3)+G4*SQR(SQR(T3));
    End;

  Procedure AFFIN;
  Const srad=206265;
   Var
      ResAL1T,ResDE1T,ResAL3T,ResDE3T,
              FX1,FY1,FZ1,FX3,FY3,FZ3,
                          A1,B1,A3,B3,
       C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,
                      C11,C12,C13,C14,
   DELTARO2,DELTAXELPRI2,DELTAYELPRI2,
                  DELTAZELPRI2,RO2OLD : Real;
   Begin
   writeln;
 {  Subroutine di AFFINAMENTO Residui in à e ë  }
     RESAL1T:= RO1 * RESAL1 / srad; RESDE1T:= RO1 * RESDE1 / srad;
     RESAL3T:= RO3 * RESAL3 / srad; RESDE3T:= RO3 * RESDE3 / srad;
 FG_LEUS;
     A1:= (3 * (1 - FF1) * (RO2 - R2COSPSI2)) / RQ; B1:= (3 * (T1 - G1) * (RO2 - R2COSPSI2)) / RQ;
     A3:= (3 * (1 - FF3) * (RO2 - R2COSPSI2)) / RQ; B3:= (3 * (T3 - G3) * (RO2 - R2COSPSI2)) / RQ;

     FX1:=FF1 * L2x + A1 * XEL2 + B1 * XELPRI2;
     FY1:=FF1 * L2y + A1 * YEL2 + B1 * YELPRI2;
     FZ1:=FF1 * L2z + A1 * ZEL2 + B1 * ZELPRI2;
     FX3:=FF3 * L2x + A3 * XEL2 + B3 * XELPRI2;
     FY3:=FF3 * L2y + A3 * YEL2 + B3 * YELPRI2;
     FZ3:=FF3 * L2z + A3 * ZEL2 + B3 * ZELPRI2;
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

 Sistema4x4;

  { Risultati del sistema lineare -> ëro2, ëx'2, ëy'2, ëz'2 }
     DELTARO2:= XP[1]; DELTAXELPRI2:= XP[2]; DELTAYELPRI2:= XP[3]; DELTAZELPRI2:= XP[4];
     RO2OLD:= RO2;
    { Parametri CORRETTI }
     RO2:= RO2 + DELTARO2;

     XELPRI2:= XELPRI2 + DELTAXELPRI2;
     YELPRI2:= YELPRI2 + DELTAYELPRI2;
     ZELPRI2:= ZELPRI2 + DELTAZELPRI2;

     XEL2:= XEL2 + L2x * DELTARO2;
     YEL2:= YEL2 + L2y * DELTARO2;
     ZEL2:= ZEL2 + L2z * DELTARO2;
        REL2:= SQRT(XEL2*XEL2 + YEL2*YEL2 + ZEL2*ZEL2);
        RELPRI2:= (XEL2 * XELPRI2 + YEL2 * YELPRI2 + ZEL2 * ZELPRI2) / REL2;
     RQ:= SQR(REL2); RC:= REL2 * RQ; { Quadrato e Cubo di r2 }

  End;

  Procedure CooElioResidui;
   Const srad=206265;
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

RO1:= L1x*Csi1+L1y*Eta1+L1z*Zit1; RO3:= L3x*Csi3+L3y*Eta3+L3z*Zit3;

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
         AFFIN; V2rpri2; FG_LEUS; CooElioResidui;
         end;
  Salta1:
    End;
END.
