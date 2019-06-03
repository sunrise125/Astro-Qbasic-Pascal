 Program PV_EQUI2;   {Vettori Posizione e Velocit… 8 Pianeti, ricavati
    (vers-Print)      a partire dalla veriabili ELLITTICHE: (a.la,k,h,q,p)

                      2 Equinozi: J2000 ed Equin.medio della DATA
                      ------------------------------------------- }
{$N+}
 Uses Crt;
 Const
       rad=pi/180;
       mu=2.95912208285591e-4;  {param.gravit.= kgaussý= 0.01720209895ý}
 Var
       XX,YY,ZZ,VXX,VYY,VZZ:          Array[0..8] of Double;
       Xd,Yd,Zd,VXd,VYd,VZd:          Array[0..8] of Double;
       i,item               : Integer;
       tau,jd,
       a,l,k,h,q,p, e, r,
       Xecl,Yecl,Zecl, Xpunto,Ypunto,Zpunto : Double;
       ft       : Text;

  Function MODULO(x,y:Double) : Double;
   Var modu : Double;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;

      Function ArcoTANG(num,den:Double): Double;
    Var atan: Double;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

  Procedure GiornoGiuliano;
 Var m,ha,mn: WORD;
        y: INTEGER;
  g,mm,yy,a,bw,kk,c,d: Double;
  dd1,dd2,dd3,dd4,dd5,tm1,tm2: STRING;

    BEGIN
    textcolor(14);textbackground(1);
     clrscr;
     gotoxy(14,1); writeln('Vettori POSIZIONE/VELOCITA'' degli 8 Pianeti Maggiori');
     gotoxy(14,2); writeln('             Teoria Planetaria: VSOP82');
   textbackground(3);
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄ EPOCA di riferimento ÄÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('³                                    ³');
gotoxy(20,7); write('³                                    ³');
gotoxy(20,8); write('ÀÄÄÄÄÄÄÄÄ Eclittica Dinamica ÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(25,5);write('  Giorno = ±±.±±±±±  TD');gotoxy(36,5); readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(25,6);write('    Mese = ±±');gotoxy(36,6); readln(m);
      until m in [1..12];
       gotoxy(25,7);write('    Anno = ±±±±');gotoxy(36,7); readln(y);
   textbackground(1);
gotoxy(32,10); write('Premere un tasto');
  readln; clrscr;
       yy:=y; mm:=m;
          IF m <= 2 then
           begin
           y:=y-1; m:=m+12;
           end;
              bw:=0; kk:=yy+mm/100+g/10000;
          IF kk > 1582.101499 then    {Data >= 15.10.1582--->Rif.Gregoriana}
           begin
          a:=int(y/100); bw:=2 - a + int(a/4);
           end;
              c:=int(365.25*y);
              d:=int(30.6001*(m+1));
      jd:=bw+c+d+g+1720994.5;

  str(int(jd):7:0,tm1);
  str(int(frac(jd)*1e5+0.5)/1e5:7:5,tm2);
  delete(tm2,1,1);
  textcolor(1);textbackground(15);
  gotoxy(15,1); writeln('(VSOP82)   JDT=',tm1+tm2,' -> 1ö8= Pianeti - Equin. J2000'); textbackground(1);
                writeln(ft,'(VSOP82)   JDT=',tm1+tm2,' -> 1ö8= Pianeti - Equin. J2000');
  textcolor(14);textbackground(1);
       { Millenni Giuliani -> tau:=(jd-2451545)/365250;}
  {Per evitare la perdita di cifre decim. in jd: 1720994.5-2451545= -730550.5}
  tau:=(bw+c+d+g-730550.5)/365250;
     end;

  Procedure PosVelo;
 Var g,gpri,variaz,F, Psi,Px,Py,Pz,Qx,Qy,Qz,
     AA,BB, Apunto,Bpunto,
     pp,efpunto, Csipunto,Etapunto    : Double;
  begin
   F:=l;
 repeat
   g:=F-k*SIN(F)+h*COS(F)-l;   {Equaz. di Keplero}
   gpri:=1-k*COS(F)-h*SIN(F);
   variaz:=-g/gpri;
   F:=F+variaz;
 until abs(variaz)<1e-9;
  Psi:=1/(1+sqrt(1-h*h-k*k));
  AA:=-a*k+a*(1-h*h*Psi)*COS(F)+a*h*k*Psi*SIN(F);
  BB:=-a*h+a*(1-k*k*Psi)*SIN(F)+a*h*k*Psi*COS(F);
 {Versori }
   Px:=(1-2*p*p); Py:=2*p*q; Pz:=-2*p*sqrt(1-p*p-q*q);
   Qx:=2*p*q; Qy:=(1-2*q*q); Qz:=2*q*sqrt(1-p*p-q*q);
 {Vettore POSIZIONE in UA: coord.rettang. eclittiche, Equin. J2000}
    Xecl:=Px*AA+Qx*BB;
    Yecl:=Py*AA+Qy*BB;
    Zecl:=Pz*AA+Qz*BB; r:=sqrt(SQR(Xecl)+SQR(Yecl)+SQR(Zecl));

  e:=sqrt(h*h+k*k);    pp:=a*(1-e*e);  {parametro orbitale}

 {Vettore VELOCITA'in UA/giorno: coord.rettang. eclittiche, Equin. J2000}
  efpunto:= sqrt(mu/pp)*sqrt(1-e*e)/r;
  Apunto:=efpunto*(-a*(1-h*h*Psi)*SIN(F)+a*h*k*Psi*COS(F));
  Bpunto:=efpunto*(a*(1-k*k*Psi)*COS(F)-a*h*k*Psi*SIN(F));

    Xpunto:=Px*Apunto+Qx*Bpunto;
    Ypunto:=Py*Apunto+Qy*Bpunto;
    Zpunto:=Pz*Apunto+Qz*Bpunto;
  end;


   Procedure VarEllitMERC;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..3] of Double;
 begin
  LA[0]:=4.4026088424; LA[1]:=26087.9031415742; LA[2]:=-9.3429e-6; LA[3]:=31e-9;
  KA[0]:=0.0446605976; KA[1]:=-0.00552114624; KA[2]:=-18.60397e-6; KA[3]:=633.32e-9;
  HA[0]:=0.20072331368; HA[1]:=0.00143750118; HA[2]:=-79.74689e-6; HA[3]:=-263.09e-9;
  QA[0]:=0.04061563384; QA[1]:=0.00065433117; QA[2]:=-10.71215e-6; QA[3]:=211.49e-9;
  PA[0]:=0.04563550461; PA[1]:=-0.00127633657; PA[2]:=-9.1335e-6; PA[3]:=180.04e-9;

    a:=0.38709830982;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=((KA[3]*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=((HA[3]*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=1;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitVENE;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..3] of Double;
 begin
  LA[0]:=3.17614669689; LA[1]:=10213.285546211; LA[2]:=2.87555e-6; LA[3]:=-30.38e-9;
  KA[0]:=-0.00449282133; KA[1]:=0.00031259019; KA[2]:=6.05913e-6; KA[3]:=-692.39e-9;
  HA[0]:=0.00506684726; HA[1]:=-0.00036121239; HA[2]:=18.39627e-6; HA[3]:=-9.71e-9;
  QA[0]:=0.00682410142; QA[1]:=0.00138133826; QA[2]:=-10.90942e-6; QA[3]:=-1859.2e-9;
  PA[0]:=0.02882285775; PA[1]:=-0.00040384791; PA[2]:=-62.32891e-6; PA[3]:=251.37e-9;

    a:=0.72332981996;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=((KA[3]*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=((HA[3]*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=2;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitTERR;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..3] of Double;
 begin
  LA[0]:=1.75347031435; LA[1]:=6283.07584918; LA[2]:=-9.9189e-6; LA[3]:=0.73e-9;
  KA[0]:=-0.00374081650; KA[1]:=-0.00082266699; KA[2]:=27.48939e-6; KA[3]:=1042.17e-9;
  HA[0]:=0.01628447663; HA[1]:=-0.00062030259; HA[2]:=-33.53888e-6; HA[3]:=711.85e-9;
  QA[0]:=0; QA[1]:=-0.00113469002; QA[2]:=12.37314e-6; QA[3]:=1270.5e-9;
  PA[0]:=0; PA[1]:=0.00010180391; PA[2]:=47.01998e-6; PA[3]:=-538.29e-9;

    a:=1.00000101778;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=((KA[3]*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=((HA[3]*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=3;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitMART;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..3] of Double;
 begin
  LA[0]:=6.20348091341; LA[1]:=3340.6124314923; LA[2]:=4.54761e-6; LA[3]:=-50.57e-9;
  KA[0]:=0.08536560252; KA[1]:=0.00376330152; KA[2]:=-246.57416e-6; KA[3]:=-3952.41e-9;
  HA[0]:=-0.03789973236; HA[1]:=0.00624657465; HA[2]:=155.27232e-6; HA[3]:=-6719.4e-9;
  QA[0]:=0.01047042574; QA[1]:=0.00017138526; QA[2]:=-40.77591e-6; QA[3]:=-1386.0e-9;
  PA[0]:=0.01228449307; PA[1]:=-0.00108020083; PA[2]:=-19.22195e-6; PA[3]:=883.73e-9;

    a:=1.52367934191+0.31e-9*tau;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=((KA[3]*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=((HA[3]*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=4;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitGIOV;
Var
       AA,LA,QA,PA:    Array[0..3] of Double;
             KA,HA:    Array[0..4] of Double;
 begin
  LA[0]:=0.59954649739; LA[1]:=529.6909650946; LA[2]:=-148.37133e-6; LA[3]:=74.82e-9;
  KA[0]:=0.04698572124; KA[1]:=0.00113010377; KA[2]:=-109.30126e-6; KA[3]:=-4.28748e-6; KA[4]:=205.39e-9;
  HA[0]:=0.01200385748; HA[1]:=0.0021714936; HA[2]:=98.58539e-6; HA[3]:=-5.13109e-6; HA[4]:=-90.07e-9;
  QA[0]:=-0.00206561098; QA[1]:=-0.00031340156; QA[2]:=-16.67392e-6; QA[3]:=769.26e-9;
  PA[0]:=0.01118377157; PA[1]:=-0.00023427562; PA[2]:=20.8676e-6; PA[3]:=507.21e-9;

    a:=5.20260319132+1.91323e-6*tau;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=(((KA[4]*tau+KA[3])*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=(((HA[4]*tau+HA[3])*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=5;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitSATU;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..4] of Double;
 begin
  LA[0]:=0.8740167565; LA[1]:=213.299095438; LA[2]:=366.59741e-6; LA[3]:=-0.3333e-6; LA[4]:=2.17e-9;
  KA[0]:=-0.00296003595; KA[1]:=-0.00529602626; KA[2]:=309.28405e-6; KA[3]:=12.96215e-6; KA[4]:=-599.59e-9;
  HA[0]:=0.05542964254; HA[1]:=-0.00375593887; HA[2]:=-319.90236e-6; HA[3]:=15.98633e-6; HA[4]:=324.51e-9;
  QA[0]:=-0.00871747436; QA[1]:=0.00080171499; QA[2]:=41.42282e-6; QA[3]:=-1.96049e-6; QA[4]:=-94.39e-9;
  PA[0]:=0.01989147301; PA[1]:=0.00059439766; PA[2]:=-52.35117e-6; PA[3]:=-1.27219e-6; PA[4]:=82.95e-9;

    a:=9.55490959574-21.38917e-6*tau;
    l:=MODULO((((LA[4]*tau+LA[3])*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=(((KA[4]*tau+KA[3])*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=(((HA[4]*tau+HA[3])*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=(((QA[4]*tau+QA[3])*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=(((PA[4]*tau+PA[3])*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=6;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitURAN;
Var
       AA,LA,QA,PA:    Array[0..3] of Double;
             KA,HA:    Array[0..4] of Double;
 begin
  LA[0]:=5.48129387159; LA[1]:=74.7815985673; LA[2]:=-8.48828e-6; LA[3]:=104.5e-9;
  KA[0]:=-0.04595132376; KA[1]:=0.0001834405; KA[2]:=-0.80849e-6; KA[3]:=-0.45396e-6; KA[4]:=21.85e-9;
  HA[0]:=0.00563791307; HA[1]:=-0.0007496435; HA[2]:=12.102e-6;  HA[3]:=-0.42088e-6; HA[4]:=-17.14e-9;
  QA[0]:=0.00185915075; QA[1]:=-0.00012449382; QA[2]:=-2.07373e-6; QA[3]:=76.21e-9;
  PA[0]:=0.00648617008; PA[1]:=-0.00011744733; PA[2]:=3.17799e-6; PA[3]:=73.17e-9;

    a:=19.21844606178-0.37163e-6*tau+97.91e-9*tau*tau;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=(((KA[4]*tau+KA[3])*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=(((HA[4]*tau+HA[3])*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=7;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure VarEllitNETT;
Var
       AA,LA,KA,HA,QA,PA:    Array[0..3] of Double;
 begin
  LA[0]:=5.31188628676; LA[1]:=38.1330356378; LA[2]:=1.02311e-6; LA[3]:=-43.4e-9;
  KA[0]:=0.00599977571; KA[1]:=8.71279e-6; KA[2]:=-1.19902e-6; KA[3]:=-40.34e-9;
  HA[0]:=0.00669242413; HA[1]:=78.24336e-6; HA[2]:=0.80801e-6; HA[3]:=-39.55e-9;
  QA[0]:=-0.01029147819; QA[1]:=-0.72727e-6; QA[2]:=-65.68e-9; QA[3]:=16.68e-9;
  PA[0]:=0.01151683985; PA[1]:=0.00002575536; PA[2]:=0.19377e-6; PA[3]:=13.31e-9;

    a:=30.11038686942-1.66346e-6*tau+68.57e-9*tau*tau;
    l:=MODULO(((LA[3]*tau+LA[2])*tau+LA[1])*tau+LA[0],2*pi);
       if l<0 then l:=2*pi+l;    {l=Long_media in RAD}
    k:=((KA[3]*tau+KA[2])*tau+KA[1])*tau+KA[0];
    h:=((HA[3]*tau+HA[2])*tau+HA[1])*tau+HA[0];
    q:=((QA[3]*tau+QA[2])*tau+QA[1])*tau+QA[0];
    p:=((PA[3]*tau+PA[2])*tau+PA[1])*tau+PA[0];
  PosVelo;
    item:=8;
  XX[item]:=Xecl;   YY[item]:=Yecl;    ZZ[item]:=Zecl;
  VXX[item]:=Xpunto;   VYY[item]:=Ypunto;    VZZ[item]:=Zpunto;
    end;

   Procedure Testata;
 begin
  textbackground(4);
  writeln('(X,Y,Z) in UA --- Coord. ECLITTICHE Eliocentriche --- (Vx,Vy,Vz) in UA/giorno');
  writeln(ft,'(X,Y,Z) in UA --- Coord. ECLITTICHE Eliocentriche --- (Vx,Vy,Vz) in UA/giorno');
  textbackground(1);
    end;

   Procedure VideoRis;
 begin
  textbackground(4);
write(item:1);  textbackground(1);
writeln(Xecl:13:8,Yecl:13:8,Zecl:13:8,'³',Xpunto:12:9,Ypunto:13:9,Zpunto:13:9);
write(ft,item:1);
writeln(ft,Xecl:13:8,Yecl:13:8,Zecl:13:8,'³',Xpunto:12:9,Ypunto:13:9,Zpunto:13:9);
    end;

   Procedure MatriceRotazione;
Var
       S11,C11,S12,C12,S13,C13,A31,A32,A33:    Array[0..6] of Double;
       S11t,C11t,S12t,C12t,S13t,C13t,A31t,A32t,A33t,
       A11t,A12t,A13t,A21t,A22t,A23t,
       Csi,Si,Co: Double;
 begin
  S11[0]:=0; S11[1]:=0; S11[2]:=-538.867722e-6; S11[3]:=-270.670e-9;
             S11[4]:=1138.205e-9; S11[5]:=8604e-12; S11[6]:=-813e-12;

  C11[0]:=1; C11[1]:=0; C11[2]:=-20.728e-9; C11[3]:=-19.147e-9;
             C11[4]:=-149.390e-9; C11[5]:=-34e-12; C11[6]:=617e-12;

  S12[0]:=-1; S12[1]:=0; S12[2]:=2.575043e-6; S12[3]:=-56.157e-9;
             S12[4]:=140.001e-9; S12[5]:=383e-12; S12[6]:=-613e-12;

  C12[0]:=0; C12[1]:=0; C12[2]:=-539.329786e-6; C12[3]:=-479.046e-9;
             C12[4]:=1144.883e-9; C12[5]:=8884e-12; C12[6]:=-830e-12;

  S13[0]:=0; S13[1]:=2.269380040e-3; S13[2]:=-24.745348e-6; S13[3]:=-2422.542e-9;
             S13[4]:=78.247e-9; S13[5]:=-468e-12; S13[6]:=-134e-12;

  C13[0]:=0; C13[1]:=-0.203607820e-3; C13[2]:=-94.040878e-6; C13[3]:=2307.025e-9;
             C13[4]:=37.729e-9; C13[5]:=-4862e-12; C13[6]:=25e-12;

  A31[0]:=0; A31[1]:=0.203607820e-3; A31[2]:=94.040878e-6; A31[3]:=-1083.606e-9;
             A31[4]:=-50.218e-9; A31[5]:=929e-12; A31[6]:=11e-12;

  A32[0]:=0; A32[1]:=2.269380040e-3; A32[2]:=-24.745348e-6; A32[3]:=-2532.307e-9;
             A32[4]:=27.473e-9; A32[5]:=643e-12; A32[6]:=-1e-12;

  A33[0]:=1; A33[1]:=0; A33[2]:=-2.595771e-6; A33[3]:=37.009e-9;
             A33[4]:=1.236e-9; A33[5]:=-13e-12; A33[6]:=0;

{
  -----------------------------------------------------------------------------
          0           1             2            3          4         5      6
  -----------------------------------------------------------------------------
  s11     0           0        -538 867 722    -270 670  1 138 205  8 604  -813
  c11    1E12         0             -20 728     -19 147   -149 390    -34   617
  s12   -1E12         0           2 575 043     -56 157    140 001    383  -613
  c12     0           0        -539 329 786    -479 046  1 144 883  8 884  -830
  s13     0     2 269 380 040   -24 745 348  -2 422 542     78 247   -468  -134
  c13     0      -203 607 820   -94 040 878   2 307 025     37 729 -4 862    25
  a31     0       203 607 820    94 040 878  -1 083 606    -50 218    929    11
  a32     0     2 269 380 040   -24 745 348  -2 532 307     27 473    643    -1
  a33    1E12         0          -2 595 771      37 009      1 236    -13     0
  -----------------------------------------------------------------------------
}

  S11t:=(((((S11[6]*tau+S11[5])*tau+S11[4])*tau+S11[3])*tau+S11[2])*tau+
        S11[1])*tau+S11[0];

  S12t:=(((((S12[6]*tau+S12[5])*tau+S12[4])*tau+S12[3])*tau+S12[2])*tau+
        S12[1])*tau+S12[0];

  S13t:=(((((S13[6]*tau+S13[5])*tau+S13[4])*tau+S13[3])*tau+S13[2])*tau+
        S13[1])*tau+S13[0];

  C11t:=(((((C11[6]*tau+C11[5])*tau+C11[4])*tau+C11[3])*tau+C11[2])*tau+
        C11[1])*tau+C11[0];

  C12t:=(((((C12[6]*tau+C12[5])*tau+C12[4])*tau+C12[3])*tau+C12[2])*tau+
        C12[1])*tau+C12[0];

  C13t:=(((((C13[6]*tau+C13[5])*tau+C13[4])*tau+C13[3])*tau+C13[2])*tau+
        C13[1])*tau+C13[0];

  A31t:=(((((A31[6]*tau+A31[5])*tau+A31[4])*tau+A31[3])*tau+A31[2])*tau+
        A31[1])*tau+A31[0];

  A32t:=(((((A32[6]*tau+A32[5])*tau+A32[4])*tau+A32[3])*tau+A32[2])*tau+
        A32[1])*tau+A32[0];

  A33t:=(((((A33[6]*tau+A33[5])*tau+A33[4])*tau+A33[3])*tau+A33[2])*tau+
        A33[1])*tau+A33[0];

   Csi:=0.243817483530*tau; Si:=SIN(Csi);  Co:=COS(Csi);

   A11t:=S11t*Si+C11t*Co;  A12t:=S12t*Si+C12t*Co;  A13t:=S13t*Si+C13t*Co;
   A21t:=C11t*Si-S11t*Co;  A22t:=C12t*Si-S12t*Co;  A23t:=C13t*Si-S13t*Co;

  { Pos/vel all'equinozio della Data }
   writeln('                             Equinozio della DATA'); textbackground(1);
   writeln(ft,'                             Equinozio della DATA'); textbackground(1);
   Testata;
     For i:=1 to 8 do begin
   Xd[i]:=A11t*XX[i]+A12t*YY[i]+A13t*ZZ[i];
   Yd[i]:=A21t*XX[i]+A22t*YY[i]+A23t*ZZ[i];
   Zd[i]:=A31t*XX[i]+A32t*YY[i]+A33t*ZZ[i];

  VXd[i]:=A11t*VXX[i]+A12t*VYY[i]+A13t*VZZ[i];
  VYd[i]:=A21t*VXX[i]+A22t*VYY[i]+A23t*VZZ[i];
  VZd[i]:=A31t*VXX[i]+A32t*VYY[i]+A33t*VZZ[i];

  textbackground(4);
write(i:1);  textbackground(1);
writeln(Xd[i]:13:8,Yd[i]:13:8,Zd[i]:13:8,'³',VXd[i]:12:9,VYd[i]:13:9,VZd[i]:13:9);
write(ft,i:1);
writeln(ft,Xd[i]:13:8,Yd[i]:13:8,Zd[i]:13:8,'³',VXd[i]:12:9,VYd[i]:13:9,VZd[i]:13:9);
     end; { ... for }

   writeln('                     ³',A11t:12:8,A12t:13:8,A13t:13:8,' ³');
   writeln('  Matrice Rotaz. A = ³',A21t:12:8,A22t:13:8,A23t:13:8,' ³');
   writeln('                     ³',A31t:12:8,A32t:13:8,A33t:13:8,' ³');

writeln(ft);
writeln(ft,'                     ³',A11t:12:8,A12t:13:8,A13t:13:8,' ³');
writeln(ft,'  Matrice Rotaz. A = ³',A21t:12:8,A22t:13:8,A23t:13:8,' ³');
writeln(ft,'                     ³',A31t:12:8,A32t:13:8,A33t:13:8,' ³');

writeln(ft);
writeln(ft,'                  Conversione da (X,Y,Z)-J2000 a (X,Y,Z)-Data');
writeln(ft,'        X,Y,Z = Coordinate eclittiche eliocentriche Equinozio J2000');
writeln(ft,'     Xd,Yd,Zd = Coordinate eclittiche eliocentriche Equinozio della Data');

writeln(ft,'                   ³ Xd ³   ³a11  a12  a13 ³ ³ X ³ ');
writeln(ft,'                   ³    ³   ³              ³ ³   ³ ');
writeln(ft,'                   ³ Yd ³ = ³a21  a22  a23 ³ù³ Y ³ ');
writeln(ft,'                   ³    ³   ³              ³ ³   ³ ');
writeln(ft,'                   ³ Zd ³   ³a31  a32  a33 ³ ³ Z ³ ');

writeln(ft);
writeln(ft,'                                          -1');
writeln(ft,'  NOTA IMP.: La Matrice di Rotaz.INVERSA A  si ricava da A invertendo righe ');
writeln(ft,'             con colonne, perchŠ A Š una Matrice UNITARIA');
 end;

 {Inizio Programma-----------------------------------------------------------}
  Begin
   ASSIGN (ft,'EQUI2.RIS');
   Rewrite(ft);
  GiornoGIULIANO;
  Testata;
  VarEllitMERC; VideoRis;
  VarEllitVENE; VideoRis;
  VarEllitTERR; VideoRis;
  VarEllitMART; VideoRis;
  VarEllitGIOV; VideoRis;
  VarEllitSATU; VideoRis;
  VarEllitURAN; VideoRis;
  VarEllitNETT; VideoRis;

   MatriceRotazione;
  Close (ft);
END.
