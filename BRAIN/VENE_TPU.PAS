 Unit VENE_TPU;    {Unit… per Calcolo Venere di PRECISIONE}

 Interface

 Uses Crt;
 Const
       Ko=0.00000001; rad=pi/180;
 Var
       i                : Integer;
       tau,so,L,L0,L1,L2,L3,L4,L5,
              B,B0,B1,B2,B3,B4,
              R,R0,R1,R2,R3,R4,
       teta,beta,Xecl,Yecl,Zecl,EPSdata,zita,Z,
       XS,YS,ZS, XS1,YS1,ZS1, XSD,YSD,ZSD  : Real;

   Procedure VEN(tau:Real;Var L,B,R: Real);
   Procedure CaricamentoDatiLONG(tau:Real; VAR L,B,R:Real);
   Procedure CaricamentoDatiLAT(tau:Real; VAR L,B,R:Real);
   Procedure CaricamentoDatiRAGGIO(tau:Real; VAR L,B,R:Real);
      Function ArcoTANG(num,den:Real): Real;
      Function ArcoSIN(valore:Real): Real;
      Function MODULO(x,y:Real) : Real;

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

      Function MODULO(x,y:Real) : Real;
   Var modu : Real;
   Begin
     modu:=(x/y-INT(x/y))*y;
     MODULO:=modu;
  End;


   Procedure CaricamentoDatiLONG(tau:Real; VAR L,B,R:Real);
Const
       Ku=0.00001;
       l0m=24; l1m=12; l2m=8; l3m=3; l4m=3; l5m=1;
Var
       L0A:           Array[1..l0m] of LongInt;
       L0B,L0C:       Array[1..l0m] of Real;
       L1A:           Array[1..l1m] of Real;
       L1B,L1C:       Array[1..l1m] of Real;
       L2A:           Array[1..l2m] of LongInt;
       L2B,L2C:       Array[1..l2m] of Real;
       L3A:           Array[1..l3m] of Integer;
       L3B,L3C:       Array[1..l3m] of Real;
       L4A:           Array[1..l4m] of Integer;
       L4B,L4C:       Array[1..l4m] of Real;
       L5A:           Array[1..l5m] of Integer;
       L5B,L5C:       Array[1..l5m] of Real;
 begin

    L0A[1]:=317614667;        L0B[1]:=0;              L0C[1]:=0;
    L0A[2]:=  1353968;        L0B[2]:=5.5931332;      L0C[2]:=10213.2855462;
    L0A[3]:=    89892;        L0B[3]:=5.30650;        L0C[3]:=20426.57109;
    L0A[4]:=     5477;        L0B[4]:=4.4163;         L0C[4]:=7860.4194;
    L0A[5]:=     3456;        L0B[5]:=2.6996;         L0C[5]:=11790.6291;
    L0A[6]:=     2372;        L0B[6]:=2.9938;         L0C[6]:=3930.2097;
    L0A[7]:=     1664;        L0B[7]:=4.2502;         L0C[7]:=1577.3435;
    L0A[8]:=     1438;        L0B[8]:=4.1575;         L0C[8]:=9683.5946;
    L0A[9]:=     1317;        L0B[9]:=5.1867;         L0C[9]:=26.2983;
   L0A[10]:=     1201;       L0B[10]:=6.1536;        L0C[10]:=30639.8566;
   L0A[11]:=      769;       L0B[11]:=0.816;         L0C[11]:=9437.763;
   L0A[12]:=      761;       L0B[12]:=1.950;         L0C[12]:=529.691;
   L0A[13]:=      708;       L0B[13]:=1.065;         L0C[13]:=775.523;
   L0A[14]:=      585;       L0B[14]:=3.998;         L0C[14]:=191.448;
   L0A[15]:=      500;       L0B[15]:=4.123;         L0C[15]:=15720.839;
   L0A[16]:=      429;       L0B[16]:=3.586;         L0C[16]:=19367.189;
   L0A[17]:=      327;       L0B[17]:=5.677;         L0C[17]:=5507.553;
   L0A[18]:=      326;       L0B[18]:=4.591;         L0C[18]:=10404.734;
   L0A[19]:=      232;       L0B[19]:=3.163;         L0C[19]:=9153.904;
   L0A[20]:=      180;       L0B[20]:=4.653;         L0C[20]:=1109.379;
   L0A[21]:=      155;       L0B[21]:=5.570;         L0C[21]:=19651.048;
   L0A[22]:=      128;       L0B[22]:=4.226;         L0C[22]:=20.775;
   L0A[23]:=      128;       L0B[23]:=0.962;         L0C[23]:=5661.332;
   L0A[24]:=      106;       L0B[24]:=1.537;         L0C[24]:=801.821;

    so:=0; FOR i:= 1 to l0m do
     begin
     so:=so+Ko*L0A[i]*(COS(L0B[i]+L0C[i]*tau));
     end;  L0:=so;
 {----------------------------------------------------------------------------}
    L1A[1]:=1021352943.053;   L1B[1]:=0;              L1C[1]:=0;
    L1A[2]:=        95.708;   L1B[2]:=2.46424;        L1C[2]:=10213.28555;
    L1A[3]:=        14.445;   L1B[3]:=0.51625;        L1C[3]:=20426.57109;
    L1A[4]:=         0.213;   L1B[4]:=1.795;          L1C[4]:=30639.857;
    L1A[5]:=         0.174;   L1B[5]:=2.655;          L1C[5]:=26.298;
    L1A[6]:=         0.152;   L1B[6]:=6.106;          L1C[6]:=1577.344;
    L1A[7]:=         0.082;   L1B[7]:=5.70;           L1C[7]:=191.45;
    L1A[8]:=         0.070;   L1B[8]:=2.68;           L1C[8]:=9437.76;
    L1A[9]:=         0.052;   L1B[9]:=3.60;           L1C[9]:=775.52;
   L1A[10]:=         0.038;  L1B[10]:=1.03;          L1C[10]:=529.69;
   L1A[11]:=         0.030;  L1B[11]:=1.25;          L1C[11]:=5507.55;
   L1A[12]:=         0.025;  L1B[12]:=6.11;          L1C[12]:=10404.73;

    so:=0; FOR i:= 1 to l1m do
     begin
     so:=so+Ku*L1A[i]*(COS(L1B[i]+L1C[i]*tau));
     end;  L1:=so;
 {----------------------------------------------------------------------------}
    L2A[1]:=54127;            L2B[1]:=0;              L2C[1]:=0;
    L2A[2]:= 3891;            L2B[2]:=0.3451;         L2C[2]:=10213.2855;
    L2A[3]:= 1338;            L2B[3]:=2.0201;         L2C[3]:=20426.5711;
    L2A[4]:=   24;            L2B[4]:=2.05;           L2C[4]:=26.30;
    L2A[5]:=   19;            L2B[5]:=3.54;           L2C[5]:=30639.86;
    L2A[6]:=   10;            L2B[6]:=3.97;           L2C[6]:=775.52;
    L2A[7]:=    7;            L2B[7]:=1.52;           L2C[7]:=1577.34;
    L2A[8]:=    6;            L2B[8]:=1.00;           L2C[8]:=191.45;

    so:=0; FOR i:= 1 to l2m do
     begin
     so:=so+Ko*L2A[i]*(COS(L2B[i]+L2C[i]*tau));
     end;  L2:=so;
 {----------------------------------------------------------------------------}
    L3A[1]:=136;              L3B[1]:=4.804;          L3C[1]:=10213.286;
    L3A[2]:= 78;              L3B[2]:=3.67;           L3C[2]:=20426.57;
    L3A[3]:= 26;              L3B[3]:=0;              L3C[3]:=0;

    so:=0; FOR i:= 1 to l3m do
     begin
     so:=so+Ko*L3A[i]*(COS(L3B[i]+L3C[i]*tau));
     end;  L3:=so;
 {----------------------------------------------------------------------------}
    L4A[1]:=114;              L4B[1]:=3.1416;         L4C[1]:=0;
    L4A[2]:=  3;              L4B[2]:=5.21;           L4C[2]:=20426.57;
    L4A[3]:=  2;              L4B[3]:=2.51;           L4C[3]:=10213.29;

    so:=0; FOR i:= 1 to l4m do
     begin
     so:=so+Ko*L4A[i]*(COS(L4B[i]+L4C[i]*tau));
     end;  L4:=so;
 {----------------------------------------------------------------------------}
    L5A[1]:=1;                L5B[1]:=3.14;           L5C[1]:=0;

    so:=0; FOR i:= 1 to l5m do
     begin
     so:=so+Ko*L5A[i]*(COS(L5B[i]+L5C[i]*tau));
     end;  L5:=so;
 {----------------------------------------------------------------------------}

    L:=((((L5*tau+L4)*tau+L3)*tau+L2)*tau+L1)*tau+L0;
    end;

   Procedure CaricamentoDatiLAT(tau:Real; VAR L,B,R:Real);
 Const
       b0m=9 ; b1m=4; b2m=4 ; b3m=4; b4m=1;
 Var
       B0A:           Array[1..B0m] of LongInt;
       B0B,B0C:       Array[1..B0m] of Real;
       B1A:           Array[1..B1m] of LongInt;
       B1B,B1C:       Array[1..B1m] of Real;
       B2A:           Array[1..B2m] of Integer;
       B2B,B2C:       Array[1..B2m] of Real;
       B3A:           Array[1..B3m] of Integer;
       B3B,B3C:       Array[1..B3m] of Real;
       B4A:           Array[1..B4m] of Integer;
       B4B,B4C:       Array[1..B4m] of Real;

  begin

    B0A[1]:=5923638;          B0B[1]:=0.2670278;      B0C[1]:=10213.2855462;
    B0A[2]:=  40108;          B0B[2]:=1.14737;        B0C[2]:=20426.57109;
    B0A[3]:=  32815;          B0B[3]:=3.14159;        B0C[3]:=0;
    B0A[4]:=   1011;          B0B[4]:=1.0895;         B0C[4]:=30639.8566;
    B0A[5]:=    149;          B0B[5]:=6.254;          B0C[5]:=18073.705;
    B0A[6]:=    138;          B0B[6]:=0.860;          B0C[6]:=1577.344;
    B0A[7]:=    130;          B0B[7]:=3.672;          B0C[7]:=9437.763;
    B0A[8]:=    120;          B0B[8]:=3.705;          B0C[8]:=2352.866;
    B0A[9]:=    108;          B0B[9]:=4.539;          B0C[9]:=22003.915;

    so:=0; FOR i:= 1 to B0m do
     begin
     so:=so+Ko*B0A[i]*(COS(B0B[i]+B0C[i]*tau));
     end;  B0:=so;
 {----------------------------------------------------------------------------}
    B1A[1]:=513348;           B1B[1]:=1.803643;       B1C[1]:=10213.285546;
    B1A[2]:=  4380;           B1B[2]:=3.3862;         B1C[2]:=20426.5711;
    B1A[3]:=   199;           B1B[3]:=0;              B1C[3]:=0;
    B1A[4]:=   197;           B1B[4]:=2.530;          B1C[4]:=30639.857;

    so:=0; FOR i:= 1 to B1m do
     begin
     so:=so+Ko*B1A[i]*(COS(B1B[i]+B1C[i]*tau));
     end;  B1:=so;
 {----------------------------------------------------------------------------}
    B2A[1]:=22378;            B2B[1]:=3.38509;        B2C[1]:=10213.28555;
    B2A[2]:=  282;            B2B[2]:=0;              B2C[2]:=0;
    B2A[3]:=  173;            B2B[3]:=5.256;          B2C[3]:=20426.571;
    B2A[4]:=   27;            B2B[4]:=3.87;           B2C[4]:=30639.857;

    so:=0; FOR i:= 1 to B2m do
     begin
     so:=so+Ko*B2A[i]*(COS(B2B[i]+B2C[i]*tau));
     end;  B2:=so;
 {----------------------------------------------------------------------------}
    B3A[1]:=647;              B3B[1]:=4.992;          B3C[1]:=10213.286;
    B3A[2]:= 20;              B3B[2]:=3.14;           B3C[2]:=0;
    B3A[3]:=  6;              B3B[3]:=0.77;           B3C[3]:=20426.57;
    B3A[4]:=  3;              B3B[4]:=5.44;           B3C[4]:=30639.86;

    so:=0; FOR i:= 1 to B3m do
     begin
     so:=so+Ko*B3A[i]*(COS(B3B[i]+B3C[i]*tau));
     end;  B3:=so;
 {----------------------------------------------------------------------------}
    B4A[1]:=14;               B4B[1]:=0.32;           B4C[1]:=10213.29;

    so:=0; FOR i:= 1 to B4m do
     begin
     so:=so+Ko*B4A[i]*(COS(B4B[i]+B4C[i]*tau));
     end;  B4:=so;
 {----------------------------------------------------------------------------}

    B:=(((B4*tau+B3)*tau+B2)*tau+B1)*tau+B0;
    end;


   Procedure CaricamentoDatiRAGGIO(tau:Real; VAR L,B,R:Real);
 Const
       r0m=12; r1m=3; r2m=3 ; r3m=1; r4m=1;
 Var
       R0A:           Array[1..R0m] of LongInt;
       R0B,R0C:       Array[1..R0m] of Real;
       R1A:           Array[1..R1m] of LongInt;
       R1B,R1C:       Array[1..R1m] of Real;
       R2A:           Array[1..R2m] of Integer;
       R2B,R2C:       Array[1..R2m] of Real;
       R3A:           Array[1..R3m] of Integer;
       R3B,R3C:       Array[1..R3m] of Real;
       R4A:           Array[1..R4m] of Integer;
       R4B,R4C:       Array[1..R4m] of Real;

    begin

    R0A[1]:=72334821;         R0B[1]:=0;              R0C[1]:=0;
    R0A[2]:=  489824;         R0B[2]:=4.021518;       R0C[2]:=10213.285546;
    R0A[3]:=    1658;         R0B[3]:=4.9021;         R0C[3]:=20426.5711;
    R0A[4]:=    1632;         R0B[4]:=2.8455;         R0C[4]:=7860.4194;
    R0A[5]:=    1378;         R0B[5]:=1.1285;         R0C[5]:=11790.6291;
    R0A[6]:=     498;         R0B[6]:=2.587;          R0C[6]:=9683.595;
    R0A[7]:=     374;         R0B[7]:=1.423;          R0C[7]:=3930.210;
    R0A[8]:=     264;         R0B[8]:=5.529;          R0C[8]:=9437.763;
    R0A[9]:=     237;         R0B[9]:=2.551;          R0C[9]:=15720.839;
   R0A[10]:=     222;        R0B[10]:=2.013;         R0C[10]:=19367.189;
   R0A[11]:=     126;        R0B[11]:=2.728;         R0C[11]:=1577.344;
   R0A[12]:=     119;        R0B[12]:=3.020;         R0C[12]:=10404.734;

    so:=0; FOR i:= 1 to R0m do
     begin
     so:=so+Ko*R0A[i]*(COS(R0B[i]+R0C[i]*tau));
     end;  R0:=so;
 {----------------------------------------------------------------------------}
    R1A[1]:=34551;            R1B[1]:=0.89199;        R1C[1]:=10213.28555;
    R1A[2]:=  234;            R1B[2]:=1.772;          R1C[2]:=20426.571;
    R1A[3]:=  234;            R1B[3]:=3.142;          R1C[3]:=0;

    so:=0; FOR i:= 1 to R1m do
     begin
     so:=so+Ko*R1A[i]*(COS(R1B[i]+R1C[i]*tau));
     end;  R1:=so;
 {----------------------------------------------------------------------------}
    R2A[1]:=1407;             R2B[1]:=5.0637;         R2C[1]:=10213.2855;
    R2A[2]:=  16;             R2B[2]:=5.47;           R2C[2]:=20426.57;
    R2A[3]:=  13;             R2B[3]:=0;              R2C[3]:=0;

    so:=0; FOR i:= 1 to R2m do
     begin
     so:=so+Ko*R2A[i]*(COS(R2B[i]+R2C[i]*tau));
     end;  R2:=so;
 {----------------------------------------------------------------------------}
    R3A[1]:=50;               R3B[1]:=3.22;           R3C[1]:=10213.29;

    so:=0; FOR i:= 1 to R3m do
     begin
     so:=so+Ko*R3A[i]*(COS(R3B[i]+R3C[i]*tau));
     end;  R3:=so;
 {----------------------------------------------------------------------------}
    R4A[1]:=1;                R4B[1]:=0.92;           R4C[1]:=10213.29;

    so:=0; FOR i:= 1 to R4m do
     begin
     so:=so+Ko*R4A[i]*(COS(R4B[i]+R4C[i]*tau));
     end;  R4:=so;
 {----------------------------------------------------------------------------}
    R:=(((R4*tau+R3)*tau+R2)*tau+R1)*tau+R0;
    end;

   Procedure VEN(tau:Real;Var L,B,R: Real);
   Begin
      CaricamentoDatiLONG(tau,L,B,R);
      CaricamentoDatiLAT(tau,L,B,R);
      CaricamentoDatiRAGGIO(tau,L,B,R);
    End;

END.
