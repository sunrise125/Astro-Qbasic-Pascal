 Unit MERC_TPU;    {Unit� per Calcolo Mercurio di PRECISIONE}

 Interface

 Uses Crt;
 Const
       Ko=0.00000001; rad=pi/180;
 Var
       i                : Integer;
       tau,so,L,L0,L1,L2,L3,L4,L5,
              B,B0,B1,B2,B3,B4,
              R,R0,R1,R2,R3,
       teta,beta,Xecl,Yecl,Zecl,EPSdata,zita,Z,
       XS,YS,ZS, XS1,YS1,ZS1, XSD,YSD,ZSD  : Real;

   Procedure MER(tau:Real;Var L,B,R: Real);
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
       l0m=38; l1m=16; l2m=10; l3m=8; l4m=6; l5m=1;
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

    L0A[1]:=440250710;        L0B[1]:=0;              L0C[1]:= 0;
    L0A[2]:= 40989415;        L0B[2]:=1.48302034;     L0C[2]:= 26087.90314157;
    L0A[3]:=  5046294;        L0B[3]:=4.4778549;      L0C[3]:= 52175.8062831;
    L0A[4]:=   855347;        L0B[4]:=1.165203;       L0C[4]:= 78263.709425;
    L0A[5]:=   165590;        L0B[5]:=4.119692;       L0C[5]:=104351.612566;
    L0A[6]:=    34562;        L0B[6]:=0.77931;        L0C[6]:=130439.51571;
    L0A[7]:=     7583;        L0B[7]:=3.7135;         L0C[7]:=156527.4188;
    L0A[8]:=     3560;        L0B[8]:=1.5120;         L0C[8]:=  1109.3786;
    L0A[9]:=     1803;        L0B[9]:=4.1033;         L0C[9]:=  5661.3320;
   L0A[10]:=     1726;       L0B[10]:=0.3583;        L0C[10]:=182615.3220;
   L0A[11]:=     1590;       L0B[11]:=2.9951;        L0C[11]:= 25028.5212;
   L0A[12]:=     1365;       L0B[12]:=4.5992;        L0C[12]:= 27197.2817;
   L0A[13]:=     1017;       L0B[13]:=0.8803;        L0C[13]:= 31749.2352;
   L0A[14]:=      714;       L0B[14]:=1.541;         L0C[14]:= 24978.525;
   L0A[15]:=      644;       L0B[15]:=5.303;         L0C[15]:= 21535.950;
   L0A[16]:=      451;       L0B[16]:=6.050;         L0C[16]:= 51116.424;
   L0A[17]:=      404;       L0B[17]:=3.282;         L0C[17]:=208703.225;
   L0A[18]:=      352;       L0B[18]:=5.242;         L0C[18]:= 20426.571;
   L0A[19]:=      345;       L0B[19]:=2.792;         L0C[19]:= 15874.618;
   L0A[20]:=      343;       L0B[20]:=5.765;         L0C[20]:=   955.600;
   L0A[21]:=      339;       L0B[21]:=5.863;         L0C[21]:= 25558.212;
   L0A[22]:=      325;       L0B[22]:=1.337;         L0C[22]:= 53285.185;
   L0A[23]:=      273;       L0B[23]:=2.495;         L0C[23]:=   529.691;
   L0A[24]:=      264;       L0B[24]:=3.917;         L0C[24]:= 57837.138;
   L0A[25]:=      260;       L0B[25]:=0.987;         L0C[25]:=  4551.953;
   L0A[26]:=      239;       L0B[26]:=0.113;         L0C[26]:=  1059.382;
   L0A[27]:=      235;       L0B[27]:=0.267;         L0C[27]:= 11322.664;
   L0A[28]:=      217;       L0B[28]:=0.660;         L0C[28]:= 13521.751;
   L0A[29]:=      209;       L0B[29]:=2.092;         L0C[29]:= 47623.853;
   L0A[30]:=      183;       L0B[30]:=2.629;         L0C[30]:= 27043.503;
   L0A[31]:=      182;       L0B[31]:=2.434;         L0C[31]:= 25661.305;
   L0A[32]:=      176;       L0B[32]:=4.536;         L0C[32]:= 51066.428;
   L0A[33]:=      173;       L0B[33]:=2.452;         L0C[33]:= 24498.830;
   L0A[34]:=      142;       L0B[34]:=3.360;         L0C[34]:= 37410.567;
   L0A[35]:=      138;       L0B[35]:=0.291;         L0C[35]:= 10213.286;
   L0A[36]:=      125;       L0B[36]:=3.721;         L0C[36]:= 39609.655;
   L0A[37]:=      118;       L0B[37]:=2.781;         L0C[37]:= 77204.327;
   L0A[38]:=      106;       L0B[38]:=4.206;         L0C[38]:= 19804.827;
    so:=0; FOR i:= 1 to l0m do
     begin
     so:=so+Ko*L0A[i]*(COS(L0B[i]+L0C[i]*tau));
     end;  L0:=so;
 {----------------------------------------------------------------------------}
    L1A[1]:=2608814706.223;   L1B[1]:=0;              L1C[1]:=     0;
    L1A[2]:=      1126.008;   L1B[2]:=6.2170397;      L1C[2]:= 26087.9031416;
    L1A[3]:=       303.471;   L1B[3]:=3.055655;       L1C[3]:= 52175.806283;
    L1A[4]:=        80.538;   L1B[4]:=6.10455;        L1C[4]:= 78263.70942;
    L1A[5]:=        21.245;   L1B[5]:=2.83532;        L1C[5]:=104351.61257;
    L1A[6]:=         5.592;   L1B[6]:=5.8268;         L1C[6]:=130439.5157;
    L1A[7]:=         1.472;   L1B[7]:=2.5185;         L1C[7]:=156527.4188;
    L1A[8]:=         0.388;   L1B[8]:=5.480;          L1C[8]:=182615.322;
    L1A[9]:=         0.352;   L1B[9]:=3.052;          L1C[9]:=  1109.379;
   L1A[10]:=         0.103;  L1B[10]:=2.149;         L1C[10]:=208703.225;
   L1A[11]:=         0.094;  L1B[11]:=6.12;          L1C[11]:= 27197.28;
   L1A[12]:=         0.091;  L1B[12]:=0.00;          L1C[12]:= 24978.52;
   L1A[13]:=         0.052;  L1B[13]:=5.62;          L1C[13]:=  5661.33;
   L1A[14]:=         0.044;  L1B[14]:=4.57;          L1C[14]:= 25028.52;
   L1A[15]:=         0.028;  L1B[15]:=3.04;          L1C[15]:= 51066.43;
   L1A[16]:=         0.027;  L1B[16]:=5.09;          L1C[16]:=234791.13;

    so:=0; FOR i:= 1 to l1m do
     begin
     so:=so+Ku*L1A[i]*(COS(L1B[i]+L1C[i]*tau));
     end;  L1:=so;
 {----------------------------------------------------------------------------}
    L2A[1]:=53050;            L2B[1]:=0;              L2C[1]:=     0;
    L2A[2]:=16904;            L2B[2]:=4.69072;        L2C[2]:= 26087.90314;
    L2A[3]:= 7397;            L2B[3]:=1.3474;         L2C[3]:= 52175.8063;
    L2A[4]:= 3018;            L2B[4]:=4.4564;         L2C[4]:= 78263.7094;
    L2A[5]:= 1107;            L2B[5]:=1.2623;         L2C[5]:=104351.6126;
    L2A[6]:=  378;            L2B[6]:=4.320;          L2C[6]:=130439.516;
    L2A[7]:=  123;            L2B[7]:=1.069;          L2C[7]:=156527.419;
    L2A[8]:=   39;            L2B[8]:=4.08;           L2C[8]:=182615.32;
    L2A[9]:=   15;            L2B[9]:=4.63;           L2C[9]:=  1109.38;
   L2A[10]:=   12;           L2B[10]:=0.79;          L2C[10]:=208703.23;

    so:=0; FOR i:= 1 to l2m do
     begin
     so:=so+Ko*L2A[i]*(COS(L2B[i]+L2C[i]*tau));
     end;  L2:=so;
 {----------------------------------------------------------------------------}
    L3A[1]:=188;              L3B[1]:=0.035;          L3C[1]:= 52175.806;
    L3A[2]:=142;              L3B[2]:=3.125;          L3C[2]:= 26087.903;
    L3A[3]:= 97;              L3B[3]:=3.00;           L3C[3]:= 78263.71;
    L3A[4]:= 44;              L3B[4]:=6.02;           L3C[4]:=104351.61;
    L3A[5]:= 35;              L3B[5]:=0;              L3C[5]:=     0;
    L3A[6]:= 18;              L3B[6]:=2.78;           L3C[6]:=130439.52;
    L3A[7]:=  7;              L3B[7]:=5.82;           L3C[7]:=156527.42;
    L3A[8]:=  3;              L3B[8]:=2.57;           L3C[8]:=182615.32;

    so:=0; FOR i:= 1 to l3m do
     begin
     so:=so+Ko*L3A[i]*(COS(L3B[i]+L3C[i]*tau));
     end;  L3:=so;
 {----------------------------------------------------------------------------}
    L4A[1]:=114;              L4B[1]:=3.1416;         L4C[1]:=     0;
    L4A[2]:=  3;              L4B[2]:=2.03;           L4C[2]:= 26087.90;
    L4A[3]:=  2;              L4B[3]:=1.42;           L4C[3]:= 78263.71;
    L4A[4]:=  2;              L4B[4]:=4.50;           L4C[4]:= 52175.81;
    L4A[5]:=  1;              L4B[5]:=4.50;           L4C[5]:=104351.61;
    L4A[6]:=  1;              L4B[6]:=1.27;           L4C[6]:=130439.52;

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
       b0m=14 ; b1m=11 ; b2m=9 ; b3m=7; b4m=2;
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

    B0A[1]:=11737529;         B0B[1]:=1.98357499;     B0C[1]:= 26087.90314157;
    B0A[2]:= 2388077;         B0B[2]:=5.0373896;      B0C[2]:= 52175.8062831;
    B0A[3]:= 1222840;         B0B[3]:=3.1415927;      B0C[3]:=     0;
    B0A[4]:=  543252;         B0B[4]:=1.796444;       B0C[4]:= 78263.709425;
    B0A[5]:=  129779;         B0B[5]:=4.832325;       B0C[5]:=104351.612566;
    B0A[6]:=   31867;         B0B[6]:=1.58088;        B0C[6]:=130439.51571;
    B0A[7]:=    7963;         B0B[7]:=4.6097;         B0C[7]:=156527.4188;
    B0A[8]:=    2014;         B0B[8]:=1.3532;         B0C[8]:=182615.3220;
    B0A[9]:=     514;         B0B[9]:=4.378;          B0C[9]:=208703.225;
   B0A[10]:=     209;        B0B[10]:=2.020;         B0C[10]:= 24978.525;
   B0A[11]:=     208;        B0B[11]:=4.918;         B0C[11]:= 27197.282;
   B0A[12]:=     132;        B0B[12]:=1.119;         B0C[12]:=234791.128;
   B0A[13]:=     121;        B0B[13]:=1.813;         B0C[13]:= 53285.185;
   B0A[14]:=     100;        B0B[14]:=5.657;         B0C[14]:= 20426.571;

    so:=0; FOR i:= 1 to B0m do
     begin
     so:=so+Ko*B0A[i]*(COS(B0B[i]+B0C[i]*tau));
     end;  B0:=so;
 {----------------------------------------------------------------------------}
    B1A[1]:=429151;           B1B[1]:=3.501698;       B1C[1]:= 26087.903142;
    B1A[2]:=146234;           B1B[2]:=3.141593;       B1C[2]:=     0;
    B1A[3]:= 22675;           B1B[3]:=0.01515;        B1C[3]:= 52175.80628;
    B1A[4]:= 10895;           B1B[4]:=0.48540;        B1C[4]:= 78263.70942;
    B1A[5]:=  6353;           B1B[5]:=3.4294;         B1C[5]:=104351.6126;
    B1A[6]:=  2496;           B1B[6]:=0.1605;         B1C[6]:=130439.5157;
    B1A[7]:=   860;           B1B[7]:=3.185;          B1C[7]:=156527.419;
    B1A[8]:=   278;           B1B[8]:=6.210;          B1C[8]:=182615.322;
    B1A[9]:=    86;           B1B[9]:=2.95;           B1C[9]:=208703.23;
   B1A[10]:=    28;          B1B[10]:=0.29;          B1C[10]:= 27197.28;
   B1A[11]:=    26;          B1B[11]:=5.98;          B1C[11]:=234791.13;

    so:=0; FOR i:= 1 to B1m do
     begin
     so:=so+Ko*B1A[i]*(COS(B1B[i]+B1C[i]*tau));
     end;  B1:=so;
 {----------------------------------------------------------------------------}
    B2A[1]:=11831;            B2B[1]:=4.79066;        B2C[1]:= 26087.90314;
    B2A[2]:= 1914;            B2B[2]:=0;              B2C[2]:=     0;
    B2A[3]:= 1045;            B2B[3]:=1.2122;         B2C[3]:= 52175.8063;
    B2A[4]:=  266;            B2B[4]:=4.434;          B2C[4]:= 78263.709;
    B2A[5]:=  170;            B2B[5]:=1.623;          B2C[5]:=104351.613;
    B2A[6]:=   96;            B2B[6]:=4.80;           B2C[6]:=130439.52;
    B2A[7]:=   45;            B2B[7]:=1.61;           B2C[7]:=156527.42;
    B2A[8]:=   18;            B2B[8]:=4.67;           B2C[8]:=182615.32;
    B2A[9]:=    7;            B2B[9]:=1.43;           B2C[9]:=208703.23;

    so:=0; FOR i:= 1 to B2m do
     begin
     so:=so+Ko*B2A[i]*(COS(B2B[i]+B2C[i]*tau));
     end;  B2:=so;
 {----------------------------------------------------------------------------}
    B3A[1]:=235;              B3B[1]:=0.354;          B3C[1]:= 26087.903;
    B3A[2]:=161;              B3B[2]:=0;              B3C[2]:=     0;
    B3A[3]:= 19;              B3B[3]:=4.36;           B3C[3]:= 52175.81;
    B3A[4]:=  6;              B3B[4]:=2.51;           B3C[4]:= 78263.71;
    B3A[5]:=  5;              B3B[5]:=6.14;           B3C[5]:=104351.61;
    B3A[6]:=  3;              B3B[6]:=3.12;           B3C[6]:=130439.52;
    B3A[7]:=  2;              B3B[7]:=6.27;           B3C[7]:=156527.42;

    so:=0; FOR i:= 1 to B3m do
     begin
     so:=so+Ko*B3A[i]*(COS(B3B[i]+B3C[i]*tau));
     end;  B3:=so;
 {----------------------------------------------------------------------------}
    B4A[1]:=4;                B4B[1]:=1.75;           B4C[1]:=26087.90;
    B4A[2]:=1;                B4B[2]:=3.14;           B4C[2]:=    0;

    so:=0; FOR i:= 1 to B4m do
     begin
     so:=so+Ko*B4A[i]*(COS(B4B[i]+B4C[i]*tau));
     end;  B4:=so;
 {----------------------------------------------------------------------------}

    B:=(((B4*tau+B3)*tau+B2)*tau+B1)*tau+B0;
    end;


   Procedure CaricamentoDatiRAGGIO(tau:Real; VAR L,B,R:Real);
 Const
       r0m=13; r1m=8; r2m=7 ; r3m=5;
 Var
       R0A:           Array[1..R0m] of LongInt;
       R0B,R0C:       Array[1..R0m] of Real;
       R1A:           Array[1..R1m] of LongInt;
       R1B,R1C:       Array[1..R1m] of Real;
       R2A:           Array[1..R2m] of Integer;
       R2B,R2C:       Array[1..R2m] of Real;
       R3A:           Array[1..R3m] of Integer;
       R3B,R3C:       Array[1..R3m] of Real;

    begin

    R0A[1]:=39528272;         R0B[1]:=0;              R0C[1]:=     0;
    R0A[2]:= 7834132;         R0B[2]:=6.1923372;      R0C[2]:= 26087.9031416;
    R0A[3]:=  795526;         R0B[3]:=2.959897;       R0C[3]:= 52175.806283;
    R0A[4]:=  121282;         R0B[4]:=6.010642;       R0C[4]:= 78263.709425;
    R0A[5]:=   21922;         R0B[5]:=2.77820;        R0C[5]:=104351.61257;
    R0A[6]:=    4354;         R0B[6]:=5.8289;         R0C[6]:=130439.5157;
    R0A[7]:=     918;         R0B[7]:=2.597;          R0C[7]:=156527.419;
    R0A[8]:=     290;         R0B[8]:=1.424;          R0C[8]:= 25028.521;
    R0A[9]:=     260;         R0B[9]:=3.028;          R0C[9]:= 27197.282;
   R0A[10]:=     202;        R0B[10]:=5.647;         R0C[10]:=182615.322;
   R0A[11]:=     201;        R0B[11]:=5.592;         R0C[11]:= 31749.235;
   R0A[12]:=     142;        R0B[12]:=6.253;         R0C[12]:= 24978.525;
   R0A[13]:=     100;        R0B[13]:=3.734;         R0C[13]:= 21535.950;

    so:=0; FOR i:= 1 to R0m do
     begin
     so:=so+Ko*R0A[i]*(COS(R0B[i]+R0C[i]*tau));
     end;  R0:=so;
 {----------------------------------------------------------------------------}
    R1A[1]:=217348;           R1B[1]:=4.656172;       R1C[1]:= 26087.903142;
    R1A[2]:= 44142;           R1B[2]:=1.42386;        R1C[2]:= 52175.80628;
    R1A[3]:= 10094;           R1B[3]:=4.47466;        R1C[3]:= 78263.70942;
    R1A[4]:=  2433;           R1B[4]:=1.2423;         R1C[4]:=104351.6126;
    R1A[5]:=  1624;           R1B[5]:=0;              R1C[5]:=     0;
    R1A[6]:=   604;           R1B[6]:=4.293;          R1C[6]:=130439.516;
    R1A[7]:=   153;           R1B[7]:=1.061;          R1C[7]:=156527.419;
    R1A[8]:=    39;           R1B[8]:=4.11;           R1C[8]:=182615.32;

    so:=0; FOR i:= 1 to R1m do
     begin
     so:=so+Ko*R1A[i]*(COS(R1B[i]+R1C[i]*tau));
     end;  R1:=so;
 {----------------------------------------------------------------------------}
    R2A[1]:=3118;             R2B[1]:=3.0823;         R2C[1]:= 26087.9031;
    R2A[2]:=1245;             R2B[2]:=6.1518;         R2C[2]:= 52175.8063;
    R2A[3]:= 425;             R2B[3]:=2.926;          R2C[3]:= 78263.709;
    R2A[4]:= 136;             R2B[4]:=5.980;          R2C[4]:=104351.613;
    R2A[5]:=  42;             R2B[5]:=2.75;           R2C[5]:=130439.52;
    R2A[6]:=  22;             R2B[6]:=3.14;           R2C[6]:=     0;
    R2A[7]:=  13;             R2B[7]:=5.80;           R2C[7]:=156527.42;

    so:=0; FOR i:= 1 to R2m do
     begin
     so:=so+Ko*R2A[i]*(COS(R2B[i]+R2C[i]*tau));
     end;  R2:=so;
 {----------------------------------------------------------------------------}
    R3A[1]:=33;               R3B[1]:=1.68;           R3C[1]:= 26087.90;
    R3A[2]:=24;               R3B[2]:=4.63;           R3C[2]:= 52175.81;
    R3A[3]:=12;               R3B[3]:=1.39;           R3C[3]:= 78263.71;
    R3A[4]:= 5;               R3B[4]:=4.44;           R3C[4]:=104351.61;
    R3A[5]:= 2;               R3B[5]:=1.21;           R3C[5]:=130439.52;

    so:=0; FOR i:= 1 to R3m do
     begin
     so:=so+Ko*R3A[i]*(COS(R3B[i]+R3C[i]*tau));
     end;  R3:=so;
 {----------------------------------------------------------------------------}
    R:=((R3*tau+R2)*tau+R1)*tau+R0;
    end;


   Procedure MER(tau:Real;Var L,B,R: Real);
   Begin
      CaricamentoDatiLONG(tau,L,B,R);
      CaricamentoDatiLAT(tau,L,B,R);
      CaricamentoDatiRAGGIO(tau,L,B,R);
    End;

END.
