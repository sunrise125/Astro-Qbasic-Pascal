 Unit SOLEXYZ;    {Unit… per Calcolo SOLE di PRECISIONE}

 Interface

 Uses Crt;
 Const
       Ko=0.00000001; JD2000=2451545;
       rad=pi/180;
 Var
       i                : Integer;
       tau,so,L,L0,L1,L2,L3,L4,
       jd,    B,B0,B1,B2,B3,B4,
              R,R0,R1,R2,R3,R4,
       teta,beta,Xecl,Yecl,Zecl,EPSdata,zita,Z,
       XS,YS,ZS, XS1,YS1,ZS1, XSD,YSD,ZSD  : Real;


   Procedure SUN(tau,jd:Real; VAR Xso,Yso,Zso: Real);
   Procedure CaricamentoDatiLONG(tau:Real; VAR L,B,R:Real);
   Procedure CaricamentoDatiLAT(tau:Real; VAR L,B,R:Real);
   Procedure CaricamentoDatiRAGGIO(tau:Real; VAR L,B,R:Real);
   Procedure CoordinateECLITTICHE;
   Procedure XYZSole(tau: Real);
   Procedure XYZPrecessione;
   Procedure CalcoloEPSData;
      Function ArcoTANG(num,den:Real): Real;
      Function ArcoSIN(valore:Real): Real;

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

   Procedure CaricamentoDatiLONG(tau:Real; VAR L,B,R:Real);
Const
       Ku=0.00001;
       l0m=64; l1m=34; l2m=20; l3m=7; l4m=2;
Var
       L0A:           Array[1..l0m] of LongInt;
       L0B,L0C:       Array[1..l0m] of Real;
       L1A:           Array[1..l1m] of Real;
       L1B,L1C:       Array[1..l1m] of Real;
       L2A:           Array[1..l2m] of Integer;
       L2B,L2C:       Array[1..l2m] of Real;
       L3A:           Array[1..l3m] of Integer;
       L3B,L3C:       Array[1..l3m] of Real;
       L4A:           Array[1..l4m] of Integer;
       L4B,L4C:       Array[1..l4m] of Real;

 begin
  L0A[1]:=175347046; L0A[2]:=3341656; L0A[3]:=34894; L0A[4]:=3497; L0A[5]:=3418;
   L0A[6]:=3136; L0A[7]:=2676; L0A[8]:=2343; L0A[9]:=1324; L0A[10]:=1273;
    L0A[11]:=1199; L0A[12]:=990; L0A[13]:=902; L0A[14]:=857; L0A[15]:=780;
     L0A[16]:=753; L0A[17]:=505; L0A[18]:=492; L0A[19]:=357; L0A[20]:=317;
      L0A[21]:=284; L0A[22]:=271; L0A[23]:=243; L0A[24]:=206; L0A[25]:=205;
       L0A[26]:=202; L0A[27]:=156; L0A[28]:=132; L0A[29]:=126; L0A[30]:=115;
        L0A[31]:=103; L0A[32]:=102; L0A[33]:=102; L0A[34]:=99; L0A[35]:=98;
        L0A[36]:=86; L0A[37]:=85; L0A[38]:=85; L0A[39]:=80; L0A[40]:=79;
       L0A[41]:=75; L0A[42]:=74; L0A[43]:=74; L0A[44]:=70; L0A[45]:=62;
      L0A[46]:=61; L0A[47]:=57; L0A[48]:=56; L0A[49]:=56; L0A[50]:=52;
     L0A[51]:=52; L0A[52]:=51; L0A[53]:=49; L0A[54]:=41; L0A[55]:=41;
    L0A[56]:=39; L0A[57]:=37; L0A[58]:=37; L0A[59]:=36; L0A[60]:=36;
    L0A[61]:=33; L0A[62]:=30; L0A[63]:=30; L0A[64]:=25;
  L0B[1]:=0; L0B[2]:=4.6692568; L0B[3]:=4.6261; L0B[4]:=2.7441; L0B[5]:=2.8289;
   L0B[6]:=3.6277; L0B[7]:=4.4181; L0B[8]:=6.1352; L0B[9]:=0.7425; L0B[10]:=2.0371;
    L0B[11]:=1.1096; L0B[12]:=5.233; L0B[13]:=2.045; L0B[14]:=3.508; L0B[15]:=1.179;
     L0B[16]:=2.533; L0B[17]:=4.583; L0B[18]:=4.205; L0B[19]:=2.920; L0B[20]:=5.849;
      L0B[21]:=1.899; L0B[22]:=0.315; L0B[23]:=0.345; L0B[24]:=4.806; L0B[25]:=1.869;
       L0B[26]:=2.458; L0B[27]:=0.833; L0B[28]:=3.411; L0B[29]:=1.083; L0B[30]:=0.645;
        L0B[31]:=0.636; L0B[32]:=0.976; L0B[33]:=4.267; L0B[34]:=6.21; L0B[35]:=0.68;
        L0B[36]:=5.98; L0B[37]:=1.30; L0B[38]:=3.67; L0B[39]:=1.81; L0B[40]:=3.04;
       L0B[41]:=1.76; L0B[42]:=3.50; L0B[43]:=4.68; L0B[44]:=0.83; L0B[45]:=3.98;
      L0B[46]:=1.82; L0B[47]:=2.78; L0B[48]:=4.39; L0B[49]:=3.47; L0B[50]:=0.19;
     L0B[51]:=1.33; L0B[52]:=0.28; L0B[53]:=0.49; L0B[54]:=5.37; L0B[55]:=2.40;
    L0B[56]:=6.17; L0B[57]:=6.04; L0B[58]:=2.57; L0B[59]:=1.71; L0B[60]:=1.78;
    L0B[61]:=0.59; L0B[62]:=0.44; L0B[63]:=2.74; L0B[64]:=3.16;
  L0C[1]:=0; L0C[2]:=6283.07585; L0C[3]:=12566.1517; L0C[4]:=5753.3849; L0C[5]:=3.5231;
   L0C[6]:=77713.7715; L0C[7]:=7860.4194; L0C[8]:=3930.2097; L0C[9]:=11506.7698; L0C[10]:=529.691;
     L0C[11]:=1577.3435; L0C[12]:=5884.927; L0C[13]:=26.298; L0C[14]:=398.149; L0C[15]:=5223.694;
      L0C[16]:=5507.553; L0C[17]:=18849.228; L0C[18]:=775.523; L0C[19]:=0.067; L0C[20]:=11790.629;
       L0C[21]:=796.298; L0C[22]:=10977.079; L0C[23]:=5486.778; L0C[24]:=2544.314; L0C[25]:=5573.143;
        L0C[26]:=6069.777; L0C[27]:=213.299; L0C[28]:=2942.463; L0C[29]:=20.775; L0C[30]:=0.980;
         L0C[31]:=4694.003; L0C[32]:=15720.839; L0C[33]:=7.114; L0C[34]:=2146.17; L0C[35]:=155.42;
         L0C[36]:=161000.69; L0C[37]:=6275.96; L0C[38]:=71430.70; L0C[39]:=17260.15; L0C[40]:=12036.46;
        L0C[41]:=5088.63; L0C[42]:=3154.69; L0C[43]:=801.82; L0C[44]:=9437.76; L0C[45]:=8827.39;
       L0C[46]:=7084.90; L0C[47]:=6286.60; L0C[48]:=14143.50; L0C[49]:=6279.55; L0C[50]:=12139.55;
      L0C[51]:=1748.02; L0C[52]:=5856.48; L0C[53]:=1194.45; L0C[54]:=8429.24; L0C[55]:=19651.05;
     L0C[56]:=10447.39; L0C[57]:=10213.29; L0C[58]:=1059.38; L0C[59]:=2352.87; L0C[60]:=6812.77;
     L0C[61]:=17789.85; L0C[62]:=83996.85; L0C[63]:=1349.87; L0C[64]:=4690.48;

    so:=0; FOR i:= 1 to l0m do
     begin
     so:=so+Ko*L0A[i]*(COS(L0B[i]+L0C[i]*tau));
     end;  L0:=so;
 {----------------------------------------------------------------------------}
  L1A[1]:=628307584.999; L1A[2]:=206.059; L1A[3]:=4.303; L1A[4]:=0.425; L1A[5]:=0.119;
   L1A[6]:=0.109; L1A[7]:=0.093; L1A[8]:=0.072; L1A[9]:=0.068; L1A[10]:=0.067;
    L1A[11]:=0.059; L1A[12]:=0.056; L1A[13]:=0.045; L1A[14]:=0.036; L1A[15]:=0.029;
     L1A[16]:=0.021; L1A[17]:=0.019; L1A[18]:=0.019; L1A[19]:=0.017; L1A[20]:=0.016;
      L1A[21]:=0.016; L1A[22]:=0.015; L1A[23]:=0.012; L1A[24]:=0.012; L1A[25]:=0.012;
       L1A[26]:=0.012; L1A[27]:=0.011; L1A[28]:=0.010; L1A[29]:=0.010; L1A[30]:=0.009;
        L1A[31]:=0.009; L1A[32]:=0.008; L1A[33]:=0.006; L1A[34]:=0.006;
  L1B[1]:=0; L1B[2]:=2.678235; L1B[3]:=2.6351; L1B[4]:=1.590; L1B[5]:=5.796;
   L1B[6]:=2.966; L1B[7]:=2.59; L1B[8]:=1.14; L1B[9]:=1.87; L1B[10]:=4.41;
    L1B[11]:=2.89; L1B[12]:=2.17; L1B[13]:=0.40; L1B[14]:=0.47; L1B[15]:=2.65;
     L1B[16]:=5.34; L1B[17]:=1.85; L1B[18]:=4.97; L1B[19]:=2.99; L1B[20]:=0.03;
      L1B[21]:=1.43; L1B[22]:=1.21; L1B[23]:=2.83; L1B[24]:=3.26; L1B[25]:=5.27;
       L1B[26]:=2.08; L1B[27]:=0.77; L1B[28]:=1.30; L1B[29]:=4.24; L1B[30]:=2.70;
        L1B[31]:=5.64; L1B[32]:=5.30; L1B[33]:=2.65; L1B[34]:=4.67;
  L1C[1]:=0; L1C[2]:=6283.07585; L1C[3]:=12566.1517; L1C[4]:=3.523; L1C[5]:=26.298;
   L1C[6]:=1577.344; L1C[7]:=18849.23; L1C[8]:=529.69; L1C[9]:=398.15; L1C[10]:=5507.55;
    L1C[11]:=5223.69; L1C[12]:=155.42; L1C[13]:=796.30; L1C[14]:=775.52; L1C[15]:=7.11;
     L1C[16]:=0.98; L1C[17]:=5486.78; L1C[18]:=213.30; L1C[19]:=6275.96; L1C[20]:=2544.31;
      L1C[21]:=2146.17; L1C[22]:=10977.08; L1C[23]:=1748.02; L1C[24]:=5088.63; L1C[25]:=1194.45;
       L1C[26]:=4694.00; L1C[27]:=553.57; L1C[28]:=6286.60; L1C[29]:=1349.87; L1C[30]:=242.73;
        L1C[31]:=951.72; L1C[32]:=2352.87; L1C[33]:=9437.76; L1C[34]:=4690.48;

    so:=0; FOR i:= 1 to l1m do
     begin
     so:=so+Ku*L1A[i]*(COS(L1B[i]+L1C[i]*tau));
     end;  L1:=so;
 {----------------------------------------------------------------------------}
  L2A[1]:=8722; L2A[2]:=991; L2A[3]:=295; L2A[4]:=27; L2A[5]:=16;
   L2A[6]:=16; L2A[7]:=9; L2A[8]:=9; L2A[9]:=7; L2A[10]:=5;
    L2A[11]:=4; L2A[12]:=4; L2A[13]:=3; L2A[14]:=3; L2A[15]:=3;
     L2A[16]:=3; L2A[17]:=3; L2A[18]:=3; L2A[19]:=2; L2A[20]:=2;
  L2B[1]:=1.0725; L2B[2]:=3.1416; L2B[3]:=0.437; L2B[4]:=0.05; L2B[5]:=5.19;
   L2B[6]:=3.69; L2B[7]:=0.30; L2B[8]:=2.06; L2B[9]:=0.83; L2B[10]:=4.66;
    L2B[11]:=1.03; L2B[12]:=3.44; L2B[13]:=5.14; L2B[14]:=6.05; L2B[15]:=1.19;
     L2B[16]:=6.12; L2B[17]:=0.30; L2B[18]:=2.28; L2B[19]:=4.38; L2B[20]:=3.75;
  L2C[1]:=6283.0758; L2C[2]:=0; L2C[3]:=12566.152; L2C[4]:=3.52; L2C[5]:=26.30;
   L2C[6]:=155.42; L2C[7]:=18849.23; L2C[8]:=77713.77; L2C[9]:=775.52; L2C[10]:=1577.34;
    L2C[11]:=7.11; L2C[12]:=5573.14; L2C[13]:=796.30; L2C[14]:=5507.55; L2C[15]:=242.73;
     L2C[16]:=529.69; L2C[17]:=398.15; L2C[18]:=553.57; L2C[19]:=5223.69; L2C[20]:=0.98;

    so:=0; FOR i:= 1 to l2m do
     begin
     so:=so+Ko*L2A[i]*(COS(L2B[i]+L2C[i]*tau));
     end;  L2:=so;
 {----------------------------------------------------------------------------}
  L3A[1]:=289; L3A[2]:=21; L3A[3]:=3; L3A[4]:=3; L3A[5]:=1;
   L3A[6]:=1; L3A[7]:=1;
  L3B[1]:=5.842; L3B[2]:=6.05; L3B[3]:=5.20; L3B[4]:=3.14; L3B[5]:=4.72;
   L3B[6]:=5.97; L3B[7]:=5.54;
  L3C[1]:=6283.076; L3C[2]:=12566.15; L3C[3]:=155.42; L3C[4]:=0; L3C[5]:=3.52;
   L3C[6]:=242.73; L3C[7]:=18849.23;

    so:=0; FOR i:= 1 to l3m do
     begin
     so:=so+Ko*L3A[i]*(COS(L3B[i]+L3C[i]*tau));
     end;  L3:=so;
 {----------------------------------------------------------------------------}

  L4A[1]:=289; L4A[2]:=21;
  L4B[1]:=5.842; L4B[2]:=6.05;
  L4C[1]:=6283.076; L4C[2]:=12566.15;

    so:=0; FOR i:= 1 to l4m do
     begin
     so:=so+Ko*L4A[i]*(COS(L4B[i]+L4C[i]*tau));
     end;  L4:=so;
 {----------------------------------------------------------------------------}

    L:=(((L4*tau+L3)*tau+L2)*tau+L1)*tau+L0;
    end;

   Procedure CaricamentoDatiLAT(tau:Real; VAR L,B,R:Real);
 Const
       b0m=5 ; b1m=7 ; b2m=4 ; b3m=3; b4m=2;
 Var
       B0A:           Array[1..B0m] of Integer;
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
  B0A[1]:=280; B0A[2]:=102; B0A[3]:=80; B0A[4]:=44; B0A[5]:=32;
  B0B[1]:=3.199; B0B[2]:=5.422; B0B[3]:=3.88; B0B[4]:=3.70; B0B[5]:=4.00;
  B0C[1]:=84334.662; B0C[2]:=5507.553; B0C[3]:=5223.69; B0C[4]:=2352.87; B0C[5]:=1577.34;
    so:=0; FOR i:= 1 to B0m do
     begin
     so:=so+Ko*B0A[i]*(COS(B0B[i]+B0C[i]*tau));
     end;  B0:=so;
 {----------------------------------------------------------------------------}
  B1A[1]:=227778; B1A[2]:=3806; B1A[3]:=3620; B1A[4]:=72; B1A[5]:=8; B1A[6]:=8; B1A[7]:=6;
  B1B[1]:=3.413766; B1B[2]:=3.3706; B1B[3]:=0; B1B[4]:=3.33; B1B[5]:=3.89; B1B[6]:=1.79; B1B[7]:=5.20;
  B1C[1]:=6283.075850; B1C[2]:=12566.1517; B1C[3]:=0; B1C[4]:=18849.23; B1C[5]:=5507.55;
   B1C[6]:=5223.69; B1C[7]:=2352.87;
    so:=0; FOR i:= 1 to B1m do
     begin
     so:=so+Ko*B1A[i]*(COS(B1B[i]+B1C[i]*tau));
     end;  B1:=so;
 {----------------------------------------------------------------------------}
  B2A[1]:=9721; B2A[2]:=233; B2A[3]:=134; B2A[4]:=7;
  B2B[1]:=5.1519; B2B[2]:=3.1416; B2B[3]:=0.644; B2B[4]:=1.07;
  B2C[1]:=6283.07585; B2C[2]:=0; B2C[3]:=12566.152; B2C[4]:=18849.23;
    so:=0; FOR i:= 1 to B2m do
     begin
     so:=so+Ko*B2A[i]*(COS(B2B[i]+B2C[i]*tau));
     end;  B2:=so;
 {----------------------------------------------------------------------------}
  B3A[1]:=276; B3A[2]:=17; B3A[3]:=4;
  B3B[1]:=0.595; B3B[2]:=3.14; B3B[3]:=0.12;
  B3C[1]:=6283.076; B3C[2]:=0; B3C[3]:=12566.15;
    so:=0; FOR i:= 1 to B3m do
     begin
     so:=so+Ko*B3A[i]*(COS(B3B[i]+B3C[i]*tau));
     end;  B3:=so;
 {----------------------------------------------------------------------------}
  B4A[1]:=6; B4A[2]:=1;
  B4B[1]:=2.27; B4B[2]:=0;
  B4C[1]:=6283.08; B4C[2]:=0;
    so:=0; FOR i:= 1 to B4m do
     begin
     so:=so+Ko*B4A[i]*(COS(B4B[i]+B4C[i]*tau));
     end;  B4:=so;
 {----------------------------------------------------------------------------}

    B:=(((B4*tau+B3)*tau+B2)*tau+B1)*tau+B0;
    end;


   Procedure CaricamentoDatiRAGGIO(tau:Real; VAR L,B,R:Real);
 Const
       r0m=40; r1m=10; r2m=6 ; r3m=2; r4m=1;
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
  R0A[1]:=100013989; R0A[2]:=1670700; R0A[3]:=13956; R0A[4]:=3084; R0A[5]:=1628;
   R0A[6]:=1576; R0A[7]:=925; R0A[8]:=542; R0A[9]:=472; R0A[10]:=346;
    R0A[11]:=329; R0A[12]:=307; R0A[13]:=243; R0A[14]:=212; R0A[15]:=186;
     R0A[16]:=175; R0A[17]:=110; R0A[18]:=98; R0A[19]:=86; R0A[20]:=86;
      R0A[21]:=65; R0A[22]:=63; R0A[23]:=57; R0A[24]:=56; R0A[25]:=49;
       R0A[26]:=47; R0A[27]:=45; R0A[28]:=43; R0A[29]:=39; R0A[30]:=38;
        R0A[31]:=37; R0A[32]:=37; R0A[33]:=36; R0A[34]:=35; R0A[35]:=33;
        R0A[36]:=32; R0A[37]:=32; R0A[38]:=28; R0A[39]:=28; R0A[40]:=26;
  R0B[1]:=0; R0B[2]:=3.0984635; R0B[3]:=3.05525; R0B[4]:=5.1985; R0B[5]:=1.1739;
   R0B[6]:=2.8469; R0B[7]:=5.453; R0B[8]:=4.564; R0B[9]:=3.661; R0B[10]:=0.964;
    R0B[11]:=5.900; R0B[12]:=0.299; R0B[13]:=4.273; R0B[14]:=5.847; R0B[15]:=5.022;
     R0B[16]:=3.012; R0B[17]:=5.055; R0B[18]:=0.89; R0B[19]:=5.69; R0B[20]:=1.27;
      R0B[21]:=0.27; R0B[22]:=0.92; R0B[23]:=2.01; R0B[24]:=5.24; R0B[25]:=3.25;
       R0B[26]:=2.58; R0B[27]:=5.54; R0B[28]:=6.01; R0B[29]:=5.36; R0B[30]:=2.39;
        R0B[31]:=0.83; R0B[32]:=4.90; R0B[33]:=1.67; R0B[34]:=1.84; R0B[35]:=0.24;
        R0B[36]:=0.18; R0B[37]:=1.78; R0B[38]:=1.21; R0B[39]:=1.90; R0B[40]:=4.59;
  R0C[1]:=0; R0C[2]:=6283.07585; R0C[3]:=12566.1517; R0C[4]:=77713.7715; R0C[5]:=5753.3849;
   R0C[6]:=7860.4194; R0C[7]:=11506.770; R0C[8]:=3930.210; R0C[9]:=5884.927; R0C[10]:=5507.553;
     R0C[11]:=5223.694; R0C[12]:=5573.143; R0C[13]:=11790.629; R0C[14]:=1577.344; R0C[15]:=10977.079;
      R0C[16]:=18849.228; R0C[17]:=5486.778; R0C[18]:=6069.78; R0C[19]:=15720.84; R0C[20]:=161000.69;
       R0C[21]:=17260.15; R0C[22]:=529.69; R0C[23]:=83996.85; R0C[24]:=71430.70; R0C[25]:=2544.31;
        R0C[26]:=775.52; R0C[27]:=9437.76; R0C[28]:=6275.96; R0C[29]:=4694.00; R0C[30]:=8827.39;
         R0C[31]:=19651.05; R0C[32]:=12139.55; R0C[33]:=12036.46; R0C[34]:=2942.46; R0C[35]:=7084.90;
         R0C[36]:=5088.63; R0C[37]:=398.15; R0C[38]:=6286.60; R0C[39]:=6279.55; R0C[40]:=10447.39;

    so:=0; FOR i:= 1 to R0m do
     begin
     so:=so+Ko*R0A[i]*(COS(R0B[i]+R0C[i]*tau));
     end;  R0:=so;
 {----------------------------------------------------------------------------}
  R1A[1]:=103019; R1A[2]:=1721; R1A[3]:=702; R1A[4]:=32; R1A[5]:=31;
   R1A[6]:=25; R1A[7]:=18; R1A[8]:=10; R1A[9]:=9; R1A[10]:=9;
  R1B[1]:=1.107490; R1B[2]:=1.0644; R1B[3]:=3.142; R1B[4]:=1.02; R1B[5]:=2.84;
   R1B[6]:=1.32; R1B[7]:=1.42; R1B[8]:=5.91; R1B[9]:=1.42; R1B[10]:=0.97;
  R1C[1]:=6283.07585; R1C[2]:=12566.1517; R1C[3]:=0; R1C[4]:=18849.23; R1C[5]:=5507.55;
   R1C[6]:=5223.69; R1C[7]:=1577.34; R1C[8]:=10977.08; R1C[9]:=6275.96; R1C[10]:=5486.78;

    so:=0; FOR i:= 1 to R1m do
     begin
     so:=so+Ko*R1A[i]*(COS(R1B[i]+R1C[i]*tau));
     end;  R1:=so;
 {----------------------------------------------------------------------------}
  R2A[1]:=4359; R2A[2]:=124; R2A[3]:=12; R2A[4]:=9; R2A[5]:=6; R2A[6]:=3;
  R2B[1]:=5.7846; R2B[2]:=5.579; R2B[3]:=3.14; R2B[4]:=3.63; R2B[5]:=1.87; R2B[6]:=5.47;
  R2C[1]:=6283.0758; R2C[2]:=12566.152; R2C[3]:=0; R2C[4]:=77713.77; R2C[5]:=5573.14;
   R2C[6]:=18849.23;

    so:=0; FOR i:= 1 to R2m do
     begin
     so:=so+Ko*R2A[i]*(COS(R2B[i]+R2C[i]*tau));
     end;  R2:=so;
 {----------------------------------------------------------------------------}
  R3A[1]:=145; R3A[2]:=7;
  R3B[1]:=4.273; R3B[2]:=3.92;
  R3C[1]:=6283.076; R3C[2]:=12566.15;

    so:=0; FOR i:= 1 to R3m do
     begin
     so:=so+Ko*R3A[i]*(COS(R3B[i]+R3C[i]*tau));
     end;  R3:=so;
 {----------------------------------------------------------------------------}
  R4A[1]:=4;
  R4B[1]:=2.56;
  R4C[1]:=6283.08;

    so:=0; FOR i:= 1 to R4m do
     begin
     so:=so+Ko*R4A[i]*(COS(R4B[i]+R4C[i]*tau));
     end;  R4:=so;
 {----------------------------------------------------------------------------}
    R:=(((R4*tau+R3)*tau+R2)*tau+R1)*tau+R0;
    end;

   Procedure CoordinateECLITTICHE;
    begin
      teta:=L+pi; beta:=-B;
        Xecl:=R*COS(beta)*COS(teta);
        Yecl:=R*COS(beta)*SIN(teta);
        Zecl:=R*SIN(beta);
  writeln;
  writeln('              Coordinate ECLITTICHE del SOLE (Equin.DINAMICO J2000.0)');
  textcolor(12);
  writeln('        X_ecl= ',Xecl:12:8,'   Y_ecl= ',Yecl:12:8,'   Z_ecl= ',Zecl:12:8);
  textcolor(14);
    end;


   Procedure XYZSole(tau: Real);
    begin
     XS:=Xecl+0.000000440360*Yecl-0.000000190919*Zecl;
     YS:=-0.000000479966*Xecl+0.917482137087*Yecl-0.397776982902*Zecl;
     ZS:=0.397776982902*Yecl+0.917482137087*Zecl;
  writeln;
  writeln('        L = ',L:15:8,' rad    ç = ',tau,' millenni');
  writeln('        B = ',B:15:8,' rad');
  writeln('        R = ',R:15:8,' UA');
  writeln('                Coordinate EQUATORIALI del SOLE (Equin. FK5 J2000.0)');
  textcolor(10);
  writeln(' Equinozio J2000.0 --> XS= ',XS:12:8,'   YS= ',YS:12:8,'   ZS= ',ZS:12:8);
  textcolor(14);

     XS1:=0.999925702634*Xecl+0.012189716217*Yecl+0.000011134016*Zecl;
     YS1:=-0.011179418036*Xecl+0.917413998946*Yecl-0.397777041885*Zecl;
     ZS1:=-0.004859003787*Xecl+0.397747363646*Yecl+0.917482111428*Zecl;
  textcolor(0);
  writeln(' Equinozio B1950.0 --> XS= ',XS1:12:8,'   YS= ',YS1:12:8,'   ZS= ',ZS1:12:8);
  textcolor(14);

    end;

   Procedure XYZPrecessione;
{
   Gli ALGORITMI dei parametri ausiliari di questo Modulo, sono stati ricavati
   da ASTROCOMPUTER Progr. 15 (Enciclop. del Gruppo Fabbri,Etas,Sonzogno-1992)

        Nota:   zita, Z, Teta ð J.Meeus p.126 Astronomical Algorithms
                    (questo procedimento Super-Ok)
  ============================================================================
}
  Var
       T0,T,zz0,
       C1,C2,C3,S1,S2,S3,
       M11,M12,M13,  M21,M22,M23,  M31,M32,M33: Real;

    Begin
   T0:=(jd-jd2000)/36525;
   T:=-T0;
          { Grandezze espresse in Arco-secondi }
  zz0:= (2306.2181+(1.39656-0.000139*T0)*T0)*T;
 zita:= zz0+(0.017998*T+(0.30188-0.000345*T0))*T*T;
    Z:= zz0+(0.018203*T+(1.09468+0.000066*T0))*T*T;
 teta:=((-0.041833*T-(0.42665+0.000217*T0))*T+(2004.3109-(0.8533+0.000217*T0)*T0))*T;
        { Trasformazione in Radianti }
 zita:= (zita/3600)*rad;
    Z:= (Z/3600)*rad;
 teta:= (teta/3600)*rad;

 C1:= COS(zita); C2:= COS(Z); C3:= COS(teta);
 S1:= SIN(zita); S2:= SIN(Z); S3:= SIN(teta);

 M11:=C1*C2*C3-S1*S2; M12:=S1*C2+C1*S2*C3; M13:=C1*S3;
 M21:=-C1*S2-S1*C2*C3; M22:=C1*C2-S1*S2*C3;  M23:=-S1*S3;
 M31:=-C2*S3;  M32:=-S2*S3; M33:=C3;

 XSD:= M11 * XS + M12 * YS + M13 * ZS;
 YSD:= M21 * XS + M22 * YS + M23 * ZS;
 ZSD:= M31 * XS + M32 * YS + M33 * ZS;
    end;

   Procedure CalcoloEPSData;
  Var
      azero,dzero, alfa,delt,
      A,B,C      :Real;

   Begin
  writeln(' Equinozio Data    -->XSD= ',XSD:12:8,'  YSD= ',YSD:12:8,'  ZSD= ',ZSD:12:8);
      EPSdata:=arctan(ZSD/YSD);
{  textcolor(10);
  writeln('                     î-data= ',EPSdata/rad:12:8,'ø approx.4^decimale');
  textcolor(14); }

  azero:=ArcoTANG(YS,XS);
  dzero:=ArcoSIN(ZS/R);
   A:=COS(dzero)*SIN(azero+zita);
   B:=COS(teta)*cos(dzero)*COS(azero+zita)-SIN(teta)*SIN(dzero);
   C:=SIN(teta)*cos(dzero)*COS(azero+zita)+COS(teta)*SIN(dzero);
     alfa:=ArcoTANG(A,B)+Z;
     delt:=ArcoSIN(C);
    EPSdata:=arctan(SIN(delt)/(SIN(alfa)*COS(delt)));
{  textcolor(12);
  writeln('                     î-data= ',EPSdata/rad:12:8,'ø  da controllare');
  textcolor(14);}
   end;

   Procedure SUN(tau,jd:Real; VAR Xso,Yso,Zso: Real);
   Begin
      CaricamentoDatiLONG(tau,L,B,R);
      CaricamentoDatiLAT(tau,L,B,R);
      CaricamentoDatiRAGGIO(tau,L,B,R);
      CoordinateECLITTICHE;
      XYZSole(tau);
      XYZPrecessione;
      CalcoloEPSData;

    End;

END.
