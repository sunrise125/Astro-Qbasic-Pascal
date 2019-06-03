Program TEF94; { --------- T.E.F. = TELAI piani agli ELEMENTI FINITI ---------
                              Aggiornamento: Giu.94
                   Passaggio dal Qbasic al Pascal: NOV.96
               Inserire nella Dir operativa il file Font: TRIP.CHR }

  Uses Crt,Graph;

    Const
         rad=pi/180;

   Var
       E,maxY,maxZ                     : Real;
       k,z,i,j,r,n,n1,n2,t,s,a,qq,h,l,f,k0,k1,k2,q,m,xc,yc,x,y    : Word;
       wo                        : String;
       ch,opz                    : Char;
       ft,fl                     : Text;

      { Dim Vett: n1=30; n2=12; t=60; y=3*30=90; qq=(90ù91/2=4095
                  N.Lib. N.vin. Travi   3*n1         3*n1*(3*n1+1)/2 }

       nu                     : Array[1..30,1..2] of Real;
       yu                     : Array[1..90] of Real;
       tu                     : Array[1..60,1..6] of Real;
    n1u,fa,ha,ma,fa1,ha1,ma1  : Array[1..30,1..2] of Real;
       au,cu,bu               : Array[1..30,1..5] of Integer;
       vu                     : Array[1..4095] of Real;
       du                     : Array[1..12,1..3] of Real;
       xnodo, ynodo           : Array[1..30] of Real;
       t0,s0,s1,ru             : Array[1..3,1..3] of Real;

  colo,color,xmed,x1,y1,x2,y2,u,v,xx,yy   : Integer;
  astr,nome,nome1,nome2,ss,blk : String;
  scelta  : String[1];
  nul     : String[15];
  strprg  : String;
  col,coli : Array[1..4] of Integer;
  labe : Array[1..4] of String[23];
  lab  : Array[1..4] of Char;

 Function Sw(x:Real): Integer;
  Begin
    if x < 0 then Sw := -1 else Sw := 1;
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

 Procedure Retina1;
  Begin
    TextBackGround(5);
    ss:='';
    for u:=41 to 80 do
     begin
      ss:=ss+chr(176);
     end;
    for u:=1 to 23 do
     begin
     gotoxy(41,u); writeln(ss);
     end;
  End;

  Procedure MainFrame;
  begin
repeat
    clrscr;
    writeln;
     textcolor(10);
    writeln('          Calcolo automatico di STRUTTURE monodimensionali PIANE agli   ');
    writeln('             ELEMENTI FINITI - Aste rettilinee a SEZIONE COSTANTE.         ');
     textcolor(14);
    writeln;
    writeln('         1) MOMENTI POSITIVI Antiorari - SFORZI Normali POSIT.=Tiranti');
    writeln('         2) FORZE Nodali POSITIVE verso il BASSO e verso DESTRA       ');
    writeln('            Ordine Numerazione NODI:si comincia dagli N1 nodi LIBERI e poi');
    writeln('            si passa agli N2 del SUOLO continuando nella progressione.    ');
    writeln('            Si effettua quindi la Numerazione delle T TRAVI.              ');
    writeln('                                 I N P U T S                             ');
    writeln('         a) N1,N2,T ed E=Modulo Elastico (E=1 aste alla Gehler)       ');
    writeln('         b) Coordin.NODI:come quelle CARTESIANE (Positive:Destra,Alto)');
    writeln('         c) Numerazione NODI alle Estremità delle TRAVI               ');
    writeln('         d) Rigidit… VINCOLI: 1e6=Rig.Infinita: 0=Rig.Nulla=Cedimento ');
    writeln('         e) AREE e MOM_INERZ: Porre per aste alla Gehler: A=1000 J=1-10');
    writeln('         f) CARICHI ESTERNI: Secondo la casistica riportata; nel caso ');
    writeln('            GENERALE introd.FORZE E COPPIE ATTIVE [Opposte alle Reaz. ');
    writeln('            e ai Mom.Inc.Perfetto con le convenzioni di segno 1) e 2)]');
    writeln('            Dimensionamento VETTORI per valori MAX: N1=30; N2=12; T=60');
    writeln;
     textcolor(12);
      write('                                 Sigla TELAIO: '); read(wo);
     textcolor(14);
      readln;
  until wo<>'';
  end;

     Procedure InserimentoDATItelaio;
  Label Sal1,Sal2,Sal3,Sal4,SAl5,Sal6,Sal7,Sal8,Sub1145;
  Var   Opz1,yes    :Word;
        qu,maa,mbb,fbb,faa,ffb,ffa,hhb,hha,fm,aj,bj,eta,tau    :Real;
   Begin
   textbackground(16);
   textcolor(10);
   gotoxy(48,4); writeln('INSERIMENTO DATI'); textcolor(14);
   gotoxy(43,6); write('N1 = Numero nodi liberi  => '); read(n1);
   gotoxy(43,7); write('N2 = Numero nodi al suolo=> '); read(n2);
   gotoxy(43,8); write(' T = Numero travi        => '); read(t);
   gotoxy(43,9); write(' E = Modulo elastico     => '); read(E);
   n:= n1+n2;   qq:= 3*n*(3*n+1) div 2;
   Retina1;
   textcolor(10);
   gotoxy (52,1); writeln('(999 per ripetere)');
   gotoxy (48,2); writeln('Coord_NODI-Vettore N');
   textcolor(14);
       FOR k:= 1 TO n do begin
       gotoxy(48,3); write('X_nodo(',k:2,')= '); read(xnodo[k]);
       gotoxy(48,4); write('Y_nodo(',k:2,')= '); read(ynodo[k]);
       gotoxy(59,3); write('        ');
       gotoxy(59,4); write('        ');
       IF ynodo[k] = 999 THEN k:= k - 2;
       end;
       FOR k:= 1 TO n do begin
       nu[k, 1]:= -ynodo[k]; nu[k, 2]:= xnodo[k];
       end;
       FOR k:= 1 TO n do begin
         FOR z:= 1 TO 2 do begin
          n1u[k, z]:= nu[k, z];
       end; end;

   textcolor(10);
   gotoxy (48,5); writeln('Nodi estremi travi-Vettore A');
   textcolor(14);
       FOR k:= 1 TO t do begin
       gotoxy(48,6); write('Nodo basso trave(',k:2,')= '); read(au[k, 1]);
       gotoxy(48,7); write('Nodo alto  trave(',k:2,')= '); read(au[k, 2]);
       gotoxy(69,6); write('    ');
       gotoxy(69,7); write('    ');
       IF au[k, 2] = 999 THEN k:= k - 2;
       end;

   textcolor(10);
   gotoxy (48,8); writeln('Rigid. Vincoli-(1e6=inf.)');
   textcolor(14);
       FOR k:= n1 + 1 TO n do begin
       gotoxy(48,9); write('Rigidit… vert. (',k:2,')= '); read(du[k-n1, 1]);
      gotoxy(48,10); write('Rigidit… oriz. (',k:2,')= '); read(du[k-n1, 2]);
      gotoxy(48,11); write('Rigidit… angol.(',k:2,')= '); read(du[k-n1, 3]);
        gotoxy(69,9); write('          ');
       gotoxy(69,10); write('          ');
       gotoxy(69,11); write('          ');
       IF du[k - n1, 3] = 999 THEN k:= k - 2;
       end;

   textcolor(10);
   gotoxy (48,12); writeln('Caratter. Elastiche - Vett. T');
   textcolor(14);
       FOR k:= 1 TO t do begin
       tu[k, 1]:= 1;
       gotoxy(48,13); write('AREA      (',k:2,')= '); read(tu[k, 2]);
       gotoxy(48,14); write('Mom-INERZ.(',k:2,')= '); read(tu[k, 3]);
       gotoxy(63,13); write('          ');
       gotoxy(63,14); write('          ');
       IF tu[k, 3] = 999 THEN k:= k - 2;
       end;

     ASSIGN(ft,wo+'.TEF');
          Rewrite(ft);
     writeln(ft,n1:4,n2:4,t:4,E:18:3);
     n:= n1 + n2; qq:= 3 * n * (3 * n + 1) div 2;

   FOR k:= 1 TO n do begin writeln(ft, nu[k, 1]:7:3, nu[k, 2]:7:3); end;
   FOR k:= 1 TO t do begin writeln(ft, au[k, 1]:3, au[k, 2]:3); end;
   FOR k:= n1 + 1 TO n do begin
   writeln(ft, du[k - n1, 1]:9:0, du[k - n1, 2]:9:0, du[k - n1, 3]:9:0); end;
   FOR k:= 1 TO t do begin tu[k, 1]:= 1; writeln(ft, tu[k, 2]:12:6, tu[k, 3]:12:8); end;

   gotoxy (50, 16); textcolor (12);
   writeln('FINE Inserimento Dati ');
   gotoxy (50, 17);
   writeln('Geometrici ed Elastici'); textcolor(14);
       ch:=readkey;
   { Inserimento Carichi Esterni }
  clrscr;
  FOR i:= 1 TO t do begin
  tu[i, 4]:= SQRT(SQR(nu[au[i, 2], 1] - nu[au[i, 1], 1]) + SQR(nu[au[i, 2], 2] - nu[au[i, 1], 2]));
  tu[i, 5]:= (nu[au[i, 1], 1] - nu[au[i, 2], 1]) / tu[i, 4];
  tu[i, 6]:= (nu[au[i, 2], 2] - nu[au[i, 1], 2]) / tu[i, 4];
   end;
          k:= 1;
 FOR i:= 1 TO t do begin
 fa[i, 1]:= 0; fa[i, 2]:= 0; ha[i, 1]:= 0; ha[i, 2]:= 0; ma[i, 1]:= 0; ma[i, 2]:= 0;
 end;

Sal1: textcolor(10); gotoxy(18,3);writeln('ASTA ',k:2); textcolor(14);
   gotoxy(10,5);writeln('1) Asta Scarica               ');
   gotoxy(10,6);writeln('2) Carico VERTICALE uniforme  ');
   gotoxy(10,7);writeln('3) Carichi VERTIC. conc. nodi ');
   gotoxy(10,8);writeln('4) Carichi ORIZZ.  conc. nodi ');
   gotoxy(10,9);writeln('5) Carico VERTIC. conc.campata');
  gotoxy(10,10);writeln('6) F+C attive (incastri perf.)');
  gotoxy(15,12);  write('Scelta opzione: '); read(Opz1);
        if Opz1=1 then goto Sal2;
        if Opz1=2 then goto Sal4;
        if Opz1=3 then goto Sal5;
        if Opz1=4 then goto Sal6;
        if Opz1=5 then goto Sal7;
        if Opz1=6 then goto Sal8;

Sal2: k:= k + 1;
IF k > t THEN goto Sal3 ELSE goto Sal1; { Asta scarica }

Sal4:
      gotoxy (53,3); writeln('ASTA ',k:2); {Carico verticale uniforme}
      gotoxy (48,5); writeln('Lunghezza Asta = ', tu[k, 4]:8:3);
      gotoxy (48,6); write('   Carico q (ton/ml) = '); read(qu);
   maa:= qu*SQR(tu[k, 4])/12 * SQR(tu[k, 6]);
   faa:= qu*tu[k, 4]/2 * ABS(tu[k, 6]);
fa[k, 1]:= fa[k, 1] + faa; fa[k, 2]:= fa[k, 2] + faa;
ma[k, 1]:= ma[k, 1] - maa * Sw(tu[k, 6]); ma[k, 2]:= ma[k, 2] + maa * Sw(tu[k, 6]);
  goto Sub1145; k:= k + 1; IF k > t THEN goto Sal3 ELSE goto Sal1;

Sal5:
      gotoxy (53,3); writeln('ASTA ',k:2); {Carichi verticali concentr. ai nodi}
      gotoxy (48,5); writeln('Lunghezza Asta = ', tu[k, 4]:8:3);
      gotoxy (48,6); write('Car. Vert.(ton) n. basso= '); read(ffb);
      gotoxy (48,7); write('Car. Vert.(ton) n. alto = '); read(ffa);
fa[k, 1]:= fa[k, 1] + ffb; fa[k, 2]:= fa[k, 2] + ffa;
goto Sub1145; k:= k + 1; IF k > t THEN goto Sal3 ELSE goto Sal1;

Sal6:
      gotoxy (53,3); writeln('ASTA ',k:2); {Carichi orizzont. concentr.ai nodi}
      gotoxy (48,5); writeln('Lunghezza Asta = ', tu[k, 4]:8:3);
      gotoxy (48,6); write('Car. Oriz.(ton) n. basso= '); read(hhb);
      gotoxy (48,7); write('Car. Oriz.(ton) n. alto = '); read(hha);
ha[k, 1]:= ha[k, 1] + hhb; ha[k, 2]:= ha[k, 2] + hha;
goto Sub1145; k:= k + 1; IF k > t THEN goto Sal3 ELSE goto Sal1;

Sal7:
      gotoxy (53,3); writeln('ASTA ',k:2); {Carichi vertic. conc. in campata}
      gotoxy (48,5); writeln('Lunghezza Asta = ', tu[k, 4]:8:3);
      gotoxy (48,6); write('Carico Concentr.(ton) = '); read(fm);
      gotoxy (48,7); write('Dist.(mt) nodo sin. = '); read(aj);

bj:= tu[k, 4] - aj; tau:= aj / tu[k, 4]; eta:= bj / tu[k, 4];
maa:= -fm * aj * SQR(eta) * tu[k, 6]; mbb:= fm * bj * SQR(tau) * tu[k, 6];
faa:= fm * eta * (1 + tau * eta - SQR(tau)); fbb:= fm * tau * (1 + tau * eta - SQR(eta));
fa[k, 1]:= fa[k, 1] + faa; fa[k, 2]:= fa[k, 2] + fbb;
ma[k, 1]:= ma[k, 1] + maa; ma[k, 2]:= ma[k, 2] + mbb;
goto Sub1145; k:= k + 1; IF k > t THEN goto Sal3 ELSE goto Sal1;

Sal8:
      gotoxy (53,3); writeln('ASTA ',k:2); {F+C attive (incastri perf.)}

     gotoxy (48,5); writeln('    Nodo basso');
     gotoxy (48,6);write('F-att. = '); read(fa[k, 1]);
     gotoxy (48,7);write('H-att. = '); read(ha[k, 1]);
     gotoxy (48,8);write('M-att. = '); read(ma[k, 1]);

     gotoxy (48,9); writeln('    Nodo alto ');
     gotoxy (48,10);write('F-att. = '); read(fa[k, 2]);
     gotoxy (48,11);write('H-att. = '); read(ha[k, 2]);
     gotoxy (48,12);write('M-att. = '); read(ma[k, 2]);

     IF ma[k, 2] = 999 THEN k:= k - 2;
goto Sub1145; k:= k + 1; IF k > t THEN goto Sal3 ELSE goto Sal1;
Sal3:

FOR k:= 1 TO t do begin writeln(ft,fa[k,1]:12:6,ha[k,1]:12:6,ma[k,1]:12:6,fa[k,2]:12:6,ha[k,2]:12:6,ma[k,2]:12:6); end;
   Close (ft);
exit;

 Sub1145:
     gotoxy (48,14); write('Altro tipo di carico? (s/ ) ');
       ch:=readkey;
       clrscr;
         if (ch='s') or (ch='S') then goto Sal1
         ELSE begin
         k:=k+1;
         goto Sal1;
         end;

    End;

     Procedure CaricaDatidaFile;
   Begin
     ASSIGN(fl,wo+'.TEF');
          Reset(fl);
          readln(fl,n1,n2,t,E);
          n:= n1 + n2; qq:= 3 * n * (3 * n + 1) div 2;
   FOR k:= 1 TO n do begin readln(fl, nu[k, 1], nu[k, 2]); end;
   FOR k:= 1 TO t do begin readln(fl, au[k, 1], au[k, 2]); end;
   FOR k:= n1 + 1 TO n do begin
   readln(fl, du[k - n1, 1], du[k - n1, 2], du[k - n1, 3]); end;
   FOR k:= 1 TO t do begin tu[k, 1]:= 1; readln(fl, tu[k, 2], tu[k, 3]); end;
FOR k:= 1 TO t do begin readln(fl,fa[k,1],ha[k,1],ma[k,1],fa[k,2],ha[k,2],ma[k,2]); end;


   Close (fl);
    End;

     Procedure Sub3330;
   Var flag :Word;
   Begin
 t0[1, 1]:= tu[cu[h, 2], 6];
 t0[1, 2]:= tu[cu[h, 2], 5];
 t0[2, 1]:= -t0[1, 2];
 t0[2, 2]:= t0[1, 1];
 t0[3, 3]:= 1;
   flag:= flag + 1; gotoxy(30,21); writeln('FLAG = ',flag:3);
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin ru[l, f]:= 0; end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin FOR k:= 1 TO 3 do begin
            ru[f, l]:= ru[f, l] + s0[f, k] * t0[k, l]; end;end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin s0[l, f]:= ru[l, f]; end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin ru[l, f]:= t0[f, l]; end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin t0[l, f]:= ru[l, f]; end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin ru[l, f]:= 0; end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin FOR k:= 1 TO 3 do begin
            ru[f, l]:= ru[f, l] + t0[f, k] * s0[k, l]; end;end;end;
 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin s0[l, f]:= ru[l, f]; end;end;
 End;

     Procedure Sub3480;
 Label Salta1,Salta2,Salta3;
 Var
      w: Real;

   Begin
      m:= 3 * n;
 gotoxy(50,21); writeln('3*N = ', m:3);
  FOR i:= 1 TO m do begin
         gotoxy(50,20); writeln('  I = ',1:3);
    FOR j:= i TO m do begin w:= 0;
      FOR q:= 1 TO i - 1 do begin
      IF i = 1 THEN goto Salta1;
        k1:= m * (q - 1) - q * (q - 1) div 2 + i;
        k2:= m * (q - 1) - q * (q - 1) div 2 + j;
             w:= w + vu[k1] * vu[k2];
      end;
Salta1: k:= m * (i - 1) - i * (i - 1) div 2 + j;
       k0:= m * (i - 1) - i * (i - 1) div 2 + i;
        IF k = k0 THEN goto Salta2;
     vu[k]:= (vu[k] - w) / vu[k0];
     goto Salta3;
Salta2:   vu[k]:= SQRT(vu[k] - w);
Salta3: end; end;
   FOR i:= 1 TO m do begin yu[i]:= yu[i]/vu[(i - 1) * m - i * (i - 1) div 2 + i];
   FOR j:= 1 TO i - 1 do begin
   yu[i]:= yu[i]-vu[(j - 1)*m-j*(j - 1) div 2+i] * yu[j] / vu[(i - 1) * m - i * (i - 1) div 2 + i];
   end;end;

  FOR i:= m DOWNTO 1 do begin
     yu[i]:= yu[i] / vu[(i - 1) * m - i * (i - 1) div 2 + i];
    FOR j:= i + 1 TO m do begin
       yu[i]:= yu[i] - vu[(i - 1) * m - i * (i - 1) div 2 + j] * yu[j] / vu[(i - 1) * m - i * (i - 1) div 2 + i];
   end;end;
 End;

     Procedure EsecuzioneCalcoli;
  Label Skip1,Skip2,Skip3,Skip4;
  Var
       som1,som2,som3,cs,sn,vv,hv,mv,
       b1,b2,b3,a1,a2,a3,v1,w1,v2,w2,g1,g2,g3,g4,
       ns,nd,ms,md,ts,td             :Real;
     rr,r0,r1,r2,r3,u1,u2,u3      :Word;
     Begin
{
 #########################################################################
     ......  ESECUZIONE Programma  (Calcoli e Risultati finali)
 #########################################################################
}
    TextBackGround(5);
FOR k:= 1 TO t do begin FOR z:= 1 TO 2 do begin
 fa1[k, z]:= fa[k, z]; ha1[k, z]:= ha[k, z]; ma1[k, z]:= ma[k, z]; end;end;

a:= 0;
FOR i:= 1 TO n do begin
  FOR j:= 1 TO 5 do begin bu[i, j]:= 0; end;end;

FOR i:= 1 TO n do begin
 s:= 1;
  FOR j:= 1 TO t do begin
IF (au[j, 1] <> i) AND (au[j, 2] <> i) THEN goto Skip1;
IF s > a THEN a:= s;
 bu[i, s]:= j; s:= s + 1;
Skip1: end;end;

  FOR i:= 1 TO t do begin
  tu[i, 4]:= SQRT(SQR(nu[au[i, 2], 1] - nu[au[i, 1], 1]) + SQR(nu[au[i, 2], 2] - nu[au[i, 1], 2]));
  tu[i, 5]:= (nu[au[i, 1], 1] - nu[au[i, 2], 1]) / tu[i, 4];
  tu[i, 6]:= (nu[au[i, 2], 2] - nu[au[i, 1], 2]) / tu[i, 4];
  end;
          FOR i:= 1 TO qq do begin  vu[i]:= 0; end;
  FOR i:= 1 TO 3 do begin FOR j:= 1 TO 3 do begin t0[i, j]:= 0; end;end;
   gotoxy(30,20); writeln('N*A  = ',n * a:2);

          FOR i:= 1 TO n do begin
  som1:= 0; som2:= 0; som3:= 0;
 r0:= 9 * (i - 1) * n - (9 * i*i - 21 * i + 12) div 2;
 r1:= r0 + 1; r2:= r0 + 3 * (n - i + 1) + 1; r3:= r0 + 6 * (n - i + 1);
    FOR l:= 1 TO a do begin FOR f:= 1 TO 2 do begin cu[l, f]:= 0; end;end;
           FOR h:= 1 TO a do begin
        IF bu[i, h] = 0 THEN goto Skip2;
      cu[h, 2]:= bu[i, h];
IF au[bu[i, h], 1] = i THEN cu[h, 1]:= au[bu[i, h], 2] ELSE cu[h, 1]:= au[bu[i, h], 1];
  IF cu[h, 1] = i THEN goto Skip2;
  IF cu[h, 1] < i THEN fa[bu[i, h], 1]:=fa[bu[i, h], 2];
  IF cu[h, 1] < i THEN ha[bu[i, h], 1]:=ha[bu[i, h], 2];
  IF cu[h, 1] < i THEN ma[bu[i, h], 1]:=ma[bu[i, h], 2];
      som1:= som1 + fa[bu[i, h], 1];
      som2:= som2 + ha[bu[i, h], 1];
      som3:= som3 + ma[bu[i, h], 1];
  yu[3 * i - 2]:= som1;
  yu[3 * i - 1]:= som2;
  yu[3 * i]:= som3;
Skip2: end;
   FOR j:= i + 1 TO n do begin
        FOR k:= 1 TO 3 do begin FOR r:= 1 TO 3 do begin s0[k, r]:= 0; end;end;
   FOR h:= 1 TO a do begin
        IF cu[h, 1] <> j THEN goto Skip3;
     s0[1, 1]:= -12 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / (tu[cu[h, 2], 4]*SQR(tu[cu[h, 2], 4]));
     s0[2, 2]:= -tu[cu[h, 2], 1] * tu[cu[h, 2], 2] / tu[cu[h, 2], 4];
     s0[3, 3]:= 2 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / tu[cu[h, 2], 4];
     s0[1, 3]:= 6 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / SQR(tu[cu[h, 2], 4]);
     s0[3, 1]:= -s0[1, 3];
 Sub3330;
       rr:= 3 * (j - i);
   u1:= r1 + rr; u2:= r2 + rr; u3:= r3 + rr;
     vu[u1]:= s0[1, 1]; vu[u1 + 1]:= s0[2, 1]; vu[u1 + 2]:= s0[3, 1];
     vu[u2 - 1]:= s0[1, 2]; vu[u2]:=s0[2, 2]; vu[u2 + 1]:= s0[3, 2];
     vu[u3 - 2]:= s0[1, 3]; vu[u3 - 1]:= s0[2, 3]; vu[u3]:= s0[3, 3];
Skip3:  end;end;

 FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin s1[l, f]:= 0; end;end;
   FOR h:= 1 TO a do begin
       FOR l:= 1 TO 3 do begin FOR f:= 1 TO 3 do begin s0[l, f]:= 0; end;end;
      IF cu[h, 1] = 0 THEN goto Skip4;
     s0[1, 1]:= 12 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / (tu[cu[h, 2], 4]*SQR(tu[cu[h, 2], 4]));
     s0[2, 2]:= tu[cu[h, 2], 1] * tu[cu[h, 2], 2] / tu[cu[h, 2], 4];
     s0[3, 3]:= 4 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / tu[cu[h, 2], 4];
     s0[1, 3]:= 6 * tu[cu[h, 2], 1] * tu[cu[h, 2], 3] / SQR(tu[cu[h, 2], 4]);
       IF cu[h, 1] > i THEN s0[1, 3]:= -s0[1, 3];
     s0[3, 1]:= s0[1, 3];
   Sub3330;
 FOR k:= 1 TO 3 do begin FOR z:= 1 TO 3 do begin s1[k, z]:= s1[k, z] + s0[k, z]; end;end;
Skip4:  end;
     vu[r1]:= s1[1, 1]; vu[r1 + 1]:= s1[2, 1]; vu[r1 + 2]:= s1[3, 1];
     vu[r2]:= s1[2, 2]; vu[r2 + 1]:= s1[3, 2]; vu[r3]:= s1[3, 3];
end;

         FOR i:= n1 + 1 TO n do begin
 r0:= 9 * (i - 1) * n - (9 * i*i - 21 * i + 12) div 2;
 r1:= r0 + 1; r2:= r0 + 3 * (n - i + 1) + 1; r3:= r0 + 6 * (n - i + 1);
   vu[r1]:= vu[r1] + du[i - n1, 1];
   vu[r2]:= vu[r2] + du[i - n1, 2];
   vu[r3]:= vu[r3] + du[i - n1, 3];
      end;
  Sub3480;
{ --------------- RISULTATI FINALI --------------- }
 clrscr;
     ASSIGN(ft,wo+'.TXT');
          Rewrite(ft);
  gotoxy(27,1); writeln('Telaio ' + wo + ' - ' + 'RISULTATI FINALI');
  writeln(ft,'                         Telaio ' + wo + ' - ' + 'RISULTATI FINALI');
  gotoxy(11,2); writeln('[MOMENTI Antiorari POSITIVI - TAGLI Sin-Alto e Des-Basso POSITIVI');
  writeln(ft,'         [MOMENTI Antiorari POSITIVI - TAGLI Sin-Alto e Des-Basso POSITIVI');
  gotoxy(16,3); writeln('SFORZI NORMALI Positivi=TIRANTI; Negativi=PUNTONI]');
  writeln(ft,'              SFORZI NORMALI Positivi=TIRANTI; Negativi=PUNTONI]');
  
   FOR i:= 1 TO t do begin
 b1:= yu[(au[i, 2] - 1) * 3 + 1]; b2:= yu[(au[i, 2] - 1) * 3 + 2]; b3:= yu[(au[i, 2] - 1) * 3 + 3];
 a1:= yu[(au[i, 1] - 1) * 3 + 1]; a2:= yu[(au[i, 1] - 1) * 3 + 2]; a3:= yu[(au[i, 1] - 1) * 3 + 3];
    v1:= a1 * tu[i, 6] + a2 * tu[i, 5];
    w1:= -a1 * tu[i, 5] + a2 * tu[i, 6];
       v2:= b1 * tu[i, 6] + b2 * tu[i, 5];
       w2:= -b1 * tu[i, 5] + b2 * tu[i, 6];

     writeln('                              TRAVE n.',i:2);
   {  writeln('                              =========='); }
     writeln('             Sinistra                        Destra');

  writeln(ft,'                              TRAVE n.',i:2);
{  writeln(ft,'                              ==========');  }
  writeln(ft,'             Sinistra                        Destra');

     g1:= (w2 - w1) * tu[i, 1] * tu[i, 2] / tu[i, 4];
     ns:= g1 - fa1[i, 1] * tu[i, 5] + ha1[i, 1] * tu[i, 6];
     nd:= g1 + fa1[i, 2] * tu[i, 5] - ha1[i, 2] * tu[i, 6];

     g2:= 6 * (v2 - v1) * tu[i, 1] * tu[i, 3] / SQR(tu[i, 4]) + (4 * a3 + 2 * b3) * tu[i, 1] * tu[i, 3] / tu[i, 4];
     ms:= g2 - ma1[i, 1];

     g3:= 6 * (v2 - v1) * tu[i, 1] * tu[i, 3] / SQR(tu[i, 4]) + (2 * a3 + 4 * b3) * tu[i, 1] * tu[i, 3] / tu[i, 4];
     md:= g3 - ma1[i, 2];

     g4:= 12*(v2-v1)*tu[i,1]*tu[i,3]/(tu[i,4]*SQR(tu[i,4]))+6*(a3+b3)*tu[i,1]*tu[i,3]/SQR(tu[i, 4]);
     ts:= g4 + fa1[i, 1] * tu[i, 6] + ha1[i, 1] * tu[i, 5];
     td:= g4 - fa1[i, 2] * tu[i, 6] - ha1[i, 2] * tu[i, 5];

writeln('               M = ',ms:8:3,'                    M = ',md:8:3);
writeln('               T = ',ts:8:3,'                    T = ',td:8:3);
writeln('               N = ',ns:8:3,'                    N = ',nd:8:3);
 
writeln(ft,'               M = ',ms:8:3,'                    M = ',md:8:3);
writeln(ft,'               T = ',ts:8:3,'                    T = ',td:8:3);
writeln(ft,'               N = ',ns:8:3,'                    N = ',nd:8:3);
 
 IF au[i, 2] > n1 THEN begin
 cs:= tu[i, 6]; sn:= tu[i, 5];
 vv:= -fa1[i, 2] * Sw(cs) + g1 * sn - g4 * cs;
 hv:= -ha1[i, 2] + g1 * cs + g4 * sn;
 mv:= g3 - ma1[i, 2];
writeln('        Reazioni vincolari:  V=',vv:8:3,'    H = ',hv:8:3,'   M = ',mv:8:3);
writeln(ft,'        Reazioni vincolari:  V=',vv:8:3,'    H = ',hv:8:3,'   M = ',mv:8:3);
end;
    end;
  Close(ft);
  repeat until keypressed;
 End;

 Procedure MessLine;
  Begin
   gotoxy(1,24);
   TextBackground(7);
   TextColor(4);
   writeln(blk);
  End;

 Procedure Retina;
  Begin
    TextBackGround(6);
    ClrScr;
    ss:='';
    for u:=1 to 79 do
     begin
      ss:=ss+chr(176);
     end;
    for u:=1 to 25 do
     begin
      writeln(ss);
     end;
  End;

 Procedure BoxB;
  Begin
   ss:='';
   xmed:=(x1+x2) div 2;
   xx:=x2-x1;
   yy:=y2-y1;

   for u:=1 to xx-2 do
    begin
     ss:=ss+chr(196);
    end;
    TextColor(4);
    TextBackGround(7);

    gotoxy(x1,y1);
    writeln(chr(218),ss,chr(191));
    for u:=1 to yy-1 do
     begin
      gotoxy(x1,wherey);
      writeln(chr(179),chr(179):(xx-1));
     end;
     gotoxy(x1,wherey);
     writeln(chr(192),ss,chr(217));
 End;

 Procedure BordEx;
  begin
   ss:='';
   xx:=x2-x1;
   yy:=y2-y1;

   for u:=1 to xx do
    begin
    ss:=ss+chr(176);
    end;
    gotoxy(x1+1,y2+1);
    TextBackGround(16);
    Textcolor(8);
    writeln(ss);
    gotoxy(x2,y1+1);
    for u:=1 to yy do
     begin
      writeln(chr(176));
      gotoxy(x2,wherey);
     end;
  End;

 Procedure Piccola;
  begin
   TextBackground(colo);
   TextColor(color);
   gotoxy(9,6+2*(i-1));
   writeln(labe[i]);
  end;

 Procedure MenuIn;
  Begin
   lab[1]:='1';
   lab[2]:='2';
   lab[3]:='3';
   lab[4]:='4';
   labe[1]:=' Inser. dati (iniziale)';
   labe[2]:=' Carica dati (da File )';
   labe[3]:=' DISEGNO del TELAIO   ';
   labe[4]:=' Uscita               ';
   for i:=1 to 4 do
    begin
     col[i]:=7;
     coli[i]:=4;
    end;

   col[1]:=9;
   coli[1]:=14;

   TextColor(9);
   Retina;
   MessLine;
   gotoxy(20,24);
   TextColor(1);
   writeln('Freccia in alto o in basso per la scelta');
   ss:='';
   x1:=3;
   y1:=4;
   x2:=37;
   y2:=14;
   BoxB;
   BordEx;
   textcolor(4);
   textbackground(7);
   gotoxy(xmed-7,y1);
   writeln('M E N U');

    for i:=1 to 4 do
     begin
      colo:=col[i];
      color:=coli[i];
      textbackground(7);
      textcolor(0);
      gotoxy(7,6+2*(i-1));
      writeln(lab[i]+'.');
      Piccola;
     end;

     gotoxy(75,23);
     i:=1;

     repeat
     ch:=readkey;
     if ch=#72 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       i:=i-1;
       if i<1 then i:=4;
       colo:=9;
       color:=14;
       Piccola;
      end;
     if ch=#80 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       i:=i+1;
       if i>4 then i:=1;
       colo:=9;
       color:=14;
       Piccola;
      end;

      gotoxy(75,23);

     until (ch=#13);

     if ch=#13 then str(i:1,scelta)

  End;


   Procedure subr;
 Var x9,y9   :Word;
   Begin
     x9:= x; y9:= y;
     x:= x9 + xc;  y:= y9 + yc;
   End;

   Procedure DisegnoTelaio;
Label vai,vai2;
   Var
        x1,y1,x2,y2,x3,y3,xw,yw,colore,ka,pas,ziz:   Word;
        dey,dez                 : Real;
        cho                     : String;
        sce                     : String[2];
        gd,gm,k                 : Integer;

   Begin
gd:=detect;initgraph(gd,gm,'');
if graphresult <> grOk then halt(1);
ka:=10; pas:=2; ziz:=ka;

vai:
  SetFillStyle(1,6);
  FloodFill(0,0,12);
  SetColor(14);

SetTextStyle(1,0,8);
SetLineStyle(0,0,1);
xw:=GetmaxX div 2; yw:=GetmaxY div 2;
Rectangle(xw-300,yw-220,xw+300,yw-155);
OutText('  T E L A I O');

 SetFillStyle(1,12);

SetTextStyle(0,0,1);
SetLineStyle(0,0,1);
Rectangle(xw-315,yw+160,xw-75,yw+208);
SetColor(10);
OutTextXY(xw-310,yw+165,' + => Ingrandisce     Grafico');
OutTextXY(xw-310,yw+180,' - => Rimpicciolisce  Grafico');
OutTextXY(xw-310,yw+195,'      Doppio INVIO per USCIRE');
SetColor(14);

 { ------------------- DISEGNO TELAIO ---------------- }
    { Centratura del telaio }
       FOR k:= 1 TO n do begin
         FOR z:= 1 TO 2 do begin
          n1u[k, z]:= nu[k, z];
          if abs(n1u[k,1]) > abs(maxY) then maxY:=abs(n1u[k,1]);
          if abs(n1u[k,2]) > abs(maxZ) then maxZ:=abs(n1u[k,2]);
       end; end;
  dey:=maxY; dez:=maxZ;

      { Inizio DISEGNO }
xc:=GetmaxX div 2-ka*round(dez/2); yc:=GetmaxY div 2+ka*round(dey/2);
      FOR i:= 1 TO t do begin
      x:= ka*round(nu[au[i, 1], 2]); y:= ka*round(nu[au[i, 1], 1]); subr; x1:= x; y1:= y;
      x:= ka*round(nu[au[i, 2], 2]); y:= ka*round(nu[au[i, 2], 1]); subr; x2:= x; y2:= y;

      x:= ka*(round(nu[au[i, 1], 2]) + round(nu[au[i, 2], 2])) div 2;
      y:= ka*(round(nu[au[i, 1], 1]) + round(nu[au[i, 2], 1])) div 2;
      subr; x3:= x; y3:= y;
      str(i:2,sce);
      SetColor(10); OutTextXY(x3-6,y3+2,sce); SetColor(14);
      LINE (x1, y1, x2, y2);
      end;

 { Disegno dei Vincoli alla base }

           FOR i:= n1 + 1 TO n do begin
      x:= ka*round(nu[i, 2]); y:= ka*round(nu[i, 1]); subr; x1:= x; y1:= y;
      IF du[i - n1, 3]=0 THEN begin
      CIRCLE (x1, y1, ziz div 4); GOTO vai2; end;
      x:= ka*round(nu[i, 2]) - ziz; y:=ka*round(nu[i, 1]); subr; x1:= x; y1:= y;
      x:= ka*round(nu[i, 2]) + ziz; y:=ka*round(nu[i, 1]); subr; x2:= x; y2:= y;
      LINE (x1, y1, x2, y2);
vai2: end;

cho:=readkey;
if cho='+' then
begin
ClearDevice;
ka:=ka+pas;
goto vai;
end;

if cho='-' then
begin
ClearDevice;
ka:=ka-pas;
goto vai;
end;

closegraph;
restoreCrtMode;
 End;
{------------------------------------------------------------------------------
  Inizio Programma }

  Begin
     ColoriSchermo;
     MainFrame;

   blk:='';
   for u:=1 to 79 do
    begin
     blk:=blk+chr(32);
    end;
   ClrScr;

 repeat
   MenuIn;
   ch:=' ';
          if scelta='1' then begin
          InserimentoDATItelaio; EsecuzioneCalcoli; end
     ELSE if scelta='2' then begin
          CaricaDatidaFile; EsecuzioneCalcoli; end
     ELSE if scelta='3' then begin
          DisegnoTelaio; end
     ELSE exit;

 ch:=readkey;
until (ch<>#13);
 end.

  { FINE Programma }
