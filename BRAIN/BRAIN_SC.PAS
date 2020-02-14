 Program BRAIN; { EFFEMERIDI Planetarie di GRANDISSIMA Precisione
                             da MERCURIO a NETTUNO:
                Teoria di Bretagnon:VSOP87 - Astronomical Algoritms di J.Meeus
                                    ----------------
                Accuratezza dei RISULTATI: Esempio di p.211-213
                  Venere al 20/12/1992 ore 0:0 TD
                    (libro)        (questo progr.)     (astron.profess.)
        à_app=   21 04 41.50        21 04 41.49           21 04 41.454
        ë_app=  -18 53 16.9        -18 53 17.0           -18 53 16.84

 Data in Ore,Min e SECONDI <- aggiunta l'Opzione dei secondi in BRAIN_SC.pas
 }
 Uses Crt,SOLEXYZ,MERC_TPU,VENE_TPU,MART_TPU,TERR_TPU,
          GIOV_TPU,SATU_TPU,URAN_TPU,NETT_TPU;

 Label  Ripeti,Ripeti2;

 Const  rad=pi/180; vluce=0.0057755183;

 Var
  iqw,colo,color,xmed,xqw1,yqw1,xqw2,yqw2,uqw,xxqw,yyqw,band,flag : Integer;
      ch,AB: Char;
  ss,blk,as,ddata  : String;
    scelta : String[1];
  col,coli : Array[1..9] of Integer;
      labe : Array[1..9] of String[27];
      lab  : Array[1..9] of Char;
           aw, L0,B0,R0,L,B,R,Xso,Yso,Zso,Xeq,Yeq,Zeq,
           Xpri,Ypri,Zpri,alfa,delt,dgeo,lam,bet,tluce,tcorr,
           ddpsi,ddeps,eps0,eps,Tg :Real;

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
    for uqw:=1 to 79 do
     begin
      ss:=ss+chr(176);
     end;
    for uqw:=1 to 25 do
     begin
      writeln(ss);
     end;
  End;

 Procedure BoxB;
  Begin
   ss:='';
   xmed:=(xqw1+xqw2) div 2;
   xxqw:=xqw2-xqw1;
   yyqw:=yqw2-yqw1;

   for uqw:=1 to xxqw-2 do
    begin
     ss:=ss+chr(196);
    end;
    TextColor(4);
    TextBackGround(7);

    gotoxy(xqw1,yqw1);
    writeln(chr(218),ss,chr(191));
    for uqw:=1 to yyqw-1 do
     begin
      gotoxy(xqw1,wherey);
      writeln(chr(179),chr(179):(xxqw-1));
     end;
     gotoxy(xqw1,wherey);
     writeln(chr(192),ss,chr(217));
 End;

 Procedure BordEx;
  begin
   ss:='';
   xxqw:=xqw2-xqw1;
   yyqw:=yqw2-yqw1;

   for uqw:=1 to xxqw do
    begin
    ss:=ss+chr(176);
    end;
    gotoxy(xqw1+1,yqw2+1);
    TextBackGround(16);
    Textcolor(8);
    writeln(ss);
    gotoxy(xqw2,yqw1+1);
    for uqw:=1 to yyqw do
     begin
      writeln(chr(176));
      gotoxy(xqw2,wherey);
     end;
  End;

 Procedure Piccola;
  begin
   TextBackground(colo);
   TextColor(color);
   gotoxy(28,7+1*(iqw-1));
   writeln(labe[iqw]);
  end;

 Procedure MenuIn;
  Begin
  lab[1]:='1';
  lab[2]:='2';
  lab[3]:='3';
  lab[4]:='4';
  lab[5]:='5';
  lab[6]:='6';
  lab[7]:='7';
  lab[8]:='8';
  lab[9]:='9';
   labe[1]:='     M E R C U R I O       ';
   labe[2]:='       V E N E R E         ';
   labe[3]:=' <<< S O L E Apparente >>> ';
   labe[4]:='         M A R T E         ';
   labe[5]:='         G I O V E         ';
   labe[6]:='       S A T U R N O       ';
   labe[7]:='         U R A N O         ';
   labe[8]:='       N E T T U N O       ';
   labe[9]:=' << Altra Data / USCITA >> ';
   for iqw:=1 to 9 do
    begin
     col[iqw]:=7;
     coli[iqw]:=4;
    end;

   col[1]:=9;
   coli[1]:=14;

   TextColor(9);
   Retina;
   MessLine;
   gotoxy(20,24);
   TextColor(1);
   writeln('Muovere le frecce o selezionare il numero');
   ss:='';
   xqw1:=23;
   yqw1:=6;
   xqw2:=57;
   yqw2:=16;
   BoxB;
   BordEx;
   textcolor(4);
   textbackground(7);
   gotoxy(xmed-15,yqw1);
   writeln(' EFFEMERIDI di Alta PRECISIONE ');
   gotoxy(xmed-14,yqw2+3);
   {writeln(' Ing.G. Matarazzo-0931/945259 ');}

    for iqw:=1 to 9 do
     begin
      colo:=col[iqw];
      color:=coli[iqw];
      textbackground(7);
      textcolor(0);
      gotoxy(26,7+1*(iqw-1));
      writeln(lab[iqw]+'.');
      Piccola;
     end;

     gotoxy(75,23);
     iqw:=1;

     repeat
     ch:=readkey;
     if ch=#72 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       iqw:=iqw-1;
       if iqw<1 then iqw:=9;
       colo:=9;
       color:=14;
       Piccola;
      end;
     if ch=#80 then
      begin
       colo:=7;
       color:=4;
       Piccola;
       iqw:=iqw+1;
       if iqw>9 then iqw:=1;
       colo:=9;
       color:=14;
       Piccola;
      end;

      gotoxy(75,23);

     until (ch=#13) or (ch=lab[1]) or (ch=lab[2]) or (ch=lab[3])
            or (ch=lab[4]) or (ch=lab[5]) or (ch=lab[6]) or (ch=lab[7])
            or (ch=lab[8]) or (ch=lab[9]);

     if ch=#13 then str(iqw:2,scelta)
     else
     begin
      for iqw:=1 to 9 do
       begin
        if ch=lab[iqw] then str(iqw:2,scelta);
       end;
     end;

  End;

  Procedure GiornoGiuliano;
 Var m,ha,mn,sc: WORD;
        y: INTEGER;
  g,mm,yy,a,bw,kk,c,d: REAL;
  dd1,dd2,dd3,dd4,dd5,tm1,tm2: STRING;

    BEGIN
    textcolor(14);textbackground(1);
     clrscr;
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄ EFFEMERIDI Planetarie ÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('³                                    ³');
gotoxy(20,7); write('³                                    ³');
gotoxy(20,8); write('³                                    ³');
gotoxy(20,9); write('³                                    ³');
gotoxy(20,10);write('³                                    ³');
gotoxy(20,11);write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(25,5);write('  Giorno = ');readln(g);
      until (g<32) and (g>=0);
repeat gotoxy(25,6);write('    Mese = ');readln(m);
      until m in [1..12];
       gotoxy(25,7);write('    Anno = ');readln(y);
repeat gotoxy(25,8);write(' Ore TDT = ');readln(ha);
      until ha in [0..23];
repeat gotoxy(25,9);write(' Min TDT = ');readln(mn);
      until mn in [0..59];
repeat gotoxy(25,10);write(' Sec TDT = ');readln(sc);
      until sc in [0..59];

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
      jd:=bw+c+d+g+1720994.5+ha/24+mn/1440+sc/86400;
  str(int(jd):7:0,tm1);
  str(int(frac(jd)*1e5+0.5)/1e5:7:5,tm2);
  delete(tm2,1,1);
gotoxy(23,12);writeln('Il Giorno Giuliano della DATA Š : ');
gotoxy(30,13);writeln(tm1+tm2);
       { Millenni Giuliani }
  tau:=(jd-2451545)/365250;
     str( g:2:0,dd1);
     str(mm:2:0,dd2);
     str(yy:4:0,dd3);
     str(ha:0,dd4);
     str(mn:0,dd5);
     ddata:=dd1+'/'+dd2+'/'+dd3+' - '+dd4+':'+dd5+' TD';
       readln;
     end;

  Procedure CancellaSchermo;
    Begin
    textcolor(14);textbackground(1); clrscr;
    End;

  Procedure ControlloFlag;
 Var
      fr                     : file of BYTE;
      byt     : Byte;
      nome     : String[12];
  Begin
     nome:='BRAIN.EXE';
     ASSIGN(fr,nome);      {legge s‚ stesso!!!}
          Reset(fr);

    seek(fr,filesize(fr)-1);
     read(fr,byt);
       if byt=15 then begin
       gotoxy (20,3);
       writeln('Necessaria l''installazione tramite INST [Invio]');
       halt;
       end;
    Close(fr);
  End;

  Procedure Nvolte;
   Var
      fr         : file of BYTE;
      NVOLTE     : Byte;
      nome       : String[12];
  Begin
     nome:='BRAIN.exe';    { 1ø ed unico intervento }

     ASSIGN(fr,nome);
          Reset(fr);

     seek(fr,filesize(fr)-2);
     read(fr,NVOLTE);
      if NVOLTE > 0 then begin
         NVOLTE:=NVOLTE-1;
         seek(fr,filesize(fr)-2);
         write(fr,NVOLTE);
         Close(fr);
         end
      else begin
         Close(fr); Halt;
      end;
  End;

   Procedure FormatData;
   Begin
       gotoxy (32,1); textcolor(10);
       writeln(ddata);textcolor(14);
   End;

   Procedure VideoRisultati;
   Begin
      teta:=L; beta:=B;
        Xecl:=R*COS(beta)*COS(teta);
        Yecl:=R*COS(beta)*SIN(teta);
        Zecl:=R*SIN(beta);
 textcolor(0);
  writeln('        L = ',L:15:8,' rad    ç = ',tau,' millenni');
  writeln('        B = ',B:15:8,' rad');
 textcolor(14);
teta:=MODULO(teta,2*pi); if teta<0 then teta:=2*pi+teta;
teta:=teta/rad;write('        LON-ecl= ',teta:9:6,'ø');
beta:=beta/rad;write('    LAT-ecl= ',beta:9:6,'ø');
  writeln('    R= ',R:11:8,' UA');
 textcolor(12);
  writeln('           Coordinate ECLITTICHE del Pianeta (Equin.Data)');
  textcolor(10);
  writeln('        X_ecl= ',Xecl:12:8,'   Y_ecl= ',Yecl:12:8,'   Z_ecl= ',Zecl:12:8);
  textcolor(14);
   End;

   Procedure  FK5Nutazione;
   Var    ome,el0,el1,Lpri,dell,debb          : Real;
   Begin
{ Correz. (Lam,Bet) per passaggio da Eq.DINAM.J2000 a Eq.FK5-J2000  }
      Tg:=10*tau; Lpri:=(-0.00031*Tg-1.397)*Tg+teta/rad; {Lpri in gradi}
      Lpri:=Lpri*rad; {Lpri in RAD}
      dell:=-0.09033+0.03916*(COS(Lpri)+SIN(Lpri))*SIN(beta)/COS(beta); {in Sec}
      debb:=0.03916*(COS(Lpri)-SIN(Lpri)); {in Sec}
      dell:=dell/3600*rad; {in RAD}
      debb:=debb/3600*rad; {in RAD}
   lam:=lam+dell; {Prima Correz. di Lam }
   bet:=bet+debb;
{ Fine Correz. Lam,Bet }
      ome:=((Tg/450000+0.0020708)*Tg-1934.136261)*Tg+125.04452;  {ome in gradi}
      ome:=MODULO(ome,360); if ome<0 then ome:=360+ome;
      el0:=Tg*36000.7698+280.4665;    {l0 in gradi}
      el0:=MODULO(el0,360); if el0<0 then el0:=360+el0;
      el1:=Tg*481267.8813+218.3165;    {l1 in gradi}
      el1:=MODULO(el1,360); if el1<0 then el1:=360+el1;

      ome:=ome*rad; {in RAD}
      el0:=el0*rad; {in RAD}
      el1:=el1*rad; {in RAD}

      ddpsi:=-17.2*SIN(ome)-1.32*SIN(2*el0)-0.23*SIN(2*el1)+0.21*SIN(2*ome); {in Sec}
      ddeps:=9.2*COS(ome)+0.57*COS(2*el0)+0.10*COS(2*el1)-0.09*COS(2*ome); {in Sec}
      ddpsi:=ddpsi/3600*rad; {in RAD}
      ddeps:=ddeps/3600*rad; {in RAD}

   lam:=lam+ddpsi;    {Seconda Correz. di Lam per effetto della Nutazione}
   eps:=eps0+ddeps;   { Correz. di î per effetto della Nutazione}

   End;

   Procedure AlfaDelta;
   Var    ialfaore,ialfamin,alfasec, idgra,idmin,deltasec,
          alfafraz,deltafraz,aux1,aux2,aux3,elon : Real;
          isegno         : Integer;
          ori: String[1];
   Begin
   dgeo:=SQRT(SQR(Xpri)+SQR(Ypri)+SQR(Zpri));
   lam:=ArcoTANG(Ypri,Xpri);
   bet:=ArcoSIN(Zpri/dgeo);
    Tg:=10*tau;

  eps0:=((Tg*0.001813-0.00059)*Tg-46.8150)*Tg+84381.448;  {eps0 in Sec}
  eps0:=eps0/3600*rad;  {eps0 in RAD}
  eps:=eps0;
       if flag=1 then FK5Nutazione;
      aux1:=SIN(lam)*COS(eps)-SIN(bet)/COS(bet)*SIN(eps);
      aux2:=COS(lam);
      aux3:=SIN(bet)*COS(eps)+COS(bet)*SIN(eps)*SIN(lam);
elon:=ArcoCOS(COS(bet)*COS(lam-teta*rad+pi))/rad;  { Elongazione in gradi }
if teta+180>360 then teta:=teta-360;
if lam-teta*rad-pi <0 then ori:='W' ELSE ori:='E'; teta:=0;
   alfa:=ArcoTANG(aux1,aux2); alfa:=alfa/(rad*15);
   delt:=ArcoSIN(aux3)/rad;
   textcolor(12);
         writeln;
 textcolor(12);
  writeln('        Coordinate ECLITTICHE GEOcentriche del Pianeta (Equin.Data)');
  textcolor(10);
  writeln('        X_geo= ',Xpri:12:8,'   Y_geo= ',Ypri:12:8,'   Z_geo= ',Zpri:12:8);
  writeln('                 Longitudine= ',lam/rad:9:5,'ø   Latitudine= ',Bet/rad:9:5,'ø');
  textcolor(14);
         writeln('                Coordinate GEOcentriche POLARI (Equin.Data)');
   textcolor(10);
         writeln('                  DIST_geoc. = ',dgeo:10:8,' UA     Elong.= ',elon:4:1,'ø'+ori);
         writeln('                    à = ',alfa:12:8,' ore;    ë = ',delt:9:7,'ø');
   textcolor(14);
{************************************************************************}
      ialfaore:= INT(alfa);
      alfafraz:= (alfa - ialfaore) * 60;
      ialfamin:= INT(alfafraz);
      alfasec:= (alfafraz - ialfamin) * 60;
      if delt >= 0 then isegno:=1
      else isegno:=-1;
      idgra:= INT(ABS(delt));
      deltafraz:= (ABS(delt) - idgra) * 60;
      idmin:= INT(deltafraz);
      deltasec:= (deltafraz - idmin) * 60;
      idgra:= idgra * isegno;
      writeln('                    à = ',ialfaore:3:0,ialfamin:3:0,alfasec:6:2,
              '         ë = ',idgra:3:0,idmin:3:0,deltasec:5:1);
{************************************************************************ }
   End;

{---------------- INIZIO PROGRAMMA -------------------}
  Begin
   blk:='';
   for uqw:=1 to 79 do
    begin
     blk:=blk+chr(32);
    end;

{    ControlloFlag;
    Nvolte;         }
 Ripeti2:
     CancellaSchermo;
     GiornoGiuliano;

 Ripeti:
   MenuIn;
   ch:=' ';

   if scelta='1' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       MER(tau,L,B,R);
       writeln('                                M E R C U R I O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       MER(tcorr,L,B,R);
       writeln('                                M E R C U R I O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='2' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       VEN(tau,L,B,R);
       writeln('                                  V E N E R E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       VEN(tcorr,L,B,R);
       writeln('                                  V E N E R E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='3' then
    Begin
     CancellaSchermo; FormatData;
      SUN(tau,jd,Xso,Yso,Zso);
       readln;
    end;

   if scelta='4' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       MAR(tau,L,B,R);
       writeln('                                   M A R T E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       MAR(tcorr,L,B,R);
       writeln('                                   M A R T E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='5' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       GIO(tau,L,B,R);
       writeln('                                   G I O V E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       GIO(tcorr,L,B,R);
       writeln('                                   G I O V E');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='6' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       SAT(tau,L,B,R);
       writeln('                                 S A T U R N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       SAT(tcorr,L,B,R);
       writeln('                                 S A T U R N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='7' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       URA(tau,L,B,R);
       writeln('                                   U R A N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       URA(tau,L,B,R);
       writeln('                                   U R A N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='8' then
    Begin
   flag:=0;
     CancellaSchermo; FormatData;
       NET(tau,L,B,R);
       writeln('                                 N E T T U N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tau,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                            SOLO Posizione geometrica');
       readln;
 tluce:=vluce*dgeo; tcorr:=tau-tluce/365250;
   flag:=1;
     CancellaSchermo; FormatData;
       NET(tcorr,L,B,R);
       writeln('                                 N E T T U N O');
       VideoRisultati;
             Xpri:=Xecl; Ypri:=Yecl; Zpri:=Zecl;
       writeln;writeln('                                   T E R R A');
       TER(tcorr,L,B,R);
       VideoRisultati;
             Xpri:=Xpri-Xecl; Ypri:=Ypri-Yecl; Zpri:=Zpri-Zecl;
       AlfaDelta; writeln;
       writeln('                           CON Aberrazione e Nutazione');
       readln;
    end;

   if scelta='9' then
    Begin
    CancellaSchermo;
gotoxy(20,4);   write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(20,5);   write('³                                    ³');
gotoxy(20,6);   write('³          a.. Altra DATA            ³');
gotoxy(20,7);   write('³                                    ³');
gotoxy(20,8);   write('³          b.. USCITA                ³');
gotoxy(20,9);   write('³                                    ³');
gotoxy(20,10);writeln('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');

repeat 
     AB:=readkey;
    if (AB = 'a') or (AB = 'A') then Goto Ripeti2;
    if (AB = 'b') or (AB = 'B') then Exit;
 until (AB<'a') and (AB>'b');

    end;

 Goto Ripeti;
  End.
