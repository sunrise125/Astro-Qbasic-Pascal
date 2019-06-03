Program ASC2BIN; {Trasferimento di un file ASCII in BINARIO applicato agli
                  elementi degli ASTEROIDI, secondo lo schema di Augusto Testa
                     PROGRAMMA GENERALE per qualsiasi tipo di Dbase:
                        a) STANDARD MPC (Minor Planet Center) - AUTOMATICO
                        b) NON-Standard, con Epoca FISSA NON specificata nel
                           Dbase e posizione dei campi da definire come INPUT
                 }
  Uses  Crt;
   Var
      fb,ft          : Text;
      ColSigla,ColM,ColPeri,ColNodo,ColIncl,
      ColEcc,Colaint,Colafrac,ColH,ColG: Integer;
      k,j,vv,xx,sce,sce1,fla : Integer;
      Epoca,aux,Anomedia,Peri,Nodo,Incl,e,a,afrac,H,G : LongInt;
      c,c1,c3,c4     : String[1];
      c2             : String[2];
      outbin         : String[12];
      s,total,k1,nomecat     : String;
      aux1           : Real;
      epB,epA,ep9,ep8,m4,m3,m2,m1 : LongInt;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(1);clrscr;
  End;

  Procedure SelezNumProv;
   Begin
            gotoxy(1,15);
writeln(' ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
writeln(' ³           1) Pianetini NUMERATI         2) Pianetini PROVVISORI            ³');
writeln(' ³                                 Scelta:                                    ³');
writeln(' ³                                                                            ³');
writeln(' ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(44,17); readln(sce1); until sce in [1..2];
         if (sce1=1) then begin
         fla:=1; outbin:='ASTNUM2.CAT';
         end;
         if (sce1=2) then begin
         fla:=2; outbin:='ASTPRO2.CAT';
         end;
  End;

  Procedure Opzioni;
   Begin
gotoxy(1,2);
writeln(' ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
writeln(' ³             Trasformazione di un DBASE di Asteroidi da file ASCII          ³');
writeln(' ³                  a file BINARIO (modello MAPPA2 di A.Testa)                ³');
writeln(' ÃÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ´');
writeln(' ³                                                                            ³');
writeln(' ³             1) STANDARD MPC (Automatico)                                   ³');
writeln(' ³                                                                            ³');
writeln(' ³             2) NON-Standard, con Epoca FISSA NON specificata nel Dbase     ³');
writeln(' ³                e posizione dei campi da definire in fase di INPUT          ³');
writeln(' ³                                                                            ³');
writeln(' ³                                 Scelta:                                    ³');
writeln(' ³                                                                            ³');
writeln(' ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
repeat gotoxy(44,12); readln(sce); until sce in [1..2];
         if sce=1 then begin
gotoxy(1,6);
writeln(' ³                                       N_Col   Lar.|          N_Col   Lar.  ³');
writeln(' ³    1) STANDARD MPC:    Sigla Corpo      1      7  |   e       71      9    ³');
writeln(' ³         Dbase del            Epoca     21      5  |  a_int    94      2    ³');
writeln(' ³       Minor Planet             M       27      9  |  a_frac   97      7    ³');
writeln(' ³          Center                w       38      9  |   H        9      5    ³');
writeln(' ³        MPCORB.DAT              ê       49      9  |   G       18      2    ³');
writeln(' ³                                i       61      8  |                        ³');

      ColSigla:=1;

           SelezNumProv;
          end;

         if sce=2 then begin
gotoxy(1,6);
writeln(' ³                                                                            ³');
writeln(' ³    1) NON Standard:                                                        ³');
writeln(' ³       (Epoca fissa)                                                        ³');
writeln(' ³                                                                            ³');
writeln(' ³                                                                            ³');
writeln(' ³                                                                            ³');
writeln(' ³                                                                            ³');
       gotoxy(28,6); write('Epoca [aaaammnn]= ');
       gotoxy(28,7); write('  Es.  19971218');
       gotoxy(46,6); readln(Epoca);

       gotoxy(60,6); write('Col_Sigla= '); readln(Colsigla);
       gotoxy(60,7); write('Col_M    = '); readln(ColM);
       gotoxy(60,8); write('Col_Peri = '); readln(ColPeri);
       gotoxy(28,9); write('Col_Nodo = '); readln(ColNodo);
       gotoxy(44,9); write('Col_Incl = '); readln(ColIncl);
       gotoxy(60,9); write('Col_Ecc  = '); readln(ColEcc);
      gotoxy(28,10); write('Col_aInt = '); readln(Colaint);
      gotoxy(44,10); write('Col_aFrac= '); readln(Colafrac);
      gotoxy(28,11); write('Col_H    = '); readln(ColH);
      gotoxy(44,11); write('Col_G    = '); readln(ColG);
           SelezNumProv;
          end;

 gotoxy(24,18); write('Nome Catalogo: '); readln(nomecat);
  End;

  Procedure SiglaEpoca;
  Begin

     c1:=copy(total,21,1);   {Sigla centinaia di anni}
          if (c1='J') then Epoca:=19*1000000;
          if (c1='K') then Epoca:=20*1000000;
  {----------------------------------------------------}
     c2:=copy(total,22,2);   {Sigla anni}
          val(c2,aux,xx);
       Epoca:=Epoca+aux*10000;
  {----------------------------------------------------}
     c3:=copy(total,24,1);   {Sigla mese}
     if (c3='1') or (c3='2') or (c3='3') or (c3='4') or (c3='5')
        or (c3='6') or (c3='7') or (c3='8') or (c3='9') then begin
          val(c3,aux,xx);
       Epoca:=Epoca+aux*100;
     end;
     if (c3='A') then Epoca:= Epoca+10*100;
     if (c3='B') then Epoca:= Epoca+11*100;
     if (c3='C') then Epoca:= Epoca+12*100;
  {----------------------------------------------------}
     c4:=copy(total,25,1);   {Sigla giorno}
     if (c4='1') or (c4='2') or (c4='3') or (c4='4') or (c4='5')
        or (c4='6') or (c4='7') or (c4='8') or (c4='9') then begin
          val(c4,aux,xx);
       Epoca:=Epoca+aux;
     end;
     if (c4='A') then Epoca:= Epoca+10;
     if (c4='B') then Epoca:= Epoca+11;
     if (c4='C') then Epoca:= Epoca+12;
     if (c4='D') then Epoca:= Epoca+13;
     if (c4='E') then Epoca:= Epoca+14;
     if (c4='F') then Epoca:= Epoca+15;
     if (c4='G') then Epoca:= Epoca+16;
     if (c4='H') then Epoca:= Epoca+17;
     if (c4='I') then Epoca:= Epoca+18;
     if (c4='J') then Epoca:= Epoca+19;
     if (c4='K') then Epoca:= Epoca+20;
     if (c4='L') then Epoca:= Epoca+21;
     if (c4='M') then Epoca:= Epoca+22;
     if (c4='N') then Epoca:= Epoca+23;
     if (c4='O') then Epoca:= Epoca+24;
     if (c4='P') then Epoca:= Epoca+25;
     if (c4='Q') then Epoca:= Epoca+26;
     if (c4='R') then Epoca:= Epoca+27;
     if (c4='S') then Epoca:= Epoca+28;
     if (c4='T') then Epoca:= Epoca+29;
     if (c4='U') then Epoca:= Epoca+30;
     if (c4='V') then Epoca:= Epoca+31;
  end;

  Procedure Trasforma(xin:LongInt; VAR n1,n2,n3,r3: LongInt);
  Var
     r1,r2: LongInt;
  Begin
   n1:=trunc(int(xin/16777216));  r1:=trunc(16777216*frac(xin/16777216));
   n2:=trunc(int(r1/65536));      r2:=trunc(65536*frac(r1/65536));
   n3:=trunc(int(r2/256));        r3:=trunc(256*frac(r1/256));
  end;

  Procedure Trasforma1(xin:LongInt; VAR n3,r3: LongInt);
  Begin
   n3:=trunc(int(xin/256));        r3:=trunc(256*frac(xin/256));
  end;

{-------------------------------------------------------------------------------
   INIZIO Programma}
 Begin
   ColoriSchermo;
   Opzioni;

   gotoxy(27,21); write('TRASFERIMENTO INIZIATO');

   ASSIGN (ft,nomecat);  {File ASCII originario}
   Reset(ft);

   ASSIGN (fb,outbin);      {File BIN finale}
   Rewrite(fb);

 while not eof(ft) do begin
          write('.');
        readln(ft,total);
     if fla=2 then begin
             For k:=ColSigla to ColSigla+6 do begin
                c:=copy(total,k,1);
                s:=s+c;
              end;
     end;

     if fla=1 then begin
                s:=s+#255#255;
             For k:=ColSigla to ColSigla+4 do begin
                c:=copy(total,k,1);
                s:=s+c;
              end;
     end;

     if sce=1 then SiglaEpoca;

 {Epoca in BINARIO}
    Trasforma(Epoca,epB,epA,ep9,ep8);
     vv:=ord(ep8);  s:=s+chr(vv);
     vv:=ord(ep9);  s:=s+chr(vv);
     vv:=ord(epA);  s:=s+chr(vv);
     vv:=ord(epB);  s:=s+chr(vv);

 {Anomedia in BINARIO}
     if sce=1 then ColM:=27;
     k1:=copy(total,ColM,9);
          val(k1,aux1,xx); Anomedia:=trunc(aux1*1e6);
    Trasforma(Anomedia,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Arg_Perielio in BINARIO}
     if sce=1 then ColPeri:=38;
     k1:=copy(total,ColPeri,9);
          val(k1,aux1,xx); Peri:=trunc(aux1*1e6);
    Trasforma(Peri,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Longit_Nodo in BINARIO}
     if sce=1 then ColNodo:=49;
     k1:=copy(total,ColNodo,9);
          val(k1,aux1,xx); Nodo:=trunc(aux1*1e6);
    Trasforma(Nodo,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Inclinaz. in BINARIO}
     if sce=1 then ColIncl:=61;
     k1:=copy(total,ColIncl,8);
          val(k1,aux1,xx); Incl:=trunc(aux1*1e6);
    Trasforma(Incl,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Eccentr. in BINARIO}
     if sce=1 then ColEcc:=71;
     k1:=copy(total,ColEcc,9);
          val(k1,aux1,xx); e:=trunc(aux1*1e8);
    Trasforma(e,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Parte intera Semid. (a) in BINARIO}
     if sce=1 then Colaint:=94;
     k1:=copy(total,Colaint,2);
          val(k1,aux1,xx); a:=trunc(aux1);
    Trasforma1(a,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);

 {Parte fraz. Semid. (a) in BINARIO}
     if sce=1 then Colafrac:=97;
     k1:=copy(total,Colafrac,7);
          val(k1,aux,xx); afrac:=aux*10;
    Trasforma(afrac,m4,m3,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);
     vv:=ord(m3);  s:=s+chr(vv);
     vv:=ord(m4);  s:=s+chr(vv);

 {Magnit. H in BINARIO}
     if sce=1 then ColH:=9;
     k1:=copy(total,ColH,4);
          val(k1,aux1,xx); H:=trunc(aux1*1e2);
    Trasforma1(H,m2,m1);
     vv:=ord(m1);  s:=s+chr(vv);
     vv:=ord(m2);  s:=s+chr(vv);

 {Magnit. G in BINARIO}
     if sce=1 then ColG:=18;
     k1:=copy(total,ColG,2);
          val(k1,aux,xx); G:=aux;
     vv:=ord(aux);  s:=s+chr(vv);

         { writeln(s);}
     write(fb,s);
   s:='';
 end; { ... while not}
     Close(fb);
    Close(ft);

   gotoxy(27,21); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
