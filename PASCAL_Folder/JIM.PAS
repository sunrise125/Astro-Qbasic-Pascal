 Program JIM; {Trasformazione Data Civile in Giorno Giuliano e Viceversa
               + Giorni tra due Date - Aggiornamento di JJ.PAS Set.98 }
 Uses Crt;

 Label  Ripeti, Skip;

 Var
      m                   : Word;
  i,colo,color,xmed,x1,y1,x2,y2,u,xx,yy,y,ncount,nr   : Integer;
  ss,blk,var1,var2 : String;
  scelta  : String[1];
  col,coli : Array[1..4] of Integer;
  labe : Array[1..4] of String[23];
  ch   : Char;
  lab  : Array[1..4] of Char;
  oggi,g,a,b,jd,c,d,e,f,alfa,z,giorno,mese,anno,difgio   : Real;

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
   gotoxy(29,8+2*(i-1));
   writeln(labe[i]);
  end;

 Procedure MenuIn;
  Begin
   lab[1]:='1';
   lab[2]:='2';
   lab[3]:='3';
   lab[4]:='4';
   labe[1]:=' Da Data Civile a G.G.';
   labe[2]:=' Da G.G a Data Civile ';
   labe[3]:=' Giorni tra 2 Date    ';
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
   writeln('Muovere le frecce o selezionare il numero');
   ss:='';
   x1:=23;
   y1:=6;
   x2:=57;
   y2:=15;
   BoxB;
   BordEx;
   textcolor(4);
   textbackground(7);
   gotoxy(xmed-8,y1);
   writeln(' GIORNO GIULIANO ');

    for i:=1 to 4 do
     begin
      colo:=col[i];
      color:=coli[i];
      textbackground(7);
      textcolor(0);
      gotoxy(27,8+2*(i-1));
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
     until (ch=#13) or (ch=lab[1]) or (ch=lab[2]) or (ch=lab[3])
                 or (ch=lab[4]);

     if ch=#13 then str(i:1,scelta)
     else
     begin
      for i:=1 to 4 do
       begin
        if ch=lab[i] then str(i:1,scelta);
       end;
     end;

  End;
{---------------- INIZIO PROGRAMMA -------------------}
  Begin
   blk:='';
   for u:=1 to 79 do
    begin
     blk:=blk+chr(32);
    end;
   ClrScr;

 Ripeti:
   MenuIn;
   ch:=' ';

   if scelta='1' then
    Begin
    textcolor(1);textbackground(7);
     clrscr;
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('³                                    ³');
gotoxy(20,7); write('³                                    ³');
gotoxy(20,8); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
   repeat
   gotoxy(25,5);write('  Giorno.dec = ');readln(g);
   until (g>0) and (g<32);
   repeat
   gotoxy(25,6);write('        Mese = ');readln(m);
   until m in [1..12];
   gotoxy(25,7);write('        Anno = ');readln(y);

          IF m < 3 then begin
           y:=y-1; m:=m+12; end;

          a:=int(y/100); b:=2 - a + int(a/4);
      jd:=int(365.25*(y+4716))+int(30.6001*(m+1))+g+b-1524.5;
      if jd<2299160.5 then jd:=jd-b;   { Data < 15.10.1582--->Rif.Gregoriana}
      oggi:=int(jd+1.5-7*int((jd+1.5)/7));

gotoxy(23,9);writeln('Il Giorno Giuliano della DATA Š : ');
    str(int(jd):7:0,var1);  str(frac(jd):7:5,var2); delete(var2,1,1);
gotoxy(30,10);writeln(var1+var2);
        textbackground(4); textcolor(14); gotoxy(33,11);
        if oggi=0 then writeln(' DOMENICA ');
        if oggi=1 then writeln(' LUNEDI ');
        if oggi=2 then writeln(' MARTEDI ');
        if oggi=3 then writeln(' MERCOLEDI ');
        if oggi=4 then writeln(' GIOVEDI ');
        if oggi=5 then writeln(' VENERDI ');
        if oggi=6 then writeln(' SABATO ');
  readln;
   end;

   if scelta='2' then
   begin
    textcolor(14);textbackground(3);
     clrscr;
gotoxy(20,4); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(20,5); write('³                                    ³');
gotoxy(20,6); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
   gotoxy(25,5);write('Giorno Giuliano=? ');readln(jd);

        z:=int(jd+0.5); f:=jd+0.5-z;
          IF z >= 2299161 then
            begin
                 alfa:=int((z-1867216.25)/36524.25);
                 a:=z+1+alfa-int(alfa/4);
            end
           else
            begin
                 a:=z;
          end;
     b:=a+1524; c:=int((b-122.1)/365.25);
     d:=int(365.25*c); e:=int((b-d)/30.6001);

  giorno:=b-d-int(30.6001*e)+f;
          IF E<13.5 then mese:=e-1 else mese:=e-13;
          IF mese>2.5 then anno:=c-4716 else anno:=c-4715;

gotoxy(25,8);writeln('La DATA Š : ',giorno:9:5,'/',mese:2:0,'/',anno:4:0);
  readln;
   end;

   if scelta='3' then    { begin exit; end; }
    Begin
    ncount:=0; clrscr;
Skip:
    textcolor(1);textbackground(7);
    if ncount=0 then nr:=0 else nr:=7;
gotoxy(20,4+nr); write('ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿');
gotoxy(20,5+nr); write('³                                    ³');
gotoxy(20,6+nr); write('³                                    ³');
gotoxy(20,7+nr); write('³                                    ³');
gotoxy(20,8+nr); write('ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ');
    ncount:=ncount+1;
    if ncount=1 then begin
        gotoxy(32,4); write(' Prima Data '); end
         else begin gotoxy(32,4+nr); write(' Seconda Data '); end;
   repeat
   gotoxy(25,5+nr);write('  Giorno.dec = ');readln(g);
   until (g>0) and (g<32);
   repeat
   gotoxy(25,6+nr);write('        Mese = ');readln(m);
   until m in [1..12];
   gotoxy(25,7+nr);write('        Anno = ');readln(y);
          IF m < 3 then begin
           y:=y-1; m:=m+12; end;

          a:=int(y/100); b:=2 - a + int(a/4);
      jd:=int(365.25*y)+int(30.6001*(m+1))+g+b+1720994.5;
      if jd<2299160.5 then jd:=jd-b;   { Data < 15.10.1582--->Rif.Gregoriana}

      if ncount=1 then begin
          difgio:=jd;
          goto Skip; end
           else
          difgio:=jd-difgio;
          gotoxy(24,17);
          writeln('Differenza giorni= ',difgio:10:5);
          readln;
   End;

   if scelta='4' then begin  exit;  end;

   goto Ripeti;
  End.
