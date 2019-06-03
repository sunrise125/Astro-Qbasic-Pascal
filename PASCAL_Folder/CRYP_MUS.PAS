Program CRYP_MUS;  { Trasferimento di un file BINARIO
                     in un altro, ma con correzione del Codice
                     x file musicale nome.NWC}

  Uses  Crt,Dos;

   Var
      wo,wo1      : String;
      fb,fw       : file of Byte;
      k,SUP,INF   : LongInt;
      byt         : Byte;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;

   gotoxy(16,2);write('Sigla FILE Originario (max 7 lettere, senza 8.NWC): '); read(wo);
   readln;
   gotoxy(16,3);write('Sigla FILE Nuovo      (max 7 lettere, senza .NWC): '); read(wo1);
   gotoxy(27,5); write('TRASFERIMENTO INIZIATO');

   ASSIGN (fb,wo+'8.NWC');     { File BINARIO originario }
   Reset(fb);

   ASSIGN (fw,wo1+'.NWC');
   Rewrite(fw);

   INF:=1;
   SUP:=filesize(fb);

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=1 to 49
     do begin
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
   end;
  k:=50; byt:=2; write(fw,byt);     gotoxy(27,7); write(k);
  k:=51; byt:=0; write(fw,byt);     gotoxy(27,7); write(k);
  k:=52; byt:=0; write(fw,byt);     gotoxy(27,7); write(k);

  k:=53; byt:=65; write(fw,byt);     gotoxy(27,7); write(k);
  k:=54; byt:=97; write(fw,byt);     gotoxy(27,7); write(k);
  k:=55; byt:=102; write(fw,byt);     gotoxy(27,7); write(k);
  k:=56; byt:=117; write(fw,byt);     gotoxy(27,7); write(k);
  k:=57; byt:=110; write(fw,byt);     gotoxy(27,7); write(k);
  k:=58; byt:=113; write(fw,byt);     gotoxy(27,7); write(k);
  k:=59; byt:=120; write(fw,byt);     gotoxy(27,7); write(k);

    For k:=60 to 68
     do begin
     byt:=0;
     write(fw,byt);
     gotoxy(27,7); write(k);
   end;

    For k:=69 to SUP
     do begin
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
   end;

     Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
