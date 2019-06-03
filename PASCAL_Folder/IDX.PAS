Program IDX; { Trasferimento di un file BINARIO   x file *.IDX
               in un altro, anche PARZIALMENTE}

  Uses  Crt,Dos;
   Label Skip;
   Var
      fb,fw       : file of Byte;
      k,SUP,INF   : LongInt;
      byt,bytt    : Byte;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;
   gotoxy(27,5); write('TRASFERIMENTO INIZIATO');

   ASSIGN (fb,'c:\mps\asteroid.idx');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'d:\asteroid.lst');
   Rewrite(fw);

  { INF:=filesize(fb)-2799;
   SUP:=filesize(fb);  }

   INF:=1;
   SUP:=filesize(fb);

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP  do begin
     if ((k mod 28 >19 ) and (k mod 28 <= 28)) or (k mod 28 <=4 )
          or (k mod 28 = 0) then begin
          goto Skip; end;
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
Skip:
        if (k mod 28 = 0) then begin
          bytt:=10;
          write(fw,bytt);
        end;
   end;

    Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
