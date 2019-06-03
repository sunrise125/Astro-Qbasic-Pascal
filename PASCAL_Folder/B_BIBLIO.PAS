Program B_BIBLIO; { Trasferimento di un file BINARIO   x file BIBLIO
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

   ASSIGN (fb,'c:\sergio\asteroid.ele');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'short.ele');
   Rewrite(fw);

   INF:=1;
   SUP:=2914632; {filesize(fb);}

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP  do begin

     if ((k mod 168 >0 ) and (k mod 168 <= 8))
        or ((k mod 168>66) and (k mod 168 <=74))
        or ((k mod 168 >82) and (k mod 168 <=140))
           then begin
          goto Skip; end;
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
Skip:
      {  if (k mod 168 = 0) then begin
          bytt:=10;
          write(fw,bytt);
        end;                     }
   end;

    Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
