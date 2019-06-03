Program BLA; { Trasferimento di un file BINARIO
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

   ASSIGN (fb,'c:\tcmicro\9142.ele');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'N9142.ele');
   Rewrite(fw);

 {  INF:=1;
   SUP:=2560;  391680; filesize(fb)-100000;}

   INF:=1;
   SUP:=filesize(fb);

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP  do begin
     if (k mod 194 >144 ) or
        ((k mod 194 >61 ) and (k mod 194 <= 90)) or
        (k mod 194 = 0 ) then begin
          goto Skip; end;
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
Skip:
   end;

    Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
