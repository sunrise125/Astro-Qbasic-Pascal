Program BCHANGE;  {Trasferimento di un file BINARIO
                   in un altro, anche PARZIALMENTE, con MODIFICA }

  Uses  Crt,Dos;
   Var
      fb,fw       : file of Byte;
      k,SUP,INF   : LongInt;
      byt         : Byte;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;
   gotoxy(27,5); write('TRASFERIMENTO INIZIATO');

   ASSIGN (fb,'ORI');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'NEW2');
   Rewrite(fw);

   INF:=1;
   SUP:=filesize(fb);

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP
     do begin
     seek(fb,k-1);
     read(fb,byt);
       if byt=66 then byt:=90;  {Modifica di tutte le B con le Z}
     write(fw,byt);
     gotoxy(27,7); write(k);
   end;
     Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
