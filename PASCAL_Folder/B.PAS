Program B;      {ex BINSLICE: Trasferimento di un file BINARIO
                in un altro, anche PARZIALMENTE}

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

   ASSIGN (fb,'AAA.EXE');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'AAA_FIN.ZIP');
   Rewrite(fw);

   INF:=570881;
   SUP:=800569; {filesize(fb)-100000;}

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP
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
