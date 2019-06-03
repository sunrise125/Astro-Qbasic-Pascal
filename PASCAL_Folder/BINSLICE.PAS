Program BINSLICE; { Trasferimento di un file BINARIO
                in un altro, anche PARZIALMENTE}

  Uses  Crt,Dos;

   Var
      fb,fw       : file of Byte;
      k           : LongInt;
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

  {Modificare a multipli di 91 Bytes gli ESTREMI di For k:=.... partendo dalla
   fine del file come indicato sotto, oppure dalla testa: For k:=1 to 546,
   oppure dal centro: For k:= 911 to 1456 etc..}

   { For k:=911 to 1456 }

    For k:=831146 to filesize(fb)-11776
     do begin
     seek(fb,k-1);
     read(fb,byt);
     write(fw,byt);
     write('.');
   end;
     Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
