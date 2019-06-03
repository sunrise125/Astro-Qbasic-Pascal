Program CRYPTIC;   { Trasferimento di un file BINARIO
                     in un altro, ma con correzione del Codice
                     x file musicale nome.NWC
                     Deve essere presente nella DIr il file CAP.NWC di
                     68 bytes}

  Uses  Crt,Dos;

   Var
      wo,wo1      : String;
      fb,fb1,fw   : file of Byte;
      k,SUP,INF   : LongInt;
      Cod3        : Integer;
      byt         : Byte;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;

   gotoxy(16,2);write('Sigla FILE Originario, senza estens. [.NWC]: '); read(wo);
   readln;
   gotoxy(16,3);write('Sigla FILE Nuovo     , senza estens. [.NWC]: '); read(wo1);
   gotoxy(27,5); write('TRASFERIMENTO INIZIATO');

   ASSIGN (fb1,'CAP.NWC');    { File di Testa BINARIO }
   Reset(fb1);

   ASSIGN (fb,wo+'.NWC');     { File BINARIO originario }
   Reset(fb);

   ASSIGN (fw,wo1+'.NWC');
   Rewrite(fw);

 k:=55; seek(fb,k-1); read(fb,byt);
      if byt=102 then Cod3:=69;
      if byt=103 then Cod3:=66;
      if byt=104 then Cod3:=67;
      if byt=121 then Cod3:=68;
      if byt=122 then Cod3:=69;
      if byt=97  then Cod3:=70;
      if byt=114 then Cod3:=65;
      if byt=115 then Cod3:=66;
      if byt=116 then Cod3:=69;
      if byt=107 then Cod3:=70;

    For k:=1 to 68
     do begin
     seek(fb1,k-1);
     read(fb1,byt);
     write(fw,byt);
     gotoxy(27,7); write(k);
   end;
     Close(fb1);

   INF:=Cod3;
   SUP:=filesize(fb);

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
