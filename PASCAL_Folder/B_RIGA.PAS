Program B_RIGA; { Trasferimento di un file BINARIO MONORIGA in un altro
                  SPEZZATO in tante righe, anche PARZIALMENTE}

  Uses  Crt,Dos;
   Const Lung=194;    {lunghezza del Record}
   Var
      fb,fw       : file of Byte;
      k,SUP,INF   : LongInt;
      byt,bytt    : Byte;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;
   gotoxy(27,5); write('TRASFERIMENTO INIZIATO');

   ASSIGN (fb,'c:\mps\1000.ele');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (fw,'dia1000.txt');
   Rewrite(fw);

 {  INF:=1;
   SUP:= filesize(fb)-100000;}

   INF:=1;
   SUP:=194000; {filesize(fb);    590;}

 textcolor(10);gotoxy(43,7); write(SUP);textcolor(14);
    For k:=INF to SUP  do begin
     seek(fb,k-1);
        if (k mod Lung = 56) then begin
        read(fb,byt);
         write(fw,byt);
        end;
        if (k mod Lung = 57) then begin
        read(fb,byt);
         write(fw,byt);
        end;
        if (k mod Lung = 58) then begin
        read(fb,byt);
         write(fw,byt);
        end;
        if (k mod Lung = 59) then begin
        read(fb,byt);
         write(fw,byt);
        end;
        if (k mod Lung = 60) then begin
        read(fb,byt);
         write(fw,byt);
        end;
     gotoxy(27,7); write(k);
        if (k mod Lung = 0) then begin
          bytt:=10;                    {Alt+10= a capo}
          write(fw,bytt);
        end;
   end;

    Close(fb);
    Close(fw);

   gotoxy(27,5); write('TRASFERIMENTO ULTIMATO');
  End.

   {FINE Programma}
