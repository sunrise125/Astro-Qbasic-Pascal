Program REC_A;       {Lettura Byte per Byte di un file BINARIO
                     e trasformazione in ASCII tramite la funz. ORD
                     Lettura file con campi da 17 bytes.

                     Passaggio da BIN ad ASCII
                     =========================
           dopo aver ... correttamente 'decriptato' le GRANDEZZE del
           file Binario

 Lettura e manipolazione dati del file esterno PPM1A.CAT
 -------------------------------------------------------
     }

  Uses  Crt,Dos;

   Var
      fb          : file of Byte;
      ft      : Text;
      k,j,SUP,INF,nro   : LongInt;
      quo,res,som: Integer;
      byt         : Byte;
      nro1,nro2,nro3,nro4 :Word;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;
   gotoxy(22,1); writeln('Lettura FILE Binario INIZIATA');

   ASSIGN (fb,'a:\PPM1A.CAT');      {File BINARIO originario}
   Reset(fb);

   ASSIGN (ft,'HASH.DAT');
   Rewrite(ft);

   INF:=1;
   SUP:=27900; {27900 n.ro elementi}

   For k:=INF to SUP do begin
    For j:=1 to 17 do begin
     seek(fb,17*(k-1)+j-1);
     read(fb,byt);

          if j=1 then nro1:=ord(byt);
          if j=2 then nro2:=ord(byt);
          if j=3 then nro3:=ord(byt);
          if j=4 then nro4:=ord(byt);

    end;
      nro:=256*256*(256*nro2+nro1)+(256*nro4+nro3);
      quo:=trunc(nro/1000); res:=nro-quo*1000;
      som:=quo+res;

  {nro stella, nro ordinale, nro chiave = somma migliaia + fraz. migliaia
                                        = es. 185325 = 185+325 = 510 = som}

      write(ft,nro:7);  write(ft,k:6);  writeln(ft,som:5);
      write(nro:7);  write(k:6);  writeln(som:5);


   end;

     Close(fb);
     Close(ft);

 writeln('                     Lettura FILE Binario ULTIMATA');
  End.

   {FINE Programma}
