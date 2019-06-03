Program TABLE; {--------- (Ordinamento File a Tabella) -----------}

  Uses  Crt;
   Var
      fb, ft      : Text;
      k,j,xxx         : Integer;
      ii         : LongInt;
      a       : String;

{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;

   ASSIGN (FB,'HASH.DAT');      {FILE ORIGINARIO}
   Reset(fb);

   ASSIGN (ft,'01010200.RIS');
   Rewrite(ft);

   FOR j:= 101 TO 200 do begin
        write (ft,j:4);
        write (j:4);

    for k:=1 to 27900 do begin
       readln(fb,a);

       val(copy(a,15,4),ii,xxx);

         if ii=j then begin
         write(copy(a, 8, 6));
         write(ft,copy(a, 8, 6));
         end;

   end;
   writeln; writeln;
   writeln(ft);

   CLOSE (fb);
   ASSIGN (FB,'HASH.DAT');     { rilettura FILE ORIGINARIO }
   Reset(fb);

   end; {for j...}

CLOSE (ft);
  End.
