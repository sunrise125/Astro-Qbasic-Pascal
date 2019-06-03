Program SLICE;     {Taglia a fette un file ASCII e ne ricostruisce un altro
                    a campi ridotti}

  Uses  Crt;
   Var
      fl,ft       : Text;
      designa           : String[7];
      epoca            : String[6];
      H,G,M,peri,nodo,ecc,incl,motog,a   :Real;
      dummy             : String[62];
      nome,aux          : String[27];

      k           : Integer;
{-------------------------------------------------------------------------------
   INIZIO Programma}
   Begin
     textcolor(14);textbackground(6);clrscr;
   ASSIGN (fl,'AST_01.TXT');
   Reset(fl);
   ASSIGN (ft,'AST_01.RIS');
   Rewrite(ft);
 k:=0; aux:='              ';
    while not eof(fl) do begin
    k:=k+1;
    readln(fl,designa,H,G,epoca,M,peri,nodo,incl,ecc,motog,a,dummy,nome);
    writeln(k:5,nome);
    writeln(ft,copy(nome+aux,2,24),H:5:2,G:6:2,epoca,M:10:5,peri:10:5,nodo:10:5,incl:9:5,ecc:10:7,motog:11:8,a:11:7);
    end;
   Close(fl);
   Close(ft);
end.
   {FINE Programma}
