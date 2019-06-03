Program PASQUA;  {Data della Pasqua cristiana: formula di Gauss
                  riveduta e ampliata, formula UNIVERSALE}

  Uses  Crt;

   Var
      y,k,p,q,f,a,b,c,d,e,M,N, Pq  : Integer;
      wo                       : String;

  {INIZIO Programma}
   Begin
     textcolor(14);textbackground(1);clrscr;
   gotoxy(32,3); write('Anno: '); readln(y);

    k:=y div 100;  a:=y mod 19;  b:=y mod 4;  c:=y mod 7;
    p:=(8*k+13) div 25; q:= k div 4;
        if y>1582 then begin
          M:=(15+k-p-q) mod 30;
          N:=(4+k-q) mod 7; end
          else begin
          M:=15; N:=6;
        end;
    d:=(19*a+M) mod 30;
    e:=(2*b+4*c+6*d+N) mod 7;
    f:=(11*d+19*e+a) div 433;

    writeln;
  {  writeln('  ',k:2,a:4,b:4,c:4,p:4,q:4,M:5,N:4,d:5,e:4,f:4); }
    wo:=' Marzo';
    Pq:=22+d+e-7*f;
       if Pq>31 then begin
          Pq:=Pq-31;
          wo:=' Aprile';
       end;

    gotoxy(26,5);
    writeln('Data Pasqua= ',Pq:2,wo);

  End.

   {FINE Programma}
