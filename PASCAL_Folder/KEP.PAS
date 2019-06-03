Program KEP;   {Soluzione dell'EQUAZIONE di KEPLERO: E = M + e * SIN(E) }

  Uses  Crt;
  Const rad=pi/180; { Trasforma gradi in radianti }
        Eps=1e-10;  { Precisione richiesta }
  Var M, Mgra, e, anomecc, variaz: Real;

   Begin
   {............ Dati del problema ........... }
   Mgra:= 149.7463873; { Anomalia media in Gradi }
   e:= 0.3978305;      { Eccentricit…, adimens. }

   M:=Mgra*rad;  { Anomalia media in Radianti }
   anomecc:= M;  {Valore iniziale di E assunto pari all'anom.media M}
   { Ciclo iterativo: si interrompe quando variaz<Eps }
 repeat
    variaz:= (M + e * SIN(anomecc) - anomecc) / (1 - e * COS(anomecc));
    anomecc:= anomecc + variaz;
 until abs(variaz) < Eps;

   writeln('     Dati:   M= ',Mgra:11:7,' Gradi   e= ',e:9:7);
   writeln('     Risultato (anomalia eccentrica E)');
   writeln('     E= ',anomecc:10:7,' Rad = ',anomecc/rad:11:7,' Gradi');
  End.
