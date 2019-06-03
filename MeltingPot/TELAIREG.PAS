Program TELAIREG; { TELAI piani REGOLARI: travi orizzontali e pilastri vertic.
                                CARICHI Verticali DISTRIBUITI nelle travi
                                e FORZE Orizzontali CONCENTRATE ai nodi,
                                pilastri incastrati alla base
   Metodo delle DEFORMAZIONI }

  Uses Crt;
    Const
         maxB=8; maxS=8;
         Ey=1;  {Modulo elastico fittizio}
   Var
       k,i,j,O,Pp,Rr,Qq,Tt   : Integer;
       Ba,St          : Integer;  { Ba=bay, campata;  St=storey, piano }
       wo                        : String;
       ft                        : Text;

       W                : Array[1..maxB] of Real;
       X,V              : Array[1..maxS] of Real;
       P                : Array[1..maxS*(maxB+2)] of Real;
       E                : Array[1..maxS*(maxB+2),1..maxS*(maxB+2)] of Real;
       L,U,Y            : Array[1..maxS,1..maxB] of Real;
       M,Z              : Array[1..maxS,1..maxB+1] of Real;
       Q,S,R,T          : Array[1..maxS+1,1..maxB+1] of Real;



  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure MainFrame;
  begin
    clrscr;
    writeln;
    writeln('   Calcolo automatico di STRUTTURE INTELAIATE a sezione costante, senza elemen-');
    writeln('   ti inclinati, con carichi verticali distribuiti in campata, orizzontali con-');
    writeln('   centrati ai nodi ed incastri al suolo  - Metodo delle DEFORMAZIONI - Dimens.');
    writeln('   massimo dei vettori:  ->  Numero Campate = 8;  Numero Piani = 8');
    textcolor(10);
    writeln('   Attn: il programma presenta un BACO nel portale semplice (1 camp. e 1 piano)');
    textcolor(14);
  end;

     Procedure DATItelaio;
   Begin
     textcolor(10);
      gotoxy(32,1);write('Sigla TELAIO: '); read(wo);
     textcolor(14);
   textcolor(12);
   gotoxy(32,7); writeln('INSERIMENTO DATI'); textcolor(14);
   gotoxy(26,8); write('B = Numero delle Campate => '); read(Ba);
   gotoxy(26,9); write('S = Numero dei   Piani   => '); read(St);
   textcolor(12);
   gotoxy (32,11); writeln('(999 per ripetere)');
   textcolor(10);
   gotoxy (28,12); writeln('Luci_Campate - Vettore W');
   textcolor(14);
       FOR k:= 1 TO Ba do begin
       gotoxy(35,13); write('Luce(',k:2,')= '); read(W[k]);
       gotoxy(44,13); write('              ');
       IF W[k] = 999 THEN k:= k - 2;
       end;
   gotoxy (28,12); writeln('                          ');
   textcolor(10);
   gotoxy (28,12); writeln('Altezze_Piani - Vettore X ');
   textcolor(14);
       FOR k:= 1 TO St do begin
       gotoxy(35,13); write('Alt.(',k:2,')= '); read(X[k]);
       gotoxy(35,13); write('              ');
       IF X[k] = 999 THEN k:= k - 2;
       end;

   gotoxy (28,12); writeln('                           ');
   textcolor(10);
   gotoxy (28,12); writeln('Mom.Inerz_Travi - Vettore L');
   textcolor(14);
       FOR i:= 1 TO St do begin
          gotoxy (28,13); write('Piano(',i:2,')');

            FOR j:= 1 TO Ba do begin
              gotoxy (38,13); write('Camp.(',j:2,')');
              gotoxy(50,13); write('J= '); read(L[i,j]);
              gotoxy(53,13); write('              ');
              IF L[i,j] = 999 THEN j:= j - 2;
            end;
       end;

   gotoxy (28,12); writeln('                           ');
   gotoxy (28,13); writeln('                           ');
   textcolor(10);
   gotoxy (25,12); writeln('Mom.Inerz_Pilastri - Vettore M');
   textcolor(14);
       FOR i:= 1 TO St do begin
          gotoxy (19,13); write('Sotto_Piano(',i:2,')');

            FOR j:= 1 TO Ba+1 do begin
              gotoxy (35,13); write('Pilastro(',j:2,')');
              gotoxy(50,13); write('J= '); read(M[i,j]);
              gotoxy(53,13); write('              ');
              IF M[i,j] = 999 THEN j:= j - 2;
            end;
       end;

   gotoxy (25,12); writeln('                           ');
   gotoxy (15,13); writeln('                                                ');
   textcolor(10);
   gotoxy (25,12); writeln('Carichi DISTRIBUITI - Vettore U');
   textcolor(14);
       FOR i:= 1 TO St do begin
          gotoxy (28,13); write('Piano(',i:2,')');

            FOR j:= 1 TO Ba do begin
              gotoxy (38,13); write('Camp.(',j:2,')');
              gotoxy(50,13); write('q= '); read(U[i,j]);
              gotoxy(53,13); write('              ');
              IF U[i,j] = 999 THEN j:= j - 2;
            end;
       end;

   gotoxy (25,12); writeln('                           ');
   gotoxy (25,13); writeln('                                            ');
   textcolor(10);
   gotoxy (22,12); writeln('Spinte CONCENTRATE Orizz. - Vettore V');
   textcolor(14);
       FOR i:= 1 TO St do begin
          gotoxy (28,13); write('Piano(',i:2,')');
              gotoxy(43,13); write('              ');
              gotoxy(40,13); write('H= '); read(V[i]);
              IF V[i] = 999 THEN i:= i - 2;
       end;
  End;


     Procedure MatriceRigidezze;

   Label Skip1,Skip2,Skip3,Skip4,Skip5,Skip6,Skip7,Skip8,
         Skip9,Skip10,Skip11,Skip12;
   Begin
   ColoriSchermo;
     { Azzeramento VETTORI }
       FOR i:= 1 TO St*(Ba+2) do begin
       FOR j:= 1 TO St*(Ba+2) do begin
            E[i,j]:=0;
       end;
       end;

       FOR i:= 1 TO St*(Ba+2) do begin
            P[i]:=0;
       end;

       FOR i:= 1 TO St+1 do begin
       FOR j:= 1 TO Ba+1 do begin
            Q[i,j]:=0; S[i,j]:=0; R[i,j]:=0;  T[i,j]:=0;
       end;
       end;

     { RIGIDEZZE travi e pilastri }

       FOR i:= 1 TO St do begin
       FOR j:= 1 TO Ba do begin
            Y[i,j]:= Ey*L[i,j]/(W[j]);     {travi}
       end;
       end;

       FOR i:= 1 TO St do begin
       FOR j:= 1 TO Ba+1 do begin
            Z[i,j]:= Ey*M[i,j]/(X[i]);     {pilastri}
       end;
       end;

     { Momenti d'incastro PERFETTO delle travi }
       FOR i:= 1 TO St do begin
       FOR j:= 1 TO Ba do begin
            T[i,j+1]:= U[i,j]*W[j]*W[j]/12;     {F.E.M.}
            R[i,j]:= -T[i,j+1];                 {momenti a destra}
       end;
       end;

     { Membri destri equaz. carico-spostamento }
       FOR i:= 1 TO St do begin
          FOR j:= 1 TO Ba+1 do begin
             O:=(i-1)*(Ba+1)+j;
             P[O]:= -(T[i,j]+R[i,j]);
          end;
       end;

       FOR i:= 1 TO St do begin
         P[St*(Ba+1)+i]:= V[i];
       end;

     { Inizio Calcoli sul TELAIO }

       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba+1 do begin
         O:=(i-1)*(Ba+1)+j;
         Pp:=O+1; Rr:=O-1; Qq:=(i-2)*(Ba+1)+j; Tt:=i*(Ba+1)+j;
           if i=1 then goto skip1;
              E[O,O]:=E[O,O]+4*Z[i-1,j];
              E[O,Qq]:=2*Z[i-1,j];
Skip1:  E[O,O]:=E[O,O]+4*Z[i,j];
           if i=St then goto skip2;
              E[O,Tt]:=2*Z[i,j];

Skip2:  if j=Ba+1 then goto skip3;
              E[O,O]:=E[O,O]+4*Y[i,j];
              E[O,Pp]:=2*Y[i,j];
           if j=1 then goto skip4;

Skip3:  E[O,O]:=E[O,O]+4*Y[i,j-1];
              E[O,Rr]:=2*Y[i,j-1];
Skip4: end;
       end;

       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba+1 do begin
               O:=(i-1)*(Ba+1)+j;
             Pp:=(Ba+1)*St+i-1;
             Qq:=Pp+1;
             Rr:=Qq+1;
             if i>1 then goto Skip5;
             E[O,Qq]:=-6*Z[i,j]/X[i];
             goto Skip6;
Skip5:   E[O,Qq]:=6*Z[i-1,j]/X[i-1]-6*Z[i,j]/X[i];
             E[O,Pp]:=-6*Z[i-1,j]/X[i-1];
             if i=St then goto Skip7;
Skip6: E[O,Rr]:=6*Z[i,j]/X[i];
Skip7: end;
       end;

       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba+1 do begin
               O:=St*(Ba+1)+i;
             Pp:=(i-1)*(Ba+1)+j;
             Qq:=Pp+Ba+1;
             Rr:=Pp-Ba-1;
             E[O,Pp]:=-6*Z[i,j]/X[i];
             if i=1 then goto Skip8;
             E[O,Pp]:=E[O,Pp]+6*Z[i-1,j]/X[i-1];
             E[O,Rr]:=6*Z[i-1,j]/X[i-1];
             if i=St then goto Skip9;
Skip8: E[O,Qq]:=-6*Z[i,j]/X[i];
Skip9: end;
       end;

       FOR i:= 1 TO St do begin
               O:=St*(Ba+1)+i;
             if i=1 then goto Skip10;
        FOR j:= 1 TO Ba+1 do begin
          E[O,O]:=E[O,O]+12*Z[i,j]/SQR(X[i])+12*Z[i-1,j]/SQR(X[i-1]);
          E[O,O-1]:=E[O,O-1]-12*Z[i-1,j]/SQR(X[i-1]);
        end;
             goto Skip11;
Skip10: FOR j:= 1 TO Ba+1 do begin
          E[O,O]:=E[O,O]+12*Z[i,j]/SQR(X[i]);
        end;

Skip11:  if i=St then goto skip12;
        FOR j:= 1 TO Ba+1 do begin
          E[O,O+1]:=E[O,O+1]-12*Z[i,j]/SQR(X[i]);
        end;
Skip12: end;

  End;

   Procedure RisoluzioneEquazCaricoSpostamento;

   Label Salta1,Salta2,Salta3,Salta4;
   Var   Hh,n,Aa,Ll: Integer;
         Kk,Gg,Xx,Zz   : Real;
   Begin
         Hh:=St*(Ba+2);
       FOR n:= 1 TO Hh do begin
           Kk:=0;
          FOR i:= n TO Hh do begin
          if abs(E[i,n]) < Kk then goto Salta1;
           Kk:=abs(E[i,n]);
           Ll:=i;
Salta1:   end;
          FOR j:= n TO Hh do begin
             Gg:=E[Ll,j];
             E[Ll,j]:=E[n,j];
             E[n,j]:=Gg;
          end;
             Xx:=P[Ll];
             P[Ll]:=P[n];
             P[n]:=Xx;
          FOR i:= n TO Hh do begin
          if E[i,n]=0 then goto Salta2;
             Zz:=E[i,n];
          FOR j:= n TO Hh do begin
             E[i,j]:=E[i,j]/Zz;
          end;
             P[i]:=P[i]/Zz;
Salta2:   end;
          if n=Hh then goto Salta3;
          FOR i:= n+1 TO Hh do begin
          if E[i,n]=0 then goto Salta4;
          FOR j:= n TO Hh do begin
             E[i,j]:=E[i,j]-E[n,j];
          end;
             P[i]:=P[i]-P[n];
Salta4:   end;
       end;

Salta3:   FOR Aa:= 2 TO Hh do begin
              i:=Hh-Aa+1;
            FOR j:= i+1 TO Hh do begin
             P[i]:=P[i]-E[i,j]*P[j];
            end;
          end;
  End;

   Procedure MomentiFinali;
   Label Jump1,Jump2,Jump3,Jump4;
   Begin
       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba do begin
          R[i,j]:=R[i,j]+4*Y[i,j]*P[(i-1)*(Ba+1)+j];
          R[i,j]:=R[i,j]+2*Y[i,j]*P[(i-1)*(Ba+1)+j+1];
            T[i,j+1]:=T[i,j+1]+2*Y[i,j]*P[(i-1)*(Ba+1)+j];
            T[i,j+1]:=T[i,j+1]+4*Y[i,j]*P[(i-1)*(Ba+1)+j+1];
        end;
       end;

       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba+1 do begin
             if i=St then goto Jump1;
          Q[i+1,j]:=4*Z[i,j]*P[i*(Ba+1)+j];
          Q[i+1,j]:=Q[i+1,j]+2*Z[i,j]*P[(i-1)*(Ba+1)+j];
          Q[i+1,j]:=Q[i+1,j]-6*Z[i,j]*P[St*(Ba+1)+i]/X[i];
          Q[i+1,j]:=Q[i+1,j]+6*Z[i,j]*P[St*(Ba+1)+i+1]/X[i];
            goto Jump2;
Jump1:  Q[i+1,j]:=2*Z[i,j]*P[(i-1)*(Ba+1)+j];
        Q[i+1,j]:=Q[i+1,j]-6*Z[i,j]*P[St*(Ba+1)+i]/X[i];
            goto Jump3;
Jump2:  S[i,j]:=2*Z[i,j]*P[i*(Ba+1)+j]+4*Z[i,j]*P[(i-1)*(Ba+1)+j];
        S[i,j]:=S[i,j]-6*Z[i,j]*(P[St*(Ba+1)+i]-P[St*(Ba+1)+i+1])/X[i];
            goto Jump4;
Jump3:  S[i,j]:=4*Z[i,j]*P[(i-1)*(Ba+1)+j]-6*Z[i,j]*P[St*(Ba+1)+i]/X[i];
Jump4:  end;
       end;
  End;


   Procedure RISULTATIFinali;
   Begin

     ASSIGN(ft,wo+'.RIS');
          Rewrite(ft);

    writeln('    RIEPILOGO Dati Telaio: '+wo+' -> Campate:',Ba:2,'  Piani:',St:2);
    write('         Luci:');
       FOR k:= 1 TO Ba do begin write(W[k]:7:3); end; writeln;
    write('      Altezze:');
       FOR k:= 1 TO St do begin write(X[k]:7:3); end; writeln;
    write('      J_travi:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba do begin write(L[i,j]:6:2); end; end; writeln;
    write('      J_pilas:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba+1 do begin write(M[i,j]:6:2); end; end; writeln;
    write(' Caric_Vert q:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba do begin write(U[i,j]:7:3); end; end; writeln;
    write('Spinte Oriz H:');
       FOR i:= 1 TO St do begin write(V[i]:7:3); end; writeln;

{---------------------------- Scrittura su file ----------------------------}
    writeln(ft,'    RIEPILOGO Dati Telaio: '+wo+' -> Campate:',Ba:2,'  Piani:',St:2);
    write(ft,'         Luci:');
       FOR k:= 1 TO Ba do begin write(ft,W[k]:7:3); end; writeln(ft);
    write(ft,'      Altezze:');
       FOR k:= 1 TO St do begin write(ft,X[k]:7:3); end; writeln(ft);
    write(ft,'      J_travi:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba do begin write(ft,L[i,j]:6:2); end; end; writeln(ft);
    write(ft,'      J_pilas:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba+1 do begin write(ft,M[i,j]:6:2); end; end; writeln(ft);
    write(ft,' Caric_Vert q:');
       FOR i:= 1 TO St do begin
            FOR j:= 1 TO Ba do begin write(ft,U[i,j]:7:3); end; end; writeln(ft);
    write(ft,'Spinte Oriz H:');
       FOR i:= 1 TO St do begin write(ft,V[i]:7:3); end; writeln(ft);

    writeln('                            Risultati FINALI');
    writeln('              MOMENTI Flettenti ai NODI (Senso Orario POSITIVO, Numeraz. NODI:)');
    writeln('              -------------------------  Sinistra-Destra, Alto-Basso)');
    writeln(ft,'                            Risultati FINALI');
    writeln(ft,'              MOMENTI Flettenti ai NODI (Senso Orario POSITIVO, Numeraz. NODI:)');
    writeln(ft,'              -------------------------  Sinistra-Destra, Alto-Basso)');
    writeln('           NODI        Sin.        Alto        Des.        Basso');
    writeln(ft,'           NODI        Sin.        Alto        Des.        Basso');
       FOR i:= 1 TO St do begin
        FOR j:= 1 TO Ba+1 do begin
      writeln('         [',i:2,',',j:2,']',t[i,j]:12:3,q[i,j]:12:3,R[i,j]:12:3,s[i,j]:12:3);
      writeln(ft,'         [',i:2,',',j:2,']',t[i,j]:12:3,q[i,j]:12:3,R[i,j]:12:3,s[i,j]:12:3);
        end;
       end;
    writeln('           Momenti d''incastro ai pilastri di base');
    writeln(ft,'           Momenti d''incastro ai pilastri di base');
       FOR i:= St+1 TO St+1 do begin
        FOR j:= 1 TO Ba+1 do begin
            writeln('             ',j:2,'       ',q[i,j]:12:3);
         writeln(ft,'             ',j:2,'       ',q[i,j]:12:3);
        end;
       end;
   Close (ft);
  End;
{------------------------------------------------------------------------------
  Inizio Programma }

  Begin
     ColoriSchermo;
     MainFrame;
     DATItelaio;
     MatriceRigidezze;
     RisoluzioneEquazCaricoSpostamento;
     MomentiFinali;
     RISULTATIFinali;
 end.

  { FINE Programma }
