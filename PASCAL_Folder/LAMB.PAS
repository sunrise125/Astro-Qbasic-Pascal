  Program LAMB;
             { PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)
             Determinazione degli elementi ORBITALI di una CONICA, NOTI:
               1   Due RAGGI Vettori r1,r2
               2 - L'angolo di trasferimento �
               3 - Il tempo di Percorrenza �=t2-t1
                   }
  Uses  Graph,Crt;

    Const
         kgauss=0.01720209895;
         uv0=29.784691694;
         rad=pi/180;
   Var
         tr       : String;
         segno    : Integer;

         r1,r2,tetagradi,teta,tau,tgiorni,
         cot,Ag, x, y, zit,z,zz, S,C,Spri,Cpri,
                          f,g,f1,g1,
         ff,gg,ff1,gg1,
         V1,V2, sigma1, fi1, cosfi,
         E1,t1menoTp,FH1,
         a,e,p,q,Triv,amin       : Real;


   Function ArcoCOS(valore:Real): Real;
     Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(abs(1-SQR(valore))));
     End;

      Function ArcoTANG(num,den:Real) : Real;
   Var atan : Real;
   Begin
     atan:=arctan(num/den);
      if den<0 then atan:=atan+pi;
      if atan<0 then atan:=atan+2*pi;
     ArcoTANG:=atan;
  End;

  Procedure ColoriSchermo;
   Begin
     textcolor(14);textbackground(6);clrscr;
  End;

  Procedure LogoDatiLAMBERT;
   Begin
     textcolor(10);
     gotoxy(7,2); writeln('PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)');
     gotoxy(7,3); writeln('---------------------------------------------------------------------');
     textcolor(14);
     gotoxy(25,5); write(' Raggio Vettore r1     [UA] = '); read(r1);
     gotoxy(25,6); write(' Raggio Vettore r2     [UA] = '); read(r2);
     gotoxy(25,7); write(' Angolo Trasfer. � [ gradi] = '); read(tetagradi);
     gotoxy(25,8); write(' Tempo di volo   � [giorni] = '); read(tgiorni);
  End;


  Procedure LAMBERT;
  Label Salta;
   Begin
       teta:=tetagradi*rad; tau:=tgiorni*kgauss;
          if tetagradi<=180 then segno:=1 ELSE segno:=-1;
     Cot:=COS(teta);
     amin:=(r1+r2+SQRT(r1*r1+r2*r2-2*r1*r2*cot))/4;
     Ag:=segno*SQRT(r1*r2*(1+Cot));
textcolor (10);

     z:=1;  { Valore iniziale di z=x�/a = Parametro adimensionale }
repeat
    {Funzioni Coniche -> S,C,S',C'}
  S:=1/6*((((((((z/342-1)*z/272+1)*z/210-1)*z/156+1)*z/110-1)*z/72+1)*z/42-1)*z/20+1);
  C:=1/2*((((((((z/306-1)*z/240+1)*z/182-1)*z/132+1)*z/90-1)*z/56+1)*z/30-1)*z/12+1);
  Spri:=(C-3*S)/(2*z);
  Cpri:=(1-z*S-2*C)/(2*z);
     { Anomalia eccentrica x e parametro y }
  y:=r1+r2-Ag*(1-z*S)/SQRT(C);
     if y<0 then
        begin
          z:=z/10;
        goto Salta
     end;

  x:=SQRT(y/C);
    F:=x*x*x*S+Ag*SQRT(y)-tau;
    F1:=x*x*x*(Spri-3*S*Cpri/(2*C))+Ag/8 *(3*S*SQRT(y)/C+Ag/x);

     zit:= z - F / F1;   { Iterazione di NEWTON }

   zz:=ABS(zit-z);
   z:=zit;
     gotoxy(25,10); write('    z = ',zit:10:8,'   x =',x:10:8);
salta:
until zz < 0.00000001;
textcolor (14);
       { Semiasse maggiore a, parametro focale p, eccentricit� e }
       a:=x*x/z;  p:=r1*r2/y*(1-cot); e:=SQRT(1-p/a); q:=a*(1-e);

         { Velocit� agli estremi dell'arco }
       V1:=SQRT(2/r1-1/a);  V2:=SQRT(2/r2-1/a);

         { Parametri di Lagrange: f,g,f',g'}
      ff:=1-y/r1;  gg:=Ag*SQRT(y);
      gg1:=1-y/r2; ff1:=(ff*gg1-1)/gg;

                  { �1 = r1�V1 = prodotto scalare
      sigma1:=(r1*r2*cot-ff*r1*r1)/gg; formula uguale a quella di sotto }
      sigma1:=(gg+r1*x*(z*S-1))/(x*x*C);

         fi1:= ArcoCOS(1/e*(p/r1-1));
         cosfi:= 1/e*(p/r1-1);
            if sigma1<0 then fi1:=2*pi-fi1;

             if e <= 0.9995 then
              begin
               E1:=ArcoCOS(1/e*(1-r1/a));
               t1menoTp:=a*SQRT(a)*(E1-sigma1/SQRT(a));
               tr:='Orbita ELLITTICA';
               Triv:=a*SQRT(a)*2*pi/kgauss;
              end
                ELSE if (e < 1.0005) and (e > 0.9995)
                 then
                  begin
                   t1menoTp:=0.5*(p*sigma1+sigma1*sigma1*sigma1/3);
                    tr:='Orbita PARABOLICA';
                  end
                 ELSE
                   begin
                   FH1:=ln(1/e*(sigma1/SQRT(-a)+(1-r1/a)));
                    t1menoTp:=-a*SQRT(-a)*(sigma1/SQRT(-a)-FH1);
                    tr:='Orbita IPERBOLICA';
                   end;

     gotoxy(32,12);write('RISULTATI FINALI');
     gotoxy(32,13);write('================');
     writeln;
     writeln('                     Semiasse magg.  a = ',a:11:8,' UA   [a_MIN = ',amin:9:6,']');
     writeln('                     Eccentricit�    e = ',e:11:8);
if e <= 0.9995 then writeln('                     Periodo Rivol.  P = ',Triv:11:6,' giorni');
     writeln('                     Param.focale    p = ',p:11:8,' UA');
     writeln('                     Dist.perielia   q = ',q:11:8,' UA');
     writeln('                     Velocit�_1     V1 = ',V1:11:8,' UA/UT0 = ',V1*uv0:6:5,' km/s');
     writeln('                     Velocit�_2     V2 = ',V2:11:8,' UA/UT0 = ',V2*uv0:6:5,' km/s');
     writeln('                     Anomalia_P1    �1 = ',fi1/rad:11:6,'�');
     writeln('                     Pass.Perie. t1-Tp = ',t1menoTp/kgauss:11:6,' giorni');
     textcolor(1);
     writeln('                                '+tr); textcolor(14);

    End;

{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiLAMBERT;
       LAMBERT;
       End.
   {FINE Programma}
