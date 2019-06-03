  Program LAMB_XYZ; { PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI
                                (Metodo BATE)
 .....  Dati in coordinate cartesiane
             Determinazione degli elementi ORBITALI di una CONICA, NOTI:
               1   Vettori posizione r1(x1,y1,z1) e r2(x2,y2,z2)
               2 - Il tempo di Percorrenza ç=t2-t1
                L'angolo di trasferimento é Š calcolato AUTOMATICAMENTE
                tramite le coordinate dei 2 punti }
  Uses  Crt;

  Const
         kgauss=0.01720209895;  uv0=29.784691694; rad=pi/180;
  Var
         segno    : Integer;
         x1,y1,z1,x2,y2,z2, Vx1,Vy1,Vz1,Vx2,Vy2,Vz2,
         r1,r2,tetagradi,teta,tau,tgiorni,
         cot,Ag, x, y, zit,z,zz, S,C,Spri,Cpri,
         f,g,f1,g1, ff,gg,ff1,gg1,
         V1,V2, sigma1, fi1, cosfi,
         E1,t1menoTp,FH1,
         a,e,p,q,Triv,amin       : Real;

   Function ArcoCOS(valore:Real): Real;
     Begin
     ArcoCOS:=pi/2-arctan(valore/SQRT(1-valore*valore));
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
     textcolor(14);textbackground(1);clrscr;
   End;

  Procedure LogoDatiLAMBERT;
   Begin
     textcolor(10);
     gotoxy(7,2); writeln('PROBLEMA di LAMBERT risolto con le VARIABILI UNIVERSALI (Metodo BATE)');
     gotoxy(7,3); writeln('------------------------ in coordinate cartesiane -------------------');
     textcolor(14);
     gotoxy(05,5); write('Vett_Pos  r1  x1 [UA] = '); read(x1);
     gotoxy(05,6); write('              y1 [UA] = '); read(y1);
     gotoxy(05,7); write('              z1 [UA] = '); read(z1);

     gotoxy(43,5); write('Vett_Pos  r2  x2 [UA] = '); read(x2);
     gotoxy(43,6); write('              y2 [UA] = '); read(y2);
     gotoxy(43,7); write('              z2 [UA] = '); read(z2);

     gotoxy(25,8); write(' Tempo di volo   ç [giorni] = '); read(tgiorni);
   End;

  Procedure LAMBERT;
  Label Salta;
   Begin
      r1:=SQRT(SQR(x1)+SQR(y1)+SQR(z1));  r2:=SQRT(SQR(x2)+SQR(y2)+SQR(z2));
      teta:=arcoCOS((x1*x2+y1*y2+z1*z2)/(r1*r2));
      if (x1*y2-x2*y1)<0 then teta:=2*pi-teta;
     gotoxy(25,9); write(' Angolo Trasfer. é [ gradi] = ',teta/rad:9:5);
    gotoxy(25,10); write(' r1= ',r1:11:8,' UA  r2= ',r2:11:8,' UA');
       tau:=tgiorni*kgauss;
          if teta<pi then segno:=1 ELSE segno:=-1;
     Cot:=COS(teta);
     amin:=(r1+r2+SQRT(r1*r1+r2*r2-2*r1*r2*cot))/4;
     Ag:=segno*SQRT(r1*r2*(1+Cot));
     textcolor (15);
     z:=20;  { Valore iniziale di z=xý/a = Parametro adimensionale }
repeat
    {Funzioni Coniche -> S,C,S',C'}
  S:=1/6*((((((((z/342-1)*z/272+1)*z/210-1)*z/156+1)*z/110-1)*z/72+1)*z/42-1)*z/20+1);
  C:=1/2*((((((((z/306-1)*z/240+1)*z/182-1)*z/132+1)*z/90-1)*z/56+1)*z/30-1)*z/12+1);
  Spri:=(C-3*S)/(2*z);
  Cpri:=(1-z*S-2*C)/(2*z);
     { Anomalia eccentrica x e parametro y }
  y:=r1+r2-Ag*(1-z*S)/SQRT(C);
 if y<0 then begin z:=z/10; goto Salta end;  {...per evitare eventuali OVERFLOW}

  x:=SQRT(y/C);
    F:=x*x*x*S+Ag*SQRT(y)-tau;
    F1:=x*x*x*(Spri-3*S*Cpri/(2*C))+Ag/8 *(3*S*SQRT(y)/C+Ag/x);

     zit:= z - F / F1;   { Iterazione di NEWTON }

   zz:=ABS(zit-z);
   z:=zit;
     gotoxy(25,4); write('    z = ',zit:10:8,'   x =',x:10:8);
salta:
until zz < 0.00000001;
textcolor (14);
       { Semiasse maggiore a, parametro focale p, eccentricit… e }
       a:=x*x/z;  p:=r1*r2/y*(1-cot); e:=SQRT(1-p/a); q:=a*(1-e);

         { Velocit… agli estremi dell'arco }
       V1:=SQRT(2/r1-1/a);  V2:=SQRT(2/r2-1/a);

         { Parametri di Lagrange: f,g,f',g'}
      ff:=1-y/r1;  gg:=Ag*SQRT(y);
      gg1:=1-y/r2; ff1:=(ff*gg1-1)/gg;

         Vx1:=(x2-ff*x1)/gg;  Vy1:=(y2-ff*y1)/gg;   Vz1:=(z2-ff*z1)/gg;
         Vx2:=(gg1*x2-x1)/gg;  Vy2:=(gg1*y2-y1)/gg;   Vz2:=(gg1*z2-z1)/gg;

                  { å1 = r1ùV1 = prodotto scalare
      sigma1:=(r1*r2*cot-ff*r1*r1)/gg; formula uguale a quella di sotto }

      sigma1:=(gg+r1*x*(z*S-1))/(x*x*C); {in coordinate cartesiane conviene
                                          usare: å1=x1ùVx1+y1ùVy1+z1ùVz1}
      {sigma1:=x1*Vx1+y1*Vy1+z1*Vz1;}

     textcolor(1); gotoxy(32,11);writeln('RISULTATI FINALI'); textcolor(14);
     writeln('                     Semiasse magg.  a = ',a:11:8,' UA   [a_MIN = ',amin:9:6,']');
     writeln('                     Eccentricit…    e = ',e:11:8);
     writeln;
     writeln('     Velocit…_1     V1 = ',V1:11:8,' UA/UT0 = ',V1*uv0:6:5,' km/s');
     writeln('                                     V1x= ',Vx1:11:8,' = ',Vx1*uv0:6:5);
     writeln('                                     V1y= ',Vy1:11:8,' = ',Vy1*uv0:6:5);
     writeln('                                     V1z= ',Vz1:11:8,' = ',Vz1*uv0:6:5);
     writeln;
     writeln('     Velocit…_2     V2 = ',V2:11:8,' UA/UT0 = ',V2*uv0:6:5,' km/s');
     writeln('                                     V2x= ',Vx2:11:8,' = ',Vx2*uv0:6:5);
     writeln('                                     V2y= ',Vy2:11:8,' = ',Vy2*uv0:6:5);
     writeln('                                     V2z= ',Vz2:11:8,' = ',Vz2*uv0:6:5);

    End;
{-----------------------------------------------------------------------------------------------
   INIZIO Programma}

   Begin
       ColoriSchermo;
       LogoDatiLAMBERT;
       LAMBERT;
       End.
   {FINE Programma}
